"""
Analysis endpoints. Execution (start/status/active/stop) is Nextflow-backed. The
read-only interactive endpoints (coverage/evidence/genome_features/consensus/annotations)
are backend-agnostic file readers — TODO: port from the old app's analysis.py, since they
just read RESULTS_DIR/<job_id>/ artifacts (the pipeline now publishes the same layout).
"""
import csv
import glob
import json
import logging
import os
import uuid
from collections import defaultdict
from typing import List, Optional

from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel

from app.config import settings
from app.models import FeatureEdit
from app import runner

router = APIRouter()
logger = logging.getLogger(__name__)


def _strip_ver(acc: str) -> str:
    """'BK006935.2' -> 'BK006935' (numeric suffix only)."""
    dot = acc.rfind(".")
    return acc[:dot] if dot > 0 and acc[dot + 1:].isdigit() else acc


class AnalysisRequest(BaseModel):
    file_ids: List[str]
    reference_genome_id: Optional[str] = None
    description: Optional[str] = None
    options: Optional[dict] = None


@router.post("/start")
async def start_analysis(request: AnalysisRequest):
    read_paths = []
    for fid in request.file_ids:
        p = runner.resolve_upload(fid)
        if not p:
            raise HTTPException(status_code=404, detail=f"File with ID {fid} not found")
        read_paths.append(p)

    if not request.reference_genome_id:
        raise HTTPException(status_code=400, detail="reference_genome_id is required")
    reference = runner.resolve_reference(request.reference_genome_id)
    if not reference:
        raise HTTPException(status_code=404, detail="Reference genome file not found")

    job_id = f"job_{request.file_ids[0]}_{str(uuid.uuid4())[:8]}"
    try:
        info = runner.start_job(job_id, read_paths, reference, request.options or {})
    except RuntimeError as e:
        raise HTTPException(status_code=500, detail=str(e))

    return {"job_id": job_id, "task_id": str(info["pid"]),
            "status": "queued",
            "message": f"Analysis job has been launched via {'Seqera Platform' if settings.LAUNCH_MODE == 'seqera' else 'Nextflow'}."}


@router.get("/status/{job_id}")
async def get_status(job_id: str):
    """Prefer results.json on disk (authoritative); fall back to the process state."""
    job_id = job_id.strip()   # tolerate stray whitespace from the "load existing job" field
    job_dir = settings.RESULTS_DIR / job_id
    results_file = job_dir / "results.json"

    if results_file.exists():
        try:
            data = json.loads(results_file.read_text() or "{}")
            if data:
                out = {"job_id": job_id, "status": data.get("status", "completed")}
                if data.get("status") == "failed" and "error" in data:
                    out["error"] = data["error"]
                if "steps" in data:
                    out["steps"] = data["steps"]
                return out
        except json.JSONDecodeError:
            pass  # being written; fall through

    proc_state = runner.job_process_state(job_id)
    if proc_state == "running":
        return {"job_id": job_id, "status": "processing"}
    if proc_state == "exited_error":
        return {"job_id": job_id, "status": "failed",
                "error": runner.tail_log(job_id) or "Nextflow exited non-zero."}
    if job_dir.exists():
        return {"job_id": job_id, "status": "processing"}
    return {"job_id": job_id, "status": "unknown"}


@router.get("/results/{job_id}")
async def get_results(job_id: str):
    job_id = job_id.strip()
    results_file = settings.RESULTS_DIR / job_id / "results.json"
    if not results_file.exists():
        raise HTTPException(status_code=404, detail="Results not found or analysis not complete")
    try:
        data = json.loads(results_file.read_text())
    except json.JSONDecodeError:
        return {"job_id": job_id, "status": "processing",
                "message": "Results file is currently being written."}

    # Normalise seq_id version suffixes (BK006935.2 -> BK006935), matching the old app.
    for m in data.get("mutations", []) or []:
        sid = m.get("seq_id", "")
        dot = sid.rfind(".")
        if dot > 0 and sid[dot + 1:].isdigit():
            m["seq_id"] = sid[:dot]
    return data


@router.get("/active_jobs")
async def active_jobs():
    return runner.list_active()


@router.post("/stop/{job_id}")
async def stop_job(job_id: str):
    if not runner.stop_job(job_id):
        raise HTTPException(status_code=404, detail="Job not found active")
    if settings.LAUNCH_MODE == "seqera":
        return {"message": f"Job {job_id}: cancellation requested in Seqera Platform."}
    return {"message": f"Job {job_id} terminated (Nextflow process + submitted Slurm jobs)."}


# ── Read-only viewers (ported from the old app; read RESULTS_DIR/<job_id>/ artifacts) ──

@router.get("/coverage/{job_id}")
async def get_coverage(job_id: str, resolution: int = 1000, chrom: str = Query(None)):
    """Downsampled per-base coverage for plotting, cached per (chrom, resolution)."""
    job_id = job_id.strip()
    job_dir = settings.RESULTS_DIR / job_id
    coverage_file = job_dir / "coverage.tsv"
    if not coverage_file.exists():
        raise HTTPException(status_code=404, detail="Coverage data not found")

    cache_file = job_dir / f"coverage_{chrom or 'all'}_{resolution}.json"
    if cache_file.exists():
        try:
            return json.loads(cache_file.read_text())
        except Exception:
            pass  # recompute

    def _downsample(points):
        step = max(1, len(points) // resolution)
        return [points[i] for i in range(0, len(points), step)]

    # chrom=None: load once, precompute global + per-chrom caches (fast chrom switching).
    if chrom is None:
        all_points, chrom_points = [], defaultdict(list)
        with open(coverage_file) as f:
            for line in f:
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                try:
                    pt = {"chrom": _strip_ver(parts[0]), "pos": int(parts[1]),
                          "depth": int(parts[2].strip())}
                except ValueError:
                    continue
                all_points.append(pt)
                chrom_points[pt["chrom"]].append(pt)

        result = {"job_id": job_id, "total_bases": len(all_points), "data": _downsample(all_points)}
        cache_file.write_text(json.dumps(result))
        for c_name, pts in chrom_points.items():
            (job_dir / f"coverage_{c_name}_{resolution}.json").write_text(
                json.dumps({"job_id": job_id, "total_bases": len(pts), "data": _downsample(pts)}))
        return result

    # Specific chromosome.
    filtered = []
    with open(coverage_file) as f:
        for line in f:
            parts = line.split("\t")
            if len(parts) >= 3 and _strip_ver(parts[0]) == chrom:
                try:
                    filtered.append({"chrom": chrom, "pos": int(parts[1]), "depth": int(parts[2].strip())})
                except ValueError:
                    continue
    result = {"job_id": job_id, "total_bases": len(filtered), "data": _downsample(filtered)}
    cache_file.write_text(json.dumps(result))
    return result


@router.get("/genome_features/{job_id}")
async def get_genome_features(job_id: str):
    """Gene features (with exon `parts`) from the reference GBK, plus any saved edits."""
    job_id = job_id.strip()
    from Bio import SeqIO  # heavy import kept local so the service starts without BioPython

    job_dir = settings.RESULTS_DIR / job_id
    results_file = job_dir / "results.json"

    # Locate a GBK reference: prefer the one named in results.json, then any GBK in UPLOAD_DIR.
    candidates = []
    if results_file.exists():
        try:
            ref_name = json.loads(results_file.read_text()).get("reference_file")
        except Exception:
            ref_name = None
        if ref_name:
            for d in (settings.REFERENCE_DIR, settings.UPLOAD_DIR):
                p = d / ref_name
                if p.exists():
                    candidates.append(p)
                if not ref_name.endswith((".gbk", ".gb")):   # FASTA used → sibling GBK
                    base = os.path.splitext(ref_name)[0]
                    for ext in (".gbk", ".gb"):
                        gp = d / (base + ext)
                        if gp.exists():
                            candidates.append(gp)
    for d in (settings.REFERENCE_DIR, settings.UPLOAD_DIR):
        if d.exists():
            candidates += [p for p in d.iterdir()
                           if p.suffix in (".gbk", ".gb") and p not in candidates]

    features, selected, total_length = [], None, 0
    for gbk in candidates:
        try:
            records = list(SeqIO.parse(str(gbk), "genbank"))
        except Exception as e:
            logger.warning("Failed to parse %s: %s", gbk, e)
            continue
        temp, length = [], 0
        for record in records:
            length += len(record.seq)
            cds_map = {}
            for feat in record.features:
                if feat.type == "CDS":
                    key = (feat.qualifiers.get("locus_tag", [""])[0]
                           or feat.qualifiers.get("gene", [""])[0])
                    if key:
                        cds_map[key] = feat.location
            for feat in record.features:
                if feat.type != "gene":
                    continue
                try:
                    start, end = int(feat.location.start), int(feat.location.end)
                    gene_name = feat.qualifiers.get("gene", [""])[0]
                    locus_tag = feat.qualifiers.get("locus_tag", [""])[0]
                    name = gene_name or locus_tag
                    loc = cds_map.get(locus_tag) or cds_map.get(gene_name) or feat.location
                    parts = ([{"start": int(p.start), "end": int(p.end)} for p in loc.parts]
                             if hasattr(loc, "parts") else [{"start": start, "end": end}])
                    temp.append({
                        "original_id": locus_tag or name, "start": start, "end": end,
                        "parts": parts, "strand": feat.location.strand, "name": name,
                        "locus_tag": locus_tag, "product": feat.qualifiers.get("product", [""])[0],
                        "seq_id": record.name,
                    })
                except Exception:
                    continue
        if temp:
            features, selected, total_length = temp, gbk, length
            break

    if not selected:
        raise HTTPException(status_code=404,
                            detail="No reference GBK with gene features found for this job")

    # Apply saved annotation edits (annotations.json), if present.
    edits_path = job_dir / "annotations.json"
    if edits_path.exists():
        try:
            edits = json.loads(edits_path.read_text())
            kept = []
            for f in features:
                e = edits.get(f.get("locus_tag", "")) or edits.get(f.get("name", ""))
                if e and e.get("action") == "delete":
                    continue
                if e and e.get("action") == "edit":
                    if e.get("new_name"):      f["name"] = e["new_name"]
                    if e.get("new_locus_tag"): f["locus_tag"] = e["new_locus_tag"]
                    if e.get("new_product"):   f["product"] = e["new_product"]
                kept.append(f)
            for key, e in edits.items():
                if e.get("action") == "add":
                    kept.append({
                        "start": int(e.get("start", 0)), "end": int(e.get("end", 0)),
                        "strand": int(e.get("strand", 1)), "name": e.get("new_name", key),
                        "locus_tag": e.get("new_locus_tag", key),
                        "product": e.get("new_product", "Added Annotation"),
                        "seq_id": e.get("seq_id", features[0]["seq_id"] if features else "unknown"),
                    })
            features = sorted(kept, key=lambda x: x.get("start", 0))
        except Exception as e:
            logger.error("Error applying annotations: %s", e)

    return {"genome_id": selected.name, "length": total_length, "features": features}


# ── Annotation edits (persisted to annotations.json; applied by genome_features) ──────────

@router.post("/annotations/{job_id}")
async def save_annotations(job_id: str, edit: FeatureEdit):
    """Persist one edit/delete/add to the job's annotations.json (keyed by feature id)."""
    job_id = job_id.strip()
    job_dir = settings.RESULTS_DIR / job_id
    if not job_dir.exists():
        raise HTTPException(status_code=404, detail="Job not found")

    annotations_file = job_dir / "annotations.json"
    current = {}
    if annotations_file.exists():
        try:
            content = annotations_file.read_text()
            if content:
                current = json.loads(content)
        except Exception as e:
            logger.warning("Could not read existing annotations: %s", e)

    key = edit.original_id
    if not key and edit.action == "add":
        key = f"add_{edit.new_locus_tag or edit.new_name or 'unnamed'}"
    if not key:
        raise HTTPException(status_code=400, detail="Missing original_id for annotation edit")

    current[key] = edit.model_dump()
    annotations_file.write_text(json.dumps(current, indent=2))
    return {"message": "Annotation saved successfully", "edit": edit.model_dump()}


@router.get("/annotations/{job_id}")
async def get_annotations(job_id: str):
    """Return saved annotation edits for a job (or an empty set)."""
    job_id = job_id.strip()
    annotations_file = settings.RESULTS_DIR / job_id / "annotations.json"
    if annotations_file.exists():
        try:
            return json.loads(annotations_file.read_text())
        except Exception:
            return {"edits": []}
    return {"edits": []}
