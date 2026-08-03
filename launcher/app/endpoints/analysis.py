"""
Analysis endpoints. Execution (start/status/active/stop) is Nextflow-backed. The
read-only interactive endpoints (coverage/evidence/genome_features/consensus/annotations)
are backend-agnostic file readers — TODO: port from the old app's analysis.py, since they
just read RESULTS_DIR/<job_id>/ artifacts (the pipeline now publishes the same layout).
"""
import json
import uuid
from typing import List, Optional

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from app.config import settings
from app import runner

router = APIRouter()


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
    reference = runner.resolve_upload(request.reference_genome_id)
    if not reference:
        raise HTTPException(status_code=404, detail="Reference genome file not found")

    job_id = f"job_{request.file_ids[0]}_{str(uuid.uuid4())[:8]}"
    info = runner.start_job(job_id, read_paths, reference, request.options or {})

    return {"job_id": job_id, "task_id": str(info["pid"]),
            "status": "queued", "message": "Analysis job has been launched via Nextflow."}


@router.get("/status/{job_id}")
async def get_status(job_id: str):
    """Prefer results.json on disk (authoritative); fall back to the process state."""
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
    return {"message": f"Job {job_id} terminated (Nextflow process + submitted Slurm jobs)."}


# ── TODO: port these read-only viewers from the old app (they read RESULTS_DIR/<job_id>/) ──
_NOT_YET = "Not yet ported to the launcher — read-only viewer coming next."


@router.get("/coverage/{job_id}")
async def coverage(job_id: str):
    raise HTTPException(status_code=501, detail=_NOT_YET)


@router.get("/genome_features/{job_id}")
async def genome_features(job_id: str):
    raise HTTPException(status_code=501, detail=_NOT_YET)
