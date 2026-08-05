"""
Consensus-genome endpoints (ported from the old app's analysis.py).

- GET /analysis/consensus/{job_id}          → FASTA with variants applied
- GET /analysis/consensus/preview/{job_id}  → proposed annotation changes (GBK modal)
- POST /analysis/consensus/download/{job_id} → GBK with approved changes + feature shifts

Reads the reference + variants from RESULTS_DIR/<job_id>/ and applies mutations above a
frequency threshold. GATK VCF is located by glob (pipeline names it <sample>.vcf, not the
old app's fixed variants.vcf).
"""
import glob
import json
import logging
import os
import re
from typing import List

from fastapi import APIRouter, HTTPException, Response

from app.config import settings
from app.models import MutationProposal

router = APIRouter()
logger = logging.getLogger(__name__)


def _parse_vcf_mutations(vcf_path: str, freq_threshold: float = 0.9, min_qual: float = 0.0) -> list:
    """Parse a VCF into breseq-GD-compatible mutation dicts {type, seq_id, pos, parts, freq}."""
    mutations = []
    with open(vcf_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, ref, alt, qual_str = cols[0], int(cols[1]), cols[3], cols[4], cols[5]

            if min_qual > 0 and qual_str != ".":
                try:
                    if float(qual_str) < min_qual:
                        continue
                except ValueError:
                    pass

            if "," in alt:
                alt = alt.split(",")[0]

            freq = 1.0
            for field in cols[7].split(";"):
                if field.startswith("AF="):
                    try:
                        freq = float(field.split("=")[1].split(",")[0])
                    except Exception:
                        pass
            if freq < freq_threshold:
                continue

            if len(ref) == 1 and len(alt) == 1:
                m_type = "SNP"
                parts = ["SNP", ".", ".", chrom, str(pos), alt]
            elif len(alt) > len(ref) and alt.startswith(ref):
                m_type = "INS"
                parts = ["INS", ".", ".", chrom, str(pos), alt[len(ref):]]
            elif len(ref) > len(alt) and ref.startswith(alt):
                m_type = "DEL"
                del_size = len(ref) - len(alt)
                pos = pos + len(alt)
                parts = ["DEL", ".", ".", chrom, str(pos), str(del_size)]
            else:
                m_type = "SUB"
                parts = ["SUB", ".", ".", chrom, str(pos), str(len(ref)), alt]

            mutations.append({"type": m_type, "seq_id": chrom, "pos": pos, "parts": parts, "freq": freq})
    return mutations


def _find_reference_gbk(job_dir: str):
    """Find the reference GBK for a job (results.json reference_file → UPLOAD_DIR), else any GBK."""
    from Bio import SeqIO
    from Bio.Seq import UndefinedSequenceError

    upload = str(settings.UPLOAD_DIR)
    seqfiles = str(settings.DATA_DIR / "seq_files")
    results_file = os.path.join(job_dir, "results.json")
    candidates = []

    if os.path.exists(results_file):
        try:
            ref_name = json.load(open(results_file)).get("reference_file")
        except Exception:
            ref_name = None
        if ref_name:
            for p in (os.path.join(upload, ref_name), os.path.join(seqfiles, ref_name)):
                if os.path.exists(p):
                    candidates.append(p)
            if not any(candidates) and not ref_name.endswith((".gbk", ".gb")):
                base = os.path.splitext(ref_name)[0]
                p_gbk = os.path.join(upload, base + ".gbk")
                if os.path.exists(p_gbk):
                    candidates.append(p_gbk)

    for d in (upload, seqfiles):
        if os.path.isdir(d):
            for fname in os.listdir(d):
                if fname.endswith((".gbk", ".gb")):
                    p = os.path.join(d, fname)
                    if p not in candidates:
                        candidates.append(p)

    for path in candidates:
        try:
            with open(path) as f:
                rec = next(SeqIO.parse(f, "genbank"))
                _ = bytes(rec.seq)
                return path
        except (StopIteration, UndefinedSequenceError, Exception):
            continue
    return None


def _find_gatk_vcf(job_dir: str):
    """Pipeline names the VCF <sample>.vcf (not the old fixed variants.vcf) — glob for it."""
    hits = sorted(glob.glob(os.path.join(job_dir, "gatk", "*.vcf")))
    return hits[0] if hits else None


@router.get("/consensus/{job_id}")
async def generate_consensus(job_id: str, freq_threshold: float = 0.9, min_qual: float = 0.0):
    """Apply variants above `freq_threshold` to the reference and return the FASTA."""
    from io import StringIO
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    job_id = job_id.strip()
    job_dir = str(settings.RESULTS_DIR / job_id)
    results_file = os.path.join(job_dir, "results.json")
    metadata = json.load(open(results_file)) if os.path.exists(results_file) else {}
    variant_caller = metadata.get("variant_caller", "breseq")

    try:
        if variant_caller == "gatk":
            gbk_path = _find_reference_gbk(job_dir)
            if not gbk_path:
                raise HTTPException(status_code=404, detail="Reference genome file not found.")
            ref_records = SeqIO.to_dict(SeqIO.parse(gbk_path, "genbank"))
            vcf_path = _find_gatk_vcf(job_dir)
            if not vcf_path:
                raise HTTPException(status_code=404, detail="VCF file not found in GATK results.")
            mutations = _parse_vcf_mutations(vcf_path, freq_threshold, min_qual)
        else:
            ref_path = os.path.join(job_dir, "breseq", "data", "reference.fasta")
            gd_path = os.path.join(job_dir, "breseq", "output", "output.gd")
            if not os.path.exists(ref_path):
                raise HTTPException(status_code=404, detail="breseq reference.fasta not found.")
            if not os.path.exists(gd_path):
                raise HTTPException(status_code=404, detail="breseq output.gd not found.")
            ref_records = SeqIO.to_dict(SeqIO.parse(ref_path, "fasta"))
            mutations = []
            with open(gd_path) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    if parts[0] not in ("SNP", "INS", "DEL", "SUB") or len(parts) < 5:
                        continue
                    try:
                        pos = int(parts[4])
                    except ValueError:
                        continue
                    freq = 1.0
                    for p in parts[5:]:
                        if p.startswith(("frequency=", "polymorphism_frequency=")):
                            try:
                                freq = float(p.split("=")[1])
                            except Exception:
                                pass
                            break
                    if freq < freq_threshold:
                        continue
                    mutations.append({"type": parts[0], "seq_id": parts[3], "pos": pos, "parts": parts, "freq": freq})

        # Apply mutations position-DESCENDING so indels don't shift unprocessed coordinates.
        list_seqs = {k: list(str(v.seq)) for k, v in ref_records.items()}
        mutations.sort(key=lambda x: x["pos"], reverse=True)
        for m in mutations:
            sid = m["seq_id"]
            if sid not in list_seqs:
                if len(list_seqs) == 1:
                    sid = next(iter(list_seqs))
                else:
                    continue
            lseq = list_seqs[sid]
            pos_0 = m["pos"] - 1
            parts = m["parts"]
            if m["type"] == "SNP" and len(parts) > 5:
                nb = parts[5]
                if len(nb) == 1 and nb not in ("=", ":"):
                    lseq[pos_0] = nb
            elif m["type"] == "DEL" and len(parts) > 5:
                try:
                    lseq[pos_0:pos_0 + int(parts[5])] = []
                except ValueError:
                    pass
            elif m["type"] == "INS" and len(parts) > 5:
                if not parts[5].startswith("frequency="):
                    lseq[pos_0 + 1:pos_0 + 1] = list(parts[5])
            elif m["type"] == "SUB" and len(parts) > 6:
                try:
                    if not parts[6].startswith("frequency="):
                        lseq[pos_0:pos_0 + int(parts[5])] = list(parts[6])
                except ValueError:
                    pass

        new_records = []
        for sid, record in ref_records.items():
            if sid in list_seqs:
                new_records.append(SeqRecord(
                    Seq("".join(list_seqs[sid])), id=f"{sid}_mutated",
                    description=f"Consensus from job {job_id} (threshold {freq_threshold})"))
            else:
                new_records.append(record)

        out = StringIO()
        SeqIO.write(new_records, out, "fasta")
        return Response(content=out.getvalue(), media_type="text/plain")

    except HTTPException:
        raise
    except Exception as e:  # noqa: BLE001
        logger.error("Error generating consensus: %s", e)
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/consensus/preview/{job_id}", response_model=List[MutationProposal])
async def get_consensus_preview(job_id: str, freq_threshold: float = 0.0, min_qual: float = 0.0):
    """Propose mutation-based annotation changes (gene hit + AA/nt suffix, nearby features,
    coverage snippet) for the consensus-GBK modal."""
    from Bio import SeqIO

    job_id = job_id.strip()
    job_dir = str(settings.RESULTS_DIR / job_id)
    results_file = os.path.join(job_dir, "results.json")

    gbk_path = _find_reference_gbk(job_dir)
    if not gbk_path:
        raise HTTPException(status_code=404, detail="No valid reference GBK found.")

    variant_caller = "breseq"
    if os.path.exists(results_file):
        try:
            variant_caller = json.load(open(results_file)).get("variant_caller", "breseq")
        except Exception:
            pass

    mutations = []
    if variant_caller == "gatk":
        vcf_path = _find_gatk_vcf(job_dir)
        if not vcf_path:
            return []
        for m in _parse_vcf_mutations(vcf_path, freq_threshold, min_qual):
            parts, m_type = m["parts"], m["type"]
            new_base = parts[5] if m_type != "SUB" else (parts[6] if len(parts) > 6 else "")
            desc = {"SNP": f"?->{new_base}", "DEL": f"DEL({parts[5]}bp)",
                    "INS": f"INS({new_base})", "SUB": f"SUB({new_base})"}.get(m_type, m_type)
            mutations.append({"id": f"{m['seq_id']}_{m['pos']}_{m_type}", "seq_id": m["seq_id"],
                              "pos": m["pos"], "type": m_type, "freq": m["freq"], "desc": desc,
                              "extra": {}, "new_base": new_base})
    else:
        gd_path = os.path.join(job_dir, "breseq", "output", "evidence", "annotated.gd")
        if not os.path.exists(gd_path):
            gd_path = os.path.join(job_dir, "breseq", "output", "output.gd")
        if not os.path.exists(gd_path):
            return []

        def parse_gd_extra(slice_):
            extra = {}
            for p in slice_:
                if "=" in p:
                    k, v = p.split("=", 1)
                    extra[k] = v
            return extra

        gd_records, mutation_ids = {}, []
        with open(gd_path) as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.split("\t")
                m_type = parts[0]
                try:
                    uid = parts[1]
                    parents = parts[2].split(",") if parts[2] and parts[2] != "." else []
                    start_idx = 7 if m_type == "SUB" else 6
                    gd_records[uid] = {"type": m_type, "uid": uid, "parents": parents,
                                       "parts": parts, "extra": parse_gd_extra(parts[start_idx:])}
                    if m_type in ("SNP", "INS", "DEL", "SUB"):
                        mutation_ids.append(uid)
                except Exception:
                    continue

        for uid in mutation_ids:
            rec = gd_records[uid]
            m_type, parts, extra = rec["type"], rec["parts"], rec["extra"]
            merged_aa = []
            for pid in rec["parents"]:
                if pid in gd_records:
                    pe = gd_records[pid]["extra"]
                    if "aa_new_seq" in pe and "aa_ref_seq" in pe:
                        cs = f"{pe['aa_ref_seq']}{pe.get('aa_position', '?')}{pe['aa_new_seq']}"
                        if pe.get("gene_name"):
                            cs = f"{pe['gene_name']}:{cs}"
                        if cs not in merged_aa:
                            merged_aa.append(cs)
                    for k in ("aa_new_seq", "aa_ref_seq", "aa_position", "codon_new_seq",
                              "codon_ref_seq", "gene_name", "snp_type"):
                        if k in pe and k not in extra:
                            extra[k] = pe[k]
            if merged_aa:
                extra["derived_aa_changes"] = merged_aa
            try:
                seq_id, pos = parts[3], int(parts[4])
                new_base = (parts[6] if len(parts) > 6 else "") if m_type == "SUB" else parts[5]
                freq = float(extra["frequency"]) if "frequency" in extra else 1.0
                desc = {"SNP": f"{extra.get('ref_seq', '?')}->{new_base}",
                        "DEL": f"DEL({extra.get('position_end', '?')})",
                        "INS": f"INS({new_base})", "SUB": f"SUB({new_base})"}.get(m_type, m_type)
                mutations.append({"id": f"{seq_id}_{pos}_{m_type}", "seq_id": seq_id, "pos": pos,
                                  "type": m_type, "freq": freq, "desc": desc, "extra": extra,
                                  "new_base": new_base})
            except Exception:
                continue

    # Coverage context (+/-20 bp) from the job's coverage.tsv
    needed = set()
    for m in mutations:
        for p in range(m["pos"] - 20, m["pos"] + 21):
            needed.add((m["seq_id"], p))
    coverage_map = {}
    cov_file = os.path.join(job_dir, "coverage.tsv")
    if os.path.exists(cov_file):
        try:
            with open(cov_file) as f:
                for line in f:
                    parts = line.split("\t")
                    if len(parts) < 3:
                        continue
                    try:
                        key = (parts[0], int(parts[1]))
                        if key in needed:
                            coverage_map[key] = int(parts[2].strip())
                    except ValueError:
                        continue
        except Exception as e:
            logger.warning("Failed to load coverage context: %s", e)

    proposals = []
    ref_records_map = {r.id: r for r in SeqIO.parse(gbk_path, "genbank")}
    for r in list(ref_records_map.values()):
        ref_records_map.setdefault(r.name, r)

    for m in mutations:
        hit_feature = gene_name = locus_tag = prod = None
        rec = ref_records_map.get(m["seq_id"])
        if not rec and len(ref_records_map) == 1:
            rec = next(iter(ref_records_map.values()))
        nearby = []
        if rec:
            for feat in rec.features:
                if feat.type in ("gene", "CDS") and feat.location.start <= (m["pos"] - 1) < feat.location.end:
                    hit_feature = feat
                    gene_name = feat.qualifiers.get("gene", [""])[0]
                    locus_tag = feat.qualifiers.get("locus_tag", [""])[0]
                    prod = feat.qualifiers.get("product", [""])[0]
                    break
            ws, we = max(0, m["pos"] - 2000), m["pos"] + 2000
            for feat in rec.features:
                if feat.type not in ("gene", "CDS", "tRNA", "rRNA"):
                    continue
                fs, fe = int(feat.location.start), int(feat.location.end)
                if fe > ws and fs < we:
                    nearby.append({"name": feat.qualifiers.get("gene", [""])[0] or feat.qualifiers.get("locus_tag", [""])[0],
                                   "start": fs, "end": fe, "strand": feat.location.strand, "type": feat.type})

        original_name = gene_name or locus_tag or "unnamed"
        extra = m.get("extra", {})
        if extra.get("derived_aa_changes"):
            suffix = "_" + "_".join(extra["derived_aa_changes"])
        elif "aa_new_seq" in extra and "aa_ref_seq" in extra and "aa_position" in extra:
            suffix = f"_{extra['aa_ref_seq']}{extra['aa_position']}{extra['aa_new_seq']}"
        elif "codon_new_seq" in extra and "codon_ref_seq" in extra:
            suffix = f"_syn_{extra['codon_ref_seq']}2{extra['codon_new_seq']}"
        else:
            suffix = "_" + m["desc"].split(" ")[0].replace("->", "")
        suffix = re.sub(r"[^A-Za-z0-9_]", "", suffix)

        if hit_feature:
            proposed, is_coding = f"{original_name}{suffix}", True
        elif "gene_name" in extra:
            proposed, is_coding = f"{extra['gene_name']}{suffix}", True
        else:
            proposed, is_coding, original_name = f"{m['seq_id']}_{m['pos']}{suffix}", False, None

        cov_snippet = [coverage_map.get((m["seq_id"], p), 0) for p in range(m["pos"] - 20, m["pos"] + 21)]
        proposals.append(MutationProposal(
            mutation_id=m["id"], seq_id=m["seq_id"], position=m["pos"], type=m["type"],
            frequency=m["freq"], gene_name=gene_name, locus_tag=locus_tag, original_product=prod,
            mutation_detail=m["desc"], proposed_name=proposed, is_coding=is_coding,
            selected=m.get("freq", 1.0) >= 0.99, coverage_context=cov_snippet, nearby_features=nearby))

    return proposals


@router.post("/consensus/download/{job_id}")
async def download_consensus_gbk(job_id: str, approved_proposals: List[MutationProposal],
                                 freq_threshold: float = 0.9, min_qual: float = 0.0):
    """Generate a GBK with variants applied, features coordinate-shifted, and approved
    proposal renames + saved manual edits merged in."""
    from io import StringIO
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, ExactPosition

    job_id = job_id.strip()
    job_dir = str(settings.RESULTS_DIR / job_id)
    results_file = os.path.join(job_dir, "results.json")

    # Build rename plan from selected proposals.
    modification_plan, new_annotations = {}, []
    for p in approved_proposals:
        if not p.selected:
            continue
        if p.is_coding and (p.locus_tag or p.gene_name):
            base = p.gene_name or p.locus_tag or "unnamed"
            for k in [x for x in (p.locus_tag, p.gene_name) if x]:
                plan = modification_plan.setdefault(k, {"suffixes": set(), "custom_name": None})
                plan["base"] = base
                if p.proposed_name.startswith(base):
                    plan["suffixes"].add(p.proposed_name[len(base):])
                else:
                    plan["custom_name"] = p.proposed_name
        elif not p.is_coding:
            new_annotations.append(p)

    gbk_path = _find_reference_gbk(job_dir)
    if not gbk_path:
        raise HTTPException(status_code=404, detail="No valid reference GBK file found.")

    variant_caller = "breseq"
    if os.path.exists(results_file):
        try:
            variant_caller = json.load(open(results_file)).get("variant_caller", "breseq")
        except Exception:
            pass

    mutations = []
    if variant_caller == "gatk":
        vcf_path = _find_gatk_vcf(job_dir)
        if vcf_path:
            mutations = _parse_vcf_mutations(vcf_path, freq_threshold, min_qual)
    else:
        gd_path = os.path.join(job_dir, "breseq", "output", "output.gd")
        if os.path.exists(gd_path):
            with open(gd_path) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    if parts[0] not in ("SNP", "INS", "DEL", "SUB") or len(parts) < 5:
                        continue
                    try:
                        pos = int(parts[4])
                    except ValueError:
                        continue
                    freq = 1.0
                    for p in parts[5:]:
                        if p.startswith("frequency="):
                            try:
                                freq = float(p.split("=")[1])
                            except Exception:
                                pass
                    if freq >= freq_threshold:
                        mutations.append({"type": parts[0], "seq_id": parts[3], "pos": pos, "parts": parts, "freq": freq})

    saved_edits = {}
    annotations_file = os.path.join(job_dir, "annotations.json")
    if os.path.exists(annotations_file):
        try:
            data = json.load(open(annotations_file))
            if isinstance(data.get("edits"), list):
                for e in data["edits"]:
                    saved_edits[e.get("original_id")] = e
            else:
                saved_edits = data
        except Exception as e:
            logger.warning("Failed to load annotations: %s", e)

    ref_records = list(SeqIO.parse(gbk_path, "genbank"))
    new_records = []
    for rec in ref_records:
        try:
            seq_list = list(str(rec.seq))
        except Exception:
            continue
        current_id = rec.id
        rec_muts = [m for m in mutations if m["seq_id"] == current_id
                    or m["seq_id"] in current_id or current_id in m["seq_id"]]
        if not rec_muts and len(ref_records) == 1 and mutations:
            rec_muts = mutations

        for m in sorted(rec_muts, key=lambda x: x["pos"], reverse=True):
            pos_0, parts = m["pos"] - 1, m["parts"]
            try:
                if m["type"] == "SNP" and len(parts) > 5 and len(parts[5]) == 1:
                    seq_list[pos_0] = parts[5]
                elif m["type"] == "DEL" and len(parts) > 5:
                    seq_list[pos_0:pos_0 + int(parts[5])] = []
                elif m["type"] == "INS" and len(parts) > 5:
                    seq_list[pos_0 + 1:pos_0 + 1] = list(parts[5])
                elif m["type"] == "SUB" and len(parts) > 6:
                    seq_list[pos_0:pos_0 + int(parts[5])] = list(parts[6])
            except Exception:
                pass
        new_seq_obj = Seq("".join(seq_list))

        adjustments = []
        for m in sorted(rec_muts, key=lambda x: x["pos"]):
            delta = 0
            if m["type"] == "DEL":
                delta = -int(m["parts"][5])
            elif m["type"] == "INS":
                delta = len(m["parts"][5])
            elif m["type"] == "SUB":
                delta = len(m["parts"][6]) - int(m["parts"][5])
            if delta:
                adjustments.append({"pos": m["pos"], "delta": delta})

        def shift_pos(p):
            original_pos = int(p)
            total = 0
            for adj in adjustments:
                mut_start_0 = adj["pos"] - 1
                if mut_start_0 < original_pos:
                    if adj["delta"] < 0:
                        del_end_0 = mut_start_0 + abs(adj["delta"])
                        if original_pos >= del_end_0:
                            total += adj["delta"]
                        elif original_pos > mut_start_0:
                            total -= (original_pos - mut_start_0)
                    else:
                        total += adj["delta"]
            return max(0, original_pos + total)

        new_features = []
        for feature in rec.features:
            locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
            gene_name = feature.qualifiers.get("gene", [""])[0]
            edit = saved_edits.get(locus_tag) or saved_edits.get(gene_name)
            if edit and edit.get("action") == "delete":
                continue
            plan = modification_plan.get(locus_tag) or modification_plan.get(gene_name)
            new_gene_val = None
            if plan:
                if plan.get("custom_name"):
                    new_gene_val = plan["custom_name"]
                elif plan.get("suffixes"):
                    base = plan.get("base") or gene_name or locus_tag or "unnamed"
                    new_gene_val = f"{base}{''.join(sorted(plan['suffixes']))}"
            try:
                loc_parts = feature.location.parts if isinstance(feature.location, CompoundLocation) else [feature.location]
                new_parts = [FeatureLocation(ExactPosition(shift_pos(pt.start)),
                                             ExactPosition(shift_pos(pt.end)), strand=pt.strand)
                             for pt in loc_parts]
                if not new_parts:
                    continue
                new_loc = new_parts[0] if len(new_parts) == 1 else CompoundLocation(new_parts)
                q = feature.qualifiers.copy()
                if new_gene_val:
                    q["gene"] = [new_gene_val]
                    q["note"] = q.get("note", []) + ["Renamed by automated consensus analysis"]
                if edit and edit.get("action") == "edit":
                    if edit.get("new_name"):
                        q["gene"] = [edit["new_name"]]
                    if edit.get("new_locus_tag"):
                        q["locus_tag"] = [edit["new_locus_tag"]]
                    if edit.get("new_product"):
                        q["product"] = [edit["new_product"]]
                new_features.append(SeqFeature(location=new_loc, type=feature.type, qualifiers=q, id=feature.id))
            except Exception:
                continue

        for p in new_annotations:
            if p.seq_id == current_id or p.seq_id in current_id:
                s = shift_pos(p.position - 1)
                new_features.append(SeqFeature(
                    location=FeatureLocation(ExactPosition(s), ExactPosition(s + 1), strand=1),
                    type="misc_feature",
                    qualifiers={"note": [p.mutation_detail, "Added by automated consensus analysis"],
                                "gene": [p.proposed_name]}))

        new_rec = SeqRecord(new_seq_obj, id=rec.id, name=rec.name,
                            description=f"{rec.description} [Modified]",
                            features=new_features, annotations=rec.annotations.copy())
        new_rec.annotations.setdefault("molecule_type", "DNA")
        new_records.append(new_rec)

    out = StringIO()
    SeqIO.write(new_records, out, "genbank")
    return Response(content=out.getvalue(), media_type="text/plain",
                    headers={"Content-Disposition": f'attachment; filename="consensus_{job_id}.gbk"'})
