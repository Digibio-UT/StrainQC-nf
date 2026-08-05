"""
Consensus-genome endpoints (ported from the old app's analysis.py).

- GET /analysis/consensus/{job_id}          → FASTA with variants applied (this file)
- GET /analysis/consensus/preview/{job_id}  → proposed annotation changes  (TODO: next increment)
- POST /analysis/consensus/download/{job_id} → GBK with approved changes    (TODO: next increment)

Reads the reference + variants from RESULTS_DIR/<job_id>/ and applies mutations above a
frequency threshold. GATK VCF is located by glob (pipeline names it <sample>.vcf, not the
old app's fixed variants.vcf).
"""
import glob
import json
import logging
import os

from fastapi import APIRouter, HTTPException, Response

from app.config import settings

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
