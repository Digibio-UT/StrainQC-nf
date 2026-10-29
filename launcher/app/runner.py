"""
Nextflow execution layer — the piece that replaces the old app's Celery worker +
slurm_runner.py (paramiko/sbatch/squeue-polling). Here a job is just a `nextflow run`
subprocess; Nextflow itself submits/monitors the Slurm jobs.
"""
import os
import shlex
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Optional

from app.config import settings
from app import seqera_runner

# In-memory job registry. Enough for a single-lab launcher; a shared instance would
# persist this (DB) — see the auth/multi-tenancy strategy question.
_JOBS: Dict[str, dict] = {}

# Option keys forwarded to the pipeline as `--<key> <value>` (names mirror nextflow.config).
_FORWARD_OPTS = {
    "variant_caller", "polymorphism_mode", "ploidy",
    "breseq_limit_fold_coverage", "breseq_base_quality_cutoff", "breseq_min_mapping_quality",
    "breseq_polymorphism_frequency_cutoff", "breseq_no_junction_prediction",
    "breseq_consensus_min_cov_each_strand", "breseq_targeted_sequencing",
    "gatk_min_base_quality", "gatk_min_mapping_quality", "gatk_max_reads_per_align_start",
    "gatk_dont_use_soft_clipped", "gatk_call_conf", "gatk_callable_depth",
}


def resolve_upload(file_id: str) -> Optional[Path]:
    """Find an uploaded file whose name starts with file_id (matches old-app behaviour)."""
    if not settings.UPLOAD_DIR.exists():
        return None
    for p in settings.UPLOAD_DIR.iterdir():
        if p.name.startswith(file_id):
            return p
    return None


def resolve_reference(ref_id: str) -> Optional[Path]:
    """Resolve a reference selection: the curated library first (friendly names), then uploads.
    Matches by exact name or prefix (the dropdown sends the full filename)."""
    for d in (settings.REFERENCE_DIR, settings.UPLOAD_DIR):
        if not d.exists():
            continue
        for p in sorted(d.iterdir()):
            if p.is_file() and p.name.startswith(ref_id):
                return p
    return None


def _seqera() -> bool:
    """True when runs go through the Seqera API (STRAINQC_LAUNCH_MODE=seqera), else local."""
    if settings.LAUNCH_MODE not in ("local", "seqera"):
        raise RuntimeError(f"Unknown STRAINQC_LAUNCH_MODE {settings.LAUNCH_MODE!r} (use 'local' or 'seqera').")
    return settings.LAUNCH_MODE == "seqera"


def _seqera_params(options: dict) -> dict:
    """The same forwarded options as the local `--key value` flags, as Launchpad parameters.
    Contamination DBs are passed as given: in Seqera mode they are object-storage paths."""
    params = {k: v for k, v in (options or {}).items() if k in _FORWARD_OPTS and v is not None}
    if settings.KRAKEN2_DB:
        params["kraken2_db"] = settings.KRAKEN2_DB
    if settings.MASH_DB:
        params["mash_db"] = settings.MASH_DB
    return params


def _write_samplesheet(job_id: str, read_paths: List[Path], launch_dir: Path) -> Path:
    sheet = launch_dir / "samplesheet.csv"
    r1 = str(read_paths[0])
    r2 = str(read_paths[1]) if len(read_paths) > 1 else ""
    sheet.write_text(f"sample,fastq_1,fastq_2\n{job_id},{r1},{r2}\n")
    return sheet


def start_job(job_id: str, read_paths: List[Path], reference: Path, options: dict) -> dict:
    """Launch `nextflow run` for a job as a detached-ish background subprocess."""
    if _seqera():
        return seqera_runner.start_job(job_id, read_paths, reference, _seqera_params(options))
    settings.ensure_dirs()
    launch_dir = settings.RESULTS_DIR / job_id / "_nf"
    launch_dir.mkdir(parents=True, exist_ok=True)
    sheet = _write_samplesheet(job_id, read_paths, launch_dir)

    cmd = [
        settings.NEXTFLOW_BIN, "run", str(settings.PIPELINE_DIR),
        "-profile", settings.NF_PROFILES,
        "-work-dir", str(settings.WORK_DIR / job_id),
        "-name", job_id,
        "--input", str(sheet),
        "--reference", str(reference),
        "--outdir", str(settings.RESULTS_DIR),
    ]
    # sample name == job_id, so publishDir ${outdir}/${meta.id} == RESULTS_DIR/<job_id>

    options = options or {}
    for key in _FORWARD_OPTS:
        if key in options and options[key] is not None:
            val = options[key]
            if isinstance(val, bool):
                val = "true" if val else "false"
            cmd += [f"--{key}", str(val)]

    # Reference DBs for contamination detection — only if configured and present on disk.
    if settings.KRAKEN2_DB and os.path.exists(settings.KRAKEN2_DB):
        cmd += ["--kraken2_db", settings.KRAKEN2_DB]
    if settings.MASH_DB and os.path.exists(settings.MASH_DB):
        cmd += ["--mash_db", settings.MASH_DB]

    # Stream run telemetry to Seqera Platform when a token is configured.
    if settings.WITH_TOWER:
        cmd += ["-with-tower"]

    if settings.NF_EXTRA_ARGS:
        cmd += shlex.split(settings.NF_EXTRA_ARGS)

    log_path = launch_dir / "nextflow.stdout.log"
    log_fh = open(log_path, "w")
    try:
        proc = subprocess.Popen(
            cmd, cwd=str(launch_dir),
            stdout=log_fh, stderr=subprocess.STDOUT,
            start_new_session=True,   # survive if the request thread goes away
        )
    except FileNotFoundError:
        log_fh.close()
        raise RuntimeError(
            f"'{settings.NEXTFLOW_BIN}' not found on PATH. Load Nextflow before starting the "
            f"launcher (e.g. `module load nextflow/25.10.2`), or set NEXTFLOW_BIN to its full path."
        )
    _JOBS[job_id] = {
        "proc": proc, "log": str(log_path), "cmd": " ".join(shlex.quote(c) for c in cmd),
        "started_at": time.time(), "launch_dir": str(launch_dir),
    }
    return {"job_id": job_id, "pid": proc.pid, "cmd": _JOBS[job_id]["cmd"]}


def job_process_state(job_id: str) -> Optional[str]:
    """'running' | 'exited_ok' | 'exited_error' | None (unknown/not tracked this session)."""
    if _seqera():
        return seqera_runner.job_process_state(job_id)
    entry = _JOBS.get(job_id)
    if not entry:
        return None
    rc = entry["proc"].poll()
    if rc is None:
        return "running"
    return "exited_ok" if rc == 0 else "exited_error"


def tail_log(job_id: str, n: int = 40) -> str:
    if _seqera():
        return seqera_runner.tail_log(job_id, n)
    entry = _JOBS.get(job_id)
    if not entry:
        return ""
    try:
        with open(entry["log"]) as fh:
            return "".join(fh.readlines()[-n:])
    except OSError:
        return ""


def list_active() -> List[dict]:
    if _seqera():
        return seqera_runner.list_active()
    out = []
    for job_id, e in _JOBS.items():
        if e["proc"].poll() is None:
            out.append({"job_id": job_id, "task_id": str(e["proc"].pid),
                        "state": "active", "started_at": e["started_at"]})
    return out


def stop_job(job_id: str) -> bool:
    """Terminate the nextflow process; it propagates `scancel` to submitted Slurm jobs."""
    if _seqera():
        return seqera_runner.stop_job(job_id)
    entry = _JOBS.get(job_id)
    if not entry or entry["proc"].poll() is not None:
        return False
    try:
        os.killpg(os.getpgid(entry["proc"].pid), 15)  # SIGTERM the process group
    except (ProcessLookupError, PermissionError):
        entry["proc"].terminate()
    return True
