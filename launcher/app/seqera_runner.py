"""
Seqera Platform execution layer (STRAINQC_LAUNCH_MODE=seqera).

Same five operations as the local runner (start / state / log / active / stop), but a job is a
run of the site's Launchpad pipeline started through the Seqera API instead of a local
`nextflow run`. The Launchpad entry supplies the compute environment, work dir, revision and
default parameters (e.g. max_cpus for the Batch VM size); per job the launcher adds the
samplesheet, reference and outdir as object-storage paths.

The launcher's data directory must be a mount of that object storage (e.g. blobfuse2), so
the files the launcher writes (uploads, samplesheets) are the files the pipeline reads, and
the results the pipeline publishes are the files the viewer endpoints read.
STRAINQC_STORAGE_PATH_MAP translates between the two path forms.
"""
import json
import re
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional

from app.config import settings

# Seqera run status -> the local runner's process states.
_RUNNING = {"SUBMITTED", "RUNNING"}
_OK = {"SUCCEEDED"}


class SeqeraError(RuntimeError):
    pass


# ── paths ──────────────────────────────────────────────────────────────────────────────

def _path_map():
    if "=" not in settings.STORAGE_PATH_MAP:
        raise SeqeraError("STRAINQC_STORAGE_PATH_MAP is not set (expected '<local prefix>=<remote prefix>').")
    local, remote = settings.STORAGE_PATH_MAP.split("=", 1)
    return Path(local).expanduser().resolve(), remote.rstrip("/")


def to_remote(path: Path) -> str:
    """/home/azureuser/StrainQC/app/uploads/x.fastq -> az://aletest/StrainQC/app/uploads/x.fastq"""
    local, remote = _path_map()
    p = Path(path).expanduser().resolve()
    try:
        rel = p.relative_to(local)
    except ValueError:
        raise SeqeraError(f"{p} is outside the storage mount {local}; the pipeline cannot read it.")
    return f"{remote}/{rel.as_posix()}" if rel.parts else remote


# ── Launchpad parameters ───────────────────────────────────────────────────────────────

def _scalar(v: str):
    v = v.strip().strip('"').strip("'")
    if v.lower() in ("true", "false"):
        return v.lower() == "true"
    if re.fullmatch(r"-?\d+", v):
        return int(v)
    if re.fullmatch(r"-?\d+\.\d+", v):
        return float(v)
    return v


def parse_params_text(text: Optional[str]) -> dict:
    """Launchpad params are YAML or JSON. Read JSON, or flat `key: value` YAML (what a
    Launchpad entry holds in practice, e.g. 'max_cpus: 4'); nested YAML is not supported."""
    text = (text or "").strip()
    if not text:
        return {}
    if text.startswith("{"):
        return json.loads(text)
    params = {}
    for line in text.splitlines():
        line = line.split(" #", 1)[0].rstrip()
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if line[0].isspace() or ":" not in line:
            raise SeqeraError(f"Launchpad parameters must be flat 'key: value' lines; cannot read: {line!r}")
        key, val = line.split(":", 1)
        params[key.strip()] = _scalar(val)
    return params


def run_name(job_id: str) -> str:
    """Seqera run names: lowercase letters, digits, '-' and '_', starting with a letter, <=80."""
    name = re.sub(r"[^a-z0-9_-]", "-", job_id.lower())
    name = re.sub(r"[-_]{2,}", "-", name).strip("-_")
    if not name or not name[0].isalpha():
        name = "strainqc-" + name
    return name[:80].rstrip("-_")


# ── API ────────────────────────────────────────────────────────────────────────────────

def _api(method: str, path: str, body: Optional[dict] = None) -> dict:
    if not (settings.SEQERA_TOKEN and settings.SEQERA_WORKSPACE_ID):
        raise SeqeraError("Seqera mode needs TOWER_ACCESS_TOKEN and SEQERA_WORKSPACE_ID.")
    sep = "&" if "?" in path else "?"
    url = f"{settings.SEQERA_API}{path}{sep}workspaceId={settings.SEQERA_WORKSPACE_ID}"
    data = json.dumps(body).encode() if body is not None else None
    req = urllib.request.Request(url, data=data, method=method, headers={
        "Authorization": f"Bearer {settings.SEQERA_TOKEN}",
        "Accept": "application/json", "Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            raw = resp.read()
    except urllib.error.HTTPError as e:
        detail = e.read().decode(errors="replace")[:500]
        raise SeqeraError(f"Seqera API {method} {path}: HTTP {e.code} {detail}")
    except urllib.error.URLError as e:
        raise SeqeraError(f"Seqera API unreachable ({settings.SEQERA_API}): {e.reason}")
    return json.loads(raw) if raw else {}


# ── job registry (job_id -> workflow id), persisted so a restart keeps it ──────────────

def _load_jobs() -> Dict[str, dict]:
    try:
        return json.loads(settings.SEQERA_JOBS_FILE.read_text())
    except (OSError, ValueError):
        return {}


def _save_job(job_id: str, entry: dict):
    jobs = _load_jobs()
    jobs[job_id] = entry
    settings.SEQERA_JOBS_FILE.parent.mkdir(parents=True, exist_ok=True)
    tmp = settings.SEQERA_JOBS_FILE.with_suffix(".tmp")
    tmp.write_text(json.dumps(jobs, indent=2))
    tmp.replace(settings.SEQERA_JOBS_FILE)


def _workflow(job_id: str) -> Optional[dict]:
    entry = _load_jobs().get(job_id)
    if not entry:
        return None
    try:
        return _api("GET", f"/workflow/{entry['workflow_id']}").get("workflow") or {}
    except SeqeraError:
        return {}          # API hiccup: report "unknown" rather than failing the status call


# ── the five runner operations ─────────────────────────────────────────────────────────

def start_job(job_id: str, read_paths: List[Path], reference: Path, params: dict) -> dict:
    """Write the samplesheet into the storage mount and launch the Launchpad pipeline."""
    if not settings.SEQERA_PIPELINE_ID:
        raise SeqeraError("Seqera mode needs SEQERA_PIPELINE_ID (the Launchpad pipeline to launch).")
    settings.ensure_dirs()
    launch_dir = settings.RESULTS_DIR / job_id / "_nf"
    launch_dir.mkdir(parents=True, exist_ok=True)
    reads = [to_remote(p) for p in read_paths]
    sheet = launch_dir / "samplesheet.csv"
    sheet.write_text(f"sample,fastq_1,fastq_2\n{job_id},{reads[0]},{reads[1] if len(reads) > 1 else ''}\n")

    launch = _api("GET", f"/pipelines/{settings.SEQERA_PIPELINE_ID}/launch").get("launch") or {}
    ce = launch.get("computeEnv") or {}
    merged = parse_params_text(launch.get("paramsText"))
    merged.update(params)
    merged.update({"input": to_remote(sheet), "reference": to_remote(reference),
                   "outdir": to_remote(settings.RESULTS_DIR)})
    # sample name == job_id, so publishDir ${outdir}/${meta.id} == RESULTS_DIR/<job_id>

    body = {"launch": {
        "computeEnvId": ce.get("id") or launch.get("computeEnvId"),
        "pipeline": launch.get("pipeline"),
        "revision": launch.get("revision"),
        "workDir": launch.get("workDir") or ce.get("config", {}).get("workDir"),
        "configProfiles": launch.get("configProfiles") or [],
        "configText": launch.get("configText"),
        "preRunScript": launch.get("preRunScript"),
        "postRunScript": launch.get("postRunScript"),
        "pullLatest": launch.get("pullLatest", False),
        "paramsText": json.dumps(merged),
        "runName": run_name(job_id),
    }}
    body["launch"] = {k: v for k, v in body["launch"].items() if v not in (None, "")}
    wf_id = _api("POST", "/workflow/launch", body).get("workflowId")
    if not wf_id:
        raise SeqeraError("Seqera accepted the launch but returned no workflowId.")
    _save_job(job_id, {"workflow_id": wf_id, "started_at": time.time(),
                       "pipeline_id": settings.SEQERA_PIPELINE_ID})
    return {"job_id": job_id, "pid": wf_id, "workflow_id": wf_id,
            "cmd": f"Seqera run {run_name(job_id)} ({wf_id})"}


def job_process_state(job_id: str) -> Optional[str]:
    """'running' | 'exited_ok' | 'exited_error' | None, mapped from the Seqera run status."""
    wf = _workflow(job_id)
    if not wf:
        return None
    status = (wf.get("status") or "").upper()
    if status in _RUNNING:
        return "running"
    if status in _OK:
        return "exited_ok"
    return "exited_error" if status else None


def tail_log(job_id: str, n: int = 40) -> str:
    wf = _workflow(job_id) or {}
    text = wf.get("errorReport") or wf.get("errorMessage") or ""
    if not text and wf.get("status"):
        text = f"Seqera run status: {wf['status']}"
    return "\n".join(text.splitlines()[-n:])


def list_active() -> List[dict]:
    out = []
    for job_id, e in _load_jobs().items():
        if job_process_state(job_id) == "running":
            out.append({"job_id": job_id, "task_id": e["workflow_id"],
                        "state": "active", "started_at": e.get("started_at")})
    return out


def stop_job(job_id: str) -> bool:
    entry = _load_jobs().get(job_id)
    if not entry or job_process_state(job_id) != "running":
        return False
    _api("POST", f"/workflow/{entry['workflow_id']}/cancel", {})
    return True
