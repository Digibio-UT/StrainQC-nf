"""Launcher Seqera mode (launcher/app/seqera_runner.py) and the runner's mode switch.

The UT HPC deployment must keep launching with a local `nextflow run` unless a site opts in
with STRAINQC_LAUNCH_MODE=seqera, so the default is pinned here. The Seqera API is replaced by
a fake that records requests; nothing here talks to Seqera.
"""

import json
import os
import subprocess
import sys

import pytest

from conftest import REPO

LAUNCHER = REPO / "launcher"
sys.path.insert(0, str(LAUNCHER))

from app import runner, seqera_runner  # noqa: E402
from app.config import settings        # noqa: E402

LAUNCHPAD = {"launch": {
    "computeEnv": {"id": "ce123", "name": "yAMP-ce-nofusion-256"},
    "pipeline": "https://github.com/Digibio-UT/StrainQC-nf",
    "revision": "v0.2.1", "workDir": "az://aletest/nf-work", "configProfiles": [],
    "pullLatest": False, "paramsText": "max_cpus: 4\nmax_memory: 30.GB\n",
}}


def test_default_launch_mode_is_local():
    env = {k: v for k, v in os.environ.items() if k != "STRAINQC_LAUNCH_MODE"}
    out = subprocess.run([sys.executable, "-c", "from app.config import settings; print(settings.LAUNCH_MODE)"],
                         cwd=LAUNCHER, env=env, capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "local"


@pytest.fixture
def site(tmp_path, monkeypatch):
    """A data dir inside a pretend blob mount, wired into the shared settings object."""
    mount = tmp_path / "mnt"
    data = mount / "app"
    for name, sub in [("DATA_DIR", ""), ("UPLOAD_DIR", "uploads"), ("RESULTS_DIR", "results"),
                      ("WORK_DIR", "work"), ("REFERENCE_DIR", "references")]:
        monkeypatch.setattr(settings, name, data / sub if sub else data)
    monkeypatch.setattr(settings, "STORAGE_PATH_MAP", f"{mount}=az://aletest/StrainQC")
    monkeypatch.setattr(settings, "SEQERA_JOBS_FILE", tmp_path / "state" / "seqera_jobs.json")
    monkeypatch.setattr(settings, "SEQERA_TOKEN", "test-token")
    monkeypatch.setattr(settings, "SEQERA_WORKSPACE_ID", "79597273081110")
    monkeypatch.setattr(settings, "SEQERA_PIPELINE_ID", "190652911420583")
    monkeypatch.setattr(settings, "KRAKEN2_DB", None)
    monkeypatch.setattr(settings, "MASH_DB", None)
    settings.ensure_dirs()
    reads = [settings.UPLOAD_DIR / "r1.fastq", settings.UPLOAD_DIR / "r2.fastq"]
    for r in reads:
        r.write_text("@r\nA\n+\nI\n")
    ref = settings.REFERENCE_DIR / "U00096.3.gbk"
    ref.write_text("LOCUS\n")
    return {"mount": mount, "reads": reads, "ref": ref}


@pytest.fixture
def fake_api(monkeypatch):
    calls, state = [], {"status": "RUNNING", "errorMessage": None}

    def api(method, path, body=None):
        calls.append((method, path, body))
        if path.startswith("/pipelines/"):
            return LAUNCHPAD
        if path == "/workflow/launch":
            return {"workflowId": "wf42"}
        if path.startswith("/workflow/") and method == "GET":
            return {"workflow": {"status": state["status"], "errorMessage": state["errorMessage"]}}
        return {}
    monkeypatch.setattr(seqera_runner, "_api", api)
    return calls, state


# ── local mode is unchanged ────────────────────────────────────────────────────────────

def test_local_mode_runs_nextflow_and_never_calls_seqera(site, fake_api, monkeypatch):
    monkeypatch.setattr(settings, "LAUNCH_MODE", "local")
    monkeypatch.setattr(settings, "WITH_TOWER", False)
    monkeypatch.setattr(settings, "NF_EXTRA_ARGS", "")
    seen = {}

    class FakeProc:
        pid = 123

        def poll(self):
            return None

    def popen(cmd, **kw):
        seen["cmd"] = cmd
        return FakeProc()
    monkeypatch.setattr(runner.subprocess, "Popen", popen)
    runner.start_job("job_x", site["reads"], site["ref"], {"variant_caller": "breseq"})
    cmd = seen["cmd"]
    assert cmd[1] == "run" and "--outdir" in cmd and str(settings.RESULTS_DIR) in cmd
    assert cmd[cmd.index("--variant_caller") + 1] == "breseq"
    assert fake_api[0] == []
    assert runner.job_process_state("job_x") == "running"


def test_unknown_launch_mode_is_an_error(monkeypatch):
    monkeypatch.setattr(settings, "LAUNCH_MODE", "slurm-direct")
    with pytest.raises(RuntimeError, match="STRAINQC_LAUNCH_MODE"):
        runner.list_active()


# ── helpers ────────────────────────────────────────────────────────────────────────────

def test_to_remote_maps_paths_under_the_mount(site):
    assert seqera_runner.to_remote(site["reads"][0]) == "az://aletest/StrainQC/app/uploads/r1.fastq"
    assert seqera_runner.to_remote(settings.RESULTS_DIR) == "az://aletest/StrainQC/app/results"


def test_to_remote_rejects_paths_outside_the_mount(site, tmp_path):
    with pytest.raises(seqera_runner.SeqeraError, match="outside the storage mount"):
        seqera_runner.to_remote(tmp_path / "elsewhere.fastq")


def test_parse_params_text():
    assert seqera_runner.parse_params_text("max_cpus: 4\nmax_memory: 30.GB\n# note\nrun_blast: false") == \
        {"max_cpus": 4, "max_memory": "30.GB", "run_blast": False}
    assert seqera_runner.parse_params_text('{"ploidy": 2}') == {"ploidy": 2}
    assert seqera_runner.parse_params_text(None) == {}
    with pytest.raises(seqera_runner.SeqeraError):
        seqera_runner.parse_params_text("outer:\n  inner: 1")


@pytest.mark.parametrize("job_id, expected", [
    ("job_07e5a9d5-d79a-4fd1_ab12cd34", "job_07e5a9d5-d79a-4fd1_ab12cd34"),
    ("Job__X..Y", "job-x-y"),
    ("123abc", "strainqc-123abc"),
])
def test_run_name_is_a_valid_seqera_name(job_id, expected):
    assert seqera_runner.run_name(job_id) == expected


# ── Seqera mode ────────────────────────────────────────────────────────────────────────

def test_seqera_start_launches_the_launchpad_pipeline(site, fake_api, monkeypatch):
    monkeypatch.setattr(settings, "LAUNCH_MODE", "seqera")
    monkeypatch.setattr(settings, "KRAKEN2_DB", "az://aletest/StrainQC/databases/kraken_db")
    info = runner.start_job("job_abc", site["reads"], site["ref"],
                            {"variant_caller": "breseq", "ploidy": 1, "not_a_param": "x"})
    assert info["workflow_id"] == "wf42"

    sheet = (settings.RESULTS_DIR / "job_abc" / "_nf" / "samplesheet.csv").read_text().splitlines()
    assert sheet[1] == ("job_abc,az://aletest/StrainQC/app/uploads/r1.fastq,"
                        "az://aletest/StrainQC/app/uploads/r2.fastq")

    calls, _ = fake_api
    method, path, body = calls[-1]
    assert (method, path) == ("POST", "/workflow/launch")
    launch = body["launch"]
    assert launch["computeEnvId"] == "ce123" and launch["revision"] == "v0.2.1"
    assert launch["workDir"] == "az://aletest/nf-work" and launch["runName"] == "job_abc"
    params = json.loads(launch["paramsText"])
    assert params == {
        "max_cpus": 4, "max_memory": "30.GB",                       # from the Launchpad entry
        "variant_caller": "breseq", "ploidy": 1,                     # forwarded options only
        "kraken2_db": "az://aletest/StrainQC/databases/kraken_db",
        "input": "az://aletest/StrainQC/app/results/job_abc/_nf/samplesheet.csv",
        "reference": "az://aletest/StrainQC/app/references/U00096.3.gbk",
        "outdir": "az://aletest/StrainQC/app/results",
    }
    assert json.loads(settings.SEQERA_JOBS_FILE.read_text())["job_abc"]["workflow_id"] == "wf42"


def test_seqera_status_log_and_stop(site, fake_api, monkeypatch):
    monkeypatch.setattr(settings, "LAUNCH_MODE", "seqera")
    runner.start_job("job_s", site["reads"], site["ref"], {})
    calls, state = fake_api

    assert runner.job_process_state("job_s") == "running"
    assert [a["job_id"] for a in runner.list_active()] == ["job_s"]
    assert runner.stop_job("job_s") is True
    assert ("POST", "/workflow/wf42/cancel", {}) in calls

    state["status"] = "SUCCEEDED"
    assert runner.job_process_state("job_s") == "exited_ok"
    assert runner.stop_job("job_s") is False

    state.update(status="FAILED", errorMessage="Process BRESEQ failed")
    assert runner.job_process_state("job_s") == "exited_error"
    assert "BRESEQ failed" in runner.tail_log("job_s")

    assert runner.job_process_state("never_started") is None


def test_seqera_start_needs_a_pipeline_id(site, fake_api, monkeypatch):
    monkeypatch.setattr(settings, "LAUNCH_MODE", "seqera")
    monkeypatch.setattr(settings, "SEQERA_PIPELINE_ID", None)
    with pytest.raises(RuntimeError, match="SEQERA_PIPELINE_ID"):
        runner.start_job("job_y", site["reads"], site["ref"], {})
