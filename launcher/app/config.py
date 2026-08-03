"""Launcher configuration (env-overridable)."""
import os
from pathlib import Path

# Repo root = two levels up from this file (launcher/app/config.py -> repo root)
REPO_ROOT = Path(__file__).resolve().parents[2]


class Settings:
    API_V1_STR = "/api/v1"
    PROJECT_NAME = "StrainQC-nf launcher"

    # Where the Nextflow pipeline lives (main.nf). Defaults to this repo.
    PIPELINE_DIR = Path(os.getenv("PIPELINE_DIR", REPO_ROOT))

    # Data directories. On the HPC these should sit on a filesystem the compute nodes see.
    DATA_DIR    = Path(os.getenv("STRAINQC_DATA_DIR", REPO_ROOT / "launcher_data"))
    UPLOAD_DIR  = Path(os.getenv("STRAINQC_UPLOAD_DIR", DATA_DIR / "uploads"))
    RESULTS_DIR = Path(os.getenv("STRAINQC_RESULTS_DIR", DATA_DIR / "results"))
    WORK_DIR    = Path(os.getenv("STRAINQC_WORK_DIR", DATA_DIR / "work"))

    # Nextflow invocation.
    NEXTFLOW_BIN = os.getenv("NEXTFLOW_BIN", "nextflow")
    # Profiles for real HPC execution (Slurm + module-provisioned tools).
    NF_PROFILES  = os.getenv("STRAINQC_NF_PROFILES", "slurm,hpc_modules")
    # Optional extra args passed verbatim to `nextflow run` (e.g. "--hpc_account xyz").
    NF_EXTRA_ARGS = os.getenv("STRAINQC_NF_EXTRA_ARGS", "")

    ALLOWED_UPLOAD_EXT = (
        ".fastq", ".fastq.gz", ".fq", ".fq.gz",
        ".fasta", ".fna", ".fa", ".gbk", ".gb",
    )

    def ensure_dirs(self):
        for d in (self.DATA_DIR, self.UPLOAD_DIR, self.RESULTS_DIR, self.WORK_DIR):
            d.mkdir(parents=True, exist_ok=True)


settings = Settings()
