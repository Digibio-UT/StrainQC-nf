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
    # Persistent curated reference genomes (friendly-named .gbk), shown in the reference
    # dropdown alongside uploads and NOT wiped when clearing uploads/. Drop standard genomes
    # here (e.g. Ecoli_K12_MG1655.gbk, Scerevisiae_S288C.gbk).
    REFERENCE_DIR = Path(os.getenv("STRAINQC_REFERENCE_DIR", DATA_DIR / "references"))

    # Nextflow invocation.
    NEXTFLOW_BIN = os.getenv("NEXTFLOW_BIN", "nextflow")
    # Profiles for real HPC execution (Slurm + module-provisioned tools).
    NF_PROFILES  = os.getenv("STRAINQC_NF_PROFILES", "slurm,hpc_modules")
    # Optional extra args passed verbatim to `nextflow run` (e.g. "--hpc_account xyz").
    NF_EXTRA_ARGS = os.getenv("STRAINQC_NF_EXTRA_ARGS", "")

    # Reference DBs for contamination detection. Site-specific → set via env (not committed).
    # UT HPC (old SlurmApp locations):
    #   STRAINQC_KRAKEN2_DB=/gpfs/helios/home/etais/hpc_marko.looke/DTU_project/SlurmApp/backend/data/kraken_db
    #   STRAINQC_MASH_DB=/gpfs/helios/home/etais/hpc_marko.looke/DTU_project/SlurmApp/backend/data/mash_db/refseq.msh
    KRAKEN2_DB = os.getenv("STRAINQC_KRAKEN2_DB")   # dir; None → Kraken2 skipped
    MASH_DB    = os.getenv("STRAINQC_MASH_DB")      # refseq.msh; None → MASH skipped

    # Seqera Platform monitoring: add `-with-tower` to `nextflow run` so runs stream to
    # Seqera. Auto-on when TOWER_ACCESS_TOKEN is in the environment (Nextflow reads the token
    # + optional TOWER_WORKSPACE_ID from env). Force with STRAINQC_WITH_TOWER=true/false.
    _with_tower = os.getenv("STRAINQC_WITH_TOWER")
    WITH_TOWER = (_with_tower.lower() in ("1", "true", "yes")
                  if _with_tower is not None else bool(os.getenv("TOWER_ACCESS_TOKEN")))

    # How runs are started. "local" (default): `nextflow run` on this machine, e.g. the UT HPC
    # login node. "seqera": launch the site's Launchpad pipeline through the Seqera Platform API
    # (e.g. DTU BRIGHT on Azure Batch); see app/seqera_runner.py.
    LAUNCH_MODE = os.getenv("STRAINQC_LAUNCH_MODE", "local").strip().lower()

    # Seqera mode. The token is TOWER_ACCESS_TOKEN (read from env; keep it in a chmod-600 file
    # outside the repo). The Launchpad pipeline supplies the compute environment, work dir,
    # revision and default parameters; the launcher adds input/reference/outdir per job.
    SEQERA_API          = os.getenv("SEQERA_API", "https://api.cloud.seqera.io").rstrip("/")
    SEQERA_TOKEN        = os.getenv("TOWER_ACCESS_TOKEN")
    SEQERA_WORKSPACE_ID = os.getenv("SEQERA_WORKSPACE_ID")
    SEQERA_PIPELINE_ID  = os.getenv("SEQERA_PIPELINE_ID")
    # Where DATA_DIR lives in object storage, as "<local prefix>=<remote prefix>", e.g.
    # "/home/azureuser/StrainQC=az://aletest/StrainQC" for a blobfuse2 mount. Every file the
    # pipeline reads or writes (uploads, references, samplesheets, outdir) must be under it.
    STORAGE_PATH_MAP = os.getenv("STRAINQC_STORAGE_PATH_MAP", "")
    # job_id -> Seqera workflow id. Keep it on local disk, not on a blob mount.
    SEQERA_JOBS_FILE = Path(os.getenv("STRAINQC_SEQERA_JOBS_FILE", DATA_DIR / "seqera_jobs.json"))

    ALLOWED_UPLOAD_EXT = (
        ".fastq", ".fastq.gz", ".fq", ".fq.gz",
        ".fasta", ".fna", ".fa", ".gbk", ".gb",
    )

    def ensure_dirs(self):
        for d in (self.DATA_DIR, self.UPLOAD_DIR, self.RESULTS_DIR, self.WORK_DIR, self.REFERENCE_DIR):
            d.mkdir(parents=True, exist_ok=True)


settings = Settings()
