# StrainQC-nf launcher

Thin FastAPI service that re-exposes the StrainQC API contract, backed by `nextflow run`
instead of the old Celery + Slurm-over-SSH worker. The existing React frontend talks to it
unchanged. Runs **on the HPC** (or any host that can launch the pipeline and reach Slurm).

## Status

Implemented — **execution core** + first read-only viewers:
`POST /upload/`, `GET /upload/files`, `POST /analysis/start` (→ `nextflow run`),
`GET /analysis/status/{job_id}`, `GET /analysis/results/{job_id}`,
`GET /analysis/active_jobs`, `POST /analysis/stop/{job_id}`, `/results_files` static mount,
`GET /analysis/coverage/{job_id}` (downsampled + cached), and
`GET /analysis/genome_features/{job_id}` (BioPython GBK parse + saved edits).

**Not yet ported:** `evidence`, `consensus*`, `annotations*` (GET/POST). Backend-agnostic
readers over `RESULTS_DIR/<job_id>/` — next increment. Note `genome_features` needs a **GBK**
reference (the synthetic test uses FASTA, so it 404s there — expected; use a GBK reference to see genes).

## Run

```bash
cd launcher
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# On the HPC, make sure `nextflow` is on PATH (module load nextflow/25.10.2) and point
# the data dirs at a filesystem the compute nodes can see:
export PIPELINE_DIR=/gpfs/helios/projects/dtu_project/StrainQC-nf
export STRAINQC_DATA_DIR=/gpfs/helios/projects/dtu_project/strainqc_launcher_data
export STRAINQC_NF_PROFILES=slurm,hpc_modules          # default
# export STRAINQC_NF_EXTRA_ARGS="--hpc_account ealloc_bioeng_h-1 --hpc_queue main"

# Contamination-detection reference DBs (else kraken/mash report "skipped"):
export STRAINQC_KRAKEN2_DB=/gpfs/helios/home/etais/hpc_marko.looke/DTU_project/SlurmApp/backend/data/kraken_db
export STRAINQC_MASH_DB=/gpfs/helios/home/etais/hpc_marko.looke/DTU_project/SlurmApp/backend/data/mash_db/refseq.msh

uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## Local frontend against the HPC launcher (dev loop)

SSH-tunnel the launcher port to your Mac, then run Vite (its proxy already targets
`127.0.0.1:8000`):

```bash
ssh -N -L 8000:localhost:8000 hpc_marko.looke@rocket.hpc.ut.ee   # on the Mac
# then, in the frontend: npm run dev  → talks to the tunneled launcher
```

## How a job runs

`/analysis/start` writes a one-line samplesheet (`sample=<job_id>`), then launches
`nextflow run $PIPELINE_DIR -profile slurm,hpc_modules --input … --reference … --outdir RESULTS_DIR`
as a background process. Nextflow submits each step as a Slurm job. Results publish to
`RESULTS_DIR/<job_id>/results.json` (+ coverage.tsv, gatk/, breseq/…), which `status`/`results`
read directly — same layout and contract as the old app.

## Config (env vars)

| Var | Default | Meaning |
|---|---|---|
| `PIPELINE_DIR` | repo root | where `main.nf` lives |
| `STRAINQC_DATA_DIR` | `launcher_data/` | base for uploads/results/work |
| `STRAINQC_NF_PROFILES` | `slurm,hpc_modules` | Nextflow `-profile` |
| `STRAINQC_NF_EXTRA_ARGS` | – | extra `nextflow run` args |
| `STRAINQC_KRAKEN2_DB` | – | Kraken2 DB dir (unset → Kraken2 skipped) |
| `STRAINQC_MASH_DB` | – | MASH `refseq.msh` (unset → MASH skipped) |
| `NEXTFLOW_BIN` | `nextflow` | nextflow executable |
