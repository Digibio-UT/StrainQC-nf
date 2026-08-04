#!/usr/bin/env bash
# Start the StrainQC-nf launcher with Nextflow + env set up. Run from launcher/ with the
# venv active:   source .venv/bin/activate && bash run.sh
# All STRAINQC_* / PIPELINE_DIR values below are UT-HPC defaults; override by exporting
# them before running.
set -euo pipefail

# Make `module` available even in a non-login shell, then load Nextflow.
if ! type module &>/dev/null; then
  for f in /etc/profile.d/lmod.sh /etc/profile.d/z00_lmod.sh /etc/profile.d/modules.sh; do
    [ -f "$f" ] && source "$f" && break
  done
fi
module load "${STRAINQC_NEXTFLOW_MODULE:-nextflow/25.10.2}"

export PIPELINE_DIR="${PIPELINE_DIR:-/gpfs/helios/projects/dtu_project/StrainQC-nf}"
export STRAINQC_DATA_DIR="${STRAINQC_DATA_DIR:-/gpfs/helios/projects/dtu_project/strainqc_launcher_data}"
export STRAINQC_KRAKEN2_DB="${STRAINQC_KRAKEN2_DB:-/gpfs/helios/home/etais/hpc_marko.looke/DTU_project/SlurmApp/backend/data/kraken_db}"
export STRAINQC_MASH_DB="${STRAINQC_MASH_DB:-/gpfs/helios/home/etais/hpc_marko.looke/DTU_project/SlurmApp/backend/data/mash_db/refseq.msh}"

echo "nextflow: $(command -v nextflow)"
echo "PIPELINE_DIR=$PIPELINE_DIR"
echo "KRAKEN2_DB=$STRAINQC_KRAKEN2_DB   MASH_DB=$STRAINQC_MASH_DB"
exec uvicorn app.main:app --host "${HOST:-0.0.0.0}" --port "${PORT:-8000}"
