#!/usr/bin/env bash
# Run the 4-sample Ottilie pilot through StrainQC-nf (GATK path + CNVkit aneuploidy).
#
# Prereqs on the HPC login node:
#   module load nextflow/25.10.2 apptainer/1.4.1
#   - reference staged (S288C R64-1-1 GenBank with Ensembl/Roman chr names)
#   - FASTQs downloaded + samplesheet written (01_data/download_pilot_fastq.sh)
#   - container image built WITH cnvkit (rebuild after this branch: docker image or SIF)
#
# Runs on Slurm (login-node compute is prohibited). Override anything via env:
#   REF=/gpfs/.../S288C_R64_ensembl_chrnames.gb \
#   CONTAINER=docker://ghcr.io/digibio-ut/strainqc:latest \
#   PROFILES=slurm,apptainer  OUTDIR=results/ottilie_pilot \
#   bash 03_run/run_pilot.sh   [-resume ...]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO="$(cd "$BENCH_DIR/../.." && pwd)"          # StrainQC-nf repo root

SHEET="${SHEET:-$BENCH_DIR/data/samplesheet_pilot.csv}"
REF="${REF:-$REPO/ottilie_test_data/S288C_reference/S288C_R64_ensembl_chrnames.gb}"
CONTAINER="${CONTAINER:-docker://ghcr.io/digibio-ut/strainqc:latest}"
PROFILES="${PROFILES:-slurm,apptainer}"
OUTDIR="${OUTDIR:-$BENCH_DIR/results/ottilie_pilot}"
PARENT="${PARENT:-NODRUG-GM2}"

for f in "$SHEET" "$REF"; do
    [[ -s "$f" ]] || { echo "ERROR: missing $f (set SHEET / REF)"; exit 1; }
done

echo "repo=$REPO"
echo "samplesheet=$SHEET"
echo "reference=$REF"
echo "container=$CONTAINER"
echo "profiles=$PROFILES  parent=$PARENT  outdir=$OUTDIR"
echo

nextflow run "$REPO" \
    -profile "$PROFILES" \
    --input "$SHEET" \
    --reference "$REF" \
    --variant_caller gatk \
    --ploidy 1 \
    --parent_sample "$PARENT" \
    --strainqc_container "$CONTAINER" \
    --outdir "$OUTDIR" \
    "$@"

echo
echo "Done. results.json per sample under: $OUTDIR/<sample>/results.json"
echo "CNV tables under: $OUTDIR/<evolved sample>/cnv/"
echo "Next: 04_comparison (compare_snvs.py, check_cnv.py, make_report.py)."
