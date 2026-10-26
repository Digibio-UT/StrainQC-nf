#!/usr/bin/env bash
# Ottilie pilot run using a conda tool env instead of the container.
#
# Why: the UT HPC apptainer/openjdk spack modules were broken (missing gcc-runtime dep), so
# neither the container path nor the module path worked. This runs the GATK path + CNVkit from
# a self-built conda env (`strainqc-tools`: bwa-mem2, gatk4, samtools, cnvkit, python, openjdk).
# Each Slurm job simply prepends that env's bin to PATH via the beforeScript (hpc_modules
# profile) — no `conda activate` (which misbehaves on the compute nodes).
#
# Adjust ENV / NFBIN / paths for your account. Then: nohup bash run_pilot_conda.sh > nf_run.log 2>&1 &
set -euo pipefail

REPO=/gpfs/helios/projects/dtu_project/StrainQC-nf
ENV=/gpfs/helios/home/etais/hpc_marko.looke/.conda/envs/strainqc-tools
NFBIN=/gpfs/helios/software/cluster_software/spack/linux-rhel9-x86_64/none-none/nextflow-25.10.2-qpujk6dk3kn76x3p3p3iysks53isng5i/bin

# Head process (login node) needs java + nextflow on PATH.
export JAVA_HOME="$ENV"
export PATH="$ENV/bin:$NFBIN:$PATH"
export NXF_ANSI_LOG=false

# Per-job beforeScript: put the tool env on PATH ($PATH stays literal → expands inside the job).
BEFORE="export PATH=$ENV/bin:\$PATH; export JAVA_HOME=$ENV; unset PYTHONPATH 2>/dev/null || true"

cd "$REPO"
exec nextflow run . -profile slurm,hpc_modules \
  --hpc_pre_command "$BEFORE" \
  --input     tests/ottilie_benchmark/data/samplesheet_pilot.csv \
  --reference tests/ottilie_benchmark/data/S288C_R64_ensembl_chrnames.gb \
  --variant_caller gatk --ploidy 1 --parent_sample NODRUG-GM2 \
  --outdir    tests/ottilie_benchmark/results/ottilie_pilot \
  -resume
