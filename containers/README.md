# StrainQC-nf container

One all-in-one image (`strainqc.def`) provides every pipeline tool, so the whole pipeline
runs under `-profile apptainer` (or `docker`) with no `module load` — the portable path to
Seqera/cloud. Under `-profile hpc_modules` the container is ignored (tools come from modules).

## Build the SIF on the HPC

```bash
cd /gpfs/helios/projects/dtu_project/StrainQC-nf
module load apptainer/1.4.1
apptainer build --fakeroot containers/strainqc.sif containers/strainqc.def
```

- `--fakeroot` lets you build unprivileged. If your site hasn't enabled it, try without it,
  or `apptainer build --sandbox` / ask the admins. The build needs internet (login node) and
  a few GB of scratch for the conda solve.
- The `%test` section verifies every tool is on PATH; a successful build prints `OK: <tool>`.

Store the SIF somewhere the compute nodes can read (gpfs), e.g.
`/gpfs/helios/projects/dtu_project/containers/strainqc.sif`.

## Run with the container instead of modules

```bash
nextflow run . -profile test,slurm,apptainer \
    --strainqc_container /gpfs/helios/projects/dtu_project/StrainQC-nf/containers/strainqc.sif
```

- `apptainer` enables the runtime; `slurm` submits jobs; `--strainqc_container <path>` points
  every process at the SIF (default is the bare name `strainqc.sif`).
- For a quick local smoke test: `-profile test,apptainer --strainqc_container <path-or-docker-uri>`.

## Later: registry image for cloud / Seqera

Convert the same recipe to a Dockerfile and push to `ghcr.io/digibio-ut/strainqc:<tag>`, then
`--strainqc_container docker://ghcr.io/digibio-ut/strainqc:<tag>`. Cloud/Seqera executors pull
by URI; no local SIF needed.
