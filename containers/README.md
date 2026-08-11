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

## Registry image (ghcr.io) — for cloud / Seqera / other sites

The same tool stack is also published as a registry image so it can run anywhere (cloud
batch, Seqera, other HPCs) without a local SIF. `containers/Dockerfile` is the registry twin
of `strainqc.def`; `.github/workflows/build-container.yml` builds and pushes it to
**`ghcr.io/digibio-ut/strainqc`** on GitHub's runners (no local Docker needed).

### One-time setup
1. Push this branch / merge to `main` (or trigger the workflow manually from the repo's
   **Actions → build-container → Run workflow**). Watch it build (~15–30 min).
2. **Make the package public** (one time): repo/org **Packages → strainqc → Package settings
   → Change visibility → Public**. Public ghcr packages are free and unlimited to pull.
   (The source repo can stay private; only the image package is public.)

Tags produced: `latest` (from `main`), `vX.Y.Z` (from git tags `v*`), and a `sha-…` per build.

### Use it
```bash
# Anywhere with a container runtime — no SIF file needed:
nextflow run Digibio-UT/StrainQC-nf -profile test,apptainer \
    --strainqc_container docker://ghcr.io/digibio-ut/strainqc:latest
```
On the HPC you can also pre-pull it into a SIF once and point at that file:
```bash
apptainer pull containers/strainqc.sif docker://ghcr.io/digibio-ut/strainqc:latest
```
Cloud/Seqera executors just take the `docker://…` URI directly. The big reference databases
(Kraken2, MASH) are **not** in the image — they're staged separately.
