# StrainQC-nf

Nextflow re-platform of the **StrainQC / Contamination Detection Platform** (DTU + UT DigiBio).
Sequencing-based strain quality control: variant calling (**breseq** or **GATK**) + contamination
detection (**Kraken2 / MASH / BLAST**). Portable across **UT HPC (Slurm), Seqera, and cloud** via
Nextflow profiles + containers — replacing the working app's bespoke `worker.py`/`slurm_runner.py`
SSH-and-sbatch executor.

> **Status: v0.2.0, validated on UT HPC.** Both callers (breseq, GATK) run end-to-end through
> the React UI with full feature parity, plus CNVkit aneuploidy calling (`--parent_sample`).
> Reproduces Ottilie et al. (2022): 96–100% SNV recall, the CRISPR-validated *PMS1* driver, and
> the chr I aneuploidy ([`tests/ottilie_benchmark`](tests/ottilie_benchmark/)). Unit, contract
> and pipeline tests run in CI on every push ([`tests/README.md`](tests/README.md)).

This is a **separate, clean repo** from the working app at `../SlurmApp` — which stays the
production source of truth until this reaches parity.

## Design

- **Contract-first.** The React frontend + FastAPI serving layer are reused unchanged, so the
  pipeline must emit `results.json` exactly as specified in [`docs/API_CONTRACT.md`](docs/API_CONTRACT.md).
  `EMIT_RESULTS` is the single process that produces it. (`manta_vcf` was removed upstream — do not emit it.)
- **Pre-compute for cloud.** Coverage downsampling and per-variant read-evidence JSON are produced
  *in-pipeline* (`SAMTOOLS_DEPTH`, `VARIANT_EVIDENCE`), so the API only serves files rather than
  re-reading BAMs on demand (the egress-cost concern on object storage).
- **Tower-compatible, not Tower-dependent.** Plain nf-core-style; runs with `nextflow run` alone,
  and also under Seqera Platform.

```
main.nf                    entry: parse samplesheet, dispatch by --variant_caller
workflows/strainqc.nf      wires the breseq | gatk paths → EMIT_RESULTS
modules/local/*.nf         breseq · gatk · qc · contamination · report  (grounded in wrappers.py)
bin/emit_results.py        results.json assembler (port of parsers.py)  ← STUB
conf/{base,test}.config    resources by label; CI test profile
nextflow.config            params, container tags, execution profiles
docs/API_CONTRACT.md       the frozen contract this pipeline must satisfy
```

## Usage

```bash
nextflow run . -profile docker \
  --input samplesheet.csv --reference ref.gbk \
  --variant_caller gatk --ploidy 1 [--parent_sample <ancestor_id>]
```

Samplesheet: `sample,fastq_1,fastq_2` (see [`assets/samplesheet.example.csv`](assets/samplesheet.example.csv)).
All parameters, with defaults and help, are in [`nextflow_schema.json`](nextflow_schema.json)
(also what Seqera Platform builds its launch form from).

Profiles (combine with commas):

| Where | Profile |
|---|---|
| UT HPC | `slurm,apptainer` (or `hpc` for both), `slurm,hpc_modules` |
| Seqera + Azure Batch | `azure`, see [`docs/deploy_seqera_azure.md`](docs/deploy_seqera_azure.md) |
| Laptop / server | `docker`, `apptainer`, `singularity`, `conda` |
| Smoke test | `test` (tiny bundled dataset) |

The default tool image is the published one for this release,
`ghcr.io/digibio-ut/strainqc:<version>`.

## Setup — reference databases (per site)

Large DBs are **never committed**. Each site stages them once and points params at them:

- `--kraken2_db <dir>` — Kraken2 database
- `--mash_db <refseq.msh>` — MASH sketch
- `--run_blast` — remote BLASTN fallback (slow; needs internet)

A dedicated `setup` step for fresh cloud sites is a **TODO** (drives Vallo's storage-cost estimate).

## Roadmap

Done: parsers ported (`bin/emit_results.py`), wiring validated, thin launcher behind the
unchanged React UI, single all-in-one image on ghcr.io, test profile + CI
([`tests/README.md`](tests/README.md)), Ottilie benchmark, CNVkit aneuploidy, Seqera launch
schema, `azure` profile.

Next:
1. **Web app on Seqera + Azure**: start runs through the Seqera API and read results from Azure
   storage instead of a local `nextflow run` + local disk.
2. **Auth / multi-tenancy** before any shared instance (current app: `allow_origins=["*"]`, no auth).
3. **Reference-database setup step** for fresh cloud sites (Kraken2 / MASH staging).
4. **Runtime parameter validation** with the nf-schema plugin (the schema is already in place).
