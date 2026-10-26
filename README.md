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

## Usage (once wired)

```bash
nextflow run . -profile docker \
  --input samplesheet.csv --reference ref.gbk \
  --variant_caller gatk --ploidy 1
```

Samplesheet: `sample,fastq_1,fastq_2` (see [`assets/samplesheet.example.csv`](assets/samplesheet.example.csv)).

Profiles: `standard` (local), `docker`, `singularity`, `conda`, `slurm` (UT HPC), `test`.

## Setup — reference databases (per site)

Large DBs are **never committed**. Each site stages them once and points params at them:

- `--kraken2_db <dir>` — Kraken2 database
- `--mash_db <refseq.msh>` — MASH sketch
- `--run_blast` — remote BLASTN fallback (slow; needs internet)

A dedicated `setup` step for fresh cloud sites is a **TODO** (drives Vallo's storage-cost estimate).

## Roadmap

1. ~~**Port parsers** → `bin/emit_results.py`~~ ✅ done (breseq HTML/GD, VCF+GBK annotation,
   kraken, mash, hc-SNP CSV, flagstat read-stats, evidence-link attach; unit-tested on synthetic
   fixtures). Fixed a double-append bug in the original BLAST parser along the way.
2. **Validate the workflow wiring** — `nextflow run . -profile test,docker -stub` (channel
   assembly in `workflows/strainqc.nf` is written but not yet run; Nextflow wasn't available
   at authoring time).
3. **`gbk→fasta` + evidence** helpers (`extract_evidence.py`) into `bin/`; resolve multi-tool
   container for `PREPARE_REFERENCE` (samtools+bwa-mem2+gatk).
4. **Pin container digests** in `nextflow.config` (tags are placeholders — ⚠️ verify).
5. **Test profile** with a tiny bundled dataset; wire CI.
6. **Thin launcher**: adapt the existing FastAPI `/analysis/start` to shell out to `nextflow run`
   (replaces Celery + Slurm SSH executor).
7. **Confirm HPC specifics** with Vallo/Hamza (partition, `--account`).
8. **Auth/multi-tenancy** before any shared instance (current app: `allow_origins=["*"]`, no auth).
```
