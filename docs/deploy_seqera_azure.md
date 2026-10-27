# Running StrainQC-nf on Seqera Platform + Azure Batch

For sites that launch pipelines from Seqera Platform with an Azure Batch compute environment
(e.g. DTU BRIGHT). The pipeline code is the same everywhere; a site only picks a profile.
Nothing in the repository is site-specific: credentials, the Batch pool and the work directory
come from the Seqera compute environment.

## 1. Before the first launch

| Need | Who | Notes |
|---|---|---|
| Seqera workspace with an **Azure Batch compute environment** | site | Batch Forge is simplest. A pool VM of 16 vCPU / 64 GB (e.g. Standard_D16s_v3) matches the defaults. |
| Azure **storage container** for the work directory and data | site | The compute environment's work dir must be `az://…`. |
| Seqera access to the **GitHub repo** | us + site | The repo is private: add a read-only GitHub credential (fine-grained token, *Contents: read*) to the workspace, or make the repo public. |
| The **container image** | done | `ghcr.io/digibio-ut/strainqc:<version>` is public; Batch nodes pull it directly. If the site only allows its own registry, mirror it there and set `strainqc_container`. |

## 2. Add the pipeline to the Launchpad

- **Pipeline to launch:** `https://github.com/Digibio-UT/StrainQC-nf`
- **Revision:** a release tag, e.g. `v0.2.1`. Pin a tag rather than `main`: each tag runs exactly
  its matching image (`strainqc_container` defaults to `ghcr.io/digibio-ut/strainqc:<version>`).
- **Config profiles:** `azure`
- **Compute environment:** the Azure Batch one from step 1.

The launch form is built from `nextflow_schema.json`: inputs, caller, CNV, contamination and
resource settings, grouped, with defaults and help text. HPC-only settings are hidden.

## 3. First run: smoke test, no data upload

Launch once with profiles **`test,azure`**. It uses the tiny dataset bundled in the repo (one
planted SNP, see `assets/test/EXPECTED.md`) and checks that the image pulls, tasks run on Batch
and `results.json` is written.

## 4. Real runs

1. Upload to the storage container: the reference (GenBank preferred, for gene annotation), the
   FASTQs, and a samplesheet `sample,fastq_1,fastq_2`. FASTQ paths can be full `az://` paths or
   relative to the samplesheet.
2. Set `input` and `reference` to their `az://` paths and choose `variant_caller`
   (`gatk` for yeast, set `ploidy`).
3. Aneuploidy: set `parent_sample` to the ancestor's sample id (GATK route only).
4. Results land in `outdir`, one folder per sample, each with a `results.json` (the format the
   StrainQC web app reads; see `docs/API_CONTRACT.md`).

### Contamination screening

Kraken2 and MASH need reference databases (several GB). Upload them **once** to the storage
account and set `kraken2_db` (the database directory) and `mash_db` (the `.msh` file) to their
`az://` paths. Each screen is skipped if its database is not given.

## 5. Running without Seqera

`nextflow run Digibio-UT/StrainQC-nf -r v0.2.1 -profile azure -w az://<container>/work …` with
the credentials in `AZURE_STORAGE_ACCOUNT_NAME` / `AZURE_STORAGE_ACCOUNT_KEY` and
`AZURE_BATCH_ACCOUNT_NAME` / `AZURE_BATCH_ACCOUNT_KEY`, plus `azure.batch.location` in a
`-c site.config`.

## Troubleshooting

| Message | Cause | Fix |
|---|---|---|
| `Not a valid Azure path: …/work` | the work directory is not on Azure storage | set the compute environment's work dir to `az://…` (or `-w az://…`) |
| `Process requirement exceeds available CPUs/memory` | the pool VM is smaller than the requested ceiling | lower `max_cpus` / `max_memory` to the VM size |
| image pull failure | registry blocked or tag missing | check `strainqc_container`; the tag must exist on ghcr.io |

## Not covered yet: the web app

The StrainQC web app currently starts runs with a local `nextflow run` and reads results from
local disk. On this setup it would start runs through the Seqera API and read `results.json` and
the per-sample files from Azure storage. That is the next integration step and needs the site's
choice of how the app is hosted and authenticated.
