# Running & testing StrainQC-nf

Two tiers: cheap wiring validation (`-stub`, no tools) and real runs (Apptainer on the HPC).

## 1. Validate the wiring — `-stub` (run this first)

Runs every process's `stub:` block (touch/echo dummy outputs) — no containers, no bio tools,
seconds. Proves the DSL parses and the channel graph connects end to end.

```bash
module load nextflow/25.10.2
cd StrainQC-nf
nextflow run . -profile test -stub
```

Expect: all processes complete, `results/testsample/results.json` exists (stub content).
If it errors, the message points at the offending process/channel — paste it back.

## 2. Real smoke test on the HPC — modules + synthetic data (current path)

Uses the tiny bundled dataset (`assets/test/`, one planted SNP at CHR1:1000) and the
proven `module load` stack (no containers yet). Multi-tool steps work because the
`hpc_modules` beforeScript loads the whole environment.

```bash
module load nextflow/25.10.2
cd StrainQC-nf
nextflow run . -profile test,slurm,hpc_modules      # one Slurm job per process
```

**Always use the `slurm` executor** — UT HPC rules prohibit running compute on login nodes.
The lightweight Nextflow orchestrator itself sits on the login node (fine); every actual
tool runs as a Slurm job. A Slurm `--account` is mandatory (default `ealloc_bioeng_h-1`);
override with `--hpc_account <acct> --hpc_queue <partition>`. For long real runs, launch the
orchestrator inside a small `sbatch`/`tmux` session so it survives disconnects.

Check the result against `assets/test/EXPECTED.md`:

```bash
cat results/testsample/results.json    # expect exactly one SNP mutation at position 1000
```

If your conda env / module names differ from the SlurmApp defaults, override the whole
activation line:

```bash
nextflow run . -profile test,hpc_modules \
    --hpc_pre_command 'module load miniconda3 ; conda activate <your_env> ; module load gatk/<ver> bwa-mem2/<ver>'
```

## 3. Containers (portable path — Seqera/cloud)

Instead of `module load`, run every process in one all-in-one Apptainer image (no HPC
modules needed). Build it once (see [`containers/README.md`](../containers/README.md)):

```bash
module load apptainer/1.4.1
apptainer build --fakeroot containers/strainqc.sif containers/strainqc.def
```

Then run with the container profile, pointing at the SIF:

```bash
nextflow run . -profile test,slurm,apptainer \
    --strainqc_container $PWD/containers/strainqc.sif
```

`-profile apptainer` enables the runtime; `--strainqc_container <path>` sets the image for
every process. This is the same pipeline that runs under `hpc_modules` — only the tool
provisioning changes — and it's what makes cloud/Seqera deployment possible (swap the SIF
path for a `docker://…` registry URI).

## 3. A real experiment

```bash
nextflow run . -profile hpc \
    --input  my_samplesheet.csv \
    --reference /path/to/genome.gbk \
    --variant_caller gatk --ploidy 1 \
    --kraken2_db /path/to/kraken_db --mash_db /path/to/refseq.msh \
    --outdir results_myrun
```

Samplesheet columns: `sample,fastq_1,fastq_2` (fastq_2 blank = single-end). Read paths may be
absolute, a URL, or relative to the samplesheet file.

## Regenerating the test dataset

```bash
python3 bin/make_testdata.py --outdir assets/test
```
