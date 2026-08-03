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

## 2. Real smoke test on the HPC — Apptainer + synthetic data

Uses the tiny bundled dataset (`assets/test/`, one planted SNP at CHR1:1000).

```bash
module load nextflow/25.10.2 apptainer/1.4.1
export NXF_APPTAINER_CACHEDIR=$HOME/.apptainer_cache   # reuse SIFs across runs
cd StrainQC-nf

# On a login node this submits nothing heavy, but prefer the Slurm executor:
nextflow run . -profile test,hpc
```

- `-profile test,apptainer` runs the executor locally (fine for this tiny job).
- `-profile test,hpc` submits each process as a Slurm job (Apptainer + `executor=slurm`).

Check the result against `assets/test/EXPECTED.md`:

```bash
cat results/testsample/results.json    # expect exactly one SNP mutation at position 1000
```

> ⚠️ The real run needs the `PREPARE_REFERENCE` container question resolved (it currently
> assumes samtools+bwa-mem2+gatk in one image; see README Roadmap #3). `-stub` does not.

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
