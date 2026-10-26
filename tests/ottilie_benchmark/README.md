# Ottilie et al. (2022) benchmark — 4-sample pilot

Reproduces published variant calls from **Ottilie et al., *Commun Biol* 5:128 (2022)**
([doi:10.1038/s42003-022-03076-7](https://doi.org/10.1038/s42003-022-03076-7)), a yeast
adaptive-laboratory-evolution study of xenobiotic resistance. Goal: confirm StrainQC-nf
reproduces a peer-reviewed analysis, not just that it runs.

The paper called variants with **GATK HaplotypeCaller on S288C R64-1-1** and parent-subtracted
— the same recipe as our GATK path — so this is an apples-to-apples comparison.

## Samples (BioProject PRJNA590203)

| sample_id | SRR | role | published truth |
|---|---|---|---|
| `NODRUG-GM2` | SRR10985539 | parent (ABC16-GM ancestor) | subtraction baseline |
| `Doxorubicin-16-R2b` | SRR10985527 | evolved | **23** SNV/INDEL (Sup 4, EAW304); incl. CRISPR-validated **PMS1 Lys724\*** + **PDR1** frameshift |
| `Carmaphycin-R9-2` | SRR10985678 | evolved | **15** SNV (Sup 4, EAW131); incl. flagship TF **YRM1 Asn660Ile** |
| `CBR110-15-R3a` | SRR10985585 | evolved | **chr I whole-chromosome duplication** (Sup 5) |

## What it checks

1. **SNV/INDEL concordance vs Supplementary Data 4** — per clone, recall (fraction of the
   paper's calls we recover) + precision, after parent-subtraction. (`compare_snvs.py`)
2. **Driver / CRISPR-validated alleles** — do we call the biologically causal mutations,
   including the CRISPR-confirmed `PMS1 Lys724*` and the flagship `YRM1` TF allele. (`compare_snvs.py`)
3. **Aneuploidy vs Supplementary Data 5** — does the new in-pipeline **CNVkit** step call the
   chr I duplication in CBR110 (and nothing spurious). (`check_cnv.py`)

## Run it

Compute runs on the **UT HPC via Slurm** (login-node compute is prohibited). The CNV step
needs the container rebuilt **with cnvkit** (this branch adds it to `containers/`).

```bash
# 0. one-time: stage the S288C reference on the HPC (already built under ottilie_test_data/)
#    rsync -a ottilie_test_data/S288C_reference <user>@rocket:/gpfs/.../ && export REF=/gpfs/.../S288C_reference/S288C_R64_ensembl_chrnames.gb

# 1. truth set (public supplementary xlsx)
bash 01_data/download_truth_set.sh

# 2. FASTQs + samplesheet (ENA direct, ~2 GB)
bash 01_data/download_pilot_fastq.sh

# 3. pipeline (Slurm + Apptainer, GATK ploidy 1, CNVkit on via --parent_sample)
module load nextflow/25.10.2 apptainer/1.4.1
bash 03_run/run_pilot.sh          # honours REF / CONTAINER / PROFILES / OUTDIR env vars

# 4. concordance + report
cd 04_comparison
python3 compare_snvs.py --results-dir ../results/ottilie_pilot --sup4 ../data/supplementary/sup4_snv_indel.xlsx
python3 check_cnv.py   --results-dir ../results/ottilie_pilot --sup5 ../data/supplementary/sup5_cnv.xlsx
python3 make_report.py                                 # -> out/SUMMARY.md
```

All comparison scripts are pure stdlib (`_xlsx.py` reads the spreadsheets — no pandas/openpyxl).

## Acceptance targets

- SNV recall ≥ 80% vs Sup 4 for the two SNV clones; documented reason for any miss
  (StrainQC-nf has no post-call filtering, so expect a few extra calls in repeat/subtelomeric
  regions — reported, not silently dropped).
- CRISPR-validated `PMS1 Lys724*` recovered.
- chr I aneuploidy detected in CBR110, no false aneuploidy elsewhere.

## Layout

```
01_data/    download_truth_set.sh, download_pilot_fastq.sh   (write into data/, gitignored)
03_run/     run_pilot.sh
04_comparison/  _xlsx.py, compare_snvs.py, check_cnv.py, make_report.py
```

Downloaded reads, truth-set xlsx, the generated samplesheet, and pipeline `results/` all land
under `data/` and `results/` (gitignored). The earlier `tests/ottilie_xenobiotic_ale/` scaffold
in the old SlurmApp repo (a Sarek-fork design) is superseded by this StrainQC-nf-native harness.
