# StrainQC-nf — Ottilie et al. (2022) benchmark

Reproduction of published variant calls for four pilot clones (*Commun Biol* 5:128), calling each sample vs S288C R64-1-1 with the GATK path (HaplotypeCaller, ploidy 1) and parent-subtracting against NODRUG-GM2.

## SNV/INDEL concordance vs Supplementary Data 4

| Sample | Clone | Truth | Recovered | Recall | Extra calls | Precision |
|---|---|--:|--:|--:|--:|--:|
| Doxorubicin-16-R2b | Doxorubicin-16--R2b | 23 | 22 | 96% | 96 | 19% |
| Carmaphycin-R9-2 | Carmaphycin--R9-2 | 15 | 15 | 100% | 41 | 27% |
| CBR110-15-R3a | CBR110-15R3a | 4 | 4 | 100% | 35 | 10% |

## Driver-allele spotlight

| Sample | Allele | CRISPR-validated | Status | Called as |
|---|---|:--:|:--:|---|
| Doxorubicin-16-R2b | PMS1 Lys724* (mismatch-repair KO; mutator) | yes | **RECOVERED** | A→T |
| Doxorubicin-16-R2b | PDR1 frameshift (multidrug-resistance TF) | — | **RECOVERED** | Δ1 bp |
| Carmaphycin-R9-2 | YRM1 Asn660Ile (flagship YRR1/YRM1 TF family) | — | **RECOVERED** | A→T |
| CBR110-15-R3a | YRR1 Thr711Asn (flagship YRR1/YRM1 TF family) | — | **RECOVERED** | G→T |

## Aneuploidy vs Supplementary Data 5 (CNVkit)

| Sample | Chromosome | Mean CN | Call | In Sup 5 truth |
|---|---|--:|---|:--:|
| CBR110-15-R3a | I | 1.285 | gain **(aneuploidy)** | ✓ |
| CBR110-15-R3a | II | 0.989 | focal |  |
| CBR110-15-R3a | III | 1.126 | focal |  |
| CBR110-15-R3a | IV | 0.998 | focal |  |
| CBR110-15-R3a | V | 0.995 | focal |  |
| CBR110-15-R3a | VI | 1.141 | neutral |  |
| CBR110-15-R3a | VII | 1.022 | focal |  |
| CBR110-15-R3a | VIII | 0.974 | neutral |  |
| CBR110-15-R3a | IX | 1.029 | focal |  |
| CBR110-15-R3a | X | 1.007 | focal |  |
| CBR110-15-R3a | XI | 0.984 | neutral |  |
| CBR110-15-R3a | XII | 0.978 | neutral |  |
| CBR110-15-R3a | XIII | 1.012 | neutral |  |
| CBR110-15-R3a | XIV | 0.982 | neutral |  |
| CBR110-15-R3a | XV | 1.019 | focal |  |
| CBR110-15-R3a | XVI | 1.001 | focal |  |
