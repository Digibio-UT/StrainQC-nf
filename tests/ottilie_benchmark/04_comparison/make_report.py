#!/usr/bin/env python3
"""Assemble SUMMARY.md from the concordance CSVs produced by compare_snvs.py + check_cnv.py."""

import argparse
import csv
import os


def read_csv(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="out")
    args = ap.parse_args()
    d = args.out_dir

    summary = read_csv(os.path.join(d, "snv_summary.csv"))
    spot = read_csv(os.path.join(d, "driver_spotlight.csv"))
    cnv = read_csv(os.path.join(d, "cnv_concordance.csv"))

    L = []
    L.append("# StrainQC-nf — Ottilie et al. (2022) benchmark\n")
    L.append("Reproduction of published variant calls for four pilot clones "
             "(*Commun Biol* 5:128), calling each sample vs S288C R64-1-1 with the GATK path "
             "(HaplotypeCaller, ploidy 1) and parent-subtracting against NODRUG-GM2.\n")

    L.append("## SNV/INDEL concordance vs Supplementary Data 4\n")
    L.append("| Sample | Clone | Truth | Recovered | Recall | Extra calls | Precision |")
    L.append("|---|---|--:|--:|--:|--:|--:|")
    for r in summary:
        recall = f"{float(r['recall'])*100:.0f}%"
        prec = f"{float(r['precision'])*100:.0f}%"
        extra = int(r["evolved_only_calls"]) - int(r["recovered"])
        L.append(f"| {r['sample']} | {r['clone']} | {r['truth']} | {r['recovered']} "
                 f"| {recall} | {extra} | {prec} |")
    L.append("")

    L.append("## Driver-allele spotlight\n")
    L.append("| Sample | Allele | CRISPR-validated | Status | Called as |")
    L.append("|---|---|:--:|:--:|---|")
    for r in spot:
        cv = "yes" if r["crispr_validated"] in ("True", "true", "1") else "—"
        L.append(f"| {r['sample']} | {r['label']} | {cv} | **{r['status']}** | {r['called_as']} |")
    L.append("")

    L.append("## Aneuploidy vs Supplementary Data 5 (CNVkit)\n")
    if cnv:
        L.append("| Sample | Chromosome | Mean CN | Call | In Sup 5 truth |")
        L.append("|---|---|--:|---|:--:|")
        for r in cnv:
            flag = "✓" if r["in_sup5_truth"] in ("True", "true", "1") else ""
            mark = " **(aneuploidy)**" if r["aneuploidy"] in ("True", "true", "1") else ""
            L.append(f"| {r['sample']} | {r['seq_id']} | {r['mean_cn']} "
                     f"| {r['call']}{mark} | {flag} |")
    else:
        L.append("_No CNV table found — run check_cnv.py first._")
    L.append("")

    out = os.path.join(d, "SUMMARY.md")
    with open(out, "w") as fh:
        fh.write("\n".join(L))
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
