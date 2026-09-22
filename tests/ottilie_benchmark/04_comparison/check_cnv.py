#!/usr/bin/env python3
"""Aneuploidy concordance: StrainQC-nf CNVkit calls vs Ottilie Supplementary Data 5.

Sup 5 lists CBR110-15R3a as a chr I "Whole Chromosome Duplication". This checks that the
new in-pipeline CNVkit step (results.json -> chromosome_copy_number / aneuploidy) recovers
that chr I gain and does not fabricate aneuploidy on other chromosomes.

Pure stdlib (uses ./_xlsx.py).
"""

import argparse
import csv
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _xlsx import read_records  # noqa: E402

# sample id (samplesheet) -> Sup-5 clone name
CNV_PILOT = {"CBR110-15-R3a": "CBR110-15R3a"}


def norm_clone(name):
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


def load_truth_cnv(sup5_path, clone_name):
    _, recs = read_records(sup5_path, header_row=1)
    key = norm_clone(clone_name)
    out = []
    for r in recs:
        if norm_clone(r.get("Clone name", "")) != key:
            continue
        out.append({
            "chrom": str(r.get("Chromosome", "")).strip(),
            "event": str(r.get("Event type", "")).strip(),
            "genes": str(r.get("Genes involved in CN", "")).strip(),
        })
    return out


def main():
    ap = argparse.ArgumentParser(description="CNVkit aneuploidy vs Ottilie Sup 5")
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--sup5", required=True)
    ap.add_argument("--out-dir", default="out")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    all_rows = []
    overall_ok = True
    for sample, clone in CNV_PILOT.items():
        with open(os.path.join(args.results_dir, sample, "results.json")) as fh:
            data = json.load(fh)
        called = {r["seq_id"]: r for r in data.get("chromosome_copy_number", [])}
        called_aneu = {a["seq_id"]: a["call"] for a in data.get("aneuploidy", [])}
        truth = load_truth_cnv(args.sup5, clone)
        truth_gain = {t["chrom"] for t in truth if "dup" in t["event"].lower()
                      or "gain" in t["event"].lower() or "amplif" in t["event"].lower()}

        print(f"== {sample}  ({clone}) ==")
        print(f"   Sup 5 truth: {[(t['chrom'], t['event']) for t in truth]}")
        print(f"   CNVkit aneuploidy calls: {called_aneu or '(none)'}")

        # recall: every truth aneuploidy chromosome must be called a gain
        missed = [c for c in truth_gain if called_aneu.get(c) != "gain"]
        # false positives: aneuploidy calls on chromosomes not in truth
        false_pos = [c for c in called_aneu if c not in truth_gain]
        recovered = [c for c in truth_gain if called_aneu.get(c) == "gain"]

        status = "PASS" if not missed and not false_pos else "CHECK"
        overall_ok = overall_ok and status == "PASS"
        print(f"   recovered gains: {recovered or '(none)'}"
              f"   missed: {missed or '(none)'}   false-positive aneuploidy: {false_pos or '(none)'}")
        print(f"   -> {status}\n")

        for chrom, row in called.items():
            all_rows.append({
                "sample": sample, "seq_id": chrom, "mean_cn": row.get("mean_cn"),
                "call": row.get("call"), "aneuploidy": row.get("aneuploidy"),
                "in_sup5_truth": chrom in truth_gain,
            })

    out_csv = os.path.join(args.out_dir, "cnv_concordance.csv")
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(all_rows[0].keys()))
        w.writeheader()
        w.writerows(all_rows)
    print(f"Wrote {out_csv}")
    sys.exit(0 if overall_ok else 1)


if __name__ == "__main__":
    main()
