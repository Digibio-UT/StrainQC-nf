#!/usr/bin/env python3
"""Summarise a CNVkit `.cns` segment file to per-chromosome copy number + aneuploidy calls.

CNVkit reports log2 copy ratios relative to its reference. Here the reference is built from
the un-evolved PARENT (see modules/local/cnv.nf), so log2 == 0 means "same dosage as the
parent". Absolute copy number is therefore:

    copy_number = baseline_ploidy * 2**log2

For the Ottilie haploid yeast benchmark baseline_ploidy == 1, so a whole-chromosome
duplication (2 copies vs the parent's 1) gives log2 ≈ 1.0 → copy_number ≈ 2.0.

Pure stdlib (matches emit_results.py — the pipeline image has only a base Python besides the
bioconda tools). Output is a small JSON that emit_results.py folds into results.json.
"""

import argparse
import json
import sys

# Whole-chromosome aneuploidy thresholds (relative to baseline ploidy).
GAIN_RATIO = 1.4    # mean CN >= 1.4 x baseline  -> gain
LOSS_RATIO = 0.6    # mean CN <= 0.6 x baseline  -> loss
# Fraction of a chromosome's length that must agree in direction to call it WHOLE-chromosome
# (vs a focal/intrachromosomal event, which we report but do not flag as aneuploidy).
WHOLE_CHROM_FRACTION = 0.7
# Organellar contigs never get an aneuploidy call (dosage is not nuclear ploidy).
NON_NUCLEAR = {"Mito", "MT", "chrMito", "chrM", "Mitochondria"}


def read_fai(path):
    """chromosome -> length, preserving file order."""
    lengths = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                lengths[parts[0]] = int(parts[1])
    return lengths


def read_cns(path):
    """Yield (chrom, start, end, log2) from a CNVkit .cns file (header-driven)."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        try:
            ci, si, ei, li = idx["chromosome"], idx["start"], idx["end"], idx["log2"]
        except KeyError:
            raise SystemExit(f"[summarize_cnv] unexpected .cns header: {header}")
        for line in fh:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            try:
                yield p[ci], int(p[si]), int(p[ei]), float(p[li])
            except (ValueError, IndexError):
                continue


def summarise(cns_path, fai_path, baseline):
    lengths = read_fai(fai_path)
    # Accumulate length-weighted copy number and altered length per chromosome.
    agg = {}  # chrom -> dict
    for chrom, start, end, log2 in read_cns(cns_path):
        seg_len = max(end - start, 0)
        if seg_len == 0:
            continue
        cn = baseline * (2.0 ** log2)
        a = agg.setdefault(chrom, {"cn_len": 0.0, "len": 0, "gain_len": 0,
                                   "loss_len": 0, "n": 0})
        a["cn_len"] += cn * seg_len
        a["len"] += seg_len
        a["n"] += 1
        if cn >= GAIN_RATIO * baseline:
            a["gain_len"] += seg_len
        elif cn <= LOSS_RATIO * baseline:
            a["loss_len"] += seg_len

    rows, aneuploidies = [], []
    # Iterate in reference (fai) order for stable, genome-ordered output.
    for chrom in lengths:
        if chrom not in agg:
            continue
        a = agg[chrom]
        covered = a["len"] or 1
        mean_cn = a["cn_len"] / covered
        gain_frac = a["gain_len"] / covered
        loss_frac = a["loss_len"] / covered
        altered_frac = max(gain_frac, loss_frac)

        call = "neutral"
        aneuploidy = False
        if chrom not in NON_NUCLEAR:
            if mean_cn >= GAIN_RATIO * baseline and gain_frac >= WHOLE_CHROM_FRACTION:
                call, aneuploidy = "gain", True
            elif mean_cn <= LOSS_RATIO * baseline and loss_frac >= WHOLE_CHROM_FRACTION:
                call, aneuploidy = "loss", True
            elif altered_frac > 0:
                call = "focal"  # a sub-chromosomal amplification/deletion, not aneuploidy

        row = {
            "seq_id": chrom,
            "mean_cn": round(mean_cn, 3),
            "length": lengths[chrom],
            "n_segments": a["n"],
            "altered_fraction": round(altered_frac, 3),
            "call": call,
            "aneuploidy": aneuploidy,
        }
        rows.append(row)
        if aneuploidy:
            aneuploidies.append({"seq_id": chrom, "call": call, "mean_cn": row["mean_cn"]})
    return rows, aneuploidies


def main():
    ap = argparse.ArgumentParser(description="CNVkit .cns -> per-chromosome CN + aneuploidy")
    ap.add_argument("--cns", required=True)
    ap.add_argument("--fai", required=True, help="reference .fai (chrom lengths + order)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--baseline-ploidy", type=float, default=1.0)
    ap.add_argument("--out", default="cnv.json")
    args = ap.parse_args()

    rows, aneuploidies = summarise(args.cns, args.fai, args.baseline_ploidy)
    out = {
        "cnv_tool": "cnvkit",
        "sample": args.sample,
        "baseline_ploidy": args.baseline_ploidy,
        "chromosome_copy_number": rows,
        "aneuploidy": aneuploidies,
    }
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"[summarize_cnv] {args.sample}: {len(rows)} chromosomes, "
          f"{len(aneuploidies)} aneuploidy call(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
