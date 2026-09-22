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

# Whole-chromosome aneuploidy thresholds, RELATIVE to the genome-wide median chromosome CN
# (the euploid baseline for this sample). Using a relative test makes calls robust to MOSAIC
# aneuploidies — an evolved population that is only partly 2n for a chromosome shows a mean CN
# of ~1.3x rather than 2x (e.g. Ottilie CBR110 chr I), which a fixed 2x cutoff would miss.
GAIN_RATIO = 1.20   # chromosome mean CN >= 1.20 x genome-median  -> gain
LOSS_RATIO = 0.75   # chromosome mean CN <= 0.75 x genome-median  -> loss
# Fraction of a chromosome's length that must agree in direction to call it WHOLE-chromosome
# (vs a focal/intrachromosomal event, which we report but do not flag as aneuploidy).
WHOLE_CHROM_FRACTION = 0.5
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


def _median(vals):
    s = sorted(vals)
    n = len(s)
    if n == 0:
        return 0.0
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2.0


def summarise(cns_path, fai_path, baseline):
    lengths = read_fai(fai_path)
    # Pass 1: length-weighted mean CN per chromosome + keep segments for the fraction test.
    agg = {}  # chrom -> {cn_len, len, n, segs:[(len,cn)]}
    for chrom, start, end, log2 in read_cns(cns_path):
        seg_len = max(end - start, 0)
        if seg_len == 0:
            continue
        cn = baseline * (2.0 ** log2)
        a = agg.setdefault(chrom, {"cn_len": 0.0, "len": 0, "n": 0, "segs": []})
        a["cn_len"] += cn * seg_len
        a["len"] += seg_len
        a["n"] += 1
        a["segs"].append((seg_len, cn))

    mean_cn = {c: a["cn_len"] / (a["len"] or 1) for c, a in agg.items()}
    # Euploid baseline for THIS sample = median of nuclear chromosome mean CNs.
    genome_median = _median([v for c, v in mean_cn.items() if c not in NON_NUCLEAR]) or baseline

    rows, aneuploidies = [], []
    for chrom in lengths:                       # reference (fai) order → stable output
        if chrom not in agg:
            continue
        a = agg[chrom]
        covered = a["len"] or 1
        mcn = mean_cn[chrom]
        gain_len = sum(l for l, cn in a["segs"] if cn >= GAIN_RATIO * genome_median)
        loss_len = sum(l for l, cn in a["segs"] if cn <= LOSS_RATIO * genome_median)
        gain_frac, loss_frac = gain_len / covered, loss_len / covered
        altered_frac = max(gain_frac, loss_frac)

        call, aneuploidy = "neutral", False
        if chrom not in NON_NUCLEAR:
            if mcn >= GAIN_RATIO * genome_median and gain_frac >= WHOLE_CHROM_FRACTION:
                call, aneuploidy = "gain", True
            elif mcn <= LOSS_RATIO * genome_median and loss_frac >= WHOLE_CHROM_FRACTION:
                call, aneuploidy = "loss", True
            elif altered_frac > 0:
                call = "focal"  # sub-chromosomal amp/del, reported but not an aneuploidy

        row = {
            "seq_id": chrom,
            "mean_cn": round(mcn, 3),
            "cn_vs_genome_median": round(mcn / genome_median, 3) if genome_median else 0,
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
