#!/usr/bin/env python3
"""
emit_results.py — assemble results.json for StrainQC-nf.

This is the ONE place that defines the results.json contract on any backend.
Port the parsing logic from the working app's backend/app/pipeline/parsers.py:
  - breseq: parse output/index.html + output.gd  -> mutations[]
  - GATK:   parse VCF                             -> mutations[] (+ evidence_links)
  - kraken.report                                 -> contaminants[]
  - mash_screen.tab                               -> mash_results[] / mash_unmapped_results[]
  - samtools depth (coverage.tsv)                 -> downsample -> coverage_file
  - samtools flagstat                             -> read_stats{}

Contract reference: docs/API_CONTRACT.md  (do NOT emit `manta_vcf` — removed upstream).

STATUS: stub. Wire real parsing before first production run.
"""
import argparse, json, sys


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--job-id", required=True)
    p.add_argument("--caller", required=True)          # breseq | gatk
    p.add_argument("--polymorphism", default="false")
    p.add_argument("--variant-src", required=True)     # breseq output/ dir OR gatk vcf
    p.add_argument("--coverage")
    p.add_argument("--flagstat")
    p.add_argument("--kraken")
    p.add_argument("--mash")
    p.add_argument("--evidence-dir")
    p.add_argument("--coverage-resolution", type=int, default=2000)
    p.add_argument("--out", default="results.json")
    args = p.parse_args()

    # TODO: replace stub with real parsers.py logic.
    results = {
        "job_id": args.job_id,
        "status": "completed",
        "variant_caller": args.caller,
        "polymorphism_mode": args.polymorphism.lower() == "true",
        "mutations": [],
        # "read_stats": {...}, "coverage_file": ..., "contaminants": [...],
        # "mash_results": [...], "unmapped_files": [...], "hc_snps_file": ...,
    }

    with open(args.out, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"[emit_results] wrote {args.out} (STUB)", file=sys.stderr)


if __name__ == "__main__":
    main()
