#!/usr/bin/env python3
"""SNV/INDEL concordance of StrainQC-nf calls vs Ottilie et al. (2022) Supplementary Data 4.

For each evolved pilot clone we:
  1. load its results.json `mutations[]` and the parent's,
  2. parent-subtract (evolved calls not shared with the parent — the paper's method),
  3. match the remaining evolved-only calls against that clone's Sup-4 truth rows,
  4. report recall / precision + a per-variant status CSV,
  5. spotlight the biologically important driver alleles (incl. the CRISPR-validated
     PMS1 Lys724* in the doxorubicin clone and the flagship YRM1 allele in carmaphycin).

Pure stdlib (uses ./_xlsx.py). Run after the pipeline produces results.json for all four
pilot samples. CNV/aneuploidy for CBR110 is handled separately by check_cnv.py.
"""

import argparse
import csv
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _xlsx import read_records  # noqa: E402

INDEL_POS_TOL = 3   # bp tolerance for indel anchor-position representation differences

# sample id (samplesheet) -> Sup-4 clone name + role. CBR110 is CNV-only (see check_cnv.py).
PILOT = {
    "NODRUG-GM2":         {"role": "parent",  "sup4_clone": None,                 "eaw": None},
    "Doxorubicin-16-R2b": {"role": "evolved", "sup4_clone": "Doxorubicin-16--R2b", "eaw": "EAW304"},
    "Carmaphycin-R9-2":   {"role": "evolved", "sup4_clone": "Carmaphycin--R9-2",   "eaw": "EAW131"},
}

# Driver alleles to spotlight: (sample, chrom, position, gene, label, crispr_validated)
SPOTLIGHT = [
    ("Doxorubicin-16-R2b", "XIV", 475560, "PMS1", "PMS1 Lys724* (mismatch-repair KO; mutator)", True),
    ("Doxorubicin-16-R2b", "VII", 470822, "PDR1", "PDR1 frameshift (multidrug-resistance TF)", False),
    ("Carmaphycin-R9-2",   "XV",  656188, "YRM1", "YRM1 Asn660Ile (flagship YRR1/YRM1 TF family)", False),
]


def norm_clone(name):
    """Collapse the inconsistent dash formatting across sheets/SRA names."""
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


def load_mutations(results_json):
    with open(results_json) as fh:
        data = json.load(fh)
    muts = []
    for m in data.get("mutations", []):
        try:
            # breseq formats positions with thousands separators ("257,908"); GATK does not.
            pos = int(str(m.get("position", "")).replace(",", "").strip())
        except ValueError:
            continue
        muts.append({
            "chrom": str(m.get("seq_id", "")).strip(),
            "pos": pos,
            "type": (m.get("type") or "").upper(),
            "mutation": m.get("mutation", ""),
            "gene": m.get("gene", "") or "",
            "freq": m.get("freq", ""),
        })
    return muts


def is_indel(t):
    return t in ("INS", "DEL", "SUB")


def snp_alt(mutation):
    """'T→A' -> ('T','A'); tolerate ascii '->'."""
    s = mutation.replace("->", "→")
    if "→" in s:
        ref, alt = s.split("→", 1)
        return ref.strip(), alt.strip()
    return None, None


def parent_subtract(evolved, parent):
    """Drop evolved calls that also appear in the parent (same chrom, ~pos, type class)."""
    p_snp = {(m["chrom"], m["pos"]) for m in parent if not is_indel(m["type"])}
    p_indel = [(m["chrom"], m["pos"]) for m in parent if is_indel(m["type"])]
    out = []
    for m in evolved:
        if is_indel(m["type"]):
            shared = any(c == m["chrom"] and abs(p - m["pos"]) <= INDEL_POS_TOL
                         for c, p in p_indel)
        else:
            shared = (m["chrom"], m["pos"]) in p_snp
        if not shared:
            out.append(m)
    return out


def load_truth(sup4_path, clone_name):
    _, recs = read_records(sup4_path, header_row=1)
    key = norm_clone(clone_name)
    rows = []
    for r in recs:
        if norm_clone(r.get("Clone Name", "")) != key:
            continue
        try:
            pos = int(str(r.get("Mutation Position", "")).strip())
        except ValueError:
            continue
        rows.append({
            "chrom": str(r.get("Chromosome", "")).strip(),
            "pos": pos,
            "ref": str(r.get("Reference Base", "")).strip(),
            "alt": str(r.get("Alternate Base", "")).strip(),
            "type": str(r.get("Type", "")).strip().upper(),  # SNP | INDEL
            "gene": str(r.get("Standard Name", "")).strip(),
            "aa": str(r.get("Amino Acid Change", "")).strip(),
        })
    return rows


def match(truth, calls):
    """Match truth rows to evolved-only calls. Returns per-truth status + unmatched calls."""
    calls = list(calls)
    used = [False] * len(calls)
    results = []
    for t in truth:
        hit = None
        t_indel = t["type"] != "SNP"
        for i, c in enumerate(calls):
            if used[i] or c["chrom"] != t["chrom"]:
                continue
            if not t_indel and not is_indel(c["type"]):
                if c["pos"] == t["pos"]:
                    _, calt = snp_alt(c["mutation"])
                    if calt is None or calt == t["alt"] or t["alt"] == "":
                        hit = i
                        break
            elif t_indel and is_indel(c["type"]):
                if abs(c["pos"] - t["pos"]) <= INDEL_POS_TOL:
                    hit = i
                    break
        if hit is not None:
            used[hit] = True
            results.append((t, "recovered", calls[hit]))
        else:
            results.append((t, "missed", None))
    fp = [calls[i] for i in range(len(calls)) if not used[i]]
    return results, fp


def main():
    ap = argparse.ArgumentParser(description="StrainQC-nf vs Ottilie Sup 4 SNV/INDEL concordance")
    ap.add_argument("--results-dir", required=True,
                    help="pipeline outdir; expects <dir>/<sample_id>/results.json")
    ap.add_argument("--sup4", required=True, help="Supplementary Data 4 .xlsx")
    ap.add_argument("--out-dir", default="out")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    def rj(sample):
        return os.path.join(args.results_dir, sample, "results.json")

    parent_muts = load_mutations(rj("NODRUG-GM2"))
    print(f"Parent NODRUG-GM2: {len(parent_muts)} raw calls vs S288C\n")

    per_variant = []
    summary = []
    for sample, cfg in PILOT.items():
        if cfg["role"] != "evolved":
            continue
        calls = load_mutations(rj(sample))
        evolved_only = parent_subtract(calls, parent_muts)
        truth = load_truth(args.sup4, cfg["sup4_clone"])
        results, fp = match(truth, evolved_only)

        tp = sum(1 for _, s, _ in results if s == "recovered")
        recall = tp / len(truth) if truth else 0.0
        precision = tp / len(evolved_only) if evolved_only else 0.0
        summary.append((sample, cfg["sup4_clone"], len(truth), tp, len(evolved_only),
                        round(recall, 3), round(precision, 3)))

        for t, status, c in results:
            per_variant.append({
                "sample": sample, "chrom": t["chrom"], "position": t["pos"],
                "type": t["type"], "gene": t["gene"], "aa_change": t["aa"],
                "ref": t["ref"], "alt": t["alt"], "status": status,
                "called_as": (c["mutation"] if c else ""), "called_pos": (c["pos"] if c else ""),
            })
        for c in fp:
            per_variant.append({
                "sample": sample, "chrom": c["chrom"], "position": c["pos"],
                "type": c["type"], "gene": c["gene"], "aa_change": "",
                "ref": "", "alt": "", "status": "extra_call_not_in_truth",
                "called_as": c["mutation"], "called_pos": c["pos"],
            })

        print(f"== {sample}  ({cfg['sup4_clone']}, {cfg['eaw']}) ==")
        print(f"   truth={len(truth)}  recovered={tp}  recall={recall:.0%}"
              f"  evolved_only_calls={len(evolved_only)}  precision={precision:.0%}")
        missed = [t for t, s, _ in results if s == "missed"]
        if missed:
            print("   missed:")
            for t in missed:
                print(f"     {t['chrom']}:{t['pos']} {t['type']} {t['gene']} {t['aa'][:24]}")
        print()

    # driver / CRISPR-validated spotlight
    print("== Driver-allele spotlight ==")
    spot_rows = []
    for sample, chrom, pos, gene, label, crispr in SPOTLIGHT:
        calls = load_mutations(rj(sample))
        evolved_only = parent_subtract(calls, parent_muts)
        found = next((c for c in evolved_only
                      if c["chrom"] == chrom and abs(c["pos"] - pos) <= INDEL_POS_TOL), None)
        status = "RECOVERED" if found else "MISSED"
        tag = " [CRISPR-validated]" if crispr else ""
        print(f"   [{status}] {sample}: {label}{tag}"
              f"  ({chrom}:{pos}{'  called '+found['mutation'] if found else ''})")
        spot_rows.append({"sample": sample, "chrom": chrom, "position": pos, "gene": gene,
                          "label": label, "crispr_validated": crispr, "status": status,
                          "called_as": found["mutation"] if found else ""})

    # write CSVs
    snv_csv = os.path.join(args.out_dir, "snv_concordance.csv")
    with open(snv_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(per_variant[0].keys()))
        w.writeheader()
        w.writerows(per_variant)
    sum_csv = os.path.join(args.out_dir, "snv_summary.csv")
    with open(sum_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "clone", "truth", "recovered", "evolved_only_calls",
                    "recall", "precision"])
        w.writerows(summary)
    spot_csv = os.path.join(args.out_dir, "driver_spotlight.csv")
    with open(spot_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(spot_rows[0].keys()))
        w.writeheader()
        w.writerows(spot_rows)
    print(f"\nWrote {snv_csv}, {sum_csv}, {spot_csv}")


if __name__ == "__main__":
    main()
