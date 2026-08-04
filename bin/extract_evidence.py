#!/usr/bin/env python3
"""
extract_evidence.py — per-variant read-evidence JSON for the GATK evidence viewer.

Port of backend/app/pipeline/wrappers.py :: extract_variant_read_evidence, made standalone
for the VARIANT_EVIDENCE Nextflow process. Parses the VCF for called variants, then for each
window pulls overlapping reads via `samtools view`, builds aligned read/ref strings from the
CIGAR, tallies alt/ref support, and writes evidence/<chrom>_<start>_<end>_<n>.json — the exact
schema the launcher's /evidence HTML renderer + emit_results.attach_evidence_links_from_dir expect.

Region/reference lookups use the RAW VCF CHROM (must match BAM/reference contig names); the
`mutation_key` uses the version-stripped seq_id so it matches results.json mutation keys.

Requires `samtools` on PATH (provided by the hpc_modules beforeScript / the process container).
Usage: extract_evidence.py --bam x.bam --ref reference.fasta --vcf x.vcf --out evidence [--window 75] [--max-reads 80]
"""
import argparse
import json
import os
import re
import subprocess
import sys
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple


def run(cmd: List[str]) -> str:
    return subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True).stdout


def _strip_ver(acc: str) -> str:
    dot = acc.rfind(".")
    return acc[:dot] if dot > 0 and acc[dot + 1:].isdigit() else acc


def _decode_cigar(cigar: str) -> List[Tuple[int, str]]:
    return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar or "")]


def _mutation_key(seq_id_stripped: str, position: int, mtype: str) -> str:
    return f"{seq_id_stripped}_{position}_{mtype}"


def parse_vcf_variants(vcf_path: str) -> List[Dict[str, Any]]:
    """Minimal VCF → variants for evidence: seq_id (raw), position, type, mutation string."""
    variants = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 8:
                continue
            chrom, pos, ref, alt = p[0], int(p[1]), p[3], p[4]
            if alt in (".", "<NON_REF>"):
                continue
            for a in alt.split(","):
                if a == "<NON_REF>":
                    continue
                if len(ref) == 1 and len(a) == 1:
                    vtype, mut = "SNP", f"{ref}→{a}"
                elif len(ref) > len(a):
                    vtype, mut = "DEL", f"Δ{len(ref) - len(a)} bp"
                elif len(a) > len(ref):
                    ins = a[len(ref):]
                    vtype = "INS"
                    mut = f"+{ins}" if len(ins) <= 5 else f"+{len(ins)} bp"
                else:
                    vtype, mut = "SUB", f"{ref}→{a}"
                variants.append({"seq_id": chrom, "position": pos, "type": vtype, "mutation": mut})
    return variants


def _parse_reference_region(ref_fasta: str, chrom: str, start: int, end: int) -> str:
    out = run(["samtools", "faidx", ref_fasta, f"{chrom}:{start}-{end}"])
    return "".join(l.strip() for l in out.splitlines() if l and not l.startswith(">")).upper()


def _build_aligned(seq: str, qual: str, cigar: str, ref_seq: str, ref_start: int) -> Dict[str, Any]:
    aref, aread, aqual = [], [], []
    mismatches, insertions, deletions = [], [], []
    si = ri = 0
    for length, op in _decode_cigar(cigar):
        if op in {"M", "=", "X"}:
            for _ in range(length):
                rb = seq[si] if si < len(seq) else "N"
                qv = ord(qual[si]) - 33 if si < len(qual) else None
                refb = ref_seq[ri] if ri < len(ref_seq) else "N"
                gpos = ref_start + ri
                aref.append(refb); aread.append(rb); aqual.append(qv)
                if rb.upper() != refb.upper():
                    mismatches.append({"position": gpos, "ref": refb, "alt": rb, "quality": qv})
                si += 1; ri += 1
        elif op == "I":
            inserted = seq[si:si + length]
            iq = [ord(c) - 33 for c in qual[si:si + length]]
            insertions.append({"position": ref_start + ri - 1, "sequence": inserted, "qualities": iq})
            aref.extend(["-"] * length); aread.extend(list(inserted)); aqual.extend(iq)
            si += length
        elif op in {"D", "N"}:
            deleted = ref_seq[ri:ri + length]
            deletions.append({"position": ref_start + ri, "sequence": deleted, "length": length})
            aref.extend(list(deleted)); aread.extend(["-"] * len(deleted)); aqual.extend([None] * len(deleted))
            ri += length
        elif op == "S":
            si += length
        elif op in {"H", "P"}:
            continue
    return {"aligned_ref": "".join(aref), "aligned_read": "".join(aread), "aligned_qual": aqual,
            "mismatches": mismatches, "insertions": insertions, "deletions": deletions}


def extract(bam: str, ref: str, variants: List[Dict[str, Any]], out_dir: str,
            window: int, max_reads: int) -> None:
    os.makedirs(out_dir, exist_ok=True)
    by_chrom: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for v in variants:
        by_chrom[v["seq_id"]].append({**v, "_pos": v["position"]})

    for chrom, muts in by_chrom.items():
        muts.sort(key=lambda m: m["_pos"])
        groups, cur = [], None
        for m in muts:
            start = max(1, m["_pos"] - window)
            end = m["_pos"] + window
            if cur is None or start > cur["end"]:
                cur = {"start": start, "end": end, "mutations": [m]}
                groups.append(cur)
            else:
                cur["end"] = max(cur["end"], end)
                cur["mutations"].append(m)

        for idx, g in enumerate(groups, start=1):
            start, end = g["start"], g["end"]
            pad = 200
            fetch_start = max(1, start - pad)
            ref_full = _parse_reference_region(ref, chrom, fetch_start, end)
            reference = ref_full[start - fetch_start:]
            rows = [l for l in run(["samtools", "view", bam, f"{chrom}:{start}-{end}"]).splitlines() if l.strip()]

            reads = []
            key_of = {id(m): _mutation_key(_strip_ver(chrom), m["_pos"], m.get("type", "VAR")) for m in g["mutations"]}
            alt_support = {key_of[id(m)]: 0 for m in g["mutations"]}
            ref_support = {key_of[id(m)]: 0 for m in g["mutations"]}

            for row in rows:
                p = row.split("\t")
                if len(p) < 11:
                    continue
                flag, pos, mapq, cigar, seq, qual = int(p[1]), int(p[3]), int(p[4]), p[5], p[9], p[10]
                ref_for = ref_full[max(0, pos - fetch_start):]
                al = _build_aligned(seq, qual, cigar, ref_for, pos)
                supported, refkeys = [], []
                mm = {i["position"]: i["alt"] for i in al["mismatches"]}
                ins_pos = {i["position"] for i in al["insertions"]}
                del_pos = set()
                for d in al["deletions"]:
                    del_pos.update(range(d["position"], d["position"] + d["length"]))
                for m in g["mutations"]:
                    key = key_of[id(m)]
                    position, mtype, ms = m["_pos"], m.get("type", ""), m.get("mutation", "")
                    if mtype in ("SNP", "SUB"):
                        alt_base = ms.split("→")[-1] if "→" in ms else None
                        rb = mm.get(position)
                        if alt_base and rb and rb.upper() == alt_base.upper():
                            supported.append(key); alt_support[key] += 1
                        elif pos <= position < pos + len(seq):
                            refkeys.append(key); ref_support[key] += 1
                    elif mtype == "INS":
                        if position in ins_pos or (position - 1) in ins_pos:
                            supported.append(key); alt_support[key] += 1
                    elif mtype == "DEL":
                        if position in del_pos:
                            supported.append(key); alt_support[key] += 1
                reads.append({
                    "name": p[0], "flag": flag, "strand": "-" if flag & 16 else "+",
                    "start": pos, "end": pos + len(al["aligned_ref"].replace("-", "")) - 1,
                    "mapq": mapq, "cigar": cigar, "sequence": seq,
                    "qualities": [ord(c) - 33 for c in qual],
                    "aligned_ref": al["aligned_ref"], "aligned_read": al["aligned_read"],
                    "aligned_qual": al["aligned_qual"], "mismatches": al["mismatches"],
                    "insertions": al["insertions"], "deletions": al["deletions"],
                    "supports": supported, "matches_ref": refkeys,
                })

            mut_positions = {m["_pos"] for m in g["mutations"]}
            reads = [r for r in reads if any(r["start"] <= x <= r["end"] for x in mut_positions)]
            truncated = len(reads) > max_reads
            reads = reads[:max_reads]

            eid = f"{chrom}_{start}_{end}_{idx}.json".replace("/", "_")
            payload = {
                "chrom": chrom, "start": start, "end": end, "reference": reference,
                "window": window, "truncated": truncated, "total_overlapping_reads": len(rows),
                "displayed_reads": len(reads),
                "mutations": [{
                    "mutation_key": key_of[id(m)], "position": m["_pos"], "type": m.get("type"),
                    "mutation": m.get("mutation"), "gene": m.get("gene"), "freq": m.get("freq"),
                    "alt_support": alt_support[key_of[id(m)]], "ref_support": ref_support[key_of[id(m)]],
                } for m in g["mutations"]],
                "reads": reads,
            }
            with open(os.path.join(out_dir, eid), "w") as fh:
                json.dump(payload, fh)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out", default="evidence")
    ap.add_argument("--window", type=int, default=75)
    ap.add_argument("--max-reads", type=int, default=80)
    args = ap.parse_args()

    variants = parse_vcf_variants(args.vcf)
    extract(args.bam, args.ref, variants, args.out, args.window, args.max_reads)
    print(f"[extract_evidence] {len(variants)} variants -> {args.out}/", file=sys.stderr)


if __name__ == "__main__":
    main()
