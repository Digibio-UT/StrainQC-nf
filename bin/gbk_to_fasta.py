#!/usr/bin/env python3
"""
gbk_to_fasta.py — convert a GenBank file to FASTA without BioPython.
Ported from backend/app/pipeline/wrappers.py :: _gbk_to_fasta.

Usage: gbk_to_fasta.py <input.gbk> <output.fasta>
"""
import re
import sys


def gbk_to_fasta(gbk_path: str, fasta_path: str) -> None:
    written = 0
    with open(gbk_path) as fin, open(fasta_path, "w") as fout:
        seq_id = None
        in_origin = False
        parts = []

        def flush():
            nonlocal written
            if seq_id and parts:
                fout.write(f">{seq_id}\n")
                seq = "".join(parts)
                for i in range(0, len(seq), 70):
                    fout.write(seq[i:i + 70] + "\n")
                written += 1

        for line in fin:
            if line.startswith("LOCUS"):
                flush()
                seq_id = line.split()[1]
                parts = []
                in_origin = False
            elif line.startswith("ORIGIN"):
                in_origin = True
            elif line.startswith("//"):
                in_origin = False
            elif in_origin:
                parts.append(re.sub(r"[\s0-9]", "", line).upper())
        flush()

    if written == 0:
        sys.exit(f"gbk_to_fasta: no records found in {gbk_path}")
    print(f"gbk_to_fasta: wrote {written} record(s) -> {fasta_path}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: gbk_to_fasta.py <input.gbk> <output.fasta>")
    gbk_to_fasta(sys.argv[1], sys.argv[2])
