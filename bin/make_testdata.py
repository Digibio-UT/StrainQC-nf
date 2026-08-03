#!/usr/bin/env python3
"""
make_testdata.py — generate a tiny, deterministic StrainQC test dataset.

Produces a ~2 kb reference and single-end reads that tile it, with ONE homozygous SNP
planted at a known position. The GATK/breseq path should therefore call exactly one SNP,
giving a cheap, self-checking smoke test that needs no external data.

Outputs (under --outdir, default assets/test/):
  reference.fasta        2 kb single-contig reference
  reads.fastq.gz         single-end reads, all reads over the SNP carry the ALT base
  samplesheet.csv        one sample pointing at reads.fastq.gz (relative name)
  EXPECTED.md            the planted variant, for asserting the run's results.json

Deterministic: fixed RNG seed, so regenerating gives identical files.
"""
import argparse
import gzip
import os
import random

REF_LEN = 2000
READ_LEN = 100
N_READS = 400          # ~20x coverage
SNP_POS = 1000         # 1-based position of the planted SNP
SEED = 42


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="assets/test")
    args = ap.parse_args()
    os.makedirs(os.path.join(args.outdir, ), exist_ok=True)

    rng = random.Random(SEED)
    bases = "ACGT"
    ref = [rng.choice(bases) for _ in range(REF_LEN)]

    ref_base = ref[SNP_POS - 1]
    alt_base = rng.choice([b for b in bases if b != ref_base])

    # reference.fasta
    ref_str = "".join(ref)
    with open(os.path.join(args.outdir, "reference.fasta"), "w") as fh:
        fh.write(">CHR1 synthetic test contig\n")
        for i in range(0, REF_LEN, 70):
            fh.write(ref_str[i:i + 70] + "\n")

    # reads: tile the genome; any read covering SNP_POS carries the ALT base (homozygous)
    with gzip.open(os.path.join(args.outdir, "reads.fastq.gz"), "wt") as fh:
        for n in range(N_READS):
            start = rng.randint(0, REF_LEN - READ_LEN)          # 0-based
            read = list(ref_str[start:start + READ_LEN])
            snp_idx = (SNP_POS - 1) - start
            if 0 <= snp_idx < READ_LEN:
                read[snp_idx] = alt_base
            seq = "".join(read)
            fh.write(f"@read_{n}_pos{start + 1}\n{seq}\n+\n{'I' * READ_LEN}\n")

    # samplesheet (bare filename; main.nf resolves relative to the samplesheet's dir)
    with open(os.path.join(args.outdir, "samplesheet.csv"), "w") as fh:
        fh.write("sample,fastq_1,fastq_2\n")
        fh.write("testsample,reads.fastq.gz,\n")

    with open(os.path.join(args.outdir, "EXPECTED.md"), "w") as fh:
        fh.write("# Expected result for the synthetic test dataset\n\n")
        fh.write(f"- Reference: CHR1, {REF_LEN} bp\n")
        fh.write(f"- Exactly **one SNP** at position **{SNP_POS}**: "
                 f"`{ref_base}→{alt_base}` (homozygous, ~{N_READS * READ_LEN // REF_LEN}x)\n")
        fh.write("- results.json should contain 1 mutation of type SNP at that position.\n")

    print(f"Wrote test data to {args.outdir}/  (SNP CHR1:{SNP_POS} {ref_base}->{alt_base})")


if __name__ == "__main__":
    main()
