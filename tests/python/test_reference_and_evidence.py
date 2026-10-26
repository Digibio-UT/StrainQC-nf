"""bin/gbk_to_fasta.py, and the agreement between the two scripts that share the evidence key.

extract_evidence.py (VARIANT_EVIDENCE) WRITES per-variant evidence files keyed by
`<seq_id>_<position>_<type>`; emit_results.py (EMIT_RESULTS) READS them back by rebuilding the
same key from its own VCF parse. The two keep separate copies of the variant-classification
logic, so if either drifts, evidence links silently vanish from the UI.
"""

import emit_results as er
import extract_evidence as ee
from gbk_to_fasta import gbk_to_fasta


def test_gbk_to_fasta_writes_every_record(fx, tmp_path):
    out = tmp_path / "ref.fasta"
    gbk_to_fasta(str(fx / "tiny.gbk"), str(out))
    lines = out.read_text().splitlines()
    headers = [l for l in lines if l.startswith(">")]
    assert headers == [">I", ">II"]
    seq_i = "".join(lines[1:lines.index(">II")])
    assert seq_i == "ACGT" * 30                                  # uppercased, digits/spaces stripped
    assert max(len(l) for l in lines if not l.startswith(">")) <= 70


def test_both_vcf_parsers_classify_identically(fx):
    writer = ee.parse_vcf_variants(str(fx / "tiny.vcf"))
    reader = er.parse_vcf(str(fx / "tiny.vcf"))
    assert [(v["type"], v["mutation"]) for v in writer] == [(m["type"], m["mutation"]) for m in reader]


def test_evidence_keys_match_between_writer_and_reader(fx):
    writer_keys = [ee._mutation_key(ee._strip_ver(v["seq_id"]), v["position"], v["type"])
                   for v in ee.parse_vcf_variants(str(fx / "tiny.vcf"))]
    reader_keys = [f'{m["seq_id"]}_{m["position"]}_{m["type"]}'
                   for m in er.parse_vcf(str(fx / "tiny.vcf"))]
    assert writer_keys == reader_keys
    assert "II_10_SUB" in reader_keys                            # versioned "II.1" stripped by both


def test_cigar_alignment_reconstruction():
    # read: 2 match, 1 inserted base, 2 match, 1 deleted ref base, 2 match, 1 soft-clip
    aln = ee._build_aligned(seq="ACTGTCAG", qual="IIIIIIII", cigar="2M1I2M1D2M1S",
                            ref_seq="ACGTTCA", ref_start=100)
    assert aln["aligned_ref"] == "AC-GTTCA"
    assert aln["aligned_read"] == "ACTGT-CA"
    assert aln["insertions"] == [{"position": 101, "sequence": "T", "qualities": [40]}]
    assert aln["deletions"] == [{"position": 104, "sequence": "T", "length": 1}]
    assert aln["mismatches"] == []
