"""Tests for bin/summarize_cnv.py — CNVkit segments → per-chromosome CN + aneuploidy calls."""

import json

import pytest

import summarize_cnv as sc

FAI = [("I", 200_000), ("II", 800_000), ("III", 300_000), ("IV", 1_500_000),
       ("V", 500_000), ("VII", 1_000_000), ("Mito", 85_000)]


@pytest.fixture
def make_cns(tmp_path):
    """Write a .cns + .fai. `segs` maps chrom -> list of (fraction_of_chrom, log2);
    chromosomes not listed get one neutral segment (log2 = 0)."""
    def _make(segs):
        fai = tmp_path / "ref.fa.fai"
        fai.write_text("".join(f"{c}\t{n}\t0\t60\t61\n" for c, n in FAI))
        lines = ["chromosome\tstart\tend\tgene\tlog2\tdepth\tprobes\tweight"]
        for chrom, length in FAI:
            start = 0
            for frac, log2 in segs.get(chrom, [(1.0, 0.0)]):
                end = start + int(length * frac)
                lines.append(f"{chrom}\t{start}\t{end}\t-\t{log2}\t100\t10\t1")
                start = end
        cns = tmp_path / "s.cns"
        cns.write_text("\n".join(lines) + "\n")
        return str(cns), str(fai)
    return _make


def _calls(rows):
    return {r["seq_id"]: r["call"] for r in rows}


def test_whole_chromosome_duplication_is_a_gain(make_cns):
    rows, aneu = sc.summarise(*make_cns({"I": [(1.0, 1.0)]}), baseline=1)
    assert aneu == [{"seq_id": "I", "call": "gain", "mean_cn": 2.0}]
    assert _calls(rows)["II"] == "neutral"


def test_mosaic_gain_is_called(make_cns):
    """Regression for Ottilie CBR110: a mosaic chr I gain sits at ~1.3x, not 2x. The old
    absolute 1.4x-of-ploidy cutoff missed it; the genome-median-relative test must not."""
    rows, aneu = sc.summarise(*make_cns({"I": [(1.0, 0.378)]}), baseline=1)   # 2**0.378 ≈ 1.30
    assert [a["seq_id"] for a in aneu] == ["I"]


def test_small_elevation_is_not_aneuploidy(make_cns):
    rows, aneu = sc.summarise(*make_cns({"I": [(1.0, 0.2)]}), baseline=1)    # ≈ 1.15x
    assert aneu == []
    assert _calls(rows)["I"] == "neutral"


def test_focal_amplification_is_reported_not_flagged(make_cns):
    rows, aneu = sc.summarise(*make_cns({"VII": [(0.3, 1.0), (0.7, 0.0)]}), baseline=1)
    assert aneu == []
    assert _calls(rows)["VII"] == "focal"


def test_whole_chromosome_loss(make_cns):
    rows, aneu = sc.summarise(*make_cns({"V": [(1.0, -1.0)]}), baseline=2)     # diploid baseline
    assert aneu == [{"seq_id": "V", "call": "loss", "mean_cn": 1.0}]


def test_mitochondrial_dosage_never_flagged(make_cns):
    rows, aneu = sc.summarise(*make_cns({"Mito": [(1.0, 3.0)]}), baseline=1)
    assert aneu == []
    assert _calls(rows)["Mito"] == "neutral"


def test_output_in_reference_order_with_median_ratio(make_cns):
    rows, _ = sc.summarise(*make_cns({"I": [(1.0, 1.0)]}), baseline=1)
    assert [r["seq_id"] for r in rows] == [c for c, _ in FAI]
    assert rows[0]["cn_vs_genome_median"] == 2.0


# ── real CNVkit output from the Ottilie pilot (golden) ────────────────────────

def test_real_cbr110_calls_only_chr1_gain(fx):
    rows, aneu = sc.summarise(str(fx / "cnv" / "CBR110-15-R3a.cns"),
                              str(fx / "cnv" / "R64-1-1.fa.fai"), baseline=1)
    assert [(a["seq_id"], a["call"]) for a in aneu] == [("I", "gain")]
    assert 1.25 < next(r for r in rows if r["seq_id"] == "I")["mean_cn"] < 1.35


def test_real_control_clone_has_no_aneuploidy(fx):
    _, aneu = sc.summarise(str(fx / "cnv" / "Carmaphycin-R9-2.cns"),
                           str(fx / "cnv" / "R64-1-1.fa.fai"), baseline=1)
    assert aneu == []


def test_cli_writes_cnv_json(fx, run_bin, load_json, tmp_path):
    run_bin("summarize_cnv.py", "--cns", fx / "cnv" / "CBR110-15-R3a.cns",
            "--fai", fx / "cnv" / "R64-1-1.fa.fai", "--sample", "CBR110",
            "--baseline-ploidy", "1", "--out", "cnv.json")
    out = load_json(tmp_path / "cnv.json")
    assert out["cnv_tool"] == "cnvkit" and out["sample"] == "CBR110"
    assert len(out["chromosome_copy_number"]) == 16          # 16 nuclear chroms in the .cns
    assert out["aneuploidy"][0]["seq_id"] == "I"
