"""Contract tests: results.json produced exactly as EMIT_RESULTS runs emit_results.py must
match the schema in docs/API_CONTRACT.md — the interface the React frontend consumes.

These guard the one thing that must not drift: a change to a parser that renames or drops a
key breaks the UI silently, while the pipeline still exits 0.
"""

import json
import shutil

import pytest

import emit_results as er

MUTATION_KEYS = {"type", "seq_id", "position", "mutation", "freq", "annotation", "gene",
                 "description", "details", "evidence_links"}
READ_STATS_KEYS = {"total_reads", "mapped_percentage", "unmapped_percentage", "unmapped_reads"}
MASH_KEYS = {"identity", "shared_hashes", "median_multiplicity", "p_value", "query_id", "comment"}
CONTAMINANT_KEYS = {"name", "rank", "percentage", "reads"}
TOP_LEVEL_ALLOWED = {
    "job_id", "status", "variant_caller", "polymorphism_mode", "reference_file", "steps",
    "mutations", "read_stats", "coverage_file", "hc_snps_file", "unmapped_files",
    "mash_results", "mash_unmapped_results", "contaminants", "blast_results",
    "snp_blast_results", "snp_query_file",
    # optional superset added with the CNVkit step (gatk path + --parent_sample)
    "chromosome_copy_number", "aneuploidy", "cnv_tool",
}


def _no_file(tmp_path, name):
    """The Nextflow placeholder idiom: assets/NO_FILE staged under a slot's conventional name."""
    p = tmp_path / name
    p.write_text("")
    return p


def _assert_contract(res):
    assert res["status"] == "completed"
    assert set(res) <= TOP_LEVEL_ALLOWED, f"undocumented keys: {set(res) - TOP_LEVEL_ALLOWED}"
    assert "manta_vcf" not in res                       # removed in 30ce748; must never return
    for m in res.get("mutations", []):
        assert MUTATION_KEYS <= set(m), f"mutation missing {MUTATION_KEYS - set(m)}"
        assert isinstance(m["details"], list)
    if "read_stats" in res:
        assert set(res["read_stats"]) == READ_STATS_KEYS
    for r in res.get("mash_results", []) + res.get("mash_unmapped_results", []):
        assert set(r) == MASH_KEYS
    for r in res.get("contaminants", []):
        assert CONTAMINANT_KEYS <= set(r)
    for cn in res.get("chromosome_copy_number", []):
        assert {"seq_id", "mean_cn", "call", "aneuploidy"} <= set(cn)


@pytest.fixture
def gatk_inputs(fx, tmp_path):
    ev = tmp_path / "evidence"
    ev.mkdir()
    (ev / "I_15.json").write_text('{"mutations":[{"mutation_key":"I_15_SNP"}]}')
    cnv = tmp_path / "cnv.json"
    cnv.write_text(json.dumps({
        "cnv_tool": "cnvkit",
        "chromosome_copy_number": [{"seq_id": "I", "mean_cn": 1.3, "call": "gain", "aneuploidy": True}],
        "aneuploidy": [{"seq_id": "I", "call": "gain", "mean_cn": 1.3}]}))
    return ev, cnv


def test_gatk_full_contract(fx, tmp_path, run_bin, load_json, gatk_inputs):
    ev, cnv = gatk_inputs
    for name in ("coverage.tsv", "flagstat.txt", "kraken.report", "mash_screen.tab"):
        shutil.copy(fx / name, tmp_path / name)
    run_bin("emit_results.py", "--job-id", "s1", "--caller", "gatk", "--polymorphism", "false",
            "--reference", fx / "tiny.gbk", "--reference-gbk", fx / "tiny.gbk",
            "--vcf", fx / "tiny.vcf", "--coverage", "coverage.tsv", "--flagstat", "flagstat.txt",
            "--kraken", "kraken.report", "--mash", "mash_screen.tab",
            "--mash-unmapped", _no_file(tmp_path, "mash_unmapped.tab"),
            "--blast", _no_file(tmp_path, "blast_results.tsv"),
            "--evidence-dir", ev, "--cnv", cnv, "--out", "results.json")
    res = load_json(tmp_path / "results.json")

    _assert_contract(res)
    assert res["variant_caller"] == "gatk" and res["polymorphism_mode"] is False
    assert res["reference_file"] == "tiny.gbk"
    assert len(res["mutations"]) == 7
    assert res["read_stats"]["total_reads"] == 1000
    assert res["coverage_file"] == "coverage.tsv"
    assert res["hc_snps_file"] == "high_confidence_snps.csv"
    assert (tmp_path / "high_confidence_snps.csv").exists()
    assert res["steps"]["haplotype_caller"] == "completed"
    assert res["steps"]["kraken"] == "completed" and res["steps"]["cnv"] == "completed"
    assert res["aneuploidy"][0]["seq_id"] == "I"
    linked = [m for m in res["mutations"] if m["evidence_links"]]
    assert [(m["seq_id"], m["position"]) for m in linked] == [("I", "15")]


def test_gatk_minimal_run_skips_optional_steps(fx, tmp_path, run_bin, load_json):
    """Every optional slot given the NO_FILE placeholder — what a run with no DBs, no
    evidence and no --parent_sample sends."""
    ph = {n: _no_file(tmp_path, n) for n in ("coverage.tsv", "flagstat.txt", "kraken.report",
          "mash_screen.tab", "mash_unmapped.tab", "blast_results.tsv", "cnv.json")}
    (tmp_path / "evidence").mkdir()
    run_bin("emit_results.py", "--job-id", "s1", "--caller", "gatk",
            "--reference", fx / "tiny.gbk", "--vcf", fx / "tiny.vcf",
            "--coverage", ph["coverage.tsv"], "--flagstat", ph["flagstat.txt"],
            "--kraken", ph["kraken.report"], "--mash", ph["mash_screen.tab"],
            "--mash-unmapped", ph["mash_unmapped.tab"], "--blast", ph["blast_results.tsv"],
            "--evidence-dir", tmp_path / "evidence", "--cnv", ph["cnv.json"],
            "--out", "results.json")
    res = load_json(tmp_path / "results.json")

    _assert_contract(res)
    for absent in ("coverage_file", "read_stats", "contaminants", "mash_results",
                   "chromosome_copy_number", "aneuploidy"):
        assert absent not in res
    assert res["steps"]["cnv"] == "skipped" and res["steps"]["mash"] == "skipped"
    assert res["unmapped_files"] == []


def test_breseq_contract(fx, tmp_path, run_bin, load_json):
    out = tmp_path / "breseq" / "output"
    out.mkdir(parents=True)
    for f in ("index.html", "summary.html"):
        shutil.copy(fx / "breseq_ecoli" / f, out / f)
    run_bin("emit_results.py", "--job-id", "b1", "--caller", "breseq",
            "--reference", "e_coli.gbk", "--breseq-dir", tmp_path / "breseq",
            "--out", "results.json")
    res = load_json(tmp_path / "results.json")

    _assert_contract(res)
    assert "polymorphism_mode" not in res               # gatk-only key
    assert res["steps"]["breseq"] == "completed"
    assert len(res["mutations"]) == 18


def test_breseq_parser_matches_golden_output(fx):
    """Golden test on a real breseq report (E. coli K-12, 18 mutations). Pins the current
    output exactly — including breseq's own formatting (thousands commas in `position`,
    non-breaking spaces, empty seq_id on single-contig references)."""
    expected = json.loads((fx / "breseq_ecoli" / "expected_results.json").read_text())
    assert er.parse_breseq_html(str(fx / "breseq_ecoli" / "index.html")) == expected["mutations"]
    assert er.parse_breseq_summary(str(fx / "breseq_ecoli" / "summary.html")) == expected["read_stats"]
