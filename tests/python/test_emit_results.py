"""Unit tests for the parsers in bin/emit_results.py (the results.json producer)."""

import emit_results as er


# ── placeholder / helper semantics ─────────────────────────────────────────────

def test_present_treats_placeholders_and_empty_as_absent(tmp_path):
    empty = tmp_path / "coverage.tsv"
    empty.write_text("")                       # stageAs-renamed NO_FILE placeholder
    named = tmp_path / "NO_FILE"
    named.write_text("x")
    real = tmp_path / "real.txt"
    real.write_text("data")
    empty_dir = tmp_path / "evidence"
    empty_dir.mkdir()

    assert not er._present("")
    assert not er._present(None)
    assert not er._present(str(tmp_path / "missing"))
    assert not er._present(str(empty))
    assert not er._present(str(named))
    assert not er._present(str(empty_dir))
    assert er._present(str(real))


def test_strip_version_only_numeric_suffix():
    assert er._strip_version("BK006935.2") == "BK006935"
    assert er._strip_version("II.1") == "II"
    assert er._strip_version("I") == "I"
    assert er._strip_version("chr.X") == "chr.X"


# ── GenBank gene table ─────────────────────────────────────────────────────────

def test_parse_gbk_genes_reads_gene_and_cds_per_record(fx):
    genes = er._parse_gbk_genes(str(fx / "tiny.gbk"))
    abc1 = [g for g in genes if g["name"] == "ABC1"]
    assert {g["seq_id"] for g in genes} == {"I", "II"}
    assert len(abc1) == 2                                   # gene + CDS features
    assert (abc1[0]["start"], abc1[0]["end"]) == (10, 40)
    comp = [g for g in genes if g["locus_tag"] == "YAL002W"]
    assert (comp[0]["start"], comp[0]["end"]) == (60, 90)   # complement(60..90)


# ── VCF → mutations ────────────────────────────────────────────────────────────

def _by_pos(muts):
    out = {}
    for m in muts:
        out.setdefault((m["seq_id"], m["position"]), []).append(m)
    return out


def test_parse_vcf_classifies_and_formats_variants(fx):
    muts = er.parse_vcf(str(fx / "tiny.vcf"), str(fx / "tiny.gbk"))
    by = _by_pos(muts)

    assert len(muts) == 7                                   # NON_REF dropped, multi-allelic split
    snp = by[("I", "15")][0]
    assert (snp["type"], snp["mutation"], snp["freq"]) == ("SNP", "A→G", "100.0%")
    assert by[("I", "20")][0]["freq"] == "75.0%"            # AD 10,30
    assert (by[("I", "25")][0]["type"], by[("I", "25")][0]["mutation"]) == ("DEL", "Δ2 bp")
    ins = by[("I", "70")][0]
    assert (ins["type"], ins["mutation"]) == ("INS", "+6 bp")   # > 5 inserted bases → length
    assert ins["description"] == "QUAL=150 LowQual"
    multi = by[("I", "100")]
    assert [m["mutation"] for m in multi] == ["G→A", "G→C"]
    assert [m["freq"] for m in multi] == ["60.0%", "30.0%"]     # per-allele AF
    sub = by[("II", "10")][0]                               # "II.1" → version stripped
    assert (sub["type"], sub["mutation"], sub["freq"]) == ("SUB", "GA→TC", "100%")


def test_parse_vcf_annotates_gene_and_product(fx):
    """Regression: gene + CDS over the same span — name comes from either, product from the
    CDS. Before the fix every GATK variant had an empty `annotation` (product)."""
    by = _by_pos(er.parse_vcf(str(fx / "tiny.vcf"), str(fx / "tiny.gbk")))
    abc1 = by[("I", "15")][0]
    assert (abc1["gene"], abc1["annotation"]) == ("ABC1", "Abc1p test protein")
    lt = by[("I", "70")][0]                                  # no /gene → locus_tag
    assert (lt["gene"], lt["annotation"]) == ("YAL002W", "hypothetical protein")
    xyz = by[("II", "10")][0]                                # gene only, no CDS/product
    assert (xyz["gene"], xyz["annotation"]) == ("XYZ1", "")
    inter = by[("I", "100")][0]                              # intergenic
    assert (inter["gene"], inter["annotation"]) == ("", "")


def test_parse_vcf_gene_lookup_respects_chromosome(fx, tmp_path):
    """A position inside ABC1's span but on chromosome II must not pick up ABC1."""
    vcf = tmp_path / "x.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                   "II\t35\t.\tA\tG\t50\tPASS\t.\n")
    m = er.parse_vcf(str(vcf), str(fx / "tiny.gbk"))[0]
    assert m["gene"] == ""


def test_parse_vcf_missing_file_returns_empty():
    assert er.parse_vcf("") == []


# ── QC / contamination parsers ─────────────────────────────────────────────────

def test_parse_flagstat_ignores_primary_mapped_line(fx):
    assert er.parse_flagstat(str(fx / "flagstat.txt")) == {
        "total_reads": 1000, "mapped_percentage": 95.0,
        "unmapped_percentage": 5.0, "unmapped_reads": 50}


def test_parse_kraken_keeps_significant_species_genus_unclassified(fx):
    res = er.parse_kraken_report(str(fx / "kraken.report"))
    assert [(r["name"], r["rank"], r["percentage"]) for r in res] == [
        ("unclassified", "U", 5.0), ("Escherichia coli", "S", 60.0), ("Escherichia", "G", 30.0)]


def test_parse_mash_filters_and_sorts(fx):
    res = er.parse_mash_screen(str(fx / "mash_screen.tab"))
    assert [r["query_id"] for r in res] == ["GCF_000005845.2.fna.gz", "GCF_salmonella.fna.gz"]
    assert res[0]["identity"] == 0.99


def test_parse_blast_aggregates_once_per_title(fx):
    """Regression: the original parser appended every hit twice."""
    res = er.parse_blast_results(str(fx / "blast_results.tsv"))
    assert len(res) == 2
    ecoli = res[0]
    assert (ecoli["title"], ecoli["count"], ecoli["identity"], ecoli["evalue"]) == (
        "Escherichia coli strain A", 2, 98.0, 1e-50)
    assert set(ecoli["queries"].split(", ")) == {"q1", "q2"}
    assert all(set(r) == {"title", "count", "identity", "evalue", "queries"} for r in res)


# ── high-confidence SNP table + evidence links ─────────────────────────────────

def test_high_confidence_snps_depth_and_freq_filter(fx, tmp_path):
    muts = er.parse_vcf(str(fx / "tiny.vcf"), str(fx / "tiny.gbk"))
    out = tmp_path / "hc.csv"
    assert er.write_high_confidence_snps(muts, str(fx / "coverage.tsv"), str(out))
    rows = out.read_text().strip().splitlines()
    assert rows[0].startswith("Position,Reference,Alternative")
    assert len(rows) == 2 and rows[1].startswith("15,A,G,100.0,30,ABC1")   # only depth>20 & freq>40


def test_high_confidence_snps_removes_empty_csv(fx, tmp_path):
    out = tmp_path / "hc.csv"
    muts = [{"type": "SNP", "position": "20", "freq": "75%", "mutation": "C→T"}]   # depth 10
    assert not er.write_high_confidence_snps(muts, str(fx / "coverage.tsv"), str(out))
    assert not out.exists()


def test_attach_evidence_links_by_mutation_key(tmp_path):
    ev = tmp_path / "evidence"
    ev.mkdir()
    (ev / "I_15.json").write_text('{"mutations":[{"mutation_key":"I_15_SNP"}]}')
    muts = [{"seq_id": "I", "position": "15", "type": "SNP"},
            {"seq_id": "I", "position": "20", "type": "SNP"}]
    er.attach_evidence_links_from_dir(muts, str(ev))
    assert muts[0]["evidence_links"] == [{"text": "reads", "url": "evidence/I_15.json"}]
    assert "evidence_links" not in muts[1]
