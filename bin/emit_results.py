#!/usr/bin/env python3
"""
emit_results.py — assemble the StrainQC results.json contract for StrainQC-nf.

This is the single point that produces results.json on any backend. It is a port of the
working app's assembly logic:
  - parsing:      backend/app/pipeline/parsers.py
  - orchestration: backend/app/pipeline/manager.py :: PipelineManager.run()

In the working app, manager.py ran the tools *and* assembled results.json. Under Nextflow
the tools run as separate processes and drop their output files; this script takes over the
assembly role — read whatever files exist, build the identical results.json.

Every input is OPTIONAL: a missing path, or the Nextflow placeholder "NO_FILE", is treated
as "that step didn't run", exactly like manager.py's os.path.exists() guards. That lets the
workflow pass placeholders for steps that were skipped (no Kraken DB, breseq path, etc.).

Contract reference: docs/API_CONTRACT.md. Note: `manta_vcf` was removed upstream — not emitted.

Container note: needs `beautifulsoup4` (breseq HTML parsing) + `biopython` is NOT required
(GBK parsing is done with a self-contained parser). See report process container in nextflow.config.
"""
import argparse
import csv
import json
import os
import re
import sys
from html.parser import HTMLParser
from typing import Any, Dict, List


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _present(path: str) -> bool:
    """A path counts as provided only if it exists, isn't a NO_* placeholder, and is
    non-empty. Empty-file detection lets Nextflow's single-placeholder idiom work even
    when the placeholder is `stageAs`-renamed to a slot's conventional name (it stays
    zero-length); real tool outputs are always non-empty."""
    if not path or os.path.basename(path).startswith("NO_") or not os.path.exists(path):
        return False
    if os.path.isdir(path):
        return any(os.scandir(path))
    return os.path.getsize(path) > 0


def _strip_version(acc: str) -> str:
    """'BK006935.2' -> 'BK006935' (only when the suffix is numeric)."""
    dot = acc.rfind(".")
    if dot > 0 and acc[dot + 1:].isdigit():
        return acc[:dot]
    return acc


# ─────────────────────────────────────────────────────────────────────────────
# Parsers (ported verbatim in behaviour from parsers.py, made self-contained)
# ─────────────────────────────────────────────────────────────────────────────

class _BreseqTableParser(HTMLParser):
    """Stdlib HTML parser (no BeautifulSoup) that captures every table's rows/cells,
    flags the table containing the breseq `mutation_header_row`, and records the index of
    the table following each <h2> (for the summary read-stats lookup). Handles nesting via
    an open-table stack."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.tables: List[Dict[str, Any]] = []   # document order; {'has_mut_header', 'rows'}
        self.h2_next_table: Dict[str, int] = {}   # h2 text -> index of next table opened
        self._stack: List[Dict[str, Any]] = []    # open tables (innermost last)
        self._row = None
        self._cell = None
        self._in_cell = False
        self._link = None
        self._in_h2 = False
        self._h2_parts: List[str] = []

    @staticmethod
    def _cls(attrs) -> List[str]:
        for k, v in attrs:
            if k == "class" and v:
                return v.split()
        return []

    @staticmethod
    def _attr(attrs, name):
        for k, v in attrs:
            if k == name:
                return v
        return None

    def handle_starttag(self, tag, attrs):
        if tag == "table":
            t = {"has_mut_header": False, "rows": []}
            self._stack.append(t)
            self.tables.append(t)
        elif tag == "tr" and self._stack:
            self._row = {"classes": self._cls(attrs), "cells": []}
            self._stack[-1]["rows"].append(self._row)
        elif tag in ("td", "th") and self._row is not None:
            self._cell = {"text": [], "links": []}
            self._in_cell = True
            if tag == "th" and "mutation_header_row" in self._cls(attrs) and self._stack:
                self._stack[-1]["has_mut_header"] = True
        elif tag == "a" and self._in_cell:
            self._link = {"href": self._attr(attrs, "href"), "text": []}
        elif tag == "h2":
            self._in_h2 = True
            self._h2_parts = []

    def handle_endtag(self, tag):
        if tag == "a" and self._link is not None:
            if self._link["href"]:
                self._cell["links"].append(("".join(self._link["text"]).strip(), self._link["href"]))
            self._link = None
        elif tag in ("td", "th") and self._in_cell:
            self._cell["text"] = "".join(self._cell["text"]).strip()
            self._row["cells"].append(self._cell)
            self._cell = None
            self._in_cell = False
        elif tag == "tr":
            self._row = None
        elif tag == "table" and self._stack:
            self._stack.pop()
        elif tag == "h2":
            self.h2_next_table["".join(self._h2_parts).strip()] = len(self.tables)
            self._in_h2 = False

    def handle_data(self, data):
        if self._in_cell:
            self._cell["text"].append(data)
        if self._link is not None:
            self._link["text"].append(data)
        if self._in_h2:
            self._h2_parts.append(data)


def parse_breseq_html(html_file_path: str) -> List[Dict[str, Any]]:
    """Parse breseq index.html 'Predicted mutations' table (richer than the GD file)."""
    mutations: List[Dict[str, Any]] = []
    if not _present(html_file_path):
        return []
    try:
        parser = _BreseqTableParser()
        with open(html_file_path, "r") as f:
            parser.feed(f.read())
        table = next((t for t in parser.tables if t["has_mut_header"]), None)
        if not table:
            return []

        valid = {"normal_table_row", "polymorphism_table_row",
                 "alternate_table_row_0", "alternate_table_row_1"}
        for row in table["rows"]:
            if not any(c in valid for c in row["classes"]):
                continue
            cells = row["cells"]
            if len(cells) < 6:
                continue
            txt = [c["text"] for c in cells]

            evidence = txt[0]
            evidence_links = [{"text": t, "url": h} for (t, h) in cells[0]["links"]]

            col1_stripped = txt[1].replace(",", "")
            if ":" in col1_stripped:
                col1_stripped = col1_stripped.split(":")[0]
            is_multi_ref = not col1_stripped.isdigit()

            seq_id = position = mutation = annotation = gene = description = ""
            freq = "100%"
            if is_multi_ref:
                seq_id, position, mutation = txt[1], txt[2], txt[3]
                if len(cells) >= 8:
                    freq, annotation, gene, description = txt[4], txt[5], txt[6], txt[7]
                elif len(cells) == 7:
                    annotation, gene, description = txt[4], txt[5], txt[6]
            else:
                position, mutation = txt[1], txt[2]
                if len(cells) >= 7:
                    freq, annotation, gene, description = txt[3], txt[4], txt[5], txt[6]
                elif len(cells) == 6:
                    annotation, gene, description = txt[3], txt[4], txt[5]

            mut_type = "SNP"
            if "→" in mutation and len(mutation) < 10:
                mut_type = "SNP"
            elif "ins" in mutation.lower() or "+" in mutation:
                mut_type = "INS"
            elif "del" in mutation.lower() or "-" in mutation or "Δ" in mutation:
                mut_type = "DEL"
            elif "IS" in mutation:
                mut_type = "MOB"

            mutations.append({
                "type": mut_type, "evidence": evidence, "evidence_links": evidence_links,
                "position": position, "mutation": mutation, "freq": freq,
                "annotation": annotation, "gene": gene, "description": description,
                "seq_id": seq_id, "details": [freq, annotation],
            })
    except Exception as e:  # noqa: BLE001 — match original defensive behaviour
        print(f"Error parsing HTML: {e}", file=sys.stderr)
        return []
    return mutations


def parse_breseq_summary(summary_file_path: str) -> Dict[str, Any]:
    """Extract read stats (total_reads, mapped/unmapped %) from breseq summary.html."""
    stats: Dict[str, Any] = {}
    if not _present(summary_file_path):
        return stats
    try:
        parser = _BreseqTableParser()
        with open(summary_file_path, "r") as f:
            parser.feed(f.read())
        idx = parser.h2_next_table.get("Read File Information")
        if idx is None or idx >= len(parser.tables):
            return stats
        for row in parser.tables[idx]["rows"]:
            cols = [c["text"] for c in row["cells"]]
            if len(cols) < 3:
                continue
            if "total" in cols[1].lower():
                try:
                    stats["total_reads"] = int(cols[2].replace(",", ""))
                    stats["mapped_percentage"] = float(cols[7].replace("%", ""))
                    stats["unmapped_percentage"] = round(100.0 - stats["mapped_percentage"], 1)
                    stats["unmapped_reads"] = int(stats["total_reads"] * (stats["unmapped_percentage"] / 100.0))
                except Exception:
                    pass
                break
    except Exception as e:  # noqa: BLE001
        print(f"Error parsing summary HTML: {e}", file=sys.stderr)
    return stats


def parse_breseq_gd(gd_file_path: str) -> List[Dict[str, Any]]:
    """Fallback: parse breseq Genome Diff file for mutation calls."""
    mutations: List[Dict[str, Any]] = []
    if not _present(gd_file_path):
        return []
    with open(gd_file_path, "r") as f:
        for row in csv.reader(f, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            record_type = row[0]
            if record_type in ["SNP", "INS", "DEL", "SUB", "MOB"]:
                mutations.append({
                    "type": record_type,
                    "seq_id": row[3] if len(row) > 3 else "",
                    "position": row[4] if len(row) > 4 else "",
                    "details": row[5:] if len(row) > 5 else [],
                })
    return mutations


def parse_kraken_report(report_file_path: str) -> List[Dict[str, Any]]:
    """Parse Kraken2 report -> significant taxa (>0.1%, ranks S/G/U)."""
    results: List[Dict[str, Any]] = []
    if not _present(report_file_path):
        return []
    with open(report_file_path, "r") as f:
        for row in csv.reader(f, delimiter="\t"):
            if not row:
                continue
            try:
                percentage = float(row[0].strip())
                num_covered = int(row[1].strip())
                rank = row[3].strip()
                tax_id = row[4].strip()
                name = row[5].strip()
                if percentage > 0.1 and rank in ["S", "G", "U"]:
                    results.append({
                        "name": name, "tax_id": tax_id, "rank": rank,
                        "percentage": percentage, "reads": num_covered,
                    })
            except (ValueError, IndexError):
                continue
    return results


def parse_blast_results(result_file_path: str) -> List[Dict[str, Any]]:
    """Parse BLAST outfmt-6 (+stitle), aggregate hits by organism title.

    NOTE: fixes a bug in the original parsers.py that appended each hit TWICE
    (once with `queries`, once without) — producing malformed contract records.
    Here each title yields exactly one BlastHit {title,count,identity,evalue,queries}.
    """
    if not _present(result_file_path):
        return []
    aggregated: Dict[str, Dict[str, Any]] = {}
    with open(result_file_path, "r") as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) < 7:
                continue
            query_id, title = row[0], row[6]
            identity, evalue = float(row[2]), float(row[4])
            if evalue > 0.01:
                continue
            agg = aggregated.setdefault(title, {
                "title": title, "count": 0, "total_identity": 0.0,
                "min_evalue": evalue, "queries": set(),
            })
            agg["count"] += 1
            agg["total_identity"] += identity
            agg["min_evalue"] = min(agg["min_evalue"], evalue)
            agg["queries"].add(query_id)

    results: List[Dict[str, Any]] = []
    for data in aggregated.values():
        q = list(data["queries"])
        q_display = ", ".join(q[:5]) + (f" (+{len(q) - 5} more)" if len(q) > 5 else "")
        results.append({
            "title": data["title"],
            "count": data["count"],
            "identity": round(data["total_identity"] / data["count"], 2),
            "evalue": data["min_evalue"],
            "queries": q_display,
        })
    results.sort(key=lambda x: x["count"], reverse=True)
    return results


def parse_mash_screen(file_path: str) -> List[Dict[str, Any]]:
    """Parse mash screen output; keep identity>0.8 & shared-hash>10, top 10 by identity."""
    results: List[Dict[str, Any]] = []
    if not _present(file_path):
        return []
    try:
        with open(file_path, "r") as f:
            for row in csv.reader(f, delimiter="\t"):
                if len(row) < 5:
                    continue
                identity = float(row[0])
                shared_hashes = row[1]
                hash_count = 0
                try:
                    if "/" in shared_hashes:
                        hash_count = int(shared_hashes.split("/")[0])
                except Exception:
                    pass
                if identity > 0.8 and hash_count > 10:
                    results.append({
                        "identity": identity,
                        "shared_hashes": shared_hashes,
                        "median_multiplicity": float(row[2]),
                        "p_value": float(row[3]),
                        "query_id": row[4],
                        "comment": row[5] if len(row) > 5 else "",
                    })
        results.sort(key=lambda x: x["identity"], reverse=True)
        results = results[:10]
    except Exception as e:  # noqa: BLE001
        print(f"Error parsing mash screen output: {e}", file=sys.stderr)
        return []
    return results


def _parse_gbk_genes(gbk_path: str) -> List[Dict[str, Any]]:
    """Extract gene/CDS features from a GenBank file without BioPython.
    Ported from wrappers.py so this script has no BioPython dependency."""
    genes: List[Dict[str, Any]] = []
    seq_id = None
    if not _present(gbk_path):
        return genes
    with open(gbk_path, "r") as f:
        cur_type = cur_loc = None
        cur_qual: Dict[str, str] = {}
        in_feature = False

        def flush():
            if in_feature and cur_type in ("gene", "CDS") and cur_loc:
                nums = re.findall(r"(\d+)", cur_loc)
                if len(nums) >= 2:
                    genes.append({
                        "seq_id": seq_id or "",
                        "start": int(nums[0]),
                        "end": int(nums[-1]),
                        "name": cur_qual.get("gene", ""),
                        "locus_tag": cur_qual.get("locus_tag", ""),
                        "product": cur_qual.get("product", ""),
                    })

        for line in f:
            if line.startswith("LOCUS"):
                seq_id = line.split()[1]
            elif line.startswith("FEATURES"):
                continue
            elif line.startswith("ORIGIN") or line.startswith("//"):
                flush()
                in_feature = False
                cur_type = None
                continue
            if len(line) > 5 and line[:5] == "     " and line[5] != " ":
                flush()
                parts = line.strip().split(None, 1)
                cur_type = parts[0]
                cur_loc = parts[1] if len(parts) > 1 else ""
                cur_qual = {}
                in_feature = True
            elif in_feature and len(line) > 21 and line[:21] == " " * 21:
                ql = line.strip()
                if ql.startswith("/"):
                    m = re.match(r'/(\w+)="?([^"]*)"?', ql)
                    if m:
                        cur_qual[m.group(1)] = m.group(2)
    return genes


def parse_vcf(vcf_file_path: str, reference_gbk_path: str = None) -> List[Dict[str, Any]]:
    """Parse a GATK VCF into the same mutation schema as parse_breseq_html.
    Optionally annotates variants with gene info from a reference GenBank file."""
    mutations: List[Dict[str, Any]] = []
    if not _present(vcf_file_path):
        return []

    gene_lookup = _parse_gbk_genes(reference_gbk_path) if reference_gbk_path else []

    def find_gene(seq_id, pos):
        for g in gene_lookup:
            if g["start"] <= pos <= g["end"]:
                if _strip_version(g["seq_id"]) == _strip_version(seq_id) or not gene_lookup:
                    return (g["name"] or g["locus_tag"] or "unknown"), g["product"]
        return "", ""

    try:
        with open(vcf_file_path, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 8:
                    continue
                chrom = _strip_version(parts[0])
                pos = int(parts[1])
                ref, alt, qual, filt = parts[3], parts[4], parts[5], parts[6]
                if alt in (".", "<NON_REF>"):
                    continue

                for alt_allele in alt.split(","):
                    if alt_allele == "<NON_REF>":
                        continue
                    ref_len, alt_len = len(ref), len(alt_allele)
                    if ref_len == 1 and alt_len == 1:
                        var_type, mutation_str = "SNP", f"{ref}→{alt_allele}"
                    elif ref_len > alt_len:
                        var_type, mutation_str = "DEL", f"Δ{len(ref[alt_len:])} bp"
                    elif alt_len > ref_len:
                        inserted = alt_allele[ref_len:]
                        var_type = "INS"
                        mutation_str = f"+{inserted}" if len(inserted) <= 5 else f"+{len(inserted)} bp"
                    else:
                        var_type, mutation_str = "SUB", f"{ref}→{alt_allele}"

                    freq_str = "100%"
                    if len(parts) >= 10:
                        fmt_map = dict(zip(parts[8].split(":"), parts[9].split(":")))
                        if "AF" in fmt_map:
                            try:
                                af_vals = fmt_map["AF"].split(",")
                                idx = alt.split(",").index(alt_allele)
                                if idx < len(af_vals):
                                    freq_str = f"{float(af_vals[idx]) * 100:.1f}%"
                            except (ValueError, IndexError):
                                pass
                        elif "AD" in fmt_map:
                            try:
                                ad_vals = [int(x) for x in fmt_map["AD"].split(",")]
                                total = sum(ad_vals)
                                idx = alt.split(",").index(alt_allele) + 1
                                if total > 0 and len(ad_vals) > 1 and idx < len(ad_vals):
                                    freq_str = f"{(ad_vals[idx] / total) * 100:.1f}%"
                            except (ValueError, IndexError):
                                pass

                    gene_name, product = find_gene(chrom, pos)
                    qual_str = f"QUAL={qual}" if qual != "." else ""
                    filter_str = filt if filt not in (".", "PASS") else ""
                    description = " ".join(x for x in [qual_str, filter_str] if x)

                    mutations.append({
                        "type": var_type,
                        "seq_id": chrom,
                        "position": str(pos),
                        "mutation": mutation_str,
                        "freq": freq_str,
                        "annotation": product or "",
                        "gene": gene_name,
                        "description": description,
                        "details": [mutation_str, f"Freq: {freq_str}"],
                        "evidence_links": [],
                    })
    except Exception as e:  # noqa: BLE001
        print(f"Error parsing VCF: {e}", file=sys.stderr)
        return []
    return mutations


def attach_evidence_links_from_dir(mutations: List[Dict[str, Any]], evidence_dir: str) -> None:
    """Rebuild the mutation-key -> evidence-json map by scanning a pre-computed evidence/
    directory (produced by the VARIANT_EVIDENCE process), then attach links in place.

    In the working app the map was returned by extract_variant_read_evidence(); here we
    reconstruct it from the JSON files' `mutations[].mutation_key` so the report step stays
    decoupled from evidence generation.
    """
    if not evidence_dir or not os.path.isdir(evidence_dir):
        return
    key_to_rel: Dict[str, str] = {}
    for fname in os.listdir(evidence_dir):
        if not fname.endswith(".json"):
            continue
        try:
            with open(os.path.join(evidence_dir, fname)) as fh:
                payload = json.load(fh)
        except Exception:
            continue
        rel = f"evidence/{fname}"
        for m in payload.get("mutations", []):
            key = m.get("mutation_key")
            if key:
                key_to_rel[key] = rel

    for m in mutations:
        seq_id = str(m.get("seq_id", "unknown"))
        position = str(m.get("position", "0")).replace(",", "")
        mtype = str(m.get("type", "VAR"))
        rel = key_to_rel.get(f"{seq_id}_{position}_{mtype}")
        if rel:
            m["evidence_links"] = [{"text": "reads", "url": rel}]


# ─────────────────────────────────────────────────────────────────────────────
# High-confidence SNP CSV (ported from manager.run())
# ─────────────────────────────────────────────────────────────────────────────

def write_high_confidence_snps(mutations, coverage_path, out_csv) -> bool:
    """Cross-reference SNP mutations with per-base depth; write filtered CSV
    (depth>20 & freq>40). Returns True if any rows were written."""
    if not _present(coverage_path) or not mutations:
        return False

    snp_positions = set()
    for m in mutations:
        if m.get("type") == "SNP":
            try:
                snp_positions.add(int(str(m["position"]).replace(",", "")))
            except Exception:
                pass

    coverage_map: Dict[int, int] = {}
    with open(coverage_path, "r") as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) >= 3:
                try:
                    p = int(row[1])
                    if p in snp_positions:
                        coverage_map[p] = int(row[2])
                except Exception:
                    pass

    count = 0
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Position", "Reference", "Alternative", "Frequency", "Depth", "Gene", "Annotation"])
        for m in mutations:
            if m.get("type") != "SNP":
                continue
            try:
                pos = int(str(m["position"]).replace(",", ""))
                depth = coverage_map.get(pos, 0)
                freq = float(str(m["freq"]).replace("%", ""))
                if depth > 20 and freq > 40:
                    ref = alt = ""
                    mut_str = m["mutation"]
                    sep = "→" if "→" in mut_str else ("->" if "->" in mut_str else None)
                    if sep:
                        ref, alt = mut_str.split(sep)[:2]
                    writer.writerow([pos, ref, alt, freq, depth, m.get("gene", ""), m.get("annotation", "")])
                    count += 1
            except Exception:
                pass
    if count == 0:
        os.remove(out_csv)
    return count > 0


# ─────────────────────────────────────────────────────────────────────────────
# Assembly (ports PipelineManager.run() ordering / step bookkeeping)
# ─────────────────────────────────────────────────────────────────────────────

def parse_flagstat(path: str) -> Dict[str, Any]:
    """samtools flagstat text -> read_stats dict (GATK path)."""
    if not _present(path):
        return {}
    total = mapped = 0
    with open(path) as f:
        for line in f:
            m = re.match(r"(\d+)\s+\+\s+\d+\s+in total", line)
            if m:
                total = int(m.group(1))
            m = re.match(r"(\d+)\s+\+\s+\d+\s+mapped\s", line)
            if m:
                mapped = int(m.group(1))
    if total <= 0:
        return {}
    mapped_pct = round((mapped / total) * 100, 1)
    return {
        "total_reads": total,
        "mapped_percentage": mapped_pct,
        "unmapped_percentage": round(100.0 - mapped_pct, 1),
        "unmapped_reads": total - mapped,
    }


def build_results(args) -> Dict[str, Any]:
    results: Dict[str, Any] = {
        "job_id": args.job_id,
        "status": "started",
        "reference_file": os.path.basename(args.reference) if args.reference else None,
        "variant_caller": args.caller,
        "steps": {},
    }
    if args.caller == "gatk":
        results["polymorphism_mode"] = str(args.polymorphism).lower() == "true"

    # ── mutations + read_stats (caller-specific) ─────────────────────────────
    if args.caller == "gatk":
        mutations = parse_vcf(args.vcf, args.reference_gbk or args.reference)
        results["steps"]["haplotype_caller" if not results.get("polymorphism_mode") else "mutect2"] = "completed"
        attach_evidence_links_from_dir(mutations, args.evidence_dir)
        if args.evidence_dir and os.path.isdir(args.evidence_dir):
            results["steps"]["extract_read_evidence"] = "completed"
        rs = parse_flagstat(args.flagstat)
        if rs:
            results["read_stats"] = rs
    else:  # breseq
        html = os.path.join(args.breseq_dir, "output", "index.html") if args.breseq_dir else ""
        if _present(html):
            mutations = parse_breseq_html(html)
        else:
            gd = os.path.join(args.breseq_dir, "output", "output.gd") if args.breseq_dir else ""
            mutations = parse_breseq_gd(gd)
        results["steps"]["breseq"] = "completed"
        rs = parse_breseq_summary(
            os.path.join(args.breseq_dir, "output", "summary.html") if args.breseq_dir else "")
        if rs:
            results["read_stats"] = rs

    results["mutations"] = mutations

    # ── coverage + high-confidence SNPs (shared) ─────────────────────────────
    if _present(args.coverage):
        results["coverage_file"] = "coverage.tsv"
        if write_high_confidence_snps(mutations, args.coverage,
                                      os.path.join(os.path.dirname(args.out) or ".", "high_confidence_snps.csv")):
            results["hc_snps_file"] = "high_confidence_snps.csv"

    # ── unmapped read files ──────────────────────────────────────────────────
    results["unmapped_files"] = [os.path.basename(f) for f in (args.unmapped or []) if _present(f)]

    # ── MASH containment (input reads + unmapped reads) ──────────────────────
    if _present(args.mash):
        results["mash_results"] = parse_mash_screen(args.mash)
        results["steps"]["mash"] = "completed"
        if _present(args.mash_unmapped):
            results["mash_unmapped_results"] = parse_mash_screen(args.mash_unmapped)
    else:
        results["steps"]["mash"] = "skipped"

    # ── contamination: Kraken2, else BLAST fallback ──────────────────────────
    if results["unmapped_files"]:
        if _present(args.kraken):
            results["contaminants"] = parse_kraken_report(args.kraken)
            results["steps"]["kraken"] = "completed"
        else:
            results["steps"]["kraken"] = "skipped_no_db"
            if _present(args.blast):
                hits = parse_blast_results(args.blast)
                results["blast_results"] = hits
                results["steps"]["blast"] = "completed" if hits else "failed_or_empty"
    else:
        results["steps"]["kraken"] = "skipped_no_unmapped"

    results["status"] = "completed"
    return results


def main():
    p = argparse.ArgumentParser(description="Assemble StrainQC results.json")
    p.add_argument("--job-id", required=True)
    p.add_argument("--caller", required=True, choices=["breseq", "gatk"])
    p.add_argument("--polymorphism", default="false")
    p.add_argument("--reference", help="reference filename (→ results.reference_file)")

    # breseq path
    p.add_argument("--breseq-dir", help="breseq output tree (contains output/index.html, summary.html, output.gd)")

    # gatk path
    p.add_argument("--vcf")
    p.add_argument("--reference-gbk", help="GenBank for VCF gene annotation (defaults to --reference)")
    p.add_argument("--evidence-dir")
    p.add_argument("--flagstat")

    # shared
    p.add_argument("--coverage")
    p.add_argument("--kraken")
    p.add_argument("--mash")
    p.add_argument("--mash-unmapped")
    p.add_argument("--blast")
    p.add_argument("--unmapped", nargs="*", default=[])
    p.add_argument("--out", default="results.json")
    args = p.parse_args()

    results = build_results(args)
    with open(args.out, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"[emit_results] wrote {args.out}: {len(results.get('mutations', []))} mutations, "
          f"status={results['status']}", file=sys.stderr)


if __name__ == "__main__":
    main()
