"""
GATK read-evidence viewer — renders the pre-computed evidence JSON
(RESULTS_DIR/<job_id>/gatk/evidence/*.json, produced by bin/extract_evidence.py) as the
pileup HTML the frontend shows in its iframe modal.

Ported verbatim from the old app's analysis.py :: get_gatk_evidence (only settings/imports
changed). The gapped-column layout keeps insertions aligned across reads.
"""
import html
import json
import os

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import HTMLResponse

from app.config import settings

router = APIRouter()


@router.get("/evidence/{job_id}", response_class=HTMLResponse)
async def get_gatk_evidence(job_id: str, path: str = Query(..., description="Relative path under gatk/, e.g. evidence/file.json")):
    job_dir = str(settings.RESULTS_DIR / job_id)
    gatk_dir = os.path.join(job_dir, "gatk")
    normalized = os.path.normpath(path).replace("\\", "/")

    if normalized.startswith("../") or normalized == ".." or normalized.startswith("/"):
        raise HTTPException(status_code=400, detail="Invalid evidence path")
    evidence_file = os.path.normpath(os.path.join(gatk_dir, normalized))
    if not evidence_file.startswith(os.path.abspath(gatk_dir) + os.sep):
        raise HTTPException(status_code=400, detail="Invalid evidence path")
    if not os.path.exists(evidence_file):
        raise HTTPException(status_code=404, detail="Evidence file not found")

    with open(evidence_file, "r") as handle:
        payload = json.load(handle)

    window_start = payload.get("start", 1)
    window_end = payload.get("end", 1)
    window_len = window_end - window_start + 1
    ref_seq = payload.get("reference", "")

    mut_positions = {}
    for m in payload.get("mutations", []):
        mut_positions[m.get("position", 0)] = m.get("type", "VAR")

    # Collect max insertion width after each ref position across all reads (gapped layout).
    ins_widths = {}
    for read in payload.get("reads", []):
        r_start = read.get("start", window_start)
        ar = read.get("aligned_ref", "")
        ref_col = r_start - window_start
        i = 0
        while i < len(ar):
            if ar[i] == "-":
                ins_abs = window_start + (ref_col - 1)
                run = 0
                while i < len(ar) and ar[i] == "-":
                    run += 1
                    i += 1
                ins_widths[ins_abs] = max(ins_widths.get(ins_abs, 0), run)
            else:
                ref_col += 1
                i += 1

    gapped_cols = []
    for idx in range(window_len):
        abs_pos = window_start + idx
        gapped_cols.append(("ref", abs_pos))
        if abs_pos in ins_widths:
            for s in range(ins_widths[abs_pos]):
                gapped_cols.append(("ins", abs_pos, s))
    total_gapped = len(gapped_cols)

    ruler_chars = [" "] * total_gapped
    gcol = 0
    for col_spec in gapped_cols:
        if col_spec[0] == "ref" and col_spec[1] % 10 == 0:
            label = str(col_spec[1])
            for j, ch in enumerate(label):
                if gcol + j < total_gapped and gapped_cols[gcol + j][0] == "ref":
                    ruler_chars[gcol + j] = ch
        gcol += 1
    ruler_text = html.escape("".join(ruler_chars))

    ref_cells = []
    for col_spec in gapped_cols:
        if col_spec[0] == "ref":
            abs_pos = col_spec[1]
            ref_idx = abs_pos - window_start
            base = ref_seq[ref_idx] if ref_idx < len(ref_seq) else "N"
            if abs_pos in mut_positions:
                ref_cells.append(f'<span class="ref-mut" title="{abs_pos}">{html.escape(base)}</span>')
            else:
                ref_cells.append(f'<span class="ref-base">{html.escape(base)}</span>')
        else:
            ref_cells.append('<span class="ins-dash">-</span>')
    ref_row_html = "".join(ref_cells)

    def build_seq_row(read):
        read_start = read.get("start", window_start)
        aligned_ref = read.get("aligned_ref", "")
        aligned_read = read.get("aligned_read", "")
        aligned_qual = read.get("aligned_qual", [])
        ref_bases, ins_bases = {}, {}
        ref_col = read_start - window_start
        ins_counter = {}
        for i in range(min(len(aligned_ref), len(aligned_read))):
            rb, rd = aligned_ref[i], aligned_read[i]
            qv = aligned_qual[i] if i < len(aligned_qual) else None
            if rb == "-":
                ins_abs = window_start + (ref_col - 1)
                sub = ins_counter.get(ins_abs, 0)
                ins_bases[(ins_abs, sub)] = (rd, qv)
                ins_counter[ins_abs] = sub + 1
            else:
                abs_pos = window_start + ref_col
                if 0 <= ref_col < window_len:
                    ref_bases[abs_pos] = (rd, qv)
                ref_col += 1
        parts = []
        for col_spec in gapped_cols:
            if col_spec[0] == "ref":
                abs_pos = col_spec[1]
                if abs_pos not in ref_bases:
                    parts.append('<span class="dot">·</span>')
                else:
                    rd, qv = ref_bases[abs_pos]
                    ref_idx = abs_pos - window_start
                    ref_base = ref_seq[ref_idx] if ref_idx < len(ref_seq) else "N"
                    if rd == "-":
                        css = "del-base"
                    elif rd.upper() != ref_base.upper():
                        css = "mm-base"
                    else:
                        css = "match-base"
                    if qv is not None and css != "del-base":
                        if qv < 10:
                            css += " q-low"
                        elif qv < 20:
                            css += " q-mid"
                    parts.append(f'<span class="{css}">{html.escape(rd if rd != "-" else chr(0x2013))}</span>')
            else:
                _, ins_abs, sub = col_spec
                key = (ins_abs, sub)
                if key not in ins_bases:
                    parts.append('<span class="dot">·</span>')
                else:
                    rd, qv = ins_bases[key]
                    css = "ins-base"
                    if qv is not None:
                        if qv < 10:
                            css += " q-low"
                        elif qv < 20:
                            css += " q-mid"
                    parts.append(f'<span class="{css}">{html.escape(rd)}</span>')
        return "".join(parts)

    table_rows = [
        f'<tr><td class="label-col"></td><td class="seq-col"><code class="ruler-code">{ruler_text}</code></td></tr>',
        f'<tr class="ref-row-top"><td class="label-col label-ref">REF</td><td class="seq-col"><code>{ref_row_html}</code></td></tr>',
    ]
    for read in payload.get("reads", []):
        strand = read.get("strand", "+")
        arrow = "→" if strand == "+" else "←"
        arrow_cls = "strand-fwd" if strand == "+" else "strand-rev"
        supports_alt = bool(read.get("supports"))
        supports_ref = bool(read.get("matches_ref"))
        if supports_alt and not supports_ref:
            row_cls = "read-row read-row-alt"
        elif supports_ref and not supports_alt:
            row_cls = "read-row read-row-ref"
        elif supports_alt and supports_ref:
            row_cls = "read-row read-row-mixed"
        else:
            row_cls = "read-row"
        label_html = (
            f'<span class="rname">{html.escape(read.get("name", "read"))}</span>'
            f'<span class="{arrow_cls}">{arrow}</span>'
            f'<span class="rmapq">MAPQ&nbsp;{read.get("mapq", 0)}</span>'
            f'<span class="rcigar">{html.escape(read.get("cigar", ""))}</span>'
        )
        table_rows.append(
            f'<tr class="{row_cls}"><td class="label-col"><div class="read-label">{label_html}</div></td>'
            f'<td class="seq-col"><code>{build_seq_row(read)}</code></td></tr>'
        )
    table_rows.append(
        f'<tr class="ref-row-bot"><td class="label-col label-ref">REF</td><td class="seq-col"><code>{ref_row_html}</code></td></tr>'
    )

    mutation_rows = []
    for m in payload.get("mutations", []):
        mutation_rows.append(
            f"<tr><td>{html.escape(str(m.get('position', '')))}</td>"
            f"<td>{html.escape(str(m.get('type', '')))}</td>"
            f"<td>{html.escape(str(m.get('mutation', '')))}</td>"
            f"<td>{html.escape(str(m.get('gene', '')))}</td>"
            f"<td>{html.escape(str(m.get('freq', '')))}</td>"
            f"<td>{m.get('alt_support', 0)}</td>"
            f"<td title='Unfiltered raw count. GATK may have excluded these reads due to low base quality or other filters.'>{m.get('ref_support', 0)}</td></tr>"
        )

    html_body = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"/><title>GATK Evidence</title>
<style>
  body {{ font-family: system-ui, sans-serif; margin: 16px; color: #1f2937; font-size: 13px; }}
  h2, h3 {{ margin: 0 0 10px; }}
  .summary {{ margin-bottom: 14px; padding: 10px 14px; background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 6px; display: flex; flex-wrap: wrap; gap: 18px; align-items: baseline; }}
  .summary strong {{ color: #1e40af; }}
  .mut-table {{ border-collapse: collapse; margin: 10px 0 16px; font-size: 0.9em; }}
  .mut-table th, .mut-table td {{ border: 1px solid #dbe2ea; padding: 6px 10px; text-align: left; }}
  .mut-table th {{ background: #f1f5f9; }}
  .stack-wrap {{ overflow-x: auto; padding: 8px 6px; background: #fafafa; }}
  .alignment-table {{ border-collapse: collapse; border-spacing: 0; }}
  .alignment-table td {{ border: none !important; padding: 0; margin: 0; }}
  .label-col {{ width: 220px; min-width: 220px; max-width: 220px; padding: 0 6px 0 2px; vertical-align: middle; font-size: 0.75em; line-height: 1.0; white-space: nowrap; overflow: hidden; }}
  .label-ref {{ color: #64748b; text-transform: uppercase; font-weight: 600; font-size: 0.75em; letter-spacing: .04em; text-align: right; }}
  .read-label {{ display: flex; align-items: center; gap: 5px; }}
  .rname {{ color: #374151; font-weight: 600; overflow: hidden; text-overflow: ellipsis; max-width: 120px; }}
  .strand-fwd {{ color: #1d4ed8; font-weight: 700; }}
  .strand-rev {{ color: #b45309; font-weight: 700; }}
  .rmapq {{ color: #6b7280; }}
  .rcigar {{ color: #94a3b8; }}
  .seq-col {{ padding-left: 6px; vertical-align: middle; }}
  .seq-col code {{ font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: 11px; white-space: pre; display: block; line-height: 1.0; }}
  .seq-col code span {{ display: inline; }}
  .ruler-code {{ color: #94a3b8; font-size: 11px; line-height: 1.0; }}
  .ref-base {{ color: #374151; }}
  .ref-mut {{ background: #fef9c3; color: #92400e; font-weight: 700; text-decoration: underline; text-decoration-color: #f59e0b; }}
  .dot {{ color: #d1d5db; }}
  .match-base {{ color: #374151; }}
  .mm-base {{ background: #fee2e2; color: #b91c1c; font-weight: 700; }}
  .del-base {{ background: #ede9fe; color: #7c3aed; font-weight: 700; }}
  .ins-base {{ background: #dcfce7; color: #15803d; font-weight: 700; }}
  .ins-dash {{ color: #94a3b8; }}
  .q-low {{ opacity: 0.4; }}
  .q-mid {{ opacity: 0.72; }}
  .read-row-alt td {{ background: #f0fdf4; }}
  .read-row-ref td {{ background: #fff7ed; }}
  .read-row-mixed td {{ background: #eff6ff; }}
  .legend {{ display: flex; flex-wrap: wrap; gap: 12px; margin-top: 14px; font-size: 0.82em; color: #475569; }}
  .leg-item {{ display: flex; align-items: center; gap: 4px; }}
  .leg-swatch {{ display: inline-block; width: 14px; height: 14px; border-radius: 2px; border: 1px solid #ccc; font-family: monospace; font-size: 10px; text-align: center; line-height: 14px; }}
</style></head>
<body>
<h2>Mutation Detailed Evidence</h2>
<div class="summary">
  <span><strong>{html.escape(str(payload.get('chrom', '')))}:{payload.get('start', 0)}&ndash;{payload.get('end', 0)}</strong></span>
  <span>Displayed reads: <strong>{payload.get('displayed_reads', 0)}</strong></span>
  <span>Total overlaps: <strong>{payload.get('total_overlapping_reads', 0)}</strong></span>
  <span>Window: <strong>{payload.get('window', 0)} bp</strong></span>
  <span>Truncated: <strong>{'yes' if payload.get('truncated') else 'no'}</strong></span>
</div>
<h3>Mutations In Window</h3>
<table class="mut-table">
  <thead><tr><th>Position</th><th>Type</th><th>Mutation</th><th>Gene</th><th>Freq (GATK)</th><th>Alt reads (raw)</th><th title='Unfiltered. GATK may have excluded these via base quality / strand bias filters.'>Ref reads (raw) &#9432;</th></tr></thead>
  <tbody>{''.join(mutation_rows)}</tbody>
</table>
<div class="stack-wrap"><table class="alignment-table">{''.join(table_rows)}</table></div>
<div class="legend">
  <span>Alignment Legend:</span>
  <span class="leg-item"><span class="leg-swatch mm-base">A</span> mismatch</span>
  <span class="leg-item"><span class="leg-swatch del-base">&ndash;</span> deletion</span>
  <span class="leg-item"><span class="leg-swatch ins-base">A</span> insertion</span>
  <span class="leg-item"><span class="leg-swatch ref-mut">A</span> mutation site</span>
  <span class="leg-item"><span class="strand-fwd">&rarr;</span>&thinsp;fwd <span class="strand-rev">&larr;</span>&thinsp;rev</span>
  <span class="leg-item"><span class="leg-swatch" style="background:#f0fdf4">&nbsp;</span> supports alt</span>
  <span class="leg-item"><span class="leg-swatch" style="background:#fff7ed">&nbsp;</span> supports ref</span>
  <span class="leg-item"><span class="leg-swatch" style="background:#eff6ff">&nbsp;</span> mixed</span>
</div>
</body></html>"""
    return HTMLResponse(content=html_body, headers={"Cache-Control": "no-store"})
