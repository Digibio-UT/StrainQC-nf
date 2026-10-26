"""Minimal dependency-free .xlsx reader (stdlib zipfile + xml only).

Just enough to read the Ottilie supplementary tables (one sheet each, string/number
cells, shared-strings and inline-string forms). Keeps the benchmark runnable without
pandas/openpyxl, matching the rest of the repo's stdlib approach (bin/emit_results.py).
"""

import re
import zipfile
from xml.etree import ElementTree as ET

_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def _col_to_idx(cell_ref):
    """'B7' -> 1 (zero-based column index)."""
    letters = re.match(r"([A-Z]+)", cell_ref).group(1)
    idx = 0
    for ch in letters:
        idx = idx * 26 + (ord(ch) - 64)
    return idx - 1


def read_rows(path):
    """Return the first worksheet as a list of rows (each a list of cell strings)."""
    z = zipfile.ZipFile(path)
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        sroot = ET.fromstring(z.read("xl/sharedStrings.xml"))
        for si in sroot.iter(_NS + "si"):
            shared.append("".join(t.text or "" for t in si.iter(_NS + "t")))
    sheet = sorted(n for n in z.namelist()
                   if re.match(r"xl/worksheets/sheet\d+\.xml", n))[0]
    root = ET.fromstring(z.read(sheet))
    rows = []
    for r in root.iter(_NS + "row"):
        cells = {}
        for c in r.findall(_NS + "c"):
            t = c.get("t")
            v = c.find(_NS + "v")
            txt = ""
            if v is not None:
                txt = shared[int(v.text)] if t == "s" else (v.text or "")
            inline = c.find(_NS + "is")
            if inline is not None:
                txt = "".join(x.text or "" for x in inline.iter(_NS + "t"))
            cells[_col_to_idx(c.get("r"))] = txt
        width = max(cells) + 1 if cells else 0
        rows.append([cells.get(i, "") for i in range(width)])
    return rows


def read_records(path, header_row=1):
    """Return (headers, list-of-dicts) using the given 0-based header row index.

    The Ottilie sheets carry a title in row 0 and column names in row 1.
    """
    rows = read_rows(path)
    headers = [h.strip() for h in rows[header_row]]
    records = []
    for r in rows[header_row + 1:]:
        if not any(str(x).strip() for x in r):
            continue
        rec = {headers[i]: (r[i] if i < len(r) else "") for i in range(len(headers))}
        records.append(rec)
    return headers, records
