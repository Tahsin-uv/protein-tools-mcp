"""Minimal mmCIF (PDBx) parser and mmCIF -> PDB-atom bridge.

Handles the constructs that matter for structural data: data_ blocks,
loop_ sections, key-value pairs, quoted values, '#' comments, and
'.'/'?' missing markers. Semicolon multiline text fields are supported in
key-value context; inside loop rows they terminate the loop gracefully
(they do not occur in _atom_site tables, which is all we need).
"""
from __future__ import annotations

from typing import Any


def _split_tokens(line: str) -> list[str]:
    """Split a CIF line into values, honoring single/double quotes."""
    out: list[str] = []
    i, n = 0, len(line)
    while i < n:
        c = line[i]
        if c in " \t":
            i += 1
            continue
        if c in "'\"":
            j = i + 1
            while True:
                j = line.find(c, j)
                if j == -1:  # unterminated quote: take rest of line
                    out.append(line[i + 1:])
                    i = n
                    break
                if j + 1 >= n or line[j + 1] in " \t":
                    out.append(line[i + 1:j])
                    i = j + 1
                    break
                j += 1
        else:
            j = line.find(" ", i)
            j2 = line.find("\t", i)
            if j2 != -1 and (j2 < j or j == -1):
                j = j2
            if j == -1:
                out.append(line[i:])
                i = n
            else:
                out.append(line[i:j])
                i = j
    return out


def _norm_val(v: str):
    """CIF missing-value markers '.' and '?' become None."""
    if v in (".", "?"):
        return None
    return v
def parse_cif(text: str) -> dict[str, Any]:
    """Parse mmCIF text into {"block", "categories": {name: {columns, rows}}}.

    Loop categories get one dict per row; key-value categories get a single
    row dict, so consumers treat both shapes identically.
    """
    lines = text.splitlines()
    categories: dict[str, dict] = {}
    block = None
    i, n = 0, len(lines)

    def add_kv(cat: str, col: str, val) -> None:
        c = categories.setdefault(cat, {"columns": [], "rows": [{}]})
        if not isinstance(c["rows"], list) or not c["rows"]:
            c["rows"] = [{}]
        if col not in c["columns"]:
            c["columns"].append(col)
        c["rows"][0][col] = val

    while i < n:
        s = lines[i].strip()
        if not s or s.startswith("#"):
            i += 1
            continue
        if s.startswith("data_"):
            block = s[5:].strip()
            i += 1
            continue
        if s.lower() == "global_" or s.startswith("save_"):
            i += 1
            continue
        if s.lower() == "loop_":
            i += 1
            headers = []
            while i < n and lines[i].strip().startswith("_"):
                headers.append(lines[i].strip())
                i += 1
            cat = headers[0].split(".", 1)[0].lstrip("_") if headers else ""
            cols = [h.split(".", 1)[1] if "." in h else h.lstrip("_")
                    for h in headers]
            rows: list[dict] = []
            terminated = False
            while i < n:
                s2 = lines[i].strip()
                if (not s2 or s2.startswith("#") or s2.startswith("_")
                        or s2.lower().startswith("loop_") or s2.startswith("data_")):
                    break
                if s2.startswith(";"):
                    i += 1
                    while i < n and not lines[i].startswith(";"):
                        i += 1
                    i += 1
                    terminated = True
                    break
                vals = [_norm_val(v) for v in _split_tokens(lines[i])]
                i += 1
                if not vals:
                    continue
                if len(vals) != len(cols):
                    vals = (vals + [None] * len(cols))[:len(cols)]
                rows.append(dict(zip(cols, vals)))
            if headers:
                c = categories.setdefault(cat, {"columns": cols, "rows": []})
                c["columns"] = cols
                c["rows"].extend(rows)
                if terminated:
                    c["multiline_truncated"] = True
            continue
        if s.startswith("_"):
            toks = _split_tokens(s)
            j = 0
            consumed_next = False
            while j < len(toks):
                name = toks[j]
                if not name.startswith("_"):
                    j += 1
                    continue
                cat, _, col = name[1:].partition(".")
                if j + 1 < len(toks):
                    add_kv(cat, col, _norm_val(toks[j + 1]))
                    j += 2
                else:
                    i += 1
                    consumed_next = True
                    if i < n and lines[i].startswith(";"):
                        buf = [lines[i][1:]]
                        i += 1
                        while i < n and not lines[i].startswith(";"):
                            buf.append(lines[i])
                            i += 1
                        i += 1
                        add_kv(cat, col, "\n".join(buf).strip())
                    elif i < n:
                        vtoks = _split_tokens(lines[i])
                        if vtoks:
                            add_kv(cat, col, _norm_val(vtoks[0]))
                    break
            if not consumed_next:
                i += 1
            continue
        i += 1  # unrecognized line: skip (guarantees progress)
    return {"block": block, "categories": categories}
def to_pdb_atoms(parsed: dict) -> dict:
    """Convert a parsed mmCIF's _atom_site table into the atom dicts that
    pdb_utils analyses consume (same shape as parse_pdb output)."""
    rows = (parsed.get("categories", {}).get("atom_site", {})
            .get("rows", []))
    atoms: list[dict] = []
    for r in rows:
        try:
            xyz = (float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"]))
        except (KeyError, TypeError, ValueError):
            continue  # skip rows without usable coordinates
        name = str(r.get("label_atom_id") or "?").strip("'\"")
        element = (r.get("type_symbol") or r.get("label_type_symbol") or "").strip()
        b = r.get("B_iso_or_equiv")
        try:
            b = float(b) if b not in (None, ".", "?") else None
        except (TypeError, ValueError):
            b = None
        atoms.append({
            "record": str(r.get("group_PDB") or "ATOM").strip().upper(),
            "name": name,
            "altloc": str(r.get("label_alt_id") or "").strip().strip("'\""),
            "resname": str(r.get("label_comp_id") or "?").strip(),
            "chain": str(r.get("auth_asym_id")
                         or r.get("label_asym_id") or "?").strip(),
            "resseq": str(r.get("auth_seq_id") or r.get("label_seq_id") or "0"),
            "icode": str(r.get("pdbx_PDB_ins_code") or "").strip().replace("?", ""),
            "element": element or (name[0] if name else "?"),
            "model": None,
            "bfactor": b,
            "xyz": xyz,
            "line": 0,
        })
    return {"atoms": atoms, "errors": [], "models": []}


