"""PROPKA .pka output parsing and protonation-state summary (stdlib only)."""
from __future__ import annotations

import re
from typing import Any

# primary determinant-table row: RESNAME num chain pKa buried% ...
_ROW = re.compile(
    r"^([A-Z+]{2,4})\s+(-?\d+)\s+([A-Za-z0-9]+)\s+(-?\d+\.\d+)\s+(\d+)\s*%")


def parse_pka(text: str) -> dict[str, Any]:
    """Parse primary rows of the PROPKA per-residue determinant table,
    plus the pI line. Continuation (interaction-detail) rows don't match."""
    rows: list[dict] = []
    pI = None
    for line in text.splitlines():
        m = _ROW.match(line)
        if m:
            rows.append({"resname": m.group(1), "resnum": int(m.group(2)),
                         "chain": m.group(3), "pka": float(m.group(4)),
                         "buried_pct": int(m.group(5))})
            continue
        if pI is None:
            mm = re.search(r"pI is\s+([\d.]+)", line)
            if mm:
                pI = float(mm.group(1))
    return {"rows": rows, "pI_folded": pI}


def summarize(text: str) -> dict[str, Any]:
    """Protonation-state summary with metalloprotein-relevant flags."""
    d = parse_pka(text)
    rows = d["rows"]
    if not rows:
        raise ValueError("no pKa rows parsed - is this a PROPKA .pka file?")
    by: dict[str, list] = {}
    for r in rows:
        by.setdefault(r["resname"], []).append(r)
    concerns = []
    for r in rows:
        res, pk = r["resname"], r["pka"]
        if res in ("ASP", "GLU") and pk > 5.5:
            concerns.append({**r, "flag": "acidic residue with raised pKa - "
                             "may be protonated near pH 7"})
        elif res == "HIS" and 5.5 <= pk <= 7.5:
            concerns.append({**r, "flag": "HIS near neutral pH - HID/HIE "
                             "choice matters; check metal coordination"})
        elif res in ("LYS", "ARG") and pk < 9.0:
            concerns.append({**r, "flag": "basic residue with lowered pKa - "
                             "may be neutral"})
        elif res == "CYS" and pk >= 99.0:
            concerns.append({**r, "flag": "CYS flagged as disulfide/buried "
                             "(99.99 marker) - will break if treated as free thiol"})
        elif res == "CYS" and pk < 9.0:
            concerns.append({**r, "flag": "CYS with low pKa - thiolate "
                             "competent; check metal binding"})
    return {"n_residues_parsed": len(rows),
            "resname_counts": {k: len(v) for k, v in sorted(by.items())},
            "pI_folded": d["pI_folded"],
            "extremes": {"highest_pka": max(rows, key=lambda r: r["pka"]),
                         "lowest_pka": min(rows, key=lambda r: r["pka"])},
            "protonation_concerns": concerns[:30],
            "n_concerns": len(concerns),
            "note": "parsed from the per-residue determinant table "
                    "(RESIDUE num chain pKa buried% ...)"}
