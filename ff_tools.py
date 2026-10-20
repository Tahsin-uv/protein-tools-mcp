"""GROMACS force-field metal parameter checks (stdlib only).

Parses [ atomtypes ] sections from .itp/.ff files and reports which metal
atomtypes exist, with their sigma/epsilon. This is the pre-flight check that
would have caught the amber99sb-ildn ZN-atomtype failure before pdb2gmx.
"""
from __future__ import annotations

import os
from typing import Any


def _strip_comment(line: str) -> str:
    return line.split(";", 1)[0].strip()


def parse_atomtypes(text: str) -> list[dict]:
    """Extract [ atomtypes ] entries from itp/ff text.

    Expected row shape: name at.num mass charge ptype sigma epsilon [s14 e14]
    """
    out: list[dict] = []
    in_section = False
    for lineno, raw in enumerate(text.splitlines(), 1):
        line = _strip_comment(raw)
        if not line:
            continue
        if line.startswith("["):
            inner = line.replace("[", " ").replace("]", " ").strip().lower()
            in_section = inner == "atomtypes"
            continue
        if not in_section:
            continue
        toks = line.split()
        if len(toks) < 7:
            continue
        try:
            mass = float(toks[2])
            charge = float(toks[3])
            sigma = float(toks[5])
            eps = float(toks[6])
        except ValueError:
            continue
        out.append({"name": toks[0].upper(), "at_num": toks[1], "mass": mass,
                    "charge": charge, "ptype": toks[4], "sigma_nm": sigma,
                    "epsilon_kj_mol": eps, "line": lineno,
                    "raw_fields": len(toks)})
    return out


def _ff_files(path: str) -> list[str]:
    if os.path.isdir(path):
        files = []
        for root, _dirs, names in os.walk(path):
            for nm in sorted(names):
                if nm.endswith((".itp", ".ff", ".top")):
                    files.append(os.path.join(root, nm))
        return files
    return [path]


def check_metal_params(ff_path: str, metals=("ZN", "FE")) -> dict[str, Any]:
    """Scan a force-field directory or single .itp/.ff file for metal
    atomtypes. Reports found entries (sigma/epsilon) and missing metals."""
    metals_u = [m.upper() for m in metals]
    files = _ff_files(ff_path)
    found: dict[str, list] = {}
    scanned = 0
    errors = []
    for f in files:
        try:
            text = open(f, errors="replace").read()
        except OSError as exc:
            errors.append({"file": f, "error": str(exc)})
            continue
        scanned += 1
        for at in parse_atomtypes(text):
            if at["name"] in metals_u:
                found.setdefault(at["name"], []).append({
                    "file": f, "line": at["line"],
                    "sigma_nm": at["sigma_nm"],
                    "epsilon_kj_mol": at["epsilon_kj_mol"],
                    "charge": at["charge"], "mass": at["mass"]})
    missing = [m for m in metals_u if m not in found]
    return {"metals_checked": metals_u,
            "scanned_files": scanned,
            "paths_scanned": files[:20],
            "read_errors": errors[:5],
            "found": {k: v[:5] for k, v in found.items()},
            "missing": missing,
            "verdict": "OK - all requested metal atomtypes present"
                       if not missing else
                       f"MISSING atomtypes: {', '.join(missing)}"
                       " (pdb2gmx WILL FAIL or need a custom .itp)"}


def preflight(source, ff_dir: str | None = None,
              metals: list[str] | None = None) -> dict[str, Any]:
    """MD pre-flight gate: structure validation + metal-site detection +
    force-field parameter check. Returns GO / NO-GO with blocking issues.

    `source` may be PDB text, a parsed structure dict, or a file path."""
    import pdb_utils  # local import keeps module import order simple

    if isinstance(source, str) and "\n" not in source:
        from pathlib import Path
        source = Path(source).read_text(errors="replace")
    v = pdb_utils.validate_pdb(source)
    p = pdb_utils._as_parsed(source)
    atoms, _ = pdb_utils._single_model(p)
    present = sorted({a["element"].strip().upper() for a in atoms}
                     & pdb_utils.METAL_ELEMENTS)
    blocking, warnings = [], []
    if v["format_errors"] or v["non_integer_resseq"]:
        blocking.append({"check": "pdb_format",
                         "detail": {"format_errors": v["format_errors"][:5],
                                    "non_integer_resseq": v["non_integer_resseq"]}})
    for key in ("chain_gaps", "unknown_residues", "altloc_present"):
        if v[key]:
            warnings.append({"check": key, "detail": v[key][:8]})
    geom = {}
    if present:
        for m in present:
            g = pdb_utils.metal_geometry(source, m)
            geom[m] = {s["site"]: s["coordination_number"]
                       for s in g["metal_sites"]}
        warnings.append({"check": "metal_sites_present",
                         "detail": "metal ions detected in the structure; "
                                   "verify each site geometry before MD"})
    param_check = None
    if ff_dir:
        param_check = check_metal_params(ff_dir, metals=metals or present
                                         or ["ZN", "FE"])
        if param_check["missing"]:
            blocking.append({"check": "ff_metal_params",
                             "detail": ("missing atomtypes: "
                                        + ", ".join(param_check["missing"])
                                        + " - supply custom .itp parameters "
                                        "(e.g. Li & Merz 12-6 single-type)")})
    return {"verdict": "GO" if not blocking else "NO-GO",
            "blocking": blocking, "warnings": warnings,
            "metals_present_in_structure": present,
            "metal_site_summary": geom, "param_check": param_check}
