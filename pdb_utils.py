"""Pure-python PDB parsing and protein-structure analysis helpers.

Deliberately stdlib-only and written against the fixed-column PDB format spec,
because the most common failure mode in MD prep pipelines is column
misalignment (e.g. atom names written with wrong padding, producing residues
like 'T A' instead of 'TA').
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict
from itertools import combinations
from typing import Any

STANDARD_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
# common protonation-state variants seen in MD prep (GROMACS/Amber naming)
VARIANT_RESIDUES = {"HID", "HIE", "HIP", "HISD", "HISE", "CYM", "CYX",
                    "ASH", "GLH", "LYN"}
SOLVENT = {"HOH", "WAT", "SOL", "TIP3", "GOL", "EDO", "SO4", "PO4", "NA", "CL"}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "HID": "H", "HIE": "H", "HIP": "H", "HISD": "H",
    "HISE": "H", "CYM": "C", "CYX": "C", "ASH": "D", "GLH": "E", "LYN": "K",
    "MSE": "M", "SEC": "U", "PYL": "O",
}

GEOMETRY_GUESS = {
    2: "linear", 3: "trigonal-planar", 4: "tetrahedral or square-planar",
    5: "square-pyramidal or trigonal-bipyramidal", 6: "octahedral",
    7: "pentagonal-bipyramidal", 8: "cubic or hexagonal-bipyramidal",
}


def parse_pdb(text: str) -> dict[str, Any]:
    """Strict fixed-column parse of ATOM/HETATM records; collects format errors.

    Tracks insertion codes (column 27) and MODEL records: every atom carries
    an "icode" and a "model" tag (None for files without MODEL records), and
    the model IDs are listed in the result. Analysis helpers restrict to the
    first model via _single_model().
    """
    atoms: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []
    model_ids: list[str] = []
    current_model: str | None = None
    for lineno, line in enumerate(text.splitlines(), 1):
        rec = line[:6].strip().upper()
        if rec == "MODEL":
            current_model = line[10:14].strip() or str(len(model_ids) + 1)
            model_ids.append(current_model)
            continue
        if rec == "ENDMDL":
            current_model = None
            continue
        if rec not in ("ATOM", "HETATM"):
            continue
        if len(line.rstrip()) < 54:
            errors.append({"line": lineno, "issue": "line too short for coordinates",
                           "text": line[:60]})
            continue
        a: dict[str, Any] = {
            "record": rec,
            "name": line[12:16].strip(),
            "altloc": line[16:17].strip(),
            "resname": line[17:20].strip(),
            "chain": line[21:22].strip(),
            "resseq": line[22:26].strip(),
            "icode": line[26:27].strip(),
            "element": line[76:78].strip(),
            "model": current_model,
            "line": lineno,
        }
        try:
            a["xyz"] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            a["bfactor"] = float(line[60:66]) if len(line.rstrip()) >= 66 else None
        except ValueError:
            errors.append({"line": lineno, "issue": "non-numeric coordinate field",
                           "text": line[:60]})
            continue
        if not a["element"]:
            a["element"] = a["name"][0] if a["name"] else "?"
        atoms.append(a)
    return {"atoms": atoms, "errors": errors, "models": model_ids}


def _as_parsed(source):
    """Accept either PDB text or an already-parsed structure dict
    (e.g. from cif_tools.to_pdb_atoms) so all analyses take both."""
    if isinstance(source, str):
        return parse_pdb(source)
    return source  # parse_pdb-shaped dict from cif_tools


def chain_key(a: dict) -> str:
    return a["chain"] or "_"


def res_key(a: dict) -> tuple[str, str, str]:
    return (a["chain"], a["resseq"], a.get("icode", ""))


METAL_ELEMENTS = {"ZN", "FE", "MN", "MG", "CU", "NI", "CO", "CD", "HG"}


def _single_model(p: dict) -> tuple[list, int]:
    """If the file contains MODEL records (NMR etc.), restrict analysis to the
    first model. Returns (atoms, n_models); n_models == 0 for plain files."""
    mods = {a["model"] for a in p["atoms"] if a["model"] is not None}
    if not mods:
        return p["atoms"], 0
    first = sorted(mods, key=lambda x: (0, int(x), "") if str(x).isdigit()
                   else (1, 0, str(x)))[0]
    kept = [a for a in p["atoms"] if a["model"] in (None, first)]
    return kept, len(mods)


def _first_altloc(atoms: list) -> list:
    """Keep one conformer per (residue, atom name): prefer blank/'A' altLoc,
    else the first seen. Standard practice for structural analyses."""
    by: dict[tuple, list] = defaultdict(list)
    order: list[tuple] = []
    for a in atoms:
        k = (a["chain"], a["resseq"], a.get("icode", ""), a["name"])
        if k not in by:
            order.append(k)
        by[k].append(a)
    out = []
    for k in order:
        group = by[k]
        pref = [a for a in group if a["altloc"] in ("", " ", "A")]
        out.append(pref[0] if pref else group[0])
    return out


def validate_pdb(source: Any) -> dict[str, Any]:
    """Validate PDB text or an mmCIF-derived structure dict (see cif_tools)."""
    p = _as_parsed(source)
    atoms, _ = _single_model(p)
    polymer = [a for a in atoms if a["resname"] in STANDARD_RESIDUES | VARIANT_RESIDUES]
    waters = [a for a in atoms if a["resname"] in ("HOH", "WAT", "SOL")]
    het = [a for a in atoms if a["record"] == "HETATM"
           and a["resname"] not in SOLVENT | STANDARD_RESIDUES | VARIANT_RESIDUES]
    residues = sorted({res_key(a) for a in polymer},
                      key=lambda r: (r[0],
                                     int(r[1]) if r[1].lstrip("-").isdigit() else 10 ** 9,
                                     r[2]))
    unknown_res = sorted({a["resname"] for a in atoms}
                         - STANDARD_RESIDUES - VARIANT_RESIDUES - SOLVENT)
    altlocs = sorted({a["altloc"] for a in atoms if a["altloc"]})
    icodes = sorted({(a["chain"], a["resseq"], a["icode"]) for a in atoms if a["icode"]})

    gaps = []
    by_chain: dict[str, list[int]] = defaultdict(list)
    for ch, rs, _ic in residues:
        if rs.lstrip("-").isdigit():
            by_chain[ch].append(int(rs))
    for ch, rsl in by_chain.items():
        nums = sorted(set(rsl))
        for a, b in zip(nums, nums[1:]):
            if b - a > 1:
                gaps.append({"chain": ch, "from": a, "to": b, "missing": b - a - 1})

    seen: Counter = Counter()
    for a in polymer:
        seen[(a["chain"], a["resseq"], a["icode"], a["name"])] += 1
    dups = sorted(":".join(k) for k, v in seen.items() if v > 1)

    out = {
        "ok": not p["errors"] and not unknown_res and not non_integer_residues(residues),
        "counts": {"atoms": len(atoms), "polymer_atoms": len(polymer),
                   "residues": len(residues), "waters": len(waters),
                   "other_het": len(het),
                   "chains": sorted({chain_key(a) for a in atoms})},
        "format_errors": p["errors"][:20],
        "unknown_residues": unknown_res,
        "altloc_present": altlocs,
        "insertion_code_residues": [f"{ch}{rs}{ic}" for ch, rs, ic in icodes[:30]],
        "chain_gaps": gaps,
        "duplicate_atom_names": dups[:20],
        "non_integer_resseq": [r for _, r, _ in residues if not r.lstrip("-").isdigit()],
    }
    if p["models"]:
        out["multi_model_note"] = (f"file contains {len(p['models'])} MODEL records; "
                                   "counts and analyses use model 1 only")
    return out


def non_integer_residues(residues) -> list[str]:
    return [r for _, r, _ in residues if not r.lstrip("-").isdigit()]
def _angle_deg(v1: tuple, v2: tuple) -> float:
    dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]
    n1 = math.sqrt(v1[0] ** 2 + v1[1] ** 2 + v1[2] ** 2)
    n2 = math.sqrt(v2[0] ** 2 + v2[1] ** 2 + v2[2] ** 2)
    if n1 == 0 or n2 == 0:
        return 0.0
    return math.degrees(math.acos(max(-1.0, min(1.0, dot / (n1 * n2)))))


def _ideal_geometries(cn: int) -> dict[str, list[float]]:
    """Ideal sorted pairwise-angle sets (degrees) for common coordination
    geometries; compared element-wise against observed sorted angles."""
    if cn == 2:
        return {"linear": [180.0]}
    if cn == 3:
        return {"trigonal_planar": [120.0] * 3}
    if cn == 4:
        return {"tetrahedral": [109.47] * 6,
                "square_planar": [90.0] * 4 + [180.0] * 2}
    if cn == 5:
        return {"trigonal_bipyramidal": [90.0] * 6 + [120.0] * 3 + [180.0],
                "square_pyramidal": [90.0] * 8 + [180.0] * 2}
    if cn == 6:
        return {"octahedral": [90.0] * 12 + [180.0] * 3}
    if cn == 7:
        # pentagonal bipyramid: 5 eq pairs at 72, 5 at 144, 10 ax-eq at 90,
        # 1 ax-ax at 180  (C(7,2) = 21 pairs)
        return {"pentagonal_bipyramidal": [90.0] * 10 + [72.0] * 5
                + [144.0] * 5 + [180.0]}
    if cn == 8:
        # cubic: 12 edges at 70.53, 12 face diagonals at 109.47, 4 at 180;
        # hexagonal bipyramid: 6 eq pairs 60/120/180-adjacent + ax pairs
        return {"cubic": [70.53] * 12 + [109.47] * 12 + [180.0] * 4,
                "hexagonal_bipyramidal": [60.0] * 6 + [120.0] * 6
                + [180.0] * 3 + [90.0] * 12 + [180.0]}
    return {}


def _classify_ligand(a: dict) -> str:
    el = a["element"].strip().upper()
    res, name = a["resname"], a["name"]
    if res in ("HOH", "WAT", "SOL"):
        return "water"
    if el == "N":
        if res in ("HIS", "HID", "HIE", "HIP", "HISD", "HISE") and name in ("ND1", "NE2"):
            return "imidazole N"
        return "N donor"
    if el == "O":
        if res in ("ASP", "ASH") and name.startswith("OD"):
            return "carboxylate O"
        if res in ("GLU", "GLH") and name.startswith("OE"):
            return "carboxylate O"
        if name in ("O", "OXT") or name.startswith("OXT"):
            return "backbone carbonyl O"
        if name in ("OG", "OG1", "OH"):
            return "hydroxyl O"
        return "O donor"
    if el == "S":
        return "thiol/thiolate S" if name == "SG" else "S donor"
    return f"metal ({el})" if el in METAL_ELEMENTS else f"{el} ({res})"


def metal_geometry(source: Any, metal_element: str = "FE",
                   cutoff: float = 3.0) -> dict[str, Any]:
    """Analyze metal sites (PDB or mmCIF): N/O/S donors only count for
    coordination; pairwise ligand-metal-ligand angles are fitted to ideal
    geometries; bidentate carboxylate chelates are detected and flagged."""
    p = _as_parsed(source)
    atoms, _ = _single_model(p)
    atoms = _first_altloc(atoms)
    target = metal_element.strip().upper()
    metals = [a for a in atoms if a["element"].strip().upper() == target]
    # spatial grid: cell size = cutoff -> all atoms within cutoff of a metal
    # lie in the 27 neighboring cells (avoids O(n_metals * n_atoms) scans)
    cs = max(cutoff, 0.1)
    grid: dict[tuple, list] = defaultdict(list)
    for a in atoms:
        grid[(int(math.floor(a["xyz"][0] / cs)),
              int(math.floor(a["xyz"][1] / cs)),
              int(math.floor(a["xyz"][2] / cs)))].append(a)
    sites = []
    for m in metals:
        mc = tuple(int(math.floor(m["xyz"][i] / cs)) for i in range(3))
        near = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    near.extend(grid.get((mc[0] + dx, mc[1] + dy, mc[2] + dz), ()))
        donors, others = [], 0
        for a in near:
            if a is m:
                continue
            d = math.dist(m["xyz"], a["xyz"])
            if d > cutoff:
                continue
            el = a["element"].strip().upper()
            c: dict[str, Any] = {
                "residue": f'{a["resname"]} {a["resseq"]}{a["icode"]}{a["chain"]}'.strip(),
                "atom": a["name"],
                "dist_A": round(d, 2),
                # store vector RELATIVE to the metal for angle computation
                "_xyz": tuple(a["xyz"][i] - m["xyz"][i] for i in range(3)),
            }
            if el in ("N", "O", "S"):
                c["type"] = _classify_ligand(a)
                donors.append(c)
            elif el == target:
                c["type"] = f"adjacent metal ({el})"
                donors.append(c)
            else:
                others += 1
        donors.sort(key=lambda c: c["dist_A"])
        cn = len(donors)

        by_res: dict[str, list] = defaultdict(list)
        for c in donors:
            by_res[c["residue"]].append(c)
        bidentate = [reslab for reslab, cs in by_res.items()
                     if len(cs) == 2
                     and all(c["atom"].startswith(("OD", "OE", "OQ")) for c in cs)
                     and any(c.get("type") == "carboxylate O" for c in cs)]

        angles = sorted(round(_angle_deg(c1["_xyz"], c2["_xyz"]), 1)
                        for c1, c2 in combinations(donors, 2))
        fit = {name_: round(math.sqrt(sum((o - i) ** 2 for o, i in zip(angles, ideal))
                                      / len(ideal)), 1)
               for name_, ideal in _ideal_geometries(cn).items()}
        best = min(fit, key=fit.get) if fit else None

        sites.append({
            "site": f'{m["resname"]} {m["resseq"]}{m["icode"]}{m["chain"]}'.strip(),
            "coordination_number": cn,
            "geometry_guess": GEOMETRY_GUESS.get(cn, "irregular / verify manually"),
            "geometry_fit": ({"best": best, "rmsd_deg": fit[best],
                              "candidates_rmsd_deg": fit} if fit else None),
            "bidentate_chelates": bidentate or None,
            "n_adjacent_metal": sum(1 for c in donors
                                    if str(c.get("type", "")).startswith("adjacent metal")),
            "ligand_contacts": [{k: v for k, v in c.items() if k != "_xyz"}
                                for c in donors],
            "other_atoms_within_cutoff": others,
        })
    return {"metal": target, "cutoff_A": cutoff, "n_sites": len(sites),
            "metal_sites": sites}


def structure_stats(source: Any) -> dict[str, Any]:
    """General PDB/mmCIF summary: atoms, residues, chains, gaps, altLocs.

    Accepts PDB text or a parsed structure dict (same as sequences()).
    Chain gaps reuse the validate_pdb gap logic (imported, not duplicated).
    """
    p = _as_parsed(source)
    atoms, n_models = _single_model(p)
    atoms = _first_altloc(atoms)
    v = validate_pdb(source) if isinstance(source, str) else None
    # residue classes
    classes = {"standard": 0, "variant_md_names": 0, "solvent": 0,
               "hetero": 0, "unknown": 0}
    seen_res: set = set()
    chains: dict[str, int] = {}
    ca_res, n_res = 0, 0
    altlocs = sum(1 for a in p["atoms"] if a.get("altloc"))
    icodes = 0
    for a in atoms:
        k = (a["chain"], a["resseq"], a.get("icode", ""))
        if a.get("icode"):
            icodes += 1
        if k in seen_res:
            continue
        seen_res.add(k)
        n_res += 1
        chains[a["chain"] or "_"] = chains.get(a["chain"] or "_", 0) + 1
        rn = a["resname"]
        if rn in STANDARD_RESIDUES:
            classes["standard"] += 1
        elif rn in VARIANT_RESIDUES:
            classes["variant_md_names"] += 1
        elif rn in SOLVENT:
            classes["solvent"] += 1
        elif rn in THREE_TO_ONE or a["record"] == "ATOM":
            classes["unknown"] += 1
        else:
            classes["hetero"] += 1
    for a in atoms:
        if a["name"] == "CA" and a["resname"] in STANDARD_RESIDUES:
            ca_res += 1
    bf = [a["bfactor"] for a in atoms if a["bfactor"] is not None]
    gaps = list(v.get("chain_gaps", [])) if v else []
    out: dict[str, Any] = {
        "n_atoms": len(atoms),
        "n_residues": n_res,
        "chains": {ch: {"length": ln} for ch, ln in chains.items()},
        "n_models": n_models,
        "residue_classes": classes,
        "ca_completeness_pct": round(100 * ca_res / n_res, 1) if n_res else 0.0,
        "chain_gaps": gaps,
        "n_altlocs": altlocs,
        "n_insertion_codes": len({k for k in seen_res if k[2]}),
        "bfactor": {"mean": round(sum(bf) / len(bf), 2),
                    "min": round(min(bf), 2),
                    "max": round(max(bf), 2)} if bf else None,
    }
    if n_models:
        out["note"] = (f"multi-model file ({n_models} MODEL records); "
                       "statistics from model 1 only")
    if not chains or all(ch == "_" for ch in chains):
        out["note"] = (out.get("note", "") + "; no chain IDs in file").strip("; ")
    return out


def sequences(source: Any) -> dict[str, Any]:
    p = _as_parsed(source)
    atoms, n_models = _single_model(p)
    atoms = _first_altloc(atoms)
    seen: dict[tuple, str] = {}
    order: dict[str, list] = {}
    for a in atoms:
        if a["resname"] not in THREE_TO_ONE:
            continue
        ch, k = chain_key(a), res_key(a)
        if (ch, k) not in seen:
            seen[(ch, k)] = THREE_TO_ONE[a["resname"]]
            order.setdefault(ch, []).append(k)
    out: dict[str, Any] = {}
    for ch, keys in order.items():
        keys.sort(key=lambda k: (int(k[1]) if k[1].lstrip("-").isdigit() else 10 ** 9, k[2]))
        seq = "".join(seen[(ch, r)] for r in keys)
        out[ch] = {"length": len(seq), "sequence": seq}
    if not out:
        nuc = {"DA", "DT", "DG", "DC", "DI", "A", "U", "G", "C", "T", "PSU"}
        found_nuc = sorted({a["resname"] for a in atoms} & nuc)
        out["_note"] = ("no standard amino-acid residues found in this file"
                        + (f"; nucleic-acid residues present: {', '.join(found_nuc[:8])}"
                           if found_nuc else
                           "; no nucleic-acid residues detected either - "
                           "run validate_pdb to see what the file contains"))
    if n_models:
        out["_note"] = (f"multi-model file ({n_models} MODEL records); "
                        "sequence extracted from model 1 only")
    return out


def plddt_summary(source: Any) -> dict[str, Any]:
    p = _as_parsed(source)
    atoms, _ = _single_model(p)
    atoms = _first_altloc(atoms)
    bf = [a["bfactor"] for a in atoms if a["bfactor"] is not None]
    if not bf:
        return {"is_alphafold_like": False,
                "note": "no B-factor values present; not an AlphaFold model?"}
    if not all(0 <= v <= 100 for v in bf):
        return {"is_alphafold_like": False,
                "mean_bfactor": round(sum(bf) / len(bf), 2),
                "note": ("B-factors outside 0-100: this is an experimental "
                         "structure (crystallographic/EM B-factors or "
                         "occupancies), not an AlphaFold model. pLDDT bands "
                         "are not applicable; run this tool on an AF-*.pdb.")}

    classes = {"very_high_gt90": 0, "confident_70_90": 0,
               "low_50_70": 0, "very_low_lt50": 0}
    for v in bf:
        if v > 90:
            classes["very_high_gt90"] += 1
        elif v >= 70:
            classes["confident_70_90"] += 1
        elif v >= 50:
            classes["low_50_70"] += 1
        else:
            classes["very_low_lt50"] += 1
    header = source[:4000].lower() if isinstance(source, str) else ""
    likely_af = "alphafold" in header
    n = len(bf)
    per_res: dict[tuple, list] = {}
    res_order: list = []
    for a in atoms:
        k = res_key(a)
        if k not in per_res:
            per_res[k] = []
            res_order.append(k)
        if a["bfactor"] is not None:
            per_res[k].append(a["bfactor"])
    low = [(k, sum(v) / len(v)) for k, v in per_res.items() if sum(v) / len(v) < 70]
    out = {
        "is_alphafold_like": likely_af or all(0 <= v <= 100 for v in bf),
        "alphafold_detected_in_header": likely_af,
        "mean_plddt": round(sum(bf) / n, 2),
        "n_atoms_scored": n,
        "bands_pct": {k: round(100 * v / n, 1) for k, v in classes.items()},
        "low_confidence_regions": [
            {"chain": ch, "resseq": rs + ic, "mean_plddt": round(v, 1)}
            for (ch, rs, ic), v in low[:50]
        ],
    }
    if not likely_af:
        out["caution"] = ("No AlphaFold marker found in the file header. If this is "
                          "an experimental structure, these B-factors are "
                          "crystallographic/EM thermal factors, NOT pLDDT confidence, "
                          "and the bands above must not be read as confidence.")
    return out

