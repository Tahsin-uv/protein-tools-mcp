"""Unit tests for correctness fixes and new analysis tools (no MCP needed)."""
import sys

import os

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pdb_utils  # noqa: E402
import xvg_tools  # noqa: E402

# --- 1. insertion codes: 100, 100A, 101 must be 3 distinct ordered residues ---
ICODE = """HEADER    test
ATOM      1  CA  ALA A 100      11.111  12.222  13.333  1.00 20.00           C
ATOM      2  CA  ALA A 100A     14.444  15.555  16.666  1.00 20.00           C
ATOM      3  CA  ALA A 101      17.777  18.888  19.999  1.00 20.00           C
END
"""
v = pdb_utils.validate_pdb(ICODE)
seqs = pdb_utils.sequences(ICODE)
assert v["counts"]["residues"] == 3, v["counts"]
assert v["insertion_code_residues"] == ["A100A"], v["insertion_code_residues"]
assert seqs["A"]["sequence"] == "AAA" and seqs["A"]["length"] == 3, seqs
print("1. insertion-code handling: OK")

# --- 2. multi-model (NMR-style): analysis restricted to model 1 ---
MODEL = """HEADER    test
MODEL        1
ATOM      1  CA  ALA A   1      11.111  12.222  13.333  1.00 20.00           C
ATOM      2  CA  GLY A   2      14.444  15.555  16.666  1.00 20.00           C
ENDMDL
MODEL        2
ATOM      3  CA  ALA A   1      21.111  22.222  23.333  1.00 20.00           C
ATOM      4  CA  GLY A   2      24.444  25.555  26.666  1.00 20.00           C
ATOM      5  CA  VAL A   3      27.777  28.888  29.999  1.00 20.00           C
ENDMDL
END
"""
v = pdb_utils.validate_pdb(MODEL)
assert v["counts"]["residues"] == 2, v["counts"]
assert "multi_model_note" in v, v
seqs = pdb_utils.sequences(MODEL)
assert seqs["A"]["sequence"] == "AG" and seqs.get("_note"), seqs
print("2. multi-model handling: OK")

# --- 3. metal geometry: angles + tetrahedral fit on the holo fixture ---
t = open(os.path.join(FIX, "A8ABE7_holo_Fe2.pdb"), errors="replace").read()
g = pdb_utils.metal_geometry(t, "FE")
assert g["n_sites"] == 2, g["n_sites"]
s = g["metal_sites"][0]
n = s["coordination_number"]
assert s["geometry_fit"] and s["geometry_fit"]["best"], s
assert len(s["ligand_contacts"]) == n
print(f"3. metal geometry: OK -> site {s['site']}: CN {n}, "
      f"fit '{s['geometry_fit']['best']}' (rmsd {s['geometry_fit']['rmsd_deg']} deg), "
      f"bidentate: {s['bidentate_chelates']}, "
      f"types: {[c['type'] for c in s['ligand_contacts']]}")

# --- 4. xvg tools on the synthetic GROMACS fixtures ---
rmsd_txt = open(os.path.join(FIX, "rmsd_core.xvg"), errors="replace").read()
r = xvg_tools.rmsd_summary(rmsd_txt)
print(f"4a. rmsd_summary: verdict = {r['verdict']} | plateau_reached = "
      f"{r['plateau_reached']} at t={r['plateau_time']} | drift = {r['drift_pct']}%")
assert "verdict" in r and "drift_pct" in r

sasa_txt = open(os.path.join(FIX, "sasa_per_residue.xvg"), errors="replace").read()
s = xvg_tools.summary(sasa_txt)
print(f"4b. xvg_summary: '{s['title']}' | columns: {list(s['columns'])} | "
      f"points: {s['n_points']}")

p = xvg_tools.rmsf_peaks(sasa_txt, top_n=5)
print(f"4c. rmsf_peaks (on SASA series): top values {[q['rmsf'] for q in p['most_flexible']]}")
assert p["most_flexible"] and p["flexible_segments"] is not None

# --- 5. error handling: malformed xvg must raise, not crash ---
try:
    xvg_tools.read_xvg("hello world\nnot numbers\n")
    raise AssertionError("expected ValueError for non-numeric xvg")
except ValueError:
    print("5. malformed xvg raises ValueError: OK")

# --- 6. altLoc: analyses keep one conformer per atom name ---
ALT = """HEADER    test
HETATM    1  ZN  ZN  A   1       0.000   0.000   0.000  1.00 20.00          ZN
ATOM      2  OD1AASP A   5       2.000   0.000   0.000  0.50 20.00           O
ATOM      3  OD1BASP A   5       2.600   0.000   0.000  0.50 20.00           O
ATOM      4  NE2AHIS A   9       0.000   2.100   0.000  1.00 20.00           N
END
"""
v = pdb_utils.validate_pdb(ALT)
assert v["altloc_present"] == ["A", "B"], v["altloc_present"]
g = pdb_utils.metal_geometry(ALT, "ZN")
s = g["metal_sites"][0]
assert s["coordination_number"] == 2, s
assert [c["atom"] for c in s["ligand_contacts"]] == ["OD1", "NE2"], s["ligand_contacts"]
print("6. altLoc handling: OK -> CN 2 (one conformer kept per atom name), "
      "validate still reports both altLocs")

# --- 7. DNA-only file: sequences() explains itself instead of returning {} ---
DNA = """HEADER    test dna only
ATOM      1  P   DG  A   1      11.111  12.222  13.333  1.00 20.00           P
ATOM      2  P   DC  A   2      14.444  15.555  16.666  1.00 20.00           P
END
"""
s = pdb_utils.sequences(DNA)
assert "_note" in s and "no standard amino-acid" in s["_note"], s
assert "DG" in s["_note"], s
print(f"7. DNA-only explanation: OK -> {s['_note']}")

print("ALL UNIT TESTS PASSED")
