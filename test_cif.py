"""mmCIF support tests: synthetic offline CIF + a synthetic Ca2+ fixture."""
import os
import sys

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cif_tools  # noqa: E402
import pdb_utils  # noqa: E402

# --- 1. synthetic mmCIF: loops, quoted key-value, missing markers ---
SYN = """data_test
_struct.title "synthetic test cif"
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.B_iso_or_equiv
ATOM N N . GLY A 1 ? 1.000 1.000 1.000 20.0
ATOM N N . ALA A 2 ? 4.000 1.000 1.000 20.0
HETATM ZN ZN . ZN A 100 ? 2.500 1.000 1.000 20.0
#
"""
parsed = cif_tools.parse_cif(SYN)
assert parsed["block"] == "test", parsed["block"]
assert parsed["categories"]["struct"]["rows"][0]["title"] == "synthetic test cif"
a = cif_tools.to_pdb_atoms(parsed)
assert len(a["atoms"]) == 3, a["atoms"]
v = pdb_utils.validate_pdb(a)
assert v["counts"]["residues"] == 2 and v["counts"]["atoms"] == 3, v["counts"]
s = pdb_utils.sequences(a)
assert s["A"]["sequence"] == "GA", s
g = pdb_utils.metal_geometry(a, "ZN")
site = g["metal_sites"][0]
assert site["coordination_number"] == 2, site
assert site["geometry_fit"]["best"] == "linear" and site["geometry_fit"]["rmsd_deg"] == 0.0, site
print("1. synthetic mmCIF parse + analyses: OK (ZN linear, rmsd 0.0 deg)")

# --- 2. synthetic Ca2+ mmCIF fixture (80 aa + 2 CA ions + pLDDT B-factors) ---
CIF = os.path.join(FIX, "fold_a0a977pl24_2ca2_model_0.cif")
parsed = cif_tools.parse_cif(open(CIF, errors="replace").read())
a = cif_tools.to_pdb_atoms(parsed)
assert a["atoms"], "no atoms parsed from the fixture cif"
v = pdb_utils.validate_pdb(a)
print("2. fixture cif:", v["counts"])
s = pdb_utils.sequences(a)
seqs = [d["sequence"] for k, d in s.items() if not k.startswith("_")]
assert seqs and all(len(x) == 80 for x in seqs), s
print("   sequence:", seqs[0][:60], f"... ({len(seqs[0])} aa)")
g = pdb_utils.metal_geometry(a, "CA")
print(f"   Ca sites: {g['n_sites']}")
assert g["n_sites"] == 2, g["n_sites"]
if g["metal_sites"]:
    ms = g["metal_sites"][0]
    fit = ms["geometry_fit"]
    print(f"   {ms['site']}: CN {ms['coordination_number']}, "
          f"fit '{fit['best']}' (rmsd {fit['rmsd_deg']} deg), "
          f"types {[c['type'] for c in ms['ligand_contacts']]}")
p = pdb_utils.plddt_summary(a)
print(f"   plddt mean: {p.get('mean_plddt')} (caution present: {'caution' in p})")
print("ALL CIF TESTS PASSED")
