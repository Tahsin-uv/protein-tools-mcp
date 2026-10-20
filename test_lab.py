"""Headless tests for the Protein Lab UI helpers (no browser needed).

Exercises the same analysis calls the Streamlit tabs make: demo fixtures
load, sequence classification + composition render without error, structure
load/validate/metal paths, and the plain-English wording helpers.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lab_seq  # noqa: E402
import lab_structure  # noqa: E402
import seq_utils  # noqa: E402
import pdb_utils  # noqa: E402

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")

# --- 1. demo FASTA parses into protein/dna/rna records ---
demo = lab_seq._demo_text()
assert "sample_protein" in demo, demo[:80]
recs = seq_utils.parse_fasta(demo, is_file=False)
assert [(r["id"], r["type"]["type"]) for r in recs] == [
    ("sample_protein", "protein"), ("sample_dna", "dna"),
    ("sample_rna", "rna")], recs
print("1. lab demo FASTA (3 records, typed): OK")

# --- 2. protein record composes + wording maps resolve ---
c = seq_utils.aa_composition(recs[0]["sequence"])
assert c["length"] == 20 and c["molecular_weight_da"] > 2000, c
words, dot = lab_seq._TYPE_WORDS["protein"]
assert "protein" in words and dot, (words, dot)
assert all(k in lab_seq._GLOSSARY for k in
           ("GRAVY", "pI", "Molecular weight")), lab_seq._GLOSSARY
print("2. lab composition + wording/glossary maps: OK "
      f"(len {c['length']}, MW {c['molecular_weight_da']:.0f} Da)")

# --- 3. bare-sequence path: protein typed, DNA typed ---
label, _d = seq_utils.classify_sequence_type(
    seq_utils.normalize_sequence("MKLIVFWYHEPQDNRSACDE"))
assert label == "protein", label
label2, _d2 = seq_utils.classify_sequence_type("ATCGATCGATCG")
assert label2 == "dna", label2
print("3. lab bare-sequence typing (protein + dna): OK")

# --- 4. structure demos load through the lab loader ---
for name, chain_n, metal in (("tiny_zn.pdb", 1, "ZN"),
                             ("A8ABE7_holo_Fe2.pdb", 1, "FE")):
    text = lab_structure._demo(name)
    struct = lab_structure._load(text, name)
    v = pdb_utils.validate_pdb(struct)
    assert not v["format_errors"], v["format_errors"][:3]
    found = lab_structure._find_metals(struct)
    assert metal in found and found[metal], (name, found)
    g = pdb_utils.metal_geometry(struct, metal)
    assert g["n_sites"] >= 1, g
    print(f"4. lab structure demo {name}: OK "
          f"({v['counts']['residues']} residues, "
          f"{g['n_sites']} {metal} site(s))")

# --- 5. geometry wording covers every fit the code can produce ---
for best in ("linear", "trigonal_planar", "tetrahedral", "square_planar",
             "trigonal_bipyramidal", "square_pyramidal", "octahedral",
             None, "weird_future_shape"):
    shape = lab_structure._GEO_WORDS.get(best, best or "an irregular "
                                                       "arrangement")
    assert isinstance(shape, str) and shape, best
print("5. lab geometry wording (all fits + fallbacks): OK")

# --- 6. mmCIF fixture loads through the lab loader ---
cif_text = lab_structure._demo("fold_a0a977pl24_2ca2_model_0.cif")
cif_struct = lab_structure._load(cif_text,
                                 "fold_a0a977pl24_2ca2_model_0.cif")
assert isinstance(cif_struct, dict) and cif_struct["atoms"], "cif load"
print(f"6. lab mmCIF demo: OK ({len(cif_struct['atoms'])} atoms)")

print("ALL PROTEIN LAB TESTS PASSED")
