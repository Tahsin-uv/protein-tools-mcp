import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import seq_utils

# 1. Known peptide ACDE: MW = 71.0788+103.1388+115.0886+129.1155+18.01528
c = seq_utils.aa_composition("ACDE")
assert c["length"] == 4, c
assert abs(c["molecular_weight_da"] - 436.44) < 0.05, c["molecular_weight_da"]
print("1. aa_composition ACDE MW: OK", c["molecular_weight_da"])

# 2. Class fractions sum to 1 over standard residues
c2 = seq_utils.aa_composition("MKLIVFWYHEPQDNRS")
tot = sum(c2["class_fractions"].values())
assert abs(tot - 1.0) < 0.01, tot
print("2. aa_composition class fractions sum: OK", tot)

# 3. GRAVY present, charge + pI present
assert isinstance(c2["gravy"], float), c2["gravy"]
assert "charge_at_pH7_estimate" in c2 and "pi_estimate" in c2
print("3. aa_composition gravy/charge/pI: OK",
      c2["gravy"], c2["charge_at_pH7_estimate"], c2["pi_estimate"])

# 4. DNA rejected with detected type named
try:
    seq_utils.aa_composition("ATCGATCGATCG")
    raise AssertionError("should have raised ValueError")
except ValueError as e:
    assert "dna" in str(e).lower(), str(e)
    print("4. aa_composition DNA rejection: OK")

# 5. Empty rejected
try:
    seq_utils.aa_composition("")
    raise AssertionError("should have raised ValueError")
except ValueError:
    print("5. aa_composition empty rejection: OK")

# 6. fasta_stats mixed file: N50 hand-check (lengths 4, 8, 12 -> total 24,
# half 12 -> N50 12)
mixed = ">p1\nACDE\n>d1\nATCGATCG\n>r1\nAUGCAUGCAUGC\n"
st = seq_utils.fasta_stats(mixed, is_file=False)
assert st["n_records"] == 3, st
assert st["n50"] == 12, st["n50"]
assert st["by_type"].get("protein") == 1 and st["by_type"].get("dna") == 1, st["by_type"]
print("6. fasta_stats mixed N50: OK", st["n50"], st["by_type"])

# 7. All-protein input: GC% null with note (never 0)
st2 = seq_utils.fasta_stats(">p1\nACDE\n>p2\nMKLIVFWY\n", is_file=False)
assert st2["gc_percent_mean"] is None, st2["gc_percent_mean"]
assert "null" in st2["gc_note"].lower() or "not applicable" in st2["gc_note"].lower()
print("7. fasta_stats all-protein GC null: OK")

print()
print("ALL aa_composition + fasta_stats TESTS PASSED")