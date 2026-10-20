import sys
import pathlib

FIX = pathlib.Path(__file__).resolve().parent / "fixtures"
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pdb_utils

# Synthetic PDB: chain A res 1-3 + insertion-code variant of res 3
# (same residue number, still length 3), chain B res 10-11,
# 1 altLoc atom, 1 insertion code.
def atom(serial, name, res, chain, seq, icode, x, alt=""):
    head = (f"ATOM  {serial:>5d} {name:^4s}{alt:1s}{res:>3s} {chain}{seq:>4d}"
            f"{icode:1s}   {x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 20.00")
    return head.ljust(76) + "  C"

lines = ["HEADER    SYNTH"]
n = 0
for seq in (1, 2, 3):
    n += 1
    lines.append(atom(n, "CA", "ALA", "A", seq, "", float(seq)))
n += 1
lines.append(atom(n, "CB", "ALA", "A", 2, "", 9.0, alt="B"))  # altLoc dup (diff atom)
n += 1
lines.append(atom(n, "CA", "GLY", "A", 3, "A", 3.0))  # insertion code on res 3
for seq in (10, 11):
    n += 1
    lines.append(atom(n, "CA", "SER", "B", seq, "", float(seq)))
# Chain B jumps 11 -> 14: a real 2-residue gap for the gap detector.
n += 1
lines.append(atom(n, "CA", "SER", "B", 14, "", 14.0))
lines.append("END")
text = "\n".join(lines) + "\n"

st = pdb_utils.structure_stats(text)
assert st["chains"] == {"A": {"length": 4}, "B": {"length": 3}}, st["chains"]
print("1. structure_stats chains: OK", st["chains"],
      "(3 seq-numbers + 1 insertion-code variant)")
assert len(st["chain_gaps"]) >= 1, st["chain_gaps"]
print("2. structure_stats gaps: OK", st["chain_gaps"])
assert st["n_altlocs"] >= 1, st["n_altlocs"]
print("3. structure_stats altlocs: OK", st["n_altlocs"])
assert st["n_insertion_codes"] == 1, st["n_insertion_codes"]
print("4. structure_stats icodes: OK")
assert st["ca_completeness_pct"] == 100.0, st["ca_completeness_pct"]
print("5. structure_stats CA completeness: OK")
assert st["residue_classes"]["standard"] == 7, st["residue_classes"]
print("6. structure_stats residue classes: OK")

# Fixture-file smoke check (130 aa chain + 2 FE + 8 donor singletons)
real = FIX / "A8ABE7_holo_Fe2.pdb"
if real.exists():
    r = pdb_utils.structure_stats(real.read_text())
    assert r["n_atoms"] == 530 and r["n_residues"] == 140, (r["n_atoms"], r["n_residues"])
    print("7. structure_stats fixture file: OK", r["n_atoms"], r["n_residues"])
else:
    print("7. structure_stats fixture file: SKIP")

print()
print("ALL structure_stats TESTS PASSED")