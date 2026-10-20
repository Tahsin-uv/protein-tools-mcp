import sys, pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import seq_utils

# Test classify_sequence - protein (needs unambiguous protein-only chars:
# E, F, H, I, L, M, P, Q, W — the T-vs-U rule otherwise claims T/U sequences)
prot = "MKTLLVAGFDEFHIKLMNPQRSWY"
r = seq_utils.classify_sequence(prot)
assert r["type"] == "protein", f"Expected protein, got {r['type']}"
print("1. classify_sequence protein: OK")

# Test classify_sequence - DNA
dna = "ATCGATCG"
r = seq_utils.classify_sequence(dna)
assert r["type"] == "dna", f"Expected dna, got {r['type']}"
print("2. classify_sequence DNA: OK")

# Test classify_sequence - RNA
rna = "AUGCAUGC"
r = seq_utils.classify_sequence(rna)
assert r["type"] == "rna", f"Expected rna, got {r['type']}"
print("3. classify_sequence RNA: OK")

# Test classify_sequence - ambiguous (both T and U)
amb = "ATCGU"
r = seq_utils.classify_sequence(amb)
assert r["type"] == "ambiguous", f"Expected ambiguous, got {r['type']}"
print("4. classify_sequence ambiguous: OK")

# Test classify_sequence - empty
r = seq_utils.classify_sequence("")
assert r["type"] == "unknown", f"Expected unknown, got {r['type']}"
print("5. classify_sequence empty: OK")

# Test parse_fasta - simple (unambiguous protein chars so the
# T-vs-U rule does not claim T/U-bearing sequences)
fasta = ">test_prot\nMKLIVFWYHEPQDNRS\n"
records = seq_utils.parse_fasta(fasta, is_file=False)
assert len(records) == 1
assert records[0]["id"] == "test_prot"
assert records[0]["sequence"] == "MKLIVFWYHEPQDNRS"
print("6. parse_fasta single: OK")

# Test parse_fasta - multi
fasta_multi = ">seq1\nMKLIVFWYHEPQDNRS\n>seq2\nATCGATCG\n"
records = seq_utils.parse_fasta(fasta_multi, is_file=False)
assert len(records) == 2
assert records[0]["id"] == "seq1"
assert records[1]["id"] == "seq2"
print("7. parse_fasta multi: OK")

# Test normalize_sequence
raw = "M K T L L V - A G F L K G *"
norm = seq_utils.normalize_sequence(raw)
assert norm == "MKTLLVAGFLKG", f"Expected MKTLLVAGFLKG, got {norm}"
print("8. normalize_sequence: OK")

# Test normalize_sequence - convert to RNA (T->U only; G/C stay)
dna_seq = "ATCGATCG"
rn = seq_utils.normalize_sequence(dna_seq, target_type="rna")
assert rn == "AUCGAUCG", f"Expected AUCGAUCG, got {rn}"
print("9. normalize_sequence to RNA: OK")

# Test _hint_uniprot_accession
assert seq_utils._hint_uniprot_accession("sp|P00533|EGFR_HUMAN") == "P00533"
assert seq_utils._hint_uniprot_accession("A0A001|NAME") == "A0A001"
assert seq_utils._hint_uniprot_accession("just_a_name") is None
print("10. _hint_uniprot_accession: OK")

print()
print("ALL seq_utils TESTS PASSED")
