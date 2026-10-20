import sys, pathlib

FIX = pathlib.Path(__file__).resolve().parent / "fixtures"
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import seq_utils

# Test read_fasta - single sequence (unambiguous protein chars so the
# T-vs-U rule does not claim it)
fasta = ">test_prot\nMKLIVFWYHEPQDNRS\n"
records = seq_utils.parse_fasta(fasta, is_file=False)
assert len(records) == 1
assert records[0]["id"] == "test_prot"
assert records[0]["sequence"] == "MKLIVFWYHEPQDNRS"
assert records[0]["type"]["type"] == "protein"
assert records[0]["length"] == 16
assert records[0]["maybe_uniprot_accession"] is None
print("1. read_fasta single: OK")

# Test read_fasta - multi sequence
fasta_multi = ">seq1 protein\nMKLIVFWYHEPQDNRS\n>seq2 dna\nATCGATCG\n>seq3 rna\nAUGCAUGC\n"
records = seq_utils.parse_fasta(fasta_multi, is_file=False)
assert len(records) == 3
assert records[0]["id"] == "seq1"
assert records[1]["id"] == "seq2"
assert records[2]["id"] == "seq3"
assert records[0]["type"]["type"] == "protein"
assert records[1]["type"]["type"] == "dna"
assert records[2]["type"]["type"] == "rna"
print("2. read_fasta multi: OK")

# Test read_fasta - from the fixture file
fasta_path = FIX / "sample.fasta"
import os
if os.path.exists(fasta_path):
    records = seq_utils.parse_fasta(fasta_path, is_file=True)
    print(f"3. read_fasta file: OK ({len(records)} records)")
else:
    print("3. read_fasta file: SKIP (file not found)")

# Test read_fasta - empty file handling
try:
    seq_utils.parse_fasta("nonexistent.fasta", is_file=True)
    assert False, "Should have raised FileNotFoundError"
except FileNotFoundError:
    print("4. read_fasta FileNotFoundError: OK")

# Test read_fasta - no records: text with no FASTA header at all
try:
    seq_utils.parse_fasta("MKTLLVAGFLKG", is_file=False)
    assert False, "Should have raised ValueError"
except ValueError as e:
    assert "No FASTA records" in str(e)
    print("5. read_fasta no records ValueError: OK")

# Test normalize_sequence used within read_fasta
seq = "M K T L L V A G F L K G"
norm = seq_utils.normalize_sequence(seq)
assert norm == "MKTLLVAGFLKG"
print("6. read_fasta normalize: OK")

print()
print("ALL read_fasta TESTS PASSED")
