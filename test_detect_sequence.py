import sys, pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import seq_utils

# Test detect_sequence_type via seq_utils (the underlying library)
# This tests the core logic that the MCP tool uses

# Test 1: detect protein — use a sequence with unambiguous protein-only chars
prot_fasta = ">test_prot\nMKTLLVAGFDEFHIKLMNPQRSWY\n"
records = seq_utils.parse_fasta(prot_fasta, is_file=False)
assert len(records) == 1
t = records[0]["type"]
assert t["type"] == "protein", f"Expected protein, got {t['type']}"
print("1. detect_sequence_type protein: OK")

# Test 2: detect DNA
dna_fasta = ">test_dna\nATCGATCG\n"
records = seq_utils.parse_fasta(dna_fasta, is_file=False)
assert len(records) == 1
t = records[0]["type"]
assert t["type"] == "dna", f"Expected dna, got {t['type']}"
print("2. detect_sequence_type DNA: OK")

# Test 3: detect RNA
rna_fasta = ">test_rna\nAUGCAUGC\n"
records = seq_utils.parse_fasta(rna_fasta, is_file=False)
assert len(records) == 1
t = records[0]["type"]
assert t["type"] == "rna", f"Expected rna, got {t['type']}"
print("3. detect_sequence_type RNA: OK")

# Test 4: detect ambiguous
amb_fasta = ">test_amb\nATCGU\n"
records = seq_utils.parse_fasta(amb_fasta, is_file=False)
assert len(records) == 1
t = records[0]["type"]
assert t["type"] == "ambiguous", f"Expected ambiguous, got {t['type']}"
print("4. detect_sequence_type ambiguous: OK")

# Test 5: detect unknown (empty)
empty_fasta = ">test\n\n"
records = seq_utils.parse_fasta(empty_fasta, is_file=False)
assert len(records) == 1
t = records[0]["type"]
assert t["type"] == "unknown", f"Expected unknown, got {t['type']}"
print("5. detect_sequence_type empty: OK")

# Test 6: detect with confidence
r = seq_utils.classify_sequence("ATCGATCG")
assert r["detail"]["confidence"] == "high", f"Expected high confidence, got {r['detail']['confidence']}"
print("6. detect_sequence_type confidence high: OK")

print()
print("ALL detect_sequence_type TESTS PASSED")
