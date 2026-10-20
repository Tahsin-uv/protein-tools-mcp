import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

# Payload-mapping for the RCSB sequence search (no network): replicate the
# seq_type_map from server.find_pdb_by_sequence.
seq_type_map = {"protein": ("protein", "blastp"),
                "dna": ("dna", "blastn"),
                "rna": ("rna", "blastn")}

import seq_utils

# 1. protein -> blastp
rec = seq_utils.parse_fasta(">p\nMKLIVFWYHEPQDNRSACDE\n", is_file=False)[0]
assert rec["type"]["type"] == "protein", rec["type"]
assert seq_type_map[rec["type"]["type"]] == ("protein", "blastp")
print("1. blast_fallback protein -> blastp: OK")

# 2. dna -> blastn
rec = seq_utils.parse_fasta(">d\nATCGATCGATCGATCG\n", is_file=False)[0]
assert rec["type"]["type"] == "dna", rec["type"]
assert seq_type_map[rec["type"]["type"]] == ("dna", "blastn")
print("2. blast_fallback dna -> blastn: OK")

# 3. rna -> blastn
rec = seq_utils.parse_fasta(">r\nAUGCAUGCAUGCAUGC\n", is_file=False)[0]
assert rec["type"]["type"] == "rna", rec["type"]
assert seq_type_map[rec["type"]["type"]] == ("rna", "blastn")
print("3. blast_fallback rna -> blastn: OK")

# 4. ambiguous rejected
import json
import server as srv
out = json.loads(srv.find_pdb_by_sequence(">a\nATCGU\n", is_file=False))
assert "error" in out and "ambiguous" in out["error"].lower(), out
print("4. blast_fallback ambiguous rejection: OK")

print()
print("ALL blast_fallback TESTS PASSED")