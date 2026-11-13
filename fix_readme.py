import pathlib

d = pathlib.Path(r"e:\Protein characterization\mcp-protein-tools\dist\protein-tools-public\README.md")
c = d.read_text()
lines = c.split("\n")

# The 4 correct tool rows (with proper backticks)
correct_lines = [
    "| `detect_sequence_type` | free | Detect whether one or more FASTA sequences are protein, DNA, RNA, ambiguous, or unknown from alphabet composition (T-vs-U logic), with confidence notes instead of silent guesses |",
    "| `read_fasta` | free | Parse one or more FASTA records (file path or raw text) into structured dicts with id, description, sequence, length, type classification, and a UniProt accession hint from the header |",
    "| `pdb_by_id` | premium (local ok) | Get info about a PDB entry by its 4-character ID: checks for a local PDB/mmCIF file and if found runs validation, sequence extraction, and metal-site analysis; also fetches RCSB metadata (requires API key for full data, gracefully degrades to local-only without) |",
    "| `find_pdb_by_sequence` | premium | Find PDB structures matching a protein sequence from a FASTA file or raw text via RCSB BLAST search. Parses the FASTA, verifies it is a protein sequence, then searches RCSB and returns matching PDB IDs with metadata. Requires API key. |",
]

print("Current lines 122-125 (0-indexed: 121-124):")
for i in range(121, 126):
    if i < len(lines):
        print(f"  {i+1}: {repr(lines[i])}")

# The corrupted lines are at indices 121, 122, 123, 124
lines[121:125] = correct_lines

print()
print("Fixed lines 122-125:")
for i in range(121, 126):
    if i < len(lines):
        print(f"  {i+1}: {lines[i]}")

d.write_text("\n".join(lines))
print()
print("README.md updated successfully!")
