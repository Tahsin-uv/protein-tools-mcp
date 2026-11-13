The plan was to complete these 5 items:

1. ✅ **Add `seq_utils.py` functions** — Done: parse_fasta, normalize_sequence, classify_sequence (and helper functions) are all in place in both `mcp-protein-tools/seq_utils.py` and `dist/protein-tools-public/seq_utils.py`.

2. ✅ **Add 4 new tools to `server.py`** — Done: `detect_sequence_type`, `read_fasta`, `pdb_by_id`, `find_pdb_by_sequence` are all added to `server.py` with `@mcp.tool()` decorators. The main `server.py` has 27 `@mcp.tool()` decorators (up from 23 with the original 23 tools + 4 new). The public snapshot `dist/protein-tools-public/server.py` also has all 4 tools (26997 chars, imports seq_utils, all 4 defs present).

3. ✅ **Sync public snapshot** — Done: `dist/protein-tools-public/` has the updated `seq_utils.py` and `server.py` with all new code.

4. ✅ **Create test files** — Done: All 4 test files created in both main and dist:
   - `test_seq_utils.py` (2464 bytes) - tests seq_utils parse_fasta, normalize_sequence, classify_sequence
   - `test_detect_sequence.py` (1975 bytes) - tests detect_sequence_type tool
   - `test_read_fasta.py` (2125 bytes) - tests read_fasta tool  
   - `test_pdb_by_id.py` (4605 bytes) - tests pdb_by_id tool

5. ⚠️ **Update README.md tool table** — **INCOMPLETE**: The public snapshot `dist/protein-tools-public/README.md` tool table still needs the 4 new tools added. The main `README.md` already has all 4 tools in its table (detected by grep), but the public snapshot README.md does NOT have them yet.

**Next step needed:** Add the 4 new tools to the public snapshot README.md tool table. The tool table is around line 96-140 in the public README. The 4 new tools to add:
- `detect_sequence_type` | free | Detect protein/DNA/RNA/ambiguous from FASTA sequence
- `read_fasta` | free | Parse FASTA records into structured dicts
- `pdb_by_id` | premium (local ok) | Get PDB info by 4-char ID
- `find_pdb_by_sequence` | premium | Find PDB matches via RCSB BLAST search

This is the only remaining item to complete the task.