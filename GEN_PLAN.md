# General-protein expansion plan (Part 2)

Status: Plan 1 (rewording) is DONE — README.md, LISTING.md, and server.py
docstring now position the server as "general protein characterization
(+ best-in-class metal-site module)". This file plans the Part 2 work:
five new general tools that make the server general-complete.

## Design principles (binding)

1. All five tools are FREE tier (local compute or existing cached lookups).
2. Stdlib-only analysis cores (`pdb_utils.py`, `seq_utils.py`), following the
   existing pattern: pure function in the `*_utils.py` module + thin
   `@mcp.tool()` wrapper in `server.py` returning `_j(...)` / `_error_response`.
3. Inputs follow existing conventions: `source: str + is_file: bool = True`
   for FASTA inputs; `pdb_path: str` for structure inputs; UniProt accessions
   reuse the premium-gated `uniprot_entry` path.
4. Every tool gets a `test_*.py` file mirroring the existing style
   (assert + "N. name: OK" prints + final "ALL ... PASSED" line).
5. README tool table + dist public snapshot + `make_public_repo.py`
   whitelist updated with each tool (see lessons from the 4-tool rollout:
   keep dist in sync at each step, don't batch it at the end).

## Tool 1: `aa_composition` (free) — protein physicochemistry

**Where:** new function `aa_composition(seq)` in `seq_utils.py` +
`@mcp.tool() def aa_composition(source, is_file=True, sequence_index=0)`
in `server.py`.

**Input:** FASTA file/text (uses existing `parse_fasta`), single record
selected by `sequence_index`; rejects non-protein sequences with a
structured error telling the caller the detected type.

**Output:** length, molecular weight (average isotopic masses, water
subtracted for peptide bonds), amino-acid counts + percentages,
class fractions (hydrophobic / polar / acidic / basic / aromatic /
Cys), net charge at pH 7 (K/R +1, D/E -1, H +0.5 approx — documented
as estimate), GRAVY hydrophobicity (Kyte-Doolittle table), and a
note that pI is an approximation (charge-curve bisection, pKa set
cited in the note).

**Edge cases:** ambiguous letters (B/Z/J/X) counted separately and
excluded from MW; empty sequence → structured error; DNA/RNA input →
error naming the detected type.

**Test:** `test_aa_composition.py` — known peptide (e.g. "ACDE" MW check),
class fractions sum to 1, DNA rejection, empty rejection.

**Effort:** ~0.5 day.

## Tool 2: `fasta_stats` (free) — multi-FASTA summary

**Where:** new function `fasta_stats(records)` in `seq_utils.py` +
`@mcp.tool() def fasta_stats(source, is_file=True)` in `server.py`.

**Input:** FASTA file/text with 1+ records (uses existing `parse_fasta`).

**Output:** n_records, length min/max/mean/median, N50 over record
lengths, total residues, per-type counts (protein/dna/rna/ambiguous/
unknown via existing `classify_sequence`), GC% computed over
dna+ambiguous records only (with n_gc_records noted; null with a note
when no dna records — never 0-by-default), per-record table
(id, length, type, GC% for nucleic records).

**Edge cases:** single-record input (N50 = that length); all-protein
input (GC% null + note); empty/garbage input → structured error.

**Test:** `test_fasta_stats.py` — mixed 3-record file (protein+dna+rna),
N50 hand-check, all-protein GC-null case.

**Effort:** ~0.5 day.


## Tool 3: `structure_stats` (free) — general PDB/mmCIF summary

**Where:** new function `structure_stats(source)` in `pdb_utils.py` +
`@mcp.tool() def structure_stats(pdb_path)` in `server.py`.

**Input:** PDB or mmCIF path (uses existing `_load_structure`, so AF3/RCSB
files work unchanged).

**Output:** n_atoms, n_residues, chains (ids + lengths), models count
(model-1-only note when multi-model), residue-class counts
(standard / variant-MD-names / solvent / hetero / unknown),
CA-only backbone completeness (% residues with CA), chain-gap count
(reuses validate_pdb gap logic — import, don't duplicate),
altLoc count, insertion-code count, B-factor mean/min/max reported raw
(the AF refusal stays inside `plddt_summary` only).

**Edge cases:** mmCIF input; NMR multi-model (model 1 + note);
nucleic-acid structures (protein chain table empty + note, reuse the
`sequences()` nucleic note pattern).

**Test:** `test_structure_stats.py` — synthetic PDB (2 chains, 1 gap,
1 altLoc, 1 insertion code) with hand-computed expectations + the
A8ABE7 fixture for a real-file smoke check.

**Effort:** 1 day.

## Tool 4: `conservation_hint` (premium lookup reuse)

**Where:** `@mcp.tool() def conservation_hint(accession_or_id)` in
`server.py`, reusing the existing premium-gated `uniprot_entry` fetch
+ `premium_tools.uniprot_summary` (extend the summarizer, see below).

**Input:** UniProt accession (e.g. "P00533").

**Output:** accession, name, organism, length, feature table grouped by
type (Metal binding / Active site / Disulfide bond / Mutagenesis /
Natural variant), each with description + position/range, plus a
"how to use" note (positions are 1-based UniProt canonical; align your
PDB sequence first — `sequence_from_pdb` + `read_fasta`).

**Premium gating:** same `_check_premium_key` + relay/local path as
`uniprot_entry` (no new auth code — call the same helper). Without a
key the tool returns the standard premium structured error.

**Summarizer change:** extend `premium_tools.uniprot_summary` to also
return `features_by_type` for the ~6 feature types above (currently
only Metal binding is extracted). Additive — existing callers fine.

**Test:** `test_conservation_hint.py` — offline summarizer test with a
canned UniProt JSON fixture (no network) + premium-gate test.

**Effort:** ~0.5 day.

## Tool 5: extend `find_pdb_by_sequence` to DNA/RNA (premium)

**Where:** extend `find_pdb_by_sequence` in `server.py` (no new tool).

**Change:** currently hard-codes `"sequence_type": "protein"`,
`"program": "blastp"`. Map the record type from `classify_sequence`:
- protein → protein/blastp (unchanged)
- dna → dna/blastn; rna → rna/blastn
- ambiguous/unknown → structured error asking caller to clarify.

**Why extend, not new:** one tool, one docstring, one test file;
premium gate, relay path, result parsing identical — two fields vary.

**Test:** `test_blast_fallback.py` — payload-mapping unit test
(no network) + ambiguous-rejection test.

**Effort:** ~1 hour.

## Rollout order + validation

1. `aa_composition` → 2. `fasta_stats` → 3. `structure_stats` →
   4. `conservation_hint` (+ summarizer extension) → 5. DNA/RNA search.
2. After each tool: run its `test_*.py` + existing suite
   (`test_server.py`, `test_units.py`, `test_cif.py`).
3. After each tool: update README tool table (main + dist), dist
   module copies, `make_public_repo.py` whitelist if new files added
   (none expected — all five extend existing modules).
4. Final: `test_server.py` tool-count bump (27 → 31), LISTING.md
   bullet for the general layer, dist final sync + verify.

## Total estimate

~3 days (0.5 + 0.5 + 1 + 0.5 + ~0) + 0.5 day rollout/verification.
No new dependencies. No new auth paths. No commercial-boundary changes.
