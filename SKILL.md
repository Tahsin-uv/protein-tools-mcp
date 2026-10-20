# protein-tools skill — how to drive the protein-tools MCP server

You are an agent with the `protein-tools` MCP server: general protein
characterization (sequence, structure, MD) plus a best-in-class metal-site
module. 31 tools. Follow this file and you will use them correctly.

## 1. Golden rules

1. **Validate before you analyze.** Any structure file (PDB/mmCIF) goes
   through `validate_pdb` first. If `ok` is false, read `format_errors`
   and fix the input — never analyze a file you haven't validated.
2. **`is_file` defaults to True.** The sequence tools (`detect_sequence_type`,
   `read_fasta`, `pdb_by_id`, `find_pdb_by_sequence`, `aa_composition`,
   `fasta_stats`) treat `source` as a *file path* unless you pass
   `is_file=False` for raw FASTA text. Forgetting this is the #1 failure:
   `FileNotFoundError: FASTA file not found: ACDEF...` means you passed
   sequence text without `is_file=False`.
3. **Type before search.** `find_pdb_by_sequence` only accepts *protein*.
   Run `detect_sequence_type` first; if the answer isn't `protein`, stop
   and say so instead of searching.
4. **Read the metal section.** `md_preflight` can return GO for a pure
   organic protein while saying nothing about metals. A GO verdict means
   "your metals are fine" ONLY if the `metal_geometry` section lists the
   expected sites. Quote the section, don't just quote the verdict.
5. **Errors are structured — obey the hint.** Every tool returns
   `{"error": ..., "hint": ...}` on failure. Read `hint`, fix the input,
   retry once. Never retry blindly.
6. **Premium tools need a key.** `uniprot_entry`, `alphafold_plddt`,
   `rcsb_entry`, `find_pdb_by_sequence` (and relay-mode index tools)
   require `PROTEIN_TOOLS_API_KEY`. Without it you get a PermissionError
   explaining how to get one — surface that to the user, don't work
   around it.
7. **Index tools need the index.** `query_metal_sites` / `index_stats`
   need the local `metal_sites.db` or the hosted relay
   (`PROTEIN_TOOLS_RELAY_URL`). Without either, the error says so.

## 2. Tool selection

| User wants | Call (in order) |
|---|---|
| Is this structure file sane? | `validate_pdb` |
| Metal coordination / geometry? | `validate_pdb` → `metal_geometry` (pick element; default cutoff 3.0 Å) |
| Sequence of a structure? | `sequence_from_pdb` |
| Is this AlphaFold model trustworthy here? | `plddt_summary` (refuses experimental B-factors — correct) |
| What is this sequence (protein/DNA/RNA)? | `detect_sequence_type` (file path OR text + `is_file=False`) |
| Parse a FASTA file? | `read_fasta` (multi-record OK; returns UniProt accession hints) |
| Protein properties (MW, pI, GRAVY)? | `aa_composition` / `fasta_stats` |
| Find PDBs for my sequence? | `detect_sequence_type` → `find_pdb_by_sequence` (protein only, premium) |
| Info on PDB 1ABC? | `pdb_by_id` (local file if present + RCSB metadata) |
| UniProt / AlphaFold / RCSB lookup? | `uniprot_entry` / `alphafold_plddt` / `rcsb_entry` (premium, 24 h cache; `cache_info` to inspect) |
| Search known metal sites (ZN tetra, ...)? | `query_metal_sites` (defaults exclude CN-0 artifacts; `min_cn=1` shows CN-1) |
| Ready for MD? | `md_preflight` (validation + metals + force-field params; READ the metal section) |
| .xvg analysis? | `xvg_summary` (any) / `rmsd_summary` (equilibration) / `rmsf_peaks` (flexible regions) |
| PROPKA pKa file? | `propka_summary` (HIS/CYS flags near working pH) |
| Run verdict / apo-vs-holo? | `md_run_report` / `compare_runs` (needs finished-mdrun markers) |
| Plots? | `plot_xvg` / `compare_plot` (needs matplotlib; graceful error otherwise) |
| Remember this run? | `log_run` → `list_runs` / `run_details` (persists in per-user app data) |

## 3. Correctness notes (do not improvise past these)

- Residue keys include PDB insertion codes (column 27): `ARG 100A` ≠ `ARG 100`.
- Multi-MODEL (NMR) files are analyzed on model 1 only.
- Coordination numbers count only N/O/S donor atoms.
- Angle-based geometry fitting compares sorted observed
  ligand-metal-ligand angles to ideal sets (linear / trigonal-planar /
  tetrahedral / square-planar / TBP / square-pyramidal / octahedral) and
  reports RMSD in degrees — quote the RMSD, not just the label.
- `plddt_summary` on an experimental structure (no AlphaFold marker)
  refuses pLDDT interpretation. That refusal is correct, not a bug.
- `query_metal_sites` defaults exclude CN-0 isolated-ion artifacts and
  CN-1. `min_cn=1` reveals unusual CN-1 chemistry; `include_artifacts=True`
  also returns CN-0 rows.
- Run-registry state (`runs.db`) and the premium cache live in per-user
  app data (`%APPDATA%\protein-tools` on Windows, `~/.protein-tools`
  elsewhere, override `PROTEIN_TOOLS_HOME`). Never assume a repo-relative
  `data/` path on the user's machine.

## 4. Worked traces (imitate these sequences)

### Trace A — "is this zinc site catalytic or structural?"

1. `validate_pdb(path)` → confirm `ok`, note chain gaps near the site.
2. `metal_geometry(path, metal="ZN")` → read CN, geometry_fit best + rmsd,
   bidentate chelates, ligand donor types.
3. Tetrahedral Cys4 (rmsd < ~15°) → structural zinc finger language.
   His2Cys2 / His3 + water → catalytic language. Say which and cite the
   numbers; if CN ≤ 2 or rmsd is large, say the site is ambiguous rather
   than forcing a label.
4. Optional: `query_metal_sites(metal="ZN", geometry="tetrahedral")` for
   PDB-wide context on how common the motif is.

### Trace B — "find structures for my sequence"

1. `detect_sequence_type(source, is_file=...)` → must be `protein`.
   If DNA/RNA/ambiguous: stop, explain, do not call the search.
2. `find_pdb_by_sequence(source, is_file=...)` → RCSB BLAST; report top
   PDB IDs. Needs `PROTEIN_TOOLS_API_KEY`.
3. `pdb_by_id(top_id)` → resolution, method, local validation if a file
   is present.

### Trace C — "is my MD run equilibrated?"

1. `md_run_report(run_dir)` → verdict, plateau, crash scan. INCOMPLETE
   with no finished-mdrun marker means the run isn't done, not that it
   failed — say so.
2. `rmsd_summary(rmsd.xvg)` for the numbers (drift %, plateau time).
3. Apo vs holo: `compare_runs(dir_a, dir_b)` — both dirs need a parsable
   `*rmsd*.xvg` or it raises; report that plainly.
4. `log_run(kind="md_report", ...)` to record the verdict for later.

