# protein-tools — MCP server for general protein characterization
# (+ best-in-class metal-site module)

A Model Context Protocol (MCP) server that gives any MCP-capable AI client
(Cline, Claude Desktop, VS Code, etc.) general protein-characterization tools:
sequence analysis (protein/DNA/RNA typing, FASTA parsing, PDB lookup),
PDB/mmCIF validation, sequence extraction, pLDDT
confidence summaries, GROMACS trajectory metrics, live database lookups
(UniProt, AlphaFold DB, RCSB) — plus a best-in-class metal-site module
(metal geometry, force-field param checks, PDB-wide metal-site index).

Built for real protein research (A8AB protein family, 150 ns GROMACS MD,
manuscript in preparation); stdlib-only analysis core
(`pdb_utils.py`, `seq_utils.py`), so the only dependency is the `mcp` SDK.
The metal-site module (geometry fitting, bidentate detection, PDB-wide
index) is the differentiator — everything else is general protein/MD tooling
any structural biologist can use.

## Install

```powershell
pip install protein-tools-mcp            # base: stdlib-only analysis core
pip install "protein-tools-mcp[plots]"   # + matplotlib for plot_xvg / compare_plot
```

Or zero-install via [uv](https://docs.astral.sh/uv/) — no venv, no paths:

```powershell
uvx protein-tools-mcp
```

(Pins to the stable 1.x SDK: `mcp<2` is load-bearing since `mcp>=2`
renamed FastMCP to MCPServer.)

## Run / test

```powershell
python test_server.py
protein-tools-setup --dry-run   # preview client registration, writes nothing
```

## Add to an MCP client (e.g. VS Code / Cline `.mcp.json`)

Easiest — auto-register (backs up configs first, `--dry-run` to preview):

```powershell
protein-tools-setup
```

Manual equivalent:

```json
{
  "mcpServers": {
    "protein-tools": {
      "command": "uvx",
      "args": ["protein-tools-mcp"],
      "env": { "PROTEIN_TOOLS_API_KEY": "your-key-here" }
    }
  }
}
```

### Source checkout (developers)

```powershell
pip install "mcp<2"          # + "matplotlib>=3.7" for the plot tools
python test_server.py
```

```json
{
  "mcpServers": {
    "protein-tools": {
      "command": "python",
      "args": ["<path-to>/protein-tools-mcp/server.py"]
    }
  }
}
```

## Tests

13 offline suites, all self-contained — fixtures are small synthetic inputs
in `fixtures/` (~96 KB, no experimental data), regenerated deterministically
by `gen_fixtures.py`:

```powershell
python gen_fixtures.py        # regenerate fixtures/ (optional)
python test_units.py          # parsing correctness, metal geometry, xvg tools
python test_cif.py            # mmCIF parse + Ca2+ fixture
python test_phase12.py        # pre-flight gate, PROPKA, registry, run reports
python test_phase3.py         # plotting, premium disk cache
python test_server.py         # live 31-tool MCP round-trip over stdio
python test_setup_client.py   # client-config registration (8 checks)
python test_lab.py            # Protein Lab web-UI helpers (headless)
# plus: test_seq_utils, test_detect_sequence, test_read_fasta,
#   test_pdb_by_id, test_aa_fasta_stats, test_blast_fallback,
#   test_conservation_hint, test_structure_stats
```

Run state (registry DB, premium cache) lives in a per-user dir, never in
the repo: `%APPDATA%\protein-tools` on Windows, `~/.protein-tools`
elsewhere, overridable via `PROTEIN_TOOLS_HOME`.

## Protein Lab (web UI for non-programmers)

Students and wet-lab biologists can skip the MCP setup entirely — a
point-and-click web app wraps the same tested analysis code:

```powershell
pip install "protein-tools-mcp[lab]"   # + streamlit
protein-lab                            # opens http://localhost:8501
```

- **🔤 Sequence tab** — paste letters, a FASTA record, or upload `.fasta`;
  get type verdicts, length/weight/charge/pI cards, a composition chart,
  and "What does this mean?" explainers + glossary.
- **🏗️ Structure tab** — upload `.pdb`/`.cif` (or one-click demo
  structures); get a health verdict, chain-gap notes, and metal-site cards
  in plain English ("the iron atom is held by 4 protein atoms in a
  pyramid" instead of `coordination_number: 4`).
- **Classroom mode** — demo datasets are the synthetic `fixtures/`, so it
  works fully offline; nothing is uploaded anywhere.

Deploy to [Hugging Face Spaces](https://huggingface.co/spaces) (Streamlit
SDK): upload the repo, set `app.py` as the entry file — dependencies come
from `spaces-requirements.txt`. The Space should be named
`protein-tools-lab` to keep the family visible in search. `test_lab.py`
covers the lab helpers headlessly (no browser needed).

## Tools

Structure tools accept **both PDB and mmCIF** (format auto-detected by
`.cif` extension or a leading `data_` record — covers AlphaFold 3 and
modern RCSB depositions):

| Tool | Tier | What it does |
|------|------|--------------|
| `validate_pdb` | free | Format check for PDB/mmCIF: counts, format errors, unknown residues, chain gaps, duplicate atom names, insertion codes, multi-model (NMR) notes |
| `metal_geometry` | free | Metal-site analysis: N/O/S donor contacts with distances + donor types, ligand-metal-ligand angles fitted to ideal geometries (octahedral, tetrahedral, ...), bidentate carboxylate chelate detection |
| `sequence_from_pdb` | free | One-letter sequence per chain (MD variants handled; model 1 for NMR) |
| `plddt_summary` | free | AlphaFold pLDDT bands + low-confidence regions; refuses to misread experimental B-factors |
| `xvg_summary` | free | Per-column stats for any GROMACS .xvg (energies, temperature, SASA, ...) |
| `rmsd_summary` | free | Equilibration assessment: drift, plateau detection, verdict |
| `rmsf_peaks` | free | Most flexible residues + contiguous flexible segments (mean + 1 stdev) |
| `md_preflight` | free | GO/NO-GO gate before MD: validation + metal sites + force-field param check |
| `check_metal_params` | free | Scan GROMACS .ff/.itp for metal atomtypes (sigma/epsilon or MISSING) |
| `propka_summary` | free | PROPKA .pka parsing: pI, extremes, protonation concerns (HIS/acid/CYS flags) |
| `md_run_report` | free | One-call run-dir verdict: RMSD plateau, crash scan, artifacts (GO/FINISHED/INCOMPLETE/CRASHED) |
| `compare_runs` | free | Side-by-side RMSD verdicts for two run dirs (apo vs holo) |
| `log_run` / `list_runs` / `run_details` | free | Persistent SQLite run registry - cross-session project memory |
| `plot_xvg` | free | Plot any .xvg to PNG (decimated for large series) |
| `compare_plot` | free | Overlay multiple .xvg files (e.g. apo vs holo) into one comparison PNG |
| `query_metal_sites` | premium (relay) | Search the hosted PDB-wide metal-site index (38,143 entries / 1,097,659 sites): filter by metal, geometry, CN range; CN-0 artifacts excluded by default. Requires `PROTEIN_TOOLS_RELAY_URL` + API key |
| `index_stats` | premium (relay) | Corpus stats for the hosted metal-site index (per-metal / per-geometry distributions, artifact volume). Requires `PROTEIN_TOOLS_RELAY_URL` + API key |
| `cache_info` | free | Inspect or clear the premium-lookup response cache |
| `uniprot_entry` | premium | Live UniProtKB lookup (name, organism, function, metal-binding sites), 24 h disk cache |
| `alphafold_plddt` | premium | Live AlphaFold DB model + confidence summary |
| `rcsb_entry` | premium | Live RCSB PDB entry metadata |
| `detect_sequence_type` | free | Detect whether a FASTA sequence is protein, DNA, RNA, ambiguous, or unknown. Accepts FASTA file path or raw FASTA text; returns type classification with confidence level and detailed notes for each sequence. |
| `read_fasta` | free | Read one or more sequences from a FASTA file or raw text. Returns full parsed records including sequence strings, type classifications, lengths, and heuristic UniProt accession hints. Handles multi-record FASTA, file-not-found errors, and empty input. |
| `pdb_by_id` | premium (local ok) | Get info about a PDB entry by its 4-character ID: checks for a local PDB/mmCIF file and if found runs validation, sequence extraction, and metal-site analysis; also fetches RCSB metadata (requires API key for full data, gracefully degrades to local-only without). |
| `find_pdb_by_sequence` | premium | Find PDB structures matching a protein sequence from a FASTA file or raw text via RCSB BLAST search. Parses the FASTA, verifies it's a protein sequence, then searches RCSB and returns matching PDB IDs with metadata. Requires API key. |

All tools return structured errors (`{"error": ..., "hint": ...}`) so AI
agents can self-correct (e.g. fix a bad path) instead of crashing.

Correctness notes: residue keys include PDB insertion codes (column 27);
multi-MODEL (NMR) files are analyzed on model 1 only; coordination numbers
count only N/O/S donor atoms; angle-based geometry fitting compares sorted
observed ligand-metal-ligand angles to ideal sets (linear / trigonal-planar /
tetrahedral / square-planar / TBP / square-pyramidal / octahedral).

Premium tools are gated by the `PROTEIN_TOOLS_API_KEY` environment variable;
without it they raise a PermissionError explaining how to get a key.
## Hosted relay (premium mode)

Set `PROTEIN_TOOLS_RELAY_URL` alongside `PROTEIN_TOOLS_API_KEY` to serve
`uniprot_entry`, `alphafold_plddt`, `rcsb_entry`, `query_metal_sites`, and
`index_stats` from the hosted relay (key validation, quotas, 24 h
server-side caching) instead of direct upstream calls and the local index
DB. Without the relay URL the server runs fully locally.

Premium tools are gated by the `PROTEIN_TOOLS_API_KEY` environment variable;
without it they return a structured error explaining how to get a key.

## Disclaimer

Free tools are research-grade (no warranty); `geometry_guess` is a first-pass
heuristic, not a substitute for manual inspection of metal coordination.
