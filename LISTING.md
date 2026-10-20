# Listing kit — protein-tools MCP server (general protein characterization
# + best-in-class metal-site module)

Everything needed to submit the server to MCP directories and to talk
about it. Update the "proof" numbers before submitting.

---

## Name

`protein-tools` — Protein Characterization Tools for AI agents

## One-liner (short description, <160 chars)

> Tested protein-analysis tools for AI agents: PDB/mmCIF validation, metal-site
> geometry, MD verdicts, UniProt/AlphaFold lookups — hardened on real research.

## Long description

protein-tools gives LLM agents reliable, tested protein-characterization
capabilities through the Model Context Protocol:

- **Structure QC for PDB *and* mmCIF** (AlphaFold 3 / RCSB): fixed-column
  format validation, chain gaps, insertion codes, altLocs, multi-model handling
- **Metal-site module (the differentiator)**: N/O/S coordination
  chemistry, ligand-metal-ligand angles fitted to ideal geometries
  (octahedral, tetrahedral, ...), bidentate carboxylate chelate detection
- **GROMACS trajectory intelligence**: .xvg summaries, RMSD equilibration
  verdicts, RMSF peak extraction, and one-call GO/NO-GO run reports
- **MD pre-flight gating**: force-field metal-parameter checks (catches
  missing ZN/FE atomtypes before pdb2gmx burns GPU-hours)
- **PROPKA parsing** with protonation-state flags (HIS near pH 7, CYS
  disulfide markers)
- **Persistent run registry** (SQLite): cross-session project memory
- **Live lookups** (UniProt / AlphaFold DB / RCSB) with a 24 h disk cache
- **PDB-wide metal-site index** (flagship): 38,143 metal-containing entries
  from the RCSB archive analyzed with uniform QC - 1,097,659 sites with
  coordination number, donor chemistry, angle-fitted geometry class,
  bidentate-chelates and adjacent-metal (cluster) flags. Query it for
  things like "all tetrahedral Cys4 zinc fingers" or "catalytic zinc with
  a labile water" in milliseconds
- **Hosted relay** for production use: key validation, per-tier quotas,
  metering, and server-side caching (`relay.py`, FastAPI on Modal)
- **Plotting**: .xvg to PNG, apo-vs-holo comparison overlays

### Provenance / why it is different

This server was built *inside* a real metalloprotein characterization study
(Zn/Fe-binding protein family, 150 ns GROMACS MD, manuscript in preparation)
and regression-tested against random RCSB entries. Every tool exists because
the author actually needed it — including the failure classes that killed real
pipeline runs (PDB column-shift corruption, missing metal atomtypes,
unconverged trajectories).

### Pricing model

Free core (all structural/trajectory/plotting tools). Premium = hosted relay
with API keys (built): hosted PDB-wide metal-site index (38,143 entries /
1,097,659 sites), cached UniProt/AlphaFold/RCSB lookups, per-tier quotas (beta 1,000 / pro 10,000 / lab 50,000 calls per
month). Phase B wires automatic key issuance to a payment webhook
(Lemon Squeezy subscriptions — merchant-of-record that onboards Bangladesh;
and/or x402 pay-per-call in USDC).

---

## Submission targets

| Directory | URL | Notes |
|-----------|-----|-------|
| mcp.so | https://mcp.so | Free submission; short + long description |
| PulseMCP | https://pulsemcp.com | Free; has categories (Science) |
| Smithery | https://smithery.ai | Needs a `smithery.yaml` config; largest registry |
| Glama | https://glama.ai/mcp/servers | Free; auto-indexes GitHub repos |
| MCP Market | https://www.mcpmarket.com | Free listing; featured placement paid |

## Submission checklist

- [ ] Public GitHub repo (README.md is already client-ready)
- [ ] Screenshots in `listing-assets/` (generate: `python test_phase3.py`)
- [ ] License file (MIT recommended)
- [ ] smithery.yaml (only for Smithery)
- [ ] Pin versions in a requirements.txt: `mcp<2` (+ matplotlib for plotting)

## Demo flow for a GIF/video

1. Drop a metalloprotein PDB into a session and ask:
   "run md_preflight on this structure against my force field"
2. Show the GO/NO-GO verdict JSON
3. Ask: "report on my last MD run" — show md_run_report verdict
4. Ask: "plot apo vs holo RMSD" — show the comparison PNG
