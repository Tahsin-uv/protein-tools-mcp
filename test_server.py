"""End-to-end smoke test for the protein-tools MCP server.

Spawns the server over stdio (exactly like a real MCP client), lists tools,
and calls the free tools against a real PDB from this workspace, plus one
premium tool with a demo key.

Set PROTEIN_TOOLS_RELAY_URL + PROTEIN_TOOLS_API_KEY in the environment to
additionally exercise relay mode (premium + index tools via the hosted relay).
"""
import asyncio
import os
import sys

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures")
PDB = os.path.join(FIX, "A8ABE7_holo_Fe2.pdb")
XVG = os.path.join(FIX, "rmsd_core.xvg")
SASA = os.path.join(FIX, "sasa_per_residue.xvg")
CIF = os.path.join(FIX, "fold_a0a977pl24_2ca2_model_0.cif")
PY = sys.executable
SERVER = os.path.join(HERE, "server.py")

RELAY = os.environ.get("PROTEIN_TOOLS_RELAY_URL", "")
RELAY_KEY = os.environ.get("PROTEIN_TOOLS_RELAY_TEST_KEY", "")


async def main() -> None:
    env = {**os.environ, "PROTEIN_TOOLS_API_KEY": RELAY_KEY or "test-key-123"}
    if RELAY:
        env["PROTEIN_TOOLS_RELAY_URL"] = RELAY
    params = StdioServerParameters(command=PY, args=[SERVER], env=env)
    async with stdio_client(params) as (read, write):
        async with ClientSession(read, write) as s:
            await s.initialize()
            tools = await s.list_tools()
            print("TOOLS:", [t.name for t in tools.tools])

            for name, args in [
                ("validate_pdb", {"pdb_path": PDB}),
                ("validate_pdb", {"pdb_path": CIF}),
                ("metal_geometry", {"pdb_path": PDB, "metal_element": "FE"}),
                ("sequence_from_pdb", {"pdb_path": PDB}),
                ("xvg_summary", {"xvg_path": XVG}),
                ("rmsd_summary", {"xvg_path": XVG}),
                ("rmsf_peaks", {"xvg_path": SASA, "top_n": 5}),
                ("md_preflight", {"pdb_path": PDB}),
                ("propka_summary", {"pka_path": os.path.join(FIX, "propka_apo.pka")}),
                ("md_run_report", {"run_dir": os.path.join(FIX, "mdrun_apo")}),
                ("list_runs", {"limit": 5}),
                ("plot_xvg", {"xvg_path": XVG}),
                ("query_metal_sites", {"metal": "ZN", "geometry": "tetrahedral",
                                       "limit": 3}),
                ("index_stats", {}),
                ("cache_info", {}),
                ("validate_pdb", {"pdb_path": os.path.join(FIX, "definitely_missing.pdb")}),
            ]:
                res = await s.call_tool(name, args)
                print(f"\n--- {name} ---")
                print(res.content[0].text[:1200])

            try:
                res = await s.call_tool("uniprot_entry", {"accession": "A8ABE7"})
                print("\n--- uniprot_entry(A8ABE7) ---")
                print(res.content[0].text[:900])
            except Exception as exc:  # network may be unavailable in sandbox
                print("\nuniprot_entry skipped:", exc)


asyncio.run(main())
