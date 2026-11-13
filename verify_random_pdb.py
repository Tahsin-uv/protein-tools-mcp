"""Verify protein-tools usability against a TOTALLY RANDOM PDB entry from RCSB.

Steps:
 1. Generate a random PDB ID and confirm it exists via the RCSB REST API.
 2. Download the actual PDB file.
 3. Run it through the MCP server (real stdio client session) with
    validate_pdb / metal_geometry / sequence_from_pdb / plddt_summary / rcsb_entry.
 4. Cross-validate the PDB-extracted sequence against RCSB's canonical
    entity sequence and, when available, the live UniProt entry.
"""
import asyncio
import json
import os
import random
import string
import urllib.error
import urllib.request

ROOT = r"e:\Protein characterization"
PDB_OUT = ROOT + r"\mcp-protein-tools\random_test.pdb"
PY = ROOT + r"\.venv\Scripts\python.exe"
SERVER = ROOT + r"\mcp-protein-tools\server.py"

import sys
sys.path.insert(0, ROOT + r"\mcp-protein-tools")
import pdb_utils  # noqa: E402  (local, for metal detection only)


def fetch(url: str, binary: bool = False):
    req = urllib.request.Request(url, headers={"User-Agent": "protein-tools-test/0.1"})
    with urllib.request.urlopen(req, timeout=30) as r:
        data = r.read()
    return data if binary else data.decode()


def random_existing_pdb_id():
    """Genuinely random: [1-9] + 3 random alphanumerics, verified to exist."""
    for _ in range(100):
        pid = random.choice(string.digits[1:]) + "".join(
            random.choices(string.ascii_lowercase + string.digits, k=3))
        try:
            d = json.loads(fetch(f"https://data.rcsb.org/rest/v1/core/entry/{pid}"))
            if d.get("rcsb_id"):
                return pid
        except urllib.error.HTTPError:
            continue
    raise RuntimeError("no random PDB ID found in 100 tries")


async def run_tools(pdb_path: str, pdb_id: str, session) -> dict:
    s = session
    res = await s.call_tool("validate_pdb", {"pdb_path": pdb_path})
    v = json.loads(res.content[0].text)
    print("--- validate_pdb ---")
    print(json.dumps(v, indent=1)[:1100])

    text = open(pdb_path, errors="replace").read()
    metals = {a["element"].strip().upper() for a in pdb_utils.parse_pdb(text)["atoms"]}
    metals &= {"FE", "ZN", "MN", "MG", "CU", "NI", "CO", "CD", "HG"}
    print(f"metals present: {sorted(metals) or 'none'}")

    res = await s.call_tool("sequence_from_pdb", {"pdb_path": pdb_path})
    seqs = json.loads(res.content[0].text)
    print("--- sequence_from_pdb ---")
    for ch, d in seqs.items():
        print(f"chain {ch}: {d['length']} aa  {d['sequence'][:60]}...")

    res = await s.call_tool("plddt_summary", {"pdb_path": pdb_path})
    print("--- plddt_summary ---")
    print(res.content[0].text[:300])

    res = await s.call_tool("rcsb_entry", {"pdb_id": pdb_id})
    print("--- rcsb_entry (premium, live) ---")
    print(res.content[0].text[:500])
    return seqs


async def hunt_metal(max_tries: int = 12):
    """Keep drawing random IDs until one contains a metal site; analyze it."""
    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    params = StdioServerParameters(
        command=PY, args=[SERVER],
        env={**os.environ, "PROTEIN_TOOLS_API_KEY": "demo-key"},
    )
    async with stdio_client(params) as (read, write):
        async with ClientSession(read, write) as s:
            await s.initialize()
            first_path = ROOT + r"\mcp-protein-tools\random_test.pdb"
            first_id = open(ROOT + r"\mcp-protein-tools\random_test.id").read().strip()
            seqs = await run_tools(first_path, first_id, s)
            print()
            cross_check(first_id, seqs)

            print(f"\n=== hunting random entries with a metal site (max {max_tries}) ===")
            for i in range(max_tries):
                pid = random_existing_pdb_id()
                path = ROOT + rf"\mcp-protein-tools\random_{pid}.pdb"
                open(path, "w", encoding="utf-8").write(
                    fetch(f"https://files.rcsb.org/download/{pid}.pdb"))
                text = open(path, errors="replace").read()
                metals = {a["element"].strip().upper() for a in pdb_utils.parse_pdb(text)["atoms"]}
                metals &= {"FE", "ZN", "MN", "MG", "CU", "NI", "CO", "CD", "HG"}
                print(f"try {i + 1}: {pid} -> metals: {sorted(metals) or 'none'}")
                if metals:
                    m = sorted(metals)[0]
                    res = await s.call_tool("metal_geometry",
                                            {"pdb_path": path, "metal_element": m})
                    g = json.loads(res.content[0].text)
                    print(f"--- metal_geometry({m}) on {pid} --- n_sites={g['n_sites']}")
                    if g["metal_sites"]:
                        print(json.dumps(g["metal_sites"][0], indent=1)[:800])
                    res = await s.call_tool("rcsb_entry", {"pdb_id": pid})
                    print("--- rcsb_entry ---")
                    print(res.content[0].text[:350])
                    return
            print("no metal-containing random entry found within the try budget")



def cross_check(pdb_id: str, seqs: dict) -> None:
    entry = json.loads(fetch(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"))
    n_ent = entry["rcsb_entry_info"]["polymer_entity_count"]
    canon, uids = [], []
    for i in range(1, n_ent + 1):
        ent = json.loads(fetch(
            f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{i}"))
        seq = ent["entity_poly"]["pdbx_seq_one_letter_code_can"].replace("\n", "")
        canon.append((i, seq))
        u = ent.get("rcsb_polymer_entity_container_identifiers", {}).get("uniprot_ids") or []
        uids.extend(u)
    print(f"--- cross-check vs RCSB ({n_ent} polymer entity/entities) ---")

    def is_subseq(s: str, t: str) -> bool:
        it = iter(t)
        return all(c in it for c in s)

    matches = 0
    for ch, d in seqs.items():
        ours = d["sequence"]
        best = max(canon, key=lambda c: (is_subseq(ours, c[1]),
                                         sum(a == b for a, b in zip(ours, c[1]))
                                         / max(len(ours), len(c[1]))))
        sub = is_subseq(ours, best[1])
        matches += sub
        cov = sum(a == b for a, b in zip(ours, best[1])) / max(len(ours), len(best[1]))
        note = ("ordered subsequence of canonical (position-wise overlap "
                f"{cov:.0%}; missing residues = expression tag / disordered "
                "residues not in the model)" if sub
                else "NOT an ordered subsequence — check extraction")
        print(f"chain {ch} ({len(ours)} aa): entity {best[0]} "
              f"({len(best[1])} aa) -> {note}")
    print(f"{matches}/{len(seqs)} extracted chains are correct ordered "
          f"subsequences of a canonical entity sequence")
    for uid in sorted(set(uids))[:2]:
        up = json.loads(fetch(f"https://rest.uniprot.org/uniprotkb/{uid}.json"))
        seq = up["sequence"]["value"]
        name = up.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value")
        in_uni = any(d["sequence"] in seq for d in seqs.values())
        print(f"UniProt {uid}: {name} ({len(seq)} aa) - "
              f"{'an extracted chain is contained in it' if in_uni else 'no chain contained'}")


if __name__ == "__main__":
    pid = random_existing_pdb_id()
    print(f"RANDOM PDB ID: {pid}")
    open(PDB_OUT, "w", encoding="utf-8").write(fetch(f"https://files.rcsb.org/download/{pid}.pdb"))
    open(ROOT + r"\mcp-protein-tools\random_test.id", "w").write(pid)
    print(f"downloaded -> {PDB_OUT} ({os.path.getsize(PDB_OUT)} bytes)\n")
    asyncio.run(hunt_metal())

