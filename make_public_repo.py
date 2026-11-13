"""Assemble the PUBLIC GitHub repository from this working tree.

The public repo ships the FREE CORE only (MIT). The commercial half - the
hosted relay, the PDB-wide metal-site index and its build pipeline - is
excluded by design (see LICENSE scope note).

What this script does:
  1. copies a whitelist of files into dist/protein-tools-public/
  2. rewrites hardcoded E:\\Protein characterization paths into portable
     __file__-based lookups, so tests run on any machine
  3. sanitizes the public README (strips private relay/monetization sections)
  4. scans the copy for secrets (API keys, tokens) and fails loudly on a leak
  5. prints per-file sizes; never touches the source tree

Usage:
  python make_public_repo.py            # build + verify
  python make_public_repo.py --out DIR  # custom output directory
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import stat
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUT = os.path.join(HERE, "dist", "protein-tools-public")

# --- what ships (free core, MIT) --------------------------------------------
PUBLIC_FILES = [
    ".gitignore",
    "LICENSE",
    "THIRD_PARTY_NOTICES.md",
    "requirements.txt",
    "pyproject.toml",
    "MANIFEST.in",
    "app_paths.py",
    "setup_client.py",
    "SKILL.md",
    "server.py",
    "pdb_utils.py",
    "cif_tools.py",
    "xvg_tools.py",
    "ff_tools.py",
    "propka_tools.py",
    "run_registry.py",
    "md_report.py",
    "plot_tools.py",
    "premium_tools.py",
    "seq_utils.py",
    "gen_fixtures.py",
    "protein_lab.py",
    "lab_app.py",
    "lab_seq.py",
    "lab_structure.py",
    "app.py",
    "spaces-requirements.txt",
    "test_lab.py",
    "test_units.py",
    "test_cif.py",
    "test_phase12.py",
    "test_phase3.py",
    "test_server.py",
    "test_seq_utils.py",
    "test_detect_sequence.py",
    "test_read_fasta.py",
    "test_pdb_by_id.py",
    "test_setup_client.py",
    "test_aa_fasta_stats.py",
    "test_blast_fallback.py",
    "test_conservation_hint.py",
    "test_structure_stats.py",
    "fixtures/sample.fasta",
    "fixtures/__init__.py",
    "fixtures/tiny_2zn.pdb",
    "fixtures/tiny_zn.pdb",
    "fixtures/small.xvg",
    "fixtures/propka_apo.pka",
    "fixtures/rmsd_core.xvg",
    "fixtures/sasa_per_residue.xvg",
    "fixtures/mdrun_apo",
    "fixtures/A8ABE7_holo_Fe2.pdb",
    "fixtures/fold_a0a977pl24_2ca2_model_0.cif",
    "README.md",
]

# --- what must never ship (commercial / bulky / sensitive) -------------------
EXCLUDED = [
    "relay.py",                 # hosted relay service (commercial)
    "builder_metal_index.py",   # index build pipeline (commercial)
    "builder_modal.py",         # detached Modal orchestration (commercial)
    "test_relay_unit.py",
    "test_relay_live.py",       # contains a live API key
    "GEN_PLAN.md",              # internal dev notes
    "BOOTSTRAP_SPEC.md",        # internal infra spec
    "data/index/metal_sites.db",     # 133 MB flagship dataset (paid asset)
    "data/pdb_index/metal_entries.json",
    "data/runs.db",
    "data/cache", "data/pdb_subset",
    "random_test.pdb", "verify_random_pdb.py",
]

ROOT_LIT = r'[eE]:\\Protein characterization'
SYS_PATH_RE = re.compile(
    r'sys\.path\.insert\(0,\s*r"' + ROOT_LIT + r'\\mcp-protein-tools"\)')
ROOT_ASSIGN_RE = re.compile(
    r'(?P<indent>^[ \t]*)ROOT = r"' + ROOT_LIT + r'"', re.MULTILINE)
# a root-prefixed literal plus any adjacent concatenated literals
ROOT_LITERAL_RE = re.compile(
    r'r"' + ROOT_LIT + r'\\(?P<rel>[^"]*)"'
    r'(?P<tail>(?:\s*r?"[^"]*")*)')

SYS_PATH_NEW = ("sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))")

ROOT_HEADER = (
    "import os\n\n"
    "# Project root holding test fixtures (PDBs/.xvg/.pka). Override with the\n"
    "# PROTEIN_TOOLS_ROOT environment variable to point at your own fixtures.\n"
    'ROOT = os.environ.get("PROTEIN_TOOLS_ROOT", '
    'r"E:\\Protein characterization")\n'
)

# Strip the private tail of the README (relay admin commands, relay test
# files, monetization internals) — the public repo documents only the free
# core. Everything from "## Hosted relay" onward is replaced with a short
# premium-pointer section.
PUBLIC_README_TAIL = """
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
"""


def _sanitize_readme(text: str) -> str:
    cut = text.find("## Hosted relay")
    if cut == -1:
        return text
    return text[:cut].rstrip() + "\n" + PUBLIC_README_TAIL.lstrip("\n")


SECRET_PATTERNS = [
    (re.compile(r"pt_live_[0-9a-f]{20,}"), "live relay API key"),
    (re.compile(r"sk-[A-Za-z0-9]{20,}"), "openai-style secret"),
    (re.compile(r"(?i)(api[_-]?key|token|secret)\s*=\s*[\"'][^\"']{16,}[\"']"),
     "hardcoded credential"),
]


def _portable(text: str) -> str:
    """Make one Python source file portable: relative imports, env-var root,
    and root-prefixed fixture literals."""
    import re as _re

    def join_lit(m: _re.Match) -> str:
        parts = [m.group("rel")]
        # strip leading separators: os.path.join treats a rooted component
        # as absolute and would collapse the path to 'E:\\...'
        parts += [p.lstrip("\\/").lstrip("\\/")
                  for p in _re.findall(r'r?"([^"]*)"', m.group("tail") or "")]
        inner = ", ".join(f'r"{p}"' for p in parts if p not in ("", "\\", "/"))
        return f"os.path.join(ROOT, {inner})" if inner else "ROOT"

    text = SYS_PATH_RE.sub(SYS_PATH_NEW, text)
    text = ROOT_ASSIGN_RE.sub(lambda _m: ROOT_HEADER.rstrip("\n"), text)
    text = ROOT_LITERAL_RE.sub(join_lit, text)
    return text


DOC_RE = re.compile(
    r'^(?:\s*#[^\n]*\n)*(?P<doc>""".*?"""|\'\'\'.*?\'\'\')\s*\n',
    re.DOTALL)


def _ensure_root(text: str) -> str:
    """Guarantee ROOT/os are defined in files that now reference them."""
    if "ROOT" not in text or "ROOT = " in text:
        return text
    m = DOC_RE.match(text)
    if m:
        return text[:m.end()] + "\n" + ROOT_HEADER + text[m.end():]
    return ROOT_HEADER + text


def _ensure_os(text: str) -> str:
    """Make sure `import os` appears BEFORE the first os.* usage (the
    sys.path rewrite introduces os usage into files that never imported it)."""
    if "os." not in text:
        return text
    first_use = text.index("os.")
    m_imp = re.search(r'^\s*import os\b', text, re.MULTILINE)
    if m_imp and m_imp.start() < first_use:
        return text
    m = DOC_RE.match(text)
    block = "import os\n"
    if m:
        return text[:m.end()] + "\n" + block + text[m.end():]
    return block + text


def assemble(out_dir: str) -> int:
    src, dst = HERE, out_dir
    if os.path.abspath(src) == os.path.abspath(dst):
        print("refusing to build over the source tree")
        return 2
    if os.path.isdir(dst):
        # git stores objects read-only, so plain rmtree raises PermissionError
        # on Windows. Force-clear the read-only bit before deleting.
        def _on_rm_error(func, path, exc_info):
            try:
                os.chmod(path, stat.S_IWRITE)
                func(path)
            except OSError:
                pass

        shutil.rmtree(dst, onerror=_on_rm_error)
    os.makedirs(dst)

    copied, missing = [], []
    for name in PUBLIC_FILES:
        s = os.path.join(src, name)
        if not os.path.exists(s):
            missing.append(name)
            continue
        d = os.path.join(dst, name)
        os.makedirs(os.path.dirname(d), exist_ok=True)
        if os.path.isdir(s):
            # fixture directories (e.g. fixtures/mdrun_apo): merge contents
            for root, _dirs, files in os.walk(s):
                for fn in files:
                    sf = os.path.join(root, fn)
                    rel = os.path.relpath(sf, s)
                    df = os.path.join(d, rel)
                    os.makedirs(os.path.dirname(df), exist_ok=True)
                    shutil.copy2(sf, df)
        else:
            shutil.copy2(s, d)
        copied.append(name)

    # make Python files portable: absolute fixture paths -> env-var lookup
    rewritten = []
    for name in copied:
        if not name.endswith(".py"):
            continue
        p = os.path.join(dst, name)
        with open(p, encoding="utf-8") as f:
            text = f.read()
        if "Protein characterization" not in text:
            continue
        new = _ensure_os(_ensure_root(_portable(text)))
        if new != text:
            with open(p, "w", encoding="utf-8") as f:
                f.write(new)
            rewritten.append(name)

    # sanitize the public README: strip private relay/monetization internals
    readme_p = os.path.join(dst, "README.md")
    if os.path.exists(readme_p):
        with open(readme_p, encoding="utf-8") as f:
            readme_text = f.read()
        sanitized = _sanitize_readme(readme_text)
        if sanitized != readme_text:
            with open(readme_p, "w", encoding="utf-8") as f:
                f.write(sanitized)
            print("  sanitized private sections out of public README.md")

    # secret scan over the whole copy
    leaks = []
    for root, _dirs, files in os.walk(dst):
        for fn in files:
            fp = os.path.join(root, fn)
            try:
                with open(fp, encoding="utf-8", errors="ignore") as f:
                    blob = f.read()
            except OSError:
                continue
            for pat, what in SECRET_PATTERNS:
                if pat.search(blob):
                    leaks.append((os.path.relpath(fp, dst), what))

    print("=" * 66)
    print("PUBLIC REPO ASSEMBLY")
    print("=" * 66)
    for name in sorted(copied):
        kb = os.path.getsize(os.path.join(dst, name)) / 1024
        print(f"  copied  {kb:8.1f} KB  {name}")
    if missing:
        print("\n  MISSING (not copied): " + ", ".join(missing))
    print("\n  path-rewritten for portability: "
          f"{', '.join(rewritten) if rewritten else 'none'}")
    print(f"  excluded from the public repo ({len(EXCLUDED)}):")
    for name in EXCLUDED:
        print(f"      {name}")

    total = sum(os.path.getsize(os.path.join(r, f))
                for r, _d, fs in os.walk(dst) for f in fs)
    print(f"\n  total public repo size: {total / 1024:.1f} KB")

    if leaks:
        print("\n  !! SECRET SCAN FAILED:")
        for path, what in leaks:
            print(f"      {path}: {what}")
        return 3

    print("  secret scan: clean")
    print("\n  next steps (manual):")
    print("   1. cd into the output dir, then 'git init' + first commit")
    print("   2. confirm the LICENSE holder/year is correct")
    return 0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=DEFAULT_OUT)
    args = ap.parse_args()
    sys.exit(assemble(args.out))


if __name__ == "__main__":
    main()
