import pathlib
import os

base = pathlib.Path(r"e:\Protein characterization\mcp-protein-tools")
dist = base / "dist" / "protein-tools-public"

print("=" * 60)
print("STATE CHECK")
print("=" * 60)
print()

# 1. Main folder test files
print("1. Test files in main folder:")
for name in ["test_seq_utils.py", "test_detect_sequence.py", "test_read_fasta.py", "test_pdb_by_id.py"]:
    p = base / name
    status = "EXISTS" if p.exists() else "MISSING"
    size = p.stat().st_size if p.exists() else 0
    print(f"   {name}: {status} ({size} bytes)")

print()

# 2. Public snapshot
print("2. Public snapshot (dist/protein-tools-public):")
print(f"   dist exists: {dist.exists()}")
if dist.exists():
    items = list(dist.iterdir())
    print(f"   total items: {len(items)}")
    
    # seq_utils.py
    seq = dist / "seq_utils.py"
    print(f"\n   seq_utils.py:")
    print(f"   exists: {seq.exists()}")
    if seq.exists():
        c = seq.read_text()
        print(f"   parse_fasta: {'def parse_fasta' in c}")
        print(f"   normalize_sequence: {'def normalize_sequence' in c}")
        print(f"   classify_sequence: {'def classify_sequence' in c}")
    
    # server.py
    srv = dist / "server.py"
    print(f"\n   server.py:")
    print(f"   exists: {srv.exists()}")
    if srv.exists():
        c = srv.read_text()
        print(f"   size: {len(c)} chars")
        print(f"   import seq_utils: {'import seq_utils' in c}")
        for t in ["detect_sequence_type", "read_fasta", "pdb_by_id", "find_pdb_by_sequence"]:
            has = ("def " + t) in c
            print(f"   def {t}: {has}")
        print(f"   @mcp.tool() count: {c.count('@mcp.tool()')}")
    
    # test files in dist
    print(f"\n   Test files in dist:")
    for name in ["test_seq_utils.py", "test_detect_sequence.py", "test_read_fasta.py", "test_pdb_by_id.py"]:
        p = dist / name
        status = "EXISTS" if p.exists() else "MISSING"
        size = p.stat().st_size if p.exists() else 0
        print(f"   {name}: {status} ({size} bytes)")
    
    # README
    rd = dist / "README.md"
    print(f"\n   README.md:")
    print(f"   exists: {rd.exists()}")
    if rd.exists():
        c = rd.read_text()
        print(f"   size: {len(c)} chars")
        for t in ["detect_sequence_type", "read_fasta", "pdb_by_id", "find_pdb_by_sequence"]:
            present = ("`" + t + "`") in c
            print(f"   `{t}` in table: {present}")
else:
    print("   DIST FOLDER DOES NOT EXIST")

print()
print("=" * 60)
print("MAIN FOLDER ALL .py FILES:")
print("=" * 60)
for f in sorted(base.glob("*.py")):
    if f.is_file():
        size = f.stat().st_size
        print(f"   {f.name}: {size} bytes")
