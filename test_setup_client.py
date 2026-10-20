"""Tests for setup_client: dry-run, merge, backup, idempotency, uninstall."""
import json
import os
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import setup_client as sc

TMP = Path(tempfile.mkdtemp(prefix="pt_setup_"))

# --- 1. dry-run on a missing file: writes nothing ---
p1 = TMP / "fresh" / ".mcp.json"
msg = sc.install_one("t", p1, None, dry_run=True)
assert "DRY-RUN would create" in msg, msg
assert not p1.exists(), "dry-run must not create files"
print("1. dry-run creates nothing: OK")

# --- 2. real install: creates file with uvx entry + placeholder key ---
msg = sc.install_one("t", p1, None, dry_run=False)
assert "created" in msg and "(new file)" in msg, msg
cfg = json.loads(p1.read_text())
e = cfg["mcpServers"]["protein-tools"]
assert e["command"] == "uvx" and e["args"] == ["protein-tools-mcp"], e
assert e["env"]["PROTEIN_TOOLS_API_KEY"] == "your-key-here", e
print("2. fresh install writes uvx entry: OK")

# --- 3. idempotent re-run: no change, no backup ---
msg = sc.install_one("t", p1, None, dry_run=False)
assert "already configured" in msg, msg
assert sorted(p1.parent.iterdir()) == [p1], "re-run must not add backups"
print("3. idempotent re-run: OK")

# --- 4. merge: other servers preserved; explicit --api-key overwrites ---
p2 = TMP / "merge" / "mcp.json"
p2.parent.mkdir(parents=True)
p2.write_text(json.dumps({
    "mcpServers": {
        "other": {"command": "node", "args": ["x.js"]},
        "protein-tools": {
            "command": "uvx", "args": ["protein-tools-mcp"],
            "env": {"PROTEIN_TOOLS_API_KEY": "pt_live_abc123"}}},
    "otherTopLevel": True}))
# no --api-key: existing entry already valid -> no change, real key kept
msg = sc.install_one("t", p2, None, dry_run=False)
assert "already configured" in msg, msg
cfg = json.loads(p2.read_text())
assert cfg["mcpServers"]["other"]["command"] == "node"
kept = cfg["mcpServers"]["protein-tools"]["env"]["PROTEIN_TOOLS_API_KEY"]
assert kept == "pt_live_abc123", kept
print("4a. existing valid entry untouched, real key kept: OK")
# explicit --api-key: entry rewritten, backup taken, others preserved
msg = sc.install_one("t", p2, "pt_live_NEW2", dry_run=False)
assert "updated" in msg and "backup" in msg, msg
cfg = json.loads(p2.read_text())
assert cfg["mcpServers"]["other"]["command"] == "node"
assert cfg["otherTopLevel"] is True
kept = cfg["mcpServers"]["protein-tools"]["env"]["PROTEIN_TOOLS_API_KEY"]
assert kept == "pt_live_NEW2", kept
baks = list(p2.parent.glob("mcp.json.bak-*"))
assert len(baks) == 1, baks
print("4b. --api-key overwrites key, preserves rest, backup written: OK")

# --- 5. JSONC tolerance: comments + trailing commas parse ---
p3 = TMP / "jsonc" / ".mcp.json"
p3.parent.mkdir(parents=True)
p3.write_text('{\n  // user comment\n  "mcpServers": {\n'
              '    "other": {"command": "x",},  /* trailing */\n  },\n}\n')
msg = sc.install_one("t", p3, None, dry_run=False)
assert "updated" in msg, msg
cfg = json.loads(p3.read_text())
assert "protein-tools" in cfg["mcpServers"] and "other" in cfg["mcpServers"]
print("5. JSONC (comments + trailing commas) tolerated: OK")

# --- 6. corrupt file: refused, untouched ---
p4 = TMP / "bad" / ".mcp.json"
p4.parent.mkdir(parents=True)
p4.write_text("{ not json {{{")
msg = sc.install_one("t", p4, None, dry_run=False)
assert msg.startswith("[t] SKIP:"), msg
assert p4.read_text() == "{ not json {{{"
print("6. corrupt config refused, untouched: OK")

# --- 7. uninstall removes only our entry ---
msg = sc.uninstall_one("t", p2, dry_run=False)
assert "removed entry" in msg, msg
cfg = json.loads(p2.read_text())
assert "protein-tools" not in cfg["mcpServers"]
assert cfg["mcpServers"]["other"]["command"] == "node"
print("7. uninstall removes only our entry: OK")

# --- 8. --api-key stores the key ---
p5 = TMP / "key" / ".mcp.json"
msg = sc.install_one("t", p5, "pt_live_NEW", dry_run=False)
cfg = json.loads(p5.read_text())
assert cfg["mcpServers"]["protein-tools"]["env"][
    "PROTEIN_TOOLS_API_KEY"] == "pt_live_NEW"
print("8. --api-key stored: OK")

print("ALL SETUP_CLIENT TESTS PASSED")
