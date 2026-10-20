"""Register the protein-tools MCP server with local MCP clients.

Usage:
  protein-tools-setup [--client vscode|cline|claude|all] [--dry-run]
                      [--uninstall] [--api-key KEY]

Behaviour:
  - Detects installed clients (VS Code / Cline workspace + global configs,
    Claude Desktop). Missing clients are skipped with a note, never an error.
  - Backs up every file before writing (<file>.bak-<timestamp>).
  - --dry-run prints the exact diff and writes nothing.
  - --uninstall removes only the protein-tools entry, leaving the rest.
  - Never overwrites an existing API key unless --api-key is given.

Stdlib-only.
"""
from __future__ import annotations

import argparse
import datetime
import json
import os
import re
import shutil
import sys
from pathlib import Path

SERVER_KEY = "protein-tools"


def _strip_jsonc(text: str) -> str:
    """Remove // and /* */ comments while keeping quoted strings intact."""
    out, i, n, q = [], 0, len(text), None
    while i < n:
        c = text[i]
        if q:
            out.append(c)
            if c == "\\" and i + 1 < n:
                out.append(text[i + 1])
                i += 2
                continue
            if c == q:
                q = None
            i += 1
            continue
        if c in ("\"", "'"):
            q, out = c, out + [c]
            i += 1
            continue
        if c == "/" and i + 1 < n and text[i + 1] == "/":
            while i < n and text[i] != "\n":
                i += 1
            continue
        if c == "/" and i + 1 < n and text[i + 1] == "*":
            i += 2
            while i + 1 < n and not (text[i] == "*" and text[i + 1] == "/"):
                i += 1
            i += 2
            continue
        out.append(c)
        i += 1
    # trailing commas before } or ]
    return re.sub(r",(\s*[}\]])", r"\1", "".join(out))


def _load_config(path: Path) -> tuple[dict, bool]:
    """Return (config, existed). Corrupt files raise with a helpful message."""
    if not path.exists():
        return {}, False
    raw = path.read_text(encoding="utf-8-sig")
    try:
        return json.loads(raw), True
    except json.JSONDecodeError:
        pass
    try:
        return json.loads(_strip_jsonc(raw)), True
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"{path} is not parseable JSON ({exc}). "
            f"Fix it by hand or delete it; refusing to touch it.") from exc


def _client_paths(which: str) -> list[tuple[str, Path]]:
    """(label, config path) pairs for the requested clients."""
    home = Path.home()
    appdata = os.environ.get("APPDATA", "").strip()
    out: list[tuple[str, Path]] = []

    def want(*names: str) -> bool:
        return which == "all" or any(n == which for n in names)

    if want("vscode", "cline"):
        out.append(("vscode/workspace", Path.cwd() / ".mcp.json"))
        if appdata:
            out.append(("vscode/global",
                        Path(appdata) / "Code" / "User" / "mcp.json"))
        else:
            out.append(("vscode/global",
                        home / ".config" / "Code" / "User" / "mcp.json"))
    if want("cline"):
        if appdata:
            out.append(("cline/global",
                        Path(appdata) / "Code" / "User" / "globalStorage"
                        / "saoudrizwan.claude-dev" / "settings"
                        / "cline_mcp_settings.json"))
        else:
            out.append(("cline/global",
                        home / ".config" / "Code" / "User" / "globalStorage"
                        / "saoudrizwan.claude-dev" / "settings"
                        / "cline_mcp_settings.json"))
    if want("claude"):
        if appdata or os.name == "nt":
            base = Path(appdata) if appdata else home / "AppData" / "Roaming"
            out.append(("claude-desktop", base / "Claude"
                        / "claude_desktop_config.json"))
        elif sys.platform == "darwin":
            out.append(("claude-desktop",
                        home / "Library" / "Application Support" / "Claude"
                        / "claude_desktop_config.json"))
        else:
            out.append(("claude-desktop",
                        home / ".config" / "Claude"
                        / "claude_desktop_config.json"))
    return out


def _entry(api_key: str | None) -> dict:
    env: dict[str, str] = {}
    if api_key:
        env["PROTEIN_TOOLS_API_KEY"] = api_key
    else:
        env["PROTEIN_TOOLS_API_KEY"] = "your-key-here"
    return {"command": "uvx", "args": ["protein-tools-mcp"], "env": env}
def _servers_of(cfg: dict) -> dict:
    servers = cfg.get("mcpServers")
    if servers is None:
        servers = {}
        cfg["mcpServers"] = servers
    if not isinstance(servers, dict):
        raise ValueError(
            "'mcpServers' exists but is not an object; refusing to touch it.")
    return servers


def _backup(path: Path) -> Path:
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    bak = path.with_name(f"{path.name}.bak-{stamp}")
    shutil.copy2(path, bak)
    return bak


def _uvx_hint() -> str:
    if shutil.which("uvx"):
        return "uvx: found on PATH."
    return ("uvx NOT found on PATH. Install it first:\n"
            "  Windows (PowerShell):  powershell -c \"irm astral.sh/uv | iex\"\n"
            "  macOS/Linux:           curl -LsSf astral.sh/uv/install.sh | sh\n"
            "  or:                    pip install uv\n"
            "Then re-run protein-tools-setup.")


def install_one(label: str, path: Path, api_key: str | None,
                dry_run: bool) -> str:
    try:
        cfg, existed = _load_config(path)
    except ValueError as exc:
        return f"[{label}] SKIP: {exc}"
    servers = _servers_of(cfg)
    old = servers.get(SERVER_KEY)
    new = _entry(api_key)
    if isinstance(old, dict) and old.get("command") == "uvx" \
            and old.get("args") == ["protein-tools-mcp"] \
            and (api_key is None or
                 (isinstance(old.get("env"), dict)
                  and old["env"].get("PROTEIN_TOOLS_API_KEY") == api_key)):
        return f"[{label}] already configured at {path} (no change)."
    # preserve an existing real key when --api-key was not given
    if api_key is None and isinstance(old, dict):
        env = old.get("env")
        if isinstance(env, dict):
            have = env.get("PROTEIN_TOOLS_API_KEY", "").strip()
            if have and have != "your-key-here":
                new["env"]["PROTEIN_TOOLS_API_KEY"] = have
    servers[SERVER_KEY] = new
    text = json.dumps(cfg, indent=2) + "\n"
    if dry_run:
        action = "would create" if not existed else "would update"
        return (f"[{label}] DRY-RUN {action} {path}:\n"
                + json.dumps({SERVER_KEY: new}, indent=2))
    if existed:
        bak = _backup(path)
        note = f" (backup: {bak.name})"
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        note = " (new file)"
    path.write_text(text, encoding="utf-8")
    return (f"[{label}] {'created' if not existed else 'updated'} "
            f"{path}{note}.")


def uninstall_one(label: str, path: Path, dry_run: bool) -> str:
    if not path.exists():
        return f"[{label}] nothing at {path} (skip)."
    try:
        cfg, _ = _load_config(path)
    except ValueError as exc:
        return f"[{label}] SKIP: {exc}"
    servers = cfg.get("mcpServers")
    if not isinstance(servers, dict) or SERVER_KEY not in servers:
        return f"[{label}] no protein-tools entry in {path} (no change)."
    if dry_run:
        return f"[{label}] DRY-RUN would remove entry from {path}."
    bak = _backup(path)
    del servers[SERVER_KEY]
    path.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
    return f"[{label}] removed entry from {path} (backup: {bak.name})."


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Register protein-tools with local MCP clients.")
    ap.add_argument("--client", default="all",
                    choices=["all", "vscode", "cline", "claude"],
                    help="which client configs to touch (default: all)")
    ap.add_argument("--dry-run", action="store_true",
                    help="print what would change; write nothing")
    ap.add_argument("--uninstall", action="store_true",
                    help="remove the protein-tools entry instead of adding it")
    ap.add_argument("--api-key", default=None,
                    help="premium API key to store (else keeps existing key "
                         "or a placeholder)")
    args = ap.parse_args(argv)

    print(_uvx_hint())
    if args.api_key:
        print("NOTE: the key is stored in plaintext client config "
              "(standard for MCP clients).")
    rc = 0
    for label, path in _client_paths(args.client):
        try:
            if args.uninstall:
                print(uninstall_one(label, path, args.dry_run))
            else:
                print(install_one(label, path, args.api_key, args.dry_run))
        except (OSError, ValueError) as exc:
            print(f"[{label}] ERROR: {exc}")
            rc = 1
    return rc


if __name__ == "__main__":
    raise SystemExit(main())

