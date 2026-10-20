"""Shared state locations for protein-tools.

Runtime state (run registry DB, premium-lookup cache) lives in a per-user
app-data directory — NOT next to the installed package. Under ``pip install``
or ``uvx`` the package directory is unpredictable or read-only, so
package-relative ``data/`` paths break.

Resolution order:
  1. ``PROTEIN_TOOLS_HOME`` env var (explicit override; tests use this).
  2. ``%APPDATA%\\protein-tools`` on Windows.
  3. ``~/.protein-tools`` everywhere else.

Stdlib-only by design (no ``platformdirs`` dependency).

Backward compatibility: if state already exists at the legacy
package-relative ``data/`` location and the new location is empty, it is
copied over once and the migration is reported (never moved/deleted, so a
downgrade keeps working).
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path


def app_dir() -> Path:
    """Per-user state directory, created on demand."""
    override = os.environ.get("PROTEIN_TOOLS_HOME", "").strip()
    if override:
        d = Path(override).expanduser()
    elif os.name == "nt":
        base = os.environ.get("APPDATA", "").strip() or str(Path.home())
        d = Path(base) / "protein-tools"
    else:
        d = Path.home() / ".protein-tools"
    d.mkdir(parents=True, exist_ok=True)
    return d


def legacy_package_dir() -> Path:
    """Old package-relative ``data/`` dir (source of a one-time migration)."""
    return Path(__file__).resolve().parent / "data"


def _migrate_dir_contents(new_dir: Path, legacy_dir: Path) -> bool:
    """Copy legacy dir contents into ``new_dir`` file-by-file.

    Returns True when at least one file was copied. Existing files are
    never overwritten.
    """
    if not legacy_dir.is_dir():
        return False
    copied = False
    for src in legacy_dir.rglob("*"):
        if not src.is_file():
            continue
        dst = new_dir / src.relative_to(legacy_dir)
        if dst.exists():
            continue
        try:
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
            copied = True
        except OSError:
            continue
    return copied


def _migrate_file(new_path: Path, legacy_path: Path) -> bool:
    """Copy one legacy state file to the new home if the new home lacks it.

    Returns True when a migration happened. Existing new-home files are
    never overwritten.
    """
    if new_path.exists() or not legacy_path.is_file():
        return False
    try:
        new_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(legacy_path, new_path)
        return True
    except OSError:
        return False  # migration is best-effort; fresh state still works


def state_dir(name: str) -> tuple[Path, bool]:
    """Per-user state subdirectory ``<app_dir>/<name>``.

    Returns ``(path, migrated)`` where ``migrated`` reports whether legacy
    package-relative state was copied over on this call.
    """
    d = app_dir() / name
    already_populated = d.exists() and any(d.iterdir())
    d.mkdir(parents=True, exist_ok=True)
    migrated = False
    if not already_populated:
        migrated = _migrate_dir_contents(d, legacy_package_dir() / name)
    return d, migrated


def state_file(name: str, filename: str) -> tuple[Path, bool]:
    """Per-user state file ``<app_dir>/<name>/<filename>``.

    Returns ``(path, migrated)`` where ``migrated`` reports whether a legacy
    package-relative file was copied over on this call.
    """
    d, migrated_dir = state_dir(name)
    p = d / filename
    migrated = migrated_dir or _migrate_file(p, legacy_package_dir() / name / filename)
    return p, migrated
