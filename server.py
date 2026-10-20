"""ProteinTools — an MCP server wrapping protein-characterization utilities.

Free tier: local structural tools (PDB validation, metal-site geometry,
sequence extraction, pLDDT summary). No account needed.

Premium tier: network lookups (UniProt, AlphaFold DB, RCSB) gated by an
API key (PROTEIN_TOOLS_API_KEY env var). In production these calls should be
proxied through your own relay service so the key does billing/quota
enforcement — see README.md for the monetization notes.
"""
from __future__ import annotations

import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

from mcp.server.fastmcp import FastMCP

import cif_tools
import ff_tools
import md_report
import pdb_utils
import plot_tools
import premium_tools
import propka_tools
import run_registry
import seq_utils
import xvg_tools
import app_paths  # per-user state dir (registry DB, premium cache)

mcp = FastMCP("protein-tools")


def _metal_index():
    """Lazy import of the local metal-site index module.

    The index builder is not part of the open-source free core, so it is
    imported on demand: without it, the two index tools explain that the
    index is available through the hosted relay instead of crashing the
    server at startup.
    """
    try:
        import builder_metal_index
        return builder_metal_index
    except ImportError as exc:
        raise RuntimeError(
            "the PDB-wide metal-site index is a premium feature served "
            "by the hosted relay - set PROTEIN_TOOLS_RELAY_URL and "
            "PROTEIN_TOOLS_API_KEY (see README 'Hosted relay'). "
            f"(original error: {exc})") from exc


PREMIUM_ENV = "PROTEIN_TOOLS_API_KEY"
RELAY_ENV = "PROTEIN_TOOLS_RELAY_URL"


def _relay_url() -> str:
    """When set, premium/index tools are served by the hosted relay
    (auth via PROTEIN_TOOLS_API_KEY) instead of direct upstream calls and
    the local index DB. This is the production/monetized mode."""
    return os.environ.get(RELAY_ENV, "").strip().rstrip("/")


def _load_text(pdb_path: str) -> str:
    return Path(pdb_path).read_text(errors="replace")


def _load_structure(path: str):
    """Load a PDB file as text, or an mmCIF file as a parsed structure dict.

    Format is detected by .cif extension or a leading 'data_' record, so
    AlphaFold 3 / modern RCSB files work with the same tools.
    """
    text = _load_text(path)
    if path.lower().endswith(".cif") or text.lstrip().startswith("data_"):
        return cif_tools.to_pdb_atoms(cif_tools.parse_cif(text))
    return text


def _j(obj) -> str:
    return json.dumps(obj, indent=2, ensure_ascii=False, default=str)


def _check_premium_key(tool: str) -> str:
    key = os.environ.get(PREMIUM_ENV, "").strip()
    if not key:
        raise PermissionError(
            f"'{tool}' is a premium tool and requires an API key. Set the "
            f"PROTEIN_TOOLS_API_KEY environment variable (see server README "
            f"for how to get a key)."
        )
    return key


def _relay_call(path: str, params: dict | None = None) -> dict:
    """GET the relay with the premium key; raises on non-200 with the
    relay's structured error message."""
    key = _check_premium_key("relay")
    url = _relay_url() + path
    if params:
        url += "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(
        url, headers={"Authorization": f"Bearer {key}",
                      "User-Agent": "protein-tools-mcp/0.1"})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return json.loads(r.read().decode())
    except urllib.error.HTTPError as exc:
        try:
            detail = json.loads(exc.read().decode())
        except Exception:
            detail = {}
        raise PermissionError(detail.get("error", f"relay HTTP {exc.code}"))


def _error_response(exc: Exception) -> str:
    """Structured error for agents: they can read the hint and self-correct."""
    hints = {
        "FileNotFoundError": "check the file path (use an absolute path)",
        "PermissionError": "this tool needs an API key - set PROTEIN_TOOLS_API_KEY",
        "ValueError": "the file content could not be parsed as requested",
        "IsADirectoryError": "the path points to a directory, not a file",
        "TimeoutError": "upstream database timed out; try again",
    }
    return _j({"error": f"{type(exc).__name__}: {exc}",
               "hint": hints.get(type(exc).__name__,
                                 "if the problem persists, report it with this message")})


def _fetch(url: str, binary: bool = False, retries: int = 2):
    last: Exception | None = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "protein-tools-mcp/0.1"})
            with urllib.request.urlopen(req, timeout=30) as r:
                data = r.read()
            return data if binary else data.decode()
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last = exc
            if attempt < retries:
                time.sleep(2 * (attempt + 1))  # brief backoff, then retry
    raise last


CACHE_DIR, _CACHE_MIGRATED = app_paths.state_dir("cache")
CACHE_TTL_HOURS = 24.0


def _cached_fetch(url: str, binary: bool = False, refresh: bool = False,
                  ttl_hours: float = CACHE_TTL_HOURS,
                  retries: int = 2):
    """Disk-cached GET for premium lookups: faster repeat calls, offline
    fallback within the TTL, and quota-friendly against upstream APIs."""
    import hashlib
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    key = hashlib.md5(url.encode("utf-8")).hexdigest()
    blob = CACHE_DIR / (key + ".cache")
    meta = CACHE_DIR / (key + ".meta")
    if not refresh and blob.exists() and meta.exists():
        try:
            if time.time() - float(meta.read_text().strip()) < ttl_hours * 3600:
                raw = blob.read_bytes()
                return raw if binary else raw.decode()
        except (OSError, ValueError):
            pass  # corrupt cache entry: refetch
    data = _fetch(url, binary=True, retries=retries)
    try:
        blob.write_bytes(data)
        meta.write_text(str(time.time()))
    except OSError:
        pass  # cache is best-effort; never block a successful fetch
    return data if binary else data.decode()


# ----------------------------- FREE TIER -----------------------------------

@mcp.tool()
def validate_pdb(pdb_path: str) -> str:
    """Validate a PDB or mmCIF structure file (AlphaFold 3 / RCSB) for
    downstream use (MD prep, docking, analysis).

    Reports atom/residue/chain counts, fixed-column format errors,
    non-standard residues, altLocs, chain gaps, duplicate atom names,
    insertion-code residues, and multi-model (NMR) notes.
    Always run this before feeding a structure into a pipeline.
    """
    try:
        return _j(pdb_utils.validate_pdb(_load_structure(pdb_path)))
    except Exception as exc:
        return _error_response(exc)


@mcp.tool()
def metal_geometry(pdb_path: str, metal_element: str = "FE",
                   cutoff_angstrom: float = 3.0) -> str:
    """Analyze metal-binding sites in a PDB file.

    Finds all atoms of the given metal element (e.g. FE, ZN, MN, MG), lists
    coordinating ligand atoms (N/O/S donors) with distances and donor types,
    computes ligand-metal-ligand angles fitted to ideal geometries
    (octahedral, tetrahedral, ...), and flags bidentate carboxylate chelates.
    """
    try:
        return _j(pdb_utils.metal_geometry(_load_structure(pdb_path), metal_element,
                                           cutoff_angstrom))
    except Exception as exc:
        return _error_response(exc)


@mcp.tool()
def structure_stats(pdb_path: str) -> str:
    """General structure summary: atoms, residues, chains, gaps, altLocs.

    Accepts PDB or mmCIF (auto-detected). Free; local compute.
    """
    try:
        return _j(pdb_utils.structure_stats(_load_structure(pdb_path)))
    except Exception as exc:
        return _error_response(exc)


@mcp.tool()
def sequence_from_pdb(pdb_path: str) -> str:
    """Extract one-letter amino-acid sequences per chain from a PDB file
    (handles MD protonation variants; for multi-model files uses model 1)."""
    try:
        return _j(pdb_utils.sequences(_load_structure(pdb_path)))
    except Exception as exc:
        return _error_response(exc)


@mcp.tool()
def plddt_summary(pdb_path: str) -> str:
    """Summarize AlphaFold-style pLDDT confidence from B-factor columns.

    Reports mean, confidence bands (>90/70-90/50-70/<50), and low-confidence
    regions. Refuses to interpret experimental-structure B-factors as pLDDT.
    """
    try:
        return _j(pdb_utils.plddt_summary(_load_structure(pdb_path)))
    except Exception as exc:
        return _error_response(exc)


@mcp.tool()
def xvg_summary(xvg_path: str) -> str:
    """Summarize any GROMACS .xvg file (RMSD, energies, temperature, SASA...):
    per-column min/max/mean/stdev/first/last plus title and axis labels."""
    try:
        return _j(xvg_tools.summary(_load_text(xvg_path)))
    except Exception as exc:
        return _error_response(exc)


@mcp.tool()
def rmsd_summary(xvg_path: str) -> str:
    """Assess equilibration from a GROMACS RMSD .xvg file: drift between the
    first and last 10% of frames, plateau detection, and an equilibrated /
    not-equilibrated verdict."""
    try:
        return _j(xvg_tools.rmsd_summary(_load_text(xvg_path)))
    except Exception as exc:
        return _error_response(exc)


@mcp.tool()
def rmsf_peaks(xvg_path: str, top_n: int = 10) -> str:
    """Extract the most flexible residues/atoms from a GROMACS RMSF .xvg
    file: top-N highest values and contiguous flexible segments
    (above mean + 1 stdev)."""
    try:
        return _j(xvg_tools.rmsf_peaks(_load_text(xvg_path), top_n))
    except Exception as exc:
        return _error_response(exc)


# ---------------------------- PREMIUM TIER ---------------------------------

@mcp.tool()
def conservation_hint(accession_or_id: str) -> str:
    """Conservation / functional hint for a UniProt entry: feature table.

    Reuses the premium UniProt fetch (same key + relay/local path as
    uniprot_entry): feature types grouped (Metal binding / Active site /
    Binding site / Disulfide bond / Mutagenesis / Natural variant) with
    1-based canonical positions. Align your PDB sequence first
    (sequence_from_pdb + read_fasta) before mapping positions.
    """
    try:
        relay = _relay_url()
        if relay:
            return _j(_relay_call("/v1/uniprot",
                                  {"accession": accession_or_id.strip()}))
        _check_premium_key("conservation_hint")
        raw = json.loads(uniprot_entry(accession_or_id))
        if "error" in raw:
            return _j(raw)
        feats = raw.get("features_by_type", {})
        n_feat = sum(len(v) for v in feats.values())
        return _j({
            "accession": raw.get("accession"),
            "name": raw.get("name"),
            "organism": raw.get("organism"),
            "length": raw.get("length"),
            "features_by_type": feats,
            "n_feature_types": len(feats),
            "n_features": n_feat,
            "how_to_use": ("positions are 1-based UniProt canonical; align "
                           "your PDB sequence first (sequence_from_pdb + "
                           "read_fasta) before mapping positions onto "
                           "structure residues."),
        })
    except Exception as exc:
        return _error_response(exc)


@mcp.tool()
def uniprot_entry(accession: str, refresh: bool = False) -> str:
    """Fetch a UniProtKB entry: protein name, organism, function, length,
    sequence, and metal-binding annotations. Premium tool (API key required).
    Results are disk-cached for 24 h; refresh=True forces a live fetch."""
    try:
        _check_premium_key("uniprot_entry")
        if _relay_url():
            return _j(_relay_call(f"/v1/uniprot/{accession}"))
        d = json.loads(_cached_fetch(
            f"https://rest.uniprot.org/uniprotkb/{accession}.json",
            refresh=refresh))
        return _j(premium_tools.uniprot_summary(d))
    except Exception as exc:
        return _error_response(exc)


@mcp.tool()
def alphafold_plddt(accession: str, refresh: bool = False) -> str:
    """Fetch the AlphaFold DB model for a UniProt accession and summarize
    its pLDDT confidence profile. Premium tool (API key required).
    Disk-cached for 24 h; refresh=True forces a live fetch."""
    try:
        _check_premium_key("alphafold_plddt")
        if _relay_url():
            return _j(_relay_call(f"/v1/alphafold/{accession}"))
        text = _cached_fetch(
            f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v4.pdb",
            refresh=refresh)
        return _j(pdb_utils.plddt_summary(text))
    except Exception as exc:
        return _error_response(exc)


@mcp.tool()
def rcsb_entry(pdb_id: str, refresh: bool = False) -> str:
    """Fetch RCSB PDB entry metadata: title, experimental method, resolution,
    deposit date. Premium tool (API key required). Disk-cached for 24 h;
    refresh=True forces a live fetch."""
    try:
        _check_premium_key("rcsb_entry")
        if _relay_url():
            return _j(_relay_call(f"/v1/rcsb/{pdb_id}"))
        d = json.loads(_cached_fetch(
            f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}",
            refresh=refresh))
        return _j(premium_tools.rcsb_summary(d))
    except Exception as exc:
        return _error_response(exc)


# ------------------------- MD LIFECYCLE TOOLS -------------------------------

@mcp.tool()
def check_metal_params(ff_path: str, metals: str = "ZN,FE") -> str:
    """Scan a GROMACS force-field directory or .itp/.ff file for metal
    atomtypes (e.g. ZN, FE) and report sigma/epsilon or MISSING entries.

    Run this BEFORE pdb2gmx on any metalloprotein - a missing metal atomtype
    will kill the topology build.
    """
    try:
        metals_list = [m.strip() for m in metals.split(",") if m.strip()]
        return _j(ff_tools.check_metal_params(ff_path, metals=metals_list))
    except Exception as exc:
        return _error_response(exc)


@mcp.tool()
def md_preflight(pdb_path: str, ff_dir: str = "") -> str:
    """MD pre-flight gate for a structure: PDB/mmCIF validation + metal-site
    detection + force-field metal parameter check. Returns GO or NO-GO with
    blocking issues. Treat NO-GO as a hard stop before any GPU/MD launch."""
    try:
        return _j(ff_tools.preflight(_load_structure(pdb_path),
                                     ff_dir=ff_dir or None))
    except Exception as exc:
        return _error_response(exc)


@mcp.tool()
def propka_summary(pka_path: str) -> str:
    """Parse a PROPKA .pka output: pI, per-residue-type counts, pKa extremes,
    and protonation-state concerns (HIS near pH 7, raised acids, CYS 99.99
    disulfide markers) - exactly the decisions pdb2gmx protonation needs."""
    try:
        return _j(propka_tools.summarize(_load_text(pka_path)))
    except Exception as exc:
        return _error_response(exc)


# ------------------------- REGISTRY / RUN REPORTS ---------------------------

@mcp.tool()
def log_run(kind: str, run_id: str = "", status: str = "logged",
            notes: str = "", metrics: str = "", paths: str = "") -> str:
    """Record a run/event in the persistent project registry (cross-session
    memory). kind: e.g. 'md_report', 'prep', 'hypothesis'. metrics/paths
    accept JSON strings or plain text. Returns the registry id."""
    try:
        def _parse(x):
            if not x:
                return None
            try:
                return json.loads(x)
            except (ValueError, TypeError):
                return x
        return _j(run_registry.log_run(kind=kind, run_id=run_id or None,
                                       status=status, notes=notes,
                                       metrics=_parse(metrics),
                                       paths=_parse(paths)))
    except Exception as exc:
        return _error_response(exc)


@mcp.tool()
def list_runs(kind: str = "", limit: int = 20) -> str:
    """List recent registry entries (newest first), optionally filtered by
    kind. Use this at session start to recall project history."""
    try:
        return _j(run_registry.list_runs(kind=kind or None, limit=limit))
    except Exception as exc:
        return _error_response(exc)


@mcp.tool()
def run_details(rid: str) -> str:
    """Full record for one registry entry, by numeric id or run_id."""
    try:
        return _j(run_registry.run_details(rid))
    except Exception as exc:
        return _error_response(exc)


@mcp.tool()
def md_run_report(run_dir: str, log: bool = True) -> str:
    """Analyze a GROMACS run directory: RMSD verdicts, other .xvg inventory,
    md.log crash scan, artifact presence, and an overall verdict (GO /
    FINISHED / INCOMPLETE / CRASHED). Optionally logs to the registry."""
    try:
        return _j(md_report.md_run_report(run_dir, log=log))
    except Exception as exc:
        return _error_response(exc)


@mcp.tool()
def compare_runs(dir_a: str, dir_b: str) -> str:
    """Side-by-side RMSD verdicts for two run directories (e.g. apo vs
    holo): plateau status, late-run means, drift, and relative difference."""
    try:
        return _j(md_report.compare_runs(dir_a, dir_b))
    except Exception as exc:
        return _error_response(exc)


@mcp.tool()
def cache_info(clear: bool = False) -> str:
    """Inspect (or with clear=True, empty) the local response cache used by
    the premium lookup tools (UniProt / AlphaFold DB / RCSB)."""
    try:
        entries, size = 0, 0
        if CACHE_DIR.exists():
            for f in CACHE_DIR.iterdir():
                entries += 1
                size += f.stat().st_size
        if clear:
            for f in CACHE_DIR.iterdir():
                f.unlink()
            entries, size = 0, 0
        return _j({"cache_dir": str(CACHE_DIR), "entries": entries,
                   "size_bytes": size, "ttl_hours": CACHE_TTL_HOURS,
                   "cleared": clear})
    except Exception as exc:
        return _error_response(exc)


# ------------------------- METAL-SITE INDEX -------------------------------

INDEX_DB_DEFAULT = str(Path(__file__).resolve().parent
                        / "data" / "index" / "metal_sites.db")


@mcp.tool()
def query_metal_sites(metal: str, geometry: str = "", min_cn: int = 2,
                      max_cn: int = 12, limit: int = 20,
                      include_artifacts: bool = False,
                      db_path: str = "") -> str:
    """Search the PDB-wide metal-site index (premium: requires the hosted
    relay — set PROTEIN_TOOLS_RELAY_URL + PROTEIN_TOOLS_API_KEY).

    Filter by metal, geometry class, coordination-number range. Defaults
    exclude CN-0 isolated-ion artifacts and CN-1; set min_cn=1 to see
    unusual CN-1 chemistry, include_artifacts=True to also return CN-0
    rows. Examples: tetrahedral Cys4 zinc fingers (metal='ZN',
    geometry='tetrahedral'), catalytic zinc (metal='ZN', min_cn=4,
    max_cn=6)."""
    try:
        if _relay_url():
            return _j(_relay_call("/v1/query_metal_sites", {
                "metal": metal, "geometry": geometry or "",
                "min_cn": min_cn, "max_cn": max_cn, "limit": limit,
                "include_artifacts": str(include_artifacts).lower()}))
        return _j(_metal_index().query(
            metal, geometry=geometry or None, min_cn=min_cn, max_cn=max_cn,
            limit=limit, include_artifacts=include_artifacts,
            db_path=db_path or INDEX_DB_DEFAULT))
    except Exception as exc:
        return _error_response(exc)


@mcp.tool()
def index_stats(db_path: str = "") -> str:
    """Corpus-level stats for the PDB-wide metal-site index (premium:
    requires the hosted relay): entry/site counts, per-metal and
    per-geometry distributions, artifact volume."""
    try:
        if _relay_url():
            return _j(_relay_call("/v1/index_stats"))
        return _j(_metal_index().status(db_path=db_path or INDEX_DB_DEFAULT))
    except Exception as exc:
        return _error_response(exc)


# --------------------------- SEQUENCE ANALYSIS ------------------------------

@mcp.tool()
def aa_composition(source: str, is_file: bool = True,
                   sequence_index: int = 0) -> str:
    """Protein physicochemistry: MW, aa composition, charge, GRAVY, pI estimate.

    Parses a FASTA file or raw text (record `sequence_index`); rejects
    non-protein sequences with a structured error naming the detected type.
    Free; local compute.
    """
    try:
        records = seq_utils.parse_fasta(source, is_file=is_file)
        if not (0 <= sequence_index < len(records)):
            raise IndexError(
                f"sequence_index {sequence_index} out of range "
                f"(0..{len(records) - 1}; {len(records)} record(s) found)")
        rec = records[sequence_index]
        comp = seq_utils.aa_composition(rec["sequence"])
        return _j({"id": rec["id"], "detected_type": rec["type"],
                   **comp})
    except Exception as exc:
        return _error_response(exc)


@mcp.tool()
def fasta_stats(source: str, is_file: bool = True) -> str:
    """Multi-FASTA summary: N50, length stats, per-type counts, GC%.

    GC% is computed over DNA records only (null with a note when none).
    Free; local compute.
    """
    try:
        return _j(seq_utils.fasta_stats(source, is_file=is_file))
    except Exception as exc:
        return _error_response(exc)


@mcp.tool()
def detect_sequence_type(source: str, is_file: bool = True, sequence_index: int = 0) -> str:
    """Detect whether a FASTA sequence is protein, DNA, RNA, ambiguous, or unknown.

    Accepts a FASTA file path (is_file=True) or raw FASTA text (is_file=False).
    Returns the type classification with confidence and detailed notes.
    """
    try:
        records = seq_utils.parse_fasta(source, is_file=is_file)
        if not records:
            raise ValueError("No FASTA records found")
        if sequence_index < 0 or sequence_index >= len(records):
            raise ValueError(f"sequence_index {sequence_index} out of range (have {len(records)} sequences)")
        rec = records[sequence_index]
        all_types = [{"id": r["id"], "type": r["type"]["type"], "confidence": r["type"]["detail"]["confidence"], "note": r["type"]["detail"]["note"]} for r in records]
        return _j({
            "source": source,
            "is_file": is_file,
            "n_sequences": len(records),
            "inspected_index": sequence_index,
            "inspected_id": rec["id"],
            "inspected_description": rec["description"],
            "type": rec["type"],
            "all_sequences_types": all_types,
        })
    except Exception as exc:
        return _error_response(exc)

@mcp.tool()
def read_fasta(source: str, is_file: bool = True, sequence_index: int = 0, normalize: bool = True) -> str:
    """Read one or more sequences from a FASTA file or text.

    Accepts a FASTA file path (is_file=True) or raw FASTA text (is_file=False).
    Returns the full parsed records including sequences, types, and metadata.
    """
    try:
        records = seq_utils.parse_fasta(source, is_file=is_file)
        if not records:
            raise ValueError("No FASTA records found")
        if sequence_index < 0 or sequence_index >= len(records):
            raise ValueError(f"sequence_index {sequence_index} out of range (have {len(records)} sequences)")
        rec = records[sequence_index]
        seq = rec["sequence"]
        if normalize:
            seq = seq_utils.normalize_sequence(seq)
        return _j({
            "source": source,
            "is_file": is_file,
            "n_sequences": len(records),
            "inspected_index": sequence_index,
            "inspected_id": rec["id"],
            "inspected_description": rec["description"],
            "type": rec["type"],
            "sequence": seq,
            "maybe_uniprot_accession": rec["maybe_uniprot_accession"],
            "length": len(seq),
        })
    except Exception as exc:
        return _error_response(exc)


@mcp.tool()
def pdb_by_id(pdb_id: str, refresh: bool = False) -> str:
    """Get info about a PDB entry by its 4-character PDB ID.

    Checks if a local PDB/mmCIF file exists for the ID, and if so returns
    local structure info (validation, sequence, metal sites). Also fetches
    RCSB metadata if an API key is set (premium tool).
    """
    try:
        pdb_id = pdb_id.strip().upper()
        if len(pdb_id) != 4 or not pdb_id.isalpha():
            raise ValueError(f"pdb_id must be a 4-letter alphabetic ID, got '{pdb_id}'")

        base_dir = Path(__file__).resolve().parent
        local_info: dict = {"local_file": None, "present": False}
        for ext in [".pdb", ".cif"]:
            candidate = base_dir / f"{pdb_id}{ext}"
            if candidate.exists():
                local_info = {
                    "local_file": str(candidate),
                    "present": True,
                    "format": ext[1:],
                }
                try:
                    atoms = _load_structure(str(candidate))
                    local_info["validate"] = pdb_utils.validate_pdb(atoms)
                    local_info["sequence"] = pdb_utils.sequences(atoms)
                    local_info["metal_sites"] = pdb_utils.metal_sites(atoms)
                except Exception as e:
                    local_info["load_error"] = str(e)
                break

        rcsb_info = None
        try:
            _check_premium_key("pdb_by_id")
            if _relay_url():
                rcsb_info = _j(_relay_call(f"/v1/rcsb/{pdb_id}"))
            else:
                d = json.loads(_cached_fetch(
                    f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}",
                    refresh=refresh))
                rcsb_info = premium_tools.rcsb_summary(d)
        except PermissionError:
            pass

        return _j({"pdb_id": pdb_id, "local": local_info, "rcsb": rcsb_info})
    except Exception as exc:
        return _error_response(exc)

@mcp.tool()
def find_pdb_by_sequence(source: str, is_file: bool = True, sequence_index: int = 0, max_results: int = 10) -> str:
    """Find PDB structures matching a protein, DNA, or RNA sequence.

    Parses the FASTA, maps the detected type to the RCSB sequence search
    (protein→blastp, dna/rna→blastn), and returns matching PDB IDs.
    Ambiguous/unknown sequences are rejected with a structured error.
    Requires an API key (premium tool).
    """
    try:
        records = seq_utils.parse_fasta(source, is_file=is_file)
        if not records:
            raise ValueError("No FASTA records found")
        if sequence_index < 0 or sequence_index >= len(records):
            raise ValueError(f"sequence_index {sequence_index} out of range (have {len(records)} sequences)")
        rec = records[sequence_index]

        if rec["type"]["type"] not in ("protein", "dna", "rna"):
            raise ValueError(f"Cannot search by {rec['type']['type']} sequence ({rec['type'].get('detail', {}).get('note', 'unresolved type')}); provide a clean protein, DNA, or RNA sequence")

        seq = rec["sequence"]
        if len(seq) < 10:
            raise ValueError(f"Sequence too short ({len(seq)} aa) for meaningful search")

        _check_premium_key("find_pdb_by_sequence")

        search_url = "https://search.rcsb.org/rcsbsearch/v2/search"
        seq_type_map = {"protein": ("protein", "blastp"),
                        "dna": ("dna", "blastn"),
                        "rna": ("rna", "blastn")}
        seq_type, program = seq_type_map[rec["type"]["type"]]
        payload = {
            "query": {
                "type": "terminal",
                "service": "sequence",
                "parameters": {
                    "evalue_cutoff": 1e-3,
                    "sequence_type": seq_type,
                    "target": "pdb",
                    "program": program,
                    "sequence": seq,
                },
            },
            "request_options": {
                "results_content_type": ["entry"],
                "pagination": {"start": 0, "rows": max_results},
            },
            "result_options": {
                "sort": [{"sort_by": "score", "direction": "desc"}],
            },
        }

        req = urllib.request.Request(
            search_url,
            data=json.dumps(payload).encode(),
            headers={
                "Content-Type": "application/json",
                "User-Agent": "protein-tools-mcp/0.1",
            })
        with urllib.request.urlopen(req, timeout=60) as r:
            search_result = json.loads(r.read().decode())

        groups = search_result.get("groups", [])
        entries = []
        for group in groups[:max_results]:
            entry_id = group.get("identifier")
            if entry_id:
                entries.append({"pdb_id": entry_id})

        return _j({
            "query_id": rec["id"],
            "query_type": rec["type"],
            "query_length": len(seq),
            "n_results": len(entries),
            "entries": entries,
            "search_url": search_url,
        })
    except PermissionError as exc:
        return _error_response(exc)
    except Exception as exc:
        return _error_response(exc)

@mcp.tool()
def plot_xvg(xvg_path: str, out_png: str = "", title: str = "",
             y_from_zero: bool = False) -> str:
    """Plot every series of a GROMACS .xvg file against column 0 into a PNG
    (large series are decimated to <=5000 points). Returns the output path."""
    try:
        return _j(plot_tools.plot_xvg(xvg_path, out_png=out_png, title=title,
                                      y_from_zero=y_from_zero))
    except Exception as exc:
        return _error_response(exc)


@mcp.tool()
def compare_plot(xvg_paths: str, out_png: str = "", labels: str = "",
                 value_column: int = 1, title: str = "",
                 y_from_zero: bool = False) -> str:
    """Overlay the value column of several .xvg files (comma-separated
    paths, e.g. apo vs holo RMSD) into one comparison PNG."""
    try:
        return _j(plot_tools.compare_plot(xvg_paths, out_png=out_png,
                                          labels=labels,
                                          value_column=value_column,
                                          title=title,
                                          y_from_zero=y_from_zero))
    except Exception as exc:
        return _error_response(exc)


def main() -> None:
    """Console entry point (``protein-tools`` script): run the MCP server."""
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
