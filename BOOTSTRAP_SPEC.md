# SPEC — PDB Metal-Site Index Bootstrap on Modal (detached, chunked)

Target: index the metal sites of 38,146 experimental PDB entries
(ZN, FE, HEM, SF4, FES, F3S map already fetched in
`data/pdb_index/metal_entries.json`) with the tested protein-tools QC
(metal_geometry v1: angle fitting, donor typing, bidentate + adjacent-metal
detection), producing `metal_sites.db` queryable by the MCP server.

Design principles:
1. **Detached** — runs on Modal; closing VS Code does not stop it.
2. **Checkpoint = volume state** — no separate state file to corrupt; every
   stage resumes by listing what exists on the volume and skipping it.
3. **Failures are per-entry, never per-run** — one bad structure cannot kill
   the batch; failed items are recorded with reasons and replayable alone.
4. **Shards before DB** — workers write JSONL shards; a single merge step
   builds SQLite (no concurrent SQLite writers over a network volume).

## 1. Volume layout (`metal-index-vol`)

```
input/metal_entries.json            # 333 KB, uploaded from local stage 1
subset/XXXX.cif.gz                  # downloaded entries (idempotent by id)
shards/build_v1/batch_000123.jsonl  # per-batch QC output, JSONL rows
shards/build_v1/batch_000123.ok     # empty marker written LAST (commit)
reports/sync_batch_000123.json      # per-batch sync report
reports/build_batch_000123.json     # per-batch build report
merged/metal_sites.db               # final SQLite (written once, at merge)
```

Shard row = one metal site:
`{entry_id, metal, site, cn, geometry, rmsd_deg, n_imidazole, n_carboxylate,
n_water, n_backbone, n_thiolate, n_donors, n_adjacent_metal, bidentate,
ligands, qc_version, batch_id}`

## 2. Modal app & functions

```python
app = modal.App("metal-index-bootstrap")
vol = modal.Volume.from_name("metal-index-vol", create_if_missing=True)

@app.function(image=img, volumes={"/vol": vol}, timeout=3600, retries=2)
def sync_batch(batch_id: str, ids: list[str]) -> dict
    # network-bound: 16 threads inside the container
    # skip files that already exist (idempotent); per-entry try/except
    # writes subset files + reports/sync_batch_<id>.json
    # backoff on 429/5xx; max_mb cap -> recorded as "too_large"

@app.function(image=img, volumes={"/vol": vol}, timeout=3600, retries=2,
              cpu=4, memory=16384)
def build_batch(batch_id: str) -> dict
    # CPU-bound: ProcessPool over ~250 files (spawn-safe module)
    # parses .cif.gz -> cif_tools -> pdb_utils.metal_geometry
    # writes shards JSONL + .ok marker atomically last
    # per-file exceptions -> error rows in the same shard (type:"error")

@app.function(image=img, volumes={"/vol": vol}, timeout=7200)
def merge_and_finalize(version: str) -> dict
    # single writer: reads all shards with .ok markers -> SQLite
    # creates indexes; runs validation queries; writes final report

@app.local_entrypoint()
def main(stage: str = "all", limit: int = 0):
    # stage: pilot | sync | build | merge | status | all
    # uploads metal_entries.json if missing; computes batches; fans out
```

Image: `debian_slim().pip_install("matplotlib")` + local source mounted
(`cif_tools`, `pdb_utils`, `builder_metal_index` core functions). QC logic
is imported, never duplicated.
## 3. Chunk math (38,146 entries)

| Stage | Batch size | # batches | Est. per batch | Fan-out | Wall est. |
|-------|-----------|-----------|----------------|---------|-----------|
| sync  | 500 ids   | 77        | 1-3 min (750 MB) | 20-40  | 15-60 min |
| build | 250 files | 153       | 2-4 min (cpu=4)  | 20-40  | 15-45 min |
| merge | -         | 1         | 2-5 min          | 1      | minutes   |

Realistic total wall time: **1-2 hours** of detached compute (vs the
optimistic-overnight single-process local run). Modal auto-scales
containers; a `pilot` mode exists for validation first.

## 4. Failure & replay matrix

| Failure | Handling | Replay cost |
|---------|----------|-------------|
| Entry download fails (network, 404, size) | recorded in batch report; batch continues | one entry |
| Batch function crashes | Modal retries (2x); timeout=1 h guards hangs | one batch |
| Parse/analyze exception | error row in shard with reason; run continues | one file |
| Merge crash | shards untouched; merge is idempotent (rebuilds DB) | minutes |
| QC version bump | shards live under `shards/build_v1/`; bump dir name | full rebuild (subset only) |

Debug loop = `status` (what's missing?) -> fix -> re-run stage -> only
unfinished work replays. No state file, no lost progress. The "guaranteed
multiple debug rounds" cost is bounded: each round replays only failures.

## 5. Merge, validation, acceptance criteria

Merge steps: collect `*.ok` shards -> executemany into `sites` (15 cols +
qc_version + batch_id) -> create indexes -> validation:
1. row count == sum of shard row counts
2. entries with zero sites are listed (some legitimately have none)
3. spot checks: 10AJ/10AK Cys4-tetrahedral present; no row with cn=0;
   geometry distribution sane per metal
4. failures report: every failure id + reason, grouped

Acceptance: >= 99% of entries synced-or-explicitly-failed; 0 unexplained
errors; spot checks pass; DB row count matches shards; final JSON report
on the volume + a local `log_run` registry entry.

## 6. Pilot before full launch

`pilot` stage: first 500 entries end-to-end (sync -> build -> merge),
validating: throughput, file-size distribution, error classes, per-batch
timing, Modal cost. Only after a clean pilot: full run.

## 7. Cost & operational notes

- Compute: ~40-80 core-hours total -> order $1-5 at typical Modal pricing
  (verify current pricing; bandwidth policy to confirm at launch)
- 38k small files on a volume is within limits but is the design point to
  watch; if slow, switch sync to tar-per-batch instead of loose files
- Tuned knobs: function timeouts (1 h) + retries (2); everything else is
  idempotency by volume layout
- After merge: download `metal_sites.db` locally, log it in the local run
  registry, expose `query_metal_sites` as an MCP tool

## 8. Out of scope (later phases)

- Incremental daily updates (list refresh -> ~40 new ids -> same pipeline)
- CA/MG/transition-metal second wave (same pipeline, new metals)
- Re-analysis on QC version bump (bump shard dir, replay)

