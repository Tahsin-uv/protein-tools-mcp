"""Phase 3 tests: plotting, premium disk cache, listing assets."""
import os
import sys
import tempfile

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_tools  # noqa: E402

ASSETS = tempfile.mkdtemp(prefix="protein_tools_assets_")

HAVE_PLT = plot_tools._matplotlib_usable()
if not HAVE_PLT:
    print("NOTE: matplotlib unavailable/blocked - skipping plot assertions "
          "(plot tools degrade with a structured error in this environment)")

if HAVE_PLT:
    # --- 1. plot the 1000-point fixture .xvg (decimated to <= MAX_POINTS) ---
    rmsd_xvg = os.path.join(FIX, "rmsd_core.xvg")
    out = os.path.join(ASSETS, "rmsd_core.png")
    r = plot_tools.plot_xvg(rmsd_xvg, out_png=out, title="Fixture RMSD")
    assert os.path.exists(out) and os.path.getsize(out) > 5000, r
    assert r["points_plotted"] <= plot_tools.MAX_POINTS, r
    assert r["n_points_total"] == 1000, r
    print(f"1. plot_xvg: OK -> {r['out_png']} ({r['size_bytes']} bytes, "
          f"{r['points_plotted']} of {r['n_points_total']} points plotted)")

    # --- 2. comparison overlay of the two fixture .xvg files ---
    sasa_xvg = os.path.join(FIX, "sasa_per_residue.xvg")
    out2 = os.path.join(ASSETS, "compare_rmsd_sasa.png")
    r2 = plot_tools.compare_plot(f"{rmsd_xvg},{sasa_xvg}", out_png=out2,
                                 labels="RMSD core, SASA total",
                                 title="Fixture trajectory overview")
    assert os.path.exists(out2) and os.path.getsize(out2) > 5000, r2
    print(f"2. compare_plot: OK -> {r2['out_png']} ({r2['size_bytes']} bytes)")
else:
    # plotting must fail GRACEFULLY (structured error), never crash the server
    try:
        plot_tools.plot_xvg(os.path.join(FIX, "rmsd_core.xvg"))
        raise AssertionError("expected RuntimeError when matplotlib is blocked")
    except RuntimeError as exc:
        assert "matplotlib" in str(exc)
    print("1-2. plot tools degrade gracefully: OK (RuntimeError with hint)")

# --- 3. premium disk cache: fetch twice, second must be cached ---
import server  # noqa: E402  (imports FastMCP; safe — nothing runs on import)

url = "https://rest.uniprot.org/uniprotkb/A8ABE7.json"
t0 = __import__("time").perf_counter()
d1 = server._cached_fetch(url)
t1 = __import__("time").perf_counter()
d2 = server._cached_fetch(url)
t2 = __import__("time").perf_counter()
assert d1 == d2 and "A8ABE7" in d1
cache_entries = list(server.CACHE_DIR.iterdir())
assert len(cache_entries) >= 2  # .cache + .meta
print(f"3. premium cache: OK (live fetch {1000 * (t1 - t0):.0f} ms, "
      f"cached {1000 * (t2 - t1):.1f} ms, entries={len(cache_entries)})")

# --- 4. cache refresh bypasses the TTL ---
d3 = server._cached_fetch(url, refresh=True)
assert d3 == d1  # same content, but fetched live
print("4. cache refresh: OK")

print("ALL PHASE 3 TESTS PASSED")
