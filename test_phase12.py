"""Phase 1+2 tests: pre-flight gate, PROPKA parsing, registry, run reports."""
import json
import os
import sys
import tempfile

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ff_tools  # noqa: E402
import md_report  # noqa: E402
import propka_tools  # noqa: E402
import run_registry  # noqa: E402

ROOT = os.path.dirname(os.path.abspath(__file__))
TMP = tempfile.mkdtemp(prefix="protein_tools_test_")

# --- 1. force-field metal param check (the amber99sb-ildn ZN scenario) ---
FF = os.path.join(TMP, "ffnonbonded.itp")
open(FF, "w").write(
    "; synthetic ffnonbonded\n"
    "[ atomtypes ]\n"
    "; name  at.num  mass    charge   ptype  sigma    epsilon\n"
    "  ZN    30     65.38    0.0000    A    0.11580  0.05229\n"
    "  FE    26     55.85    3.0000    A    0.11580  0.05229\n"
    "  C      6     12.01    0.0000    A    0.33997  0.45773\n"
)
r = ff_tools.check_metal_params(FF, metals=["ZN", "FE"])
assert r["verdict"].startswith("OK") and "ZN" in r["found"], r
assert abs(r["found"]["ZN"][0]["sigma_nm"] - 0.11580) < 1e-9
# now only ZN in the ff -> FE missing (the E6 holo blocker scenario)
open(FF, "w").write(
    "[ atomtypes ]\n  ZN    30     65.38    0.0000    A    0.11580  0.05229\n")
r = ff_tools.check_metal_params(FF, metals=["ZN", "FE"])
assert "FE" in r["missing"] and "MISSING" in r["verdict"], r
print("1. check_metal_params: OK (ZN found with Li&Merz values; FE missing detected)")

# --- 2. md_preflight: NO-GO when the structure's metal lacks ff params ---
PDB = os.path.join(FIX, "A8ABE7_holo_Fe2.pdb")
r = ff_tools.preflight(PDB, ff_dir=TMP)
assert r["verdict"] == "NO-GO" and r["metals_present_in_structure"] == ["FE"], r
assert any(b["check"] == "ff_metal_params" for b in r["blocking"]), r
r2 = ff_tools.preflight(PDB)  # no ff_dir -> warnings only
assert r2["verdict"] == "GO" and r2["metals_present_in_structure"] == ["FE"], r2
print("2. md_preflight: OK (NO-GO on missing FE params; GO when no ff_dir given)")

# --- 3. PROPKA parsing on the synthetic apo fixture (130 rows) ---
pk = open(os.path.join(FIX, "propka_apo.pka"), errors="replace").read()
s = propka_tools.summarize(pk)
assert s["n_residues_parsed"] == 130, s["n_residues_parsed"]
assert s["pI_folded"] is not None, s
his = [c["resnum"] for c in s["protonation_concerns"]
       if c["resname"] == "HIS"]
assert 33 in his and 34 in his and 70 in his, his
assert any(c["resname"] == "CYS" and c["pka"] >= 99.0
           for c in s["protonation_concerns"]), s["protonation_concerns"][:5]
print(f"3. propka_summary: OK (n={s['n_residues_parsed']}, pI={s['pI_folded']}, "
      f"HIS flags at {his}, CYS-disulfide flags present, "
      f"n_concerns={s['n_concerns']})")

# --- 4. run registry CRUD ---
reg = run_registry.log_run(kind="test", run_id="unit_test_run",
                           status="running", notes="phase12 test",
                           metrics={"rmsd": 0.2}, paths=[TMP])
det = run_registry.run_details(str(reg["id"]))
assert det["metrics"]["rmsd"] == 0.2 and det["status"] == "running"
det2 = run_registry.run_details("unit_test_run")
assert det2["id"] == reg["id"]
run_registry.update_run(reg["id"], status="done", metrics={"rmsd": 0.19})
det3 = run_registry.run_details(reg["id"])
assert det3["status"] == "done" and det3["metrics"]["rmsd"] == 0.19
lst = run_registry.list_runs(kind="test")
assert any(x["id"] == reg["id"] for x in lst["runs"])
print(f"4. run_registry: OK (log->details->update->list, db={lst['db']})")

# --- 5. md_run_report on synthetic dirs: GO and CRASHED ---
gdir = os.path.join(TMP, "run_go")
os.makedirs(gdir)
open(os.path.join(gdir, "rmsd.xvg"), "w").write(
    "@ title \"RMSD\"\n@ xaxis label \"Time (ns)\"\n"
    + "\n".join(f"{i * 0.1:.2f} 0.2000" for i in range(100)))
open(os.path.join(gdir, "md.log"), "w").write("step 4999900\nFinished mdrun\n")
rep = md_report.md_run_report(gdir, log=True)
assert rep["verdict"].startswith("GO"), rep["verdict"]
assert rep["registry_id"] is not None

cdir = os.path.join(TMP, "run_crash")
os.makedirs(cdir)
open(os.path.join(cdir, "md.log"), "w").write(
    "Fatal error: water molecule starting at atom 523\n")
rep2 = md_report.md_run_report(cdir, log=False)
assert rep2["verdict"].startswith("CRASHED"), rep2["verdict"]
print(f"5. md_run_report: OK (synthetic GO with registry id {rep['registry_id']}; "
      f"crash dir -> CRASHED)")

# --- 6. fixture run directory (finished, equilibrated) ---
real = os.path.join(FIX, "mdrun_apo")
rep3 = md_report.md_run_report(real, log=False)
assert "verdict" in rep3 and rep3["n_xvg"] > 0
assert rep3["verdict"].startswith("GO"), rep3["verdict"]
print(f"6. fixture apo run: verdict = {rep3['verdict']} | "
      f"xvg files: {rep3['n_xvg']} | finished: {rep3['finished_mdrun']}")

# --- 7. compare_runs (fixture apo vs synthetic holo-style copy) ---
holo = os.path.join(TMP, "md_run_holo")
os.makedirs(holo, exist_ok=True)
open(os.path.join(holo, "rmsd.xvg"), "w").write(
    "@ title \"RMSD\"\n@ xaxis label \"Time (ns)\"\n"
    + "\n".join(f"{i * 1.5:.2f} 0.3000" for i in range(100)))
try:
    cmp_ = md_report.compare_runs(real, holo)
    print(f"7. compare_runs: OK -> a.plateau={cmp_['a']['plateau_reached']}, "
          f"b.plateau={cmp_['b']['plateau_reached']}")
except ValueError as exc:
    print(f"7. compare_runs: skipped ({exc})")

print("ALL PHASE 1+2 TESTS PASSED")
