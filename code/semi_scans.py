"""Semiempirical Eu-O rigid scans (Sparkle/PM7, Sparkle/PM7*, RM1) on the shared ligand geometries (Table 2, Table S3,
Figure 3). Writes the MOPAC inputs/outputs to ./scan_runs and the heats of formation to ./scans_semiempirical.json.

Usage: python semi_scans.py            (MOPAC 23.2.5 executable from the environment variable MOPAC_EXE, or 'mopac')
Re-run on 2026-10-13 with MOPAC 23.2.5 (Linux): identical to data/scans_semiempirical.json for water, methanol,
methoxide and acetate; ether and ester within 0.001 kcal/mol (their MMFF ligand geometries depend on the RDKit version).
"""
import os, re, json, subprocess, numpy as np, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ligands import ligands, place_eu, R_GRID

M = os.environ.get("MOPAC_EXE", "mopac")
ENV = dict(os.environ)
if os.path.isabs(M):   # Linux release: libmopac.so lives in ../lib next to bin/mopac
    ENV["LD_LIBRARY_PATH"] = os.path.join(os.path.dirname(os.path.dirname(M)), "lib") + os.pathsep + ENV.get("LD_LIBRARY_PATH", "")
W = os.path.abspath("scan_runs")
os.makedirs(W, exist_ok=True)
open(os.path.join(W, "eu_fix.txt"), "w").write("XFAC_Eu O 2.207500\nALPB_Eu O 2.500000\n")
METHODS = {"Sparkle/PM7": "PM7", "Sparkle/PM7*": "PM7 EXTERNAL=eu_fix.txt", "RM1": "RM1"}


def heat(name, keys, el, xyz):
    fn = os.path.join(W, name + ".mop")
    with open(fn, "w") as f:
        f.write(f"{keys} 1SCF\n{name}\n\n")
        for e, p in zip(el, xyz):
            f.write(f"{e} {p[0]:.6f} 0 {p[1]:.6f} 0 {p[2]:.6f} 0\n")
    subprocess.run([M, name + ".mop"], cwd=W, env=ENV, capture_output=True, timeout=300)
    txt = open(os.path.join(W, name + ".out"), errors="ignore").read()
    m = re.search(r"FINAL HEAT OF FORMATION =\s+(-?\d+\.\d+) KCAL", txt)
    return float(m.group(1)) if m else None


L = ligands()
out = {"r": R_GRID}
for mname, keys in METHODS.items():
    tag = mname.replace("/", "_").replace("*", "s")
    d = {}
    d["Eu3+"] = heat(f"{tag}_Eu", f"{keys} CHARGE=3", ["Eu"], np.zeros((1, 3)))
    for lname, lig in L.items():
        d[lname + "_ligand"] = heat(f"{tag}_{lname}_lig", f"{keys} CHARGE={lig['charge']}", lig["el"], lig["xyz"])
        curve = []
        for r in R_GRID:
            el, xyz = place_eu(lig, r)
            curve.append(heat(f"{tag}_{lname}_{r:.1f}", f"{keys} CHARGE={lig['charge'] + 3}", el, xyz))
        d[lname] = curve
        be = [None if (c is None or d["Eu3+"] is None or d[lname + "_ligand"] is None) else c - d["Eu3+"] - d[lname + "_ligand"] for c in curve]
        ok = [(r, b) for r, b in zip(R_GRID, be) if b is not None]
        rmin, bmin = min(ok, key=lambda t: t[1]) if ok else (None, None)
        print(f"{mname:13s} {lname:10s} min at r={rmin} BE={bmin}", flush=True)
    out[mname] = d
json.dump(out, open("scans_semiempirical.json", "w"), indent=1)
print("done")
