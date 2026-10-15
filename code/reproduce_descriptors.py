"""Rebuild the descriptor table (Table S1) and the final-geometry file from the raw MOPAC outputs.

This is the pipeline that produced data/results_final_v2.csv and data/final_geoms_final.json:

  1. analyze.py            first-coordination-sphere descriptors (cutoff 3.0 A) and termination status of every
                           original run (<job>/<job>.aux, .out) and of every rerun (<job>/<job>_rerun.aux, .out)
  2. merge                 a rerun (repetition of a run that failed to reach self-consistence, from the same
                           starting geometry, with SHIFT=15 PULAY CAMP ITRY=999 added) replaces the original run
                           unless it failed again (status "ended"); column "rerun" = 1 when the rerun was used;
                           ncycles is the cycle count of the selected run
  3. final geometries      geometry, Mulliken charges, heat of formation and heat-of-formation trajectory of
                           the selected run, read from its .aux file (trajectory kept every 5th point for the
                           original runs and in full for the reruns, as archived)
  4. residue_descriptors   residue-level descriptors with the strict carboxylate definition (kind A only)

Status values: converged (gradient norm < 1.0 kcal/mol/A), cycles (3000-cycle budget exhausted), ended (SCF failure).

Usage:
    python reproduce_descriptors.py [--runs ../mopac_optimizations] [--out ./descriptors_out] [--compare ../data]

With --compare, the rebuilt files are compared with the archived ones (byte-for-byte and value-by-value).
Requires: Python 3.9+, numpy, pandas. Runtime: about 1 min (reads ~1.7 GB of .aux files).
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import analyze                                   # noqa: E402
from mopio import read_aux                       # noqa: E402
from residue_descriptors import residue_row      # noqa: E402

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--runs", default=os.path.join(HERE, "..", "mopac_optimizations"))
ap.add_argument("--out", default=os.path.join(os.getcwd(), "descriptors_out"))
ap.add_argument("--compare", default=None, help="folder with the archived results_final_v2.csv and final_geoms_final.json")
args = ap.parse_args()
RUNS = os.path.abspath(args.runs)
os.makedirs(args.out, exist_ok=True)
man = json.load(open(os.path.join(RUNS, "manifest.json")))


def analyze_stem(rec, stem):
    """analyze.analyze_job applied to the files <job>/<stem>.* (the original run when stem == job)."""
    original = analyze.last_stage
    analyze.last_stage = lambda jd, job, max_stage=6: (stem, 0)
    try:
        return analyze.analyze_job(RUNS, rec)
    finally:
        analyze.last_stage = original


def same_start(job):
    """True when the rerun starts from the same geometry as the original input (geometry block of the .mop)."""
    geo = lambda stem: open(os.path.join(RUNS, job, stem + ".mop")).read().splitlines()[3:]
    return geo(job) == geo(job + "_rerun")


def final_geometry(job, stem, traj_step, with_status=True):
    """Same fields as the archived final_geoms_final.json (code of the original extract_final.py; the entries of
    the rerun jobs were written without the 'stage' and 'status' fields, which are reproduced here as archived)."""
    jd = os.path.join(RUNS, job)
    aux = os.path.join(jd, stem + ".aux")
    if not os.path.exists(aux):
        return None
    a = read_aux(aux)
    if "xyz" not in a:
        return None
    g = dict(stage=0, status=analyze.out_status(os.path.join(jd, stem + ".out"))) if with_status else {}
    g.update(el=a["el"], xyz=np.round(a["xyz"], 4).tolist(), heat=a.get("heat"), gnorm=a.get("gnorm"),
             charges=(np.round(a["charges"], 4).tolist() if "charges" in a else None),
             ncycles=len(a.get("heat_traj", [])), heat_traj=[round(x, 3) for x in a.get("heat_traj", [])][::traj_step])
    return g


# ---------------------------------------------------------------- 1. descriptors of original runs and reruns
orig = pd.DataFrame([analyze_stem(rec, rec["job"]) for rec in man])
rerun_jobs = [rec["job"] for rec in man if os.path.exists(os.path.join(RUNS, rec["job"], rec["job"] + "_rerun.aux"))]
rer = pd.DataFrame([analyze_stem(rec, rec["job"] + "_rerun") for rec in man if rec["job"] in rerun_jobs]).set_index("job")
orig_csv = os.path.join(args.out, "results_original_runs.csv")
orig.to_csv(orig_csv, index=False)
orig = pd.read_csv(orig_csv)          # as in the original pipeline, the original-run table is re-read from its CSV
print(f"analysed {len(orig)} original runs and {len(rer)} reruns")

# ---------------------------------------------------------------- 2. merge
merged = orig.copy()
merged["rerun"] = 0
use_rerun = []
for i, job in enumerate(merged.job):
    if job in rer.index and rer.loc[job, "status"] not in ("missing", "running", "ended"):
        r = rer.loc[job]
        for c in orig.columns:
            if c == "job":
                continue
            merged.at[i, c] = r[c]
        if not same_start(job):          # a restart from the last geometry would count the cycles of both runs
            merged.at[i, "ncycles"] = orig.loc[i, "ncycles"] + r["ncycles"]
        merged.at[i, "rerun"] = 1
        use_rerun.append(job)
merged["scf_fail"] = 0
merged_csv = os.path.join(args.out, "results_merged.csv")
merged.to_csv(merged_csv, index=False)
print(f"merge: rerun used for {len(use_rerun)} jobs")

# ---------------------------------------------------------------- 3. final geometries
G = {}
for rec in man:
    job = rec["job"]
    g = final_geometry(job, job + "_rerun", 1, with_status=False) if job in use_rerun else final_geometry(job, job, 5)
    if g is not None:
        G[job] = g
geoms_json = os.path.join(args.out, "final_geoms_final.json")
json.dump(G, open(geoms_json, "w"))

# ---------------------------------------------------------------- 4. residue-level descriptors (residue_descriptors.py)
df = pd.read_csv(merged_csv)
mrec = {m["job"]: m for m in man}
new = {job: residue_row(G[job]["el"], G[job]["xyz"], mrec[job]["residues"]) for job in df.job if job in G and mrec[job].get("residues")}
nd = pd.DataFrame.from_dict(new, orient="index")
for col in nd.columns:
    df[col] = df.job.map(nd[col])
df = df.drop(columns=[c for c in ("chelate_ring_O", "hydroxyl_bound") if c in df.columns])
final_csv = os.path.join(args.out, "results_final_v2.csv")
df.to_csv(final_csv, index=False)
# multi-structure XYZ with the same final geometries (MOPAC atom order), as distributed in data/ (code of make_SI.py)
xyz_all = os.path.join(args.out, "final_geometries_all.xyz")
with open(xyz_all, "w") as f:
    for rec in man:
        job = rec["job"]
        if job not in G:
            continue
        g = G[job]; row = df[df.job == job].iloc[0]
        f.write(f"{len(g['el'])}\n{job} model={row.model} pattern={row.pattern} conf={row.conf} start={row['mode']} status={row.status} heat_kcal={g.get('heat')}\n")
        for e, p in zip(g["el"], g["xyz"]):
            f.write(f"{e:2s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}\n")
print("written", final_csv, df.shape, ",", geoms_json, f"({len(G)} structures) and", xyz_all)

# ---------------------------------------------------------------- comparison with the archived files
if args.compare:
    for name in ("results_final_v2.csv", "final_geoms_final.json", "final_geometries_all.xyz"):
        a, b = os.path.join(args.out, name), os.path.join(args.compare, name)
        same_bytes = open(a, "rb").read() == open(b, "rb").read()
        if name.endswith(".csv"):
            x, y = pd.read_csv(a), pd.read_csv(b)
            num = x.select_dtypes("number").columns
            same_values = (list(x.columns) == list(y.columns) and x.shape == y.shape
                           and np.allclose(x[num].astype(float), y[num].astype(float), rtol=0, atol=1e-12, equal_nan=True)
                           and x.drop(columns=num).fillna("").astype(str).equals(y.drop(columns=num).fillna("").astype(str)))
        elif name.endswith(".json"):
            same_values = json.load(open(a)) == json.load(open(b))
        else:
            same_values = open(a).read().split() == open(b).read().split()
        print(f"{name}: identical bytes = {same_bytes}; identical values = {same_values}")
