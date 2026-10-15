"""Rebuild the descriptor table of the revised article (data/results_v2_all_models.csv, Table S1) and check it against the
MOPAC outputs of the package.

Two links are verified, for the 388 optimizations of the four models (382 after removing the six superseded ester starts):

  A. MOPAC outputs -> archived values.  For every optimization the selected run is located (<job>.out, or <job>_rerun.out
     when a repetition of a run that failed to reach self-consistence exists and did not fail again) and its termination
     status, final heat of formation, gradient norm, number of cycles, net atomic charge of Eu and final Cartesian
     geometry are read from the .out file and compared with the archived table and geometry file.
     The .out file of a run that did not converge does not contain its last geometry; for those runs (28) the geometry
     of data/final_geoms_v2_all.json, extracted from the .aux files that are not part of this package, cannot be checked.
  B. Geometries -> descriptors.  All first-sphere descriptors (cutoff 3.0 A), the residue-level descriptors, the chain
     deformation and the chemical-integrity criterion are recomputed from data/final_geoms_v2_all.json, the manifest
     (data/manifest_v2_all.json) and the starting geometries (<job>_start.xyz), and the summaries of the vibrational
     analysis (data/vibrational_summary_v2.csv) and of the displacement test (data/displacement_test_results.csv) are
     merged in.  The result must be identical to data/results_v2_all_models.csv.
  C. Displacement test.  Status, energy change, heavy-atom RMSD from the parent structure and first-sphere changes of the
     186 re-optimizations are recomputed from their .out files and compared with data/displacement_test_results.csv.

Status values: converged (gradient norm < 1.0 kcal/mol/A), cycles (3000-cycle budget exhausted), ended (SCF failure).
analysis_set = converged, heavy-atom connectivity intact and not superseded.

Usage (from code/):  python reproduce_descriptors_v2.py [--sm ..] [--out descriptors_out] [--compare]
Writes results_v2_all_models.csv, TableS9_integrity.csv and check_outputs.csv to --out.  Needs numpy and pandas.
"""
import argparse, json, os, re, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import analyze                                    # noqa: E402  (out_status, kabsch_rmsd, dihedral, type lists)
import sphere                                     # noqa: E402
from residue_descriptors import residue_row       # noqa: E402
from integrity_check import bond_matrix, read_xyz  # noqa: E402

PECT, WAT = analyze.PECTIN_TYPES, analyze.WATER_TYPES
OPT_FOLDERS = ["mopac_optimizations", os.path.join("revision", "R2_sparkle_pm7_std"), os.path.join("revision", "R5_sparkle_rm1"),
               os.path.join("revision", "R6_M6_ester_restart")]


def read_out(fn):
    """Final geometry, heat of formation, gradient norm, net atomic charges and number of cycles from a MOPAC .out file."""
    txt = open(fn, errors="ignore").read(); o = {}
    idx = [m.start() for m in re.finditer(r"CARTESIAN COORDINATES", txt)]
    if idx:
        el = []; xyz = []
        for l in txt[idx[-1]:].splitlines()[1:]:
            t = l.split()
            if len(t) == 5 and t[0].isdigit(): el.append(t[1]); xyz.append([float(v) for v in t[2:5]])
            elif el: break
        o["el"] = el; o["xyz"] = np.array(xyz)
    m = re.findall(r"FINAL HEAT OF FORMATION\s*=\s*([-\d.]+)\s*KCAL/MOL", txt)
    if m: o["heat"] = float(m[-1])
    m = re.findall(r"GRADIENT NORM\s*=\s*([\d.]+)", txt)
    if m: o["gnorm"] = float(m[-1])
    k = txt.rfind("NET ATOMIC CHARGES")
    if k >= 0:
        ch = []
        for l in txt[k:].splitlines()[1:]:
            t = l.split()
            if len(t) >= 3 and t[0].isdigit(): ch.append(float(t[2]))
            elif ch: break
        o["charges"] = np.array(ch)
    o["ncycles"] = len(re.findall(r"^\s*CYCLE:\s*\d+", txt, flags=re.M))
    return o


def descriptors(rec, el, xyz, el0=None, xyz0=None, cutoff=3.0):
    """First-sphere descriptors of one geometry (definitions of analyze.analyze_job and residue_descriptors.residue_row)."""
    row = {}; el = list(el); xyz = np.asarray(xyz, float)
    if "Eu" not in el: return row
    eu = el.index("Eu"); types = sphere.classify_oxygens(el, xyz); dist = {i: float(np.linalg.norm(xyz[i] - xyz[eu])) for i in types}
    fs = sorted((dist[i], i, types[i]) for i in types if dist[i] < cutoff)
    n_p = sum(1 for d, i, t in fs if t in PECT); n_w = sum(1 for d, i, t in fs if t in WAT); CN = len(fs)
    row.update(CN=CN, n_pectin=n_p, n_water=n_w, f_pectin=(n_p / CN if CN else np.nan), f_water=(n_w / CN if CN else np.nan))
    for t in ("carboxylate", "hydroxyl", "ring_or_glycosidic", "carbonyl_ester", "ether_ester", "acid_OH", "carbonyl_acid", "alkoxide", "hydroxide", "oxonium"):
        row["n_" + t] = sum(1 for d, i, tt in fs if tt == t)
    dp = [d for d, i, t in fs if t in PECT]; dw = [d for d, i, t in fs if t in WAT]; dc = [d for d, i, t in fs if t == "carboxylate"]
    row.update(d_pectin_mean=(np.mean(dp) if dp else np.nan), d_pectin_min=(min(dp) if dp else np.nan), d_water_mean=(np.mean(dw) if dw else np.nan),
               d_water_min=(min(dw) if dw else np.nan), d_water_max=(max(dw) if dw else np.nan), d_carbox_mean=(np.mean(dc) if dc else np.nan),
               d_all_mean=(np.mean([d for d, i, t in fs]) if fs else np.nan), d_all_std=(np.std([d for d, i, t in fs]) if fs else np.nan),
               d_min_any=min(dist.values()) if dist else np.nan)
    row["n_pectin_3.6"] = sum(1 for i in types if types[i] in PECT and dist[i] < 3.6)
    row["n_water_3.6"] = sum(1 for i in types if types[i] in WAT and dist[i] < 3.6)
    row["n_water_defects"] = sum(1 for i in types if types[i] in ("hydroxide", "oxonium", "oxide"))
    row["n_acidOH"] = sum(1 for i in types if types[i] == "acid_OH"); row["collapsed"] = int(any(d < 2.0 for d, i, t in fs))
    res = rec.get("residues"); same = el0 is not None and len(el0) == len(el)
    if res:
        row.update(residue_row(el, xyz, res, cutoff))
        if same:
            xyz0 = np.asarray(xyz0, float); ci = [i for i in range(len(el)) if i < eu and el[i] != "H"]
            row["chain_rmsd"] = analyze.kabsch_rmsd(xyz0[ci], xyz[ci]); D = analyze.dihedral; tors = []
            for i in range(len(res) - 1):
                A, B = res[i], res[i + 1]
                tors.append((round(D(xyz[A["O5"]], xyz[A["C1"]], xyz[B["O4"]], xyz[B["C4"]]), 1), round(D(xyz[A["C1"]], xyz[B["O4"]], xyz[B["C4"]], xyz[B["C5"]]), 1),
                             round(D(xyz0[A["O5"]], xyz0[A["C1"]], xyz0[B["O4"]], xyz0[B["C4"]]), 1), round(D(xyz0[A["C1"]], xyz0[B["O4"]], xyz0[B["C4"]], xyz0[B["C5"]]), 1)))
            row["torsions_final"] = ";".join(f"{p:.0f}/{s:.0f}" for p, s, p0, s0 in tors)
            dphi = [abs(((p - p0 + 180) % 360) - 180) for p, s, p0, s0 in tors]; dpsi = [abs(((s - s0 + 180) % 360) - 180) for p, s, p0, s0 in tors]
            row["max_dtorsion"] = max(dphi + dpsi) if (dphi + dpsi) else np.nan
    if same:
        x0 = np.asarray(xyz0, float)
        row["eu_min_chainO_start"] = (float(min(np.linalg.norm(x0[i] - x0[eu]) for i in types if types[i] in PECT))
                                      if any(types[i] in PECT for i in types) else np.nan)
    return row


def integrity(rec, el, xyz, sel, sxyz, q):
    """Heavy-atom connectivity of the final structure against the starting structure (integrity_check.bond_matrix)."""
    xyz = np.asarray(xyz, float); Ms, _ = bond_matrix(sel, sxyz); Mf, Df = bond_matrix(el, xyz); heavy = np.array(el) != "H"; n = len(el)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n) if heavy[i] and heavy[j]]
    broken = [(i, j) for i, j in pairs if Ms[i, j] and not Mf[i, j]]; new = [(i, j) for i, j in pairs if Mf[i, j] and not Ms[i, j]]
    ieu = el.index("Eu"); has_q = q is not None and len(q) == n
    d = dict(n_broken_bonds=len(broken), n_new_bonds=len(new),
             broken_bonds=";".join(f"{el[i]}{i + 1}-{el[j]}{j + 1}:{Df[i, j]:.2f}" for i, j in broken),
             new_bonds=";".join(f"{el[i]}{i + 1}-{el[j]}{j + 1}:{Df[i, j]:.2f}" for i, j in new),
             q_Eu=(round(float(q[ieu]), 3) if has_q else np.nan))
    lost = []
    for k, rd in enumerate(rec.get("residues") or []):
        c5, c6, oa, ob = rd["C5"], rd["C6"], rd["O6a"], rd["O6b"]
        if np.linalg.norm(xyz[c5] - xyz[c6]) > 1.9:
            u, v = xyz[oa] - xyz[c6], xyz[ob] - xyz[c6]
            lost.append(dict(res=f"{rd['kind']}{k}", d_C5C6=np.linalg.norm(xyz[c5] - xyz[c6]),
                             angle=np.degrees(np.arccos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))),
                             d_CO=(np.linalg.norm(u) + np.linalg.norm(v)) / 2, q=(float(q[c6] + q[oa] + q[ob]) if has_q else np.nan),
                             d_Eu=min(np.linalg.norm(xyz[oa] - xyz[ieu]), np.linalg.norm(xyz[ob] - xyz[ieu])),
                             d_Eu_start=min(np.linalg.norm(sxyz[oa] - sxyz[ieu]), np.linalg.norm(sxyz[ob] - sxyz[ieu]))))
    d.update(n_CO2_lost=len(lost), CO2_residue=";".join(l["res"] for l in lost), CO2_d_C5C6=";".join(f"{l['d_C5C6']:.2f}" for l in lost),
             CO2_angle_OCO=";".join(f"{l['angle']:.1f}" for l in lost), CO2_d_CO=";".join(f"{l['d_CO']:.3f}" for l in lost),
             CO2_charge=";".join(f"{l['q']:.3f}" for l in lost), CO2_min_d_EuO=";".join(f"{l['d_Eu']:.2f}" for l in lost),
             CO2_min_d_EuO_start=";".join(f"{l['d_Eu_start']:.2f}" for l in lost), intact=(len(broken) == 0 and len(new) == 0))
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sm", default=os.path.join(HERE, "..")); ap.add_argument("--out", default=os.path.join(os.getcwd(), "descriptors_out"))
    ap.add_argument("--compare", action="store_true", help="compare the rebuilt table with data/results_v2_all_models.csv")
    a = ap.parse_args(); SM = os.path.abspath(a.sm); DATA = os.path.join(SM, "data"); os.makedirs(a.out, exist_ok=True)
    man = json.load(open(os.path.join(DATA, "manifest_v2_all.json"))); G = json.load(open(os.path.join(DATA, "final_geoms_v2_all.json")))
    ARCH = pd.read_csv(os.path.join(DATA, "results_v2_all_models.csv")).set_index("job")

    def folder(job):
        for f in OPT_FOLDERS:
            p = os.path.join(SM, f, job)
            if os.path.isdir(p): return p
        return None

    def start_geometry(rec):
        for j in (rec["job"], rec.get("base")):
            p = folder(j) if j else None
            if p and os.path.exists(os.path.join(p, j + "_start.xyz")): return read_xyz(os.path.join(p, j + "_start.xyz"))
        return None, None

    rows = []; integ = []; chk = []
    for rec in man:
        job = rec["job"]; jd = folder(job); g = G[job]; el = g["el"]; xyz = np.array(g["xyz"], float)
        # ---- A. selected run and values of its .out file
        stem = job
        rer = os.path.join(jd, job + "_rerun.out") if jd else None
        if rer and os.path.exists(rer) and analyze.out_status(rer) not in ("missing", "running", "ended"): stem = job + "_rerun"
        outf = os.path.join(jd, stem + ".out") if jd else None; status = analyze.out_status(outf) if outf else "missing"
        o = read_out(outf) if outf and os.path.exists(outf) else {}
        c = dict(job=job, run=stem, status=status)
        if status == "converged":
            c.update(heat_out=o.get("heat"), gnorm_out=o.get("gnorm"), ncycles_out=o.get("ncycles"))
            if "xyz" in o and o["el"] == list(el): c["geometry_max_deviation"] = float(np.abs(o["xyz"] - xyz).max())
            if "charges" in o and len(o["charges"]) == len(el) and "Eu" in el: c["q_Eu_out"] = round(float(o["charges"][el.index("Eu")]), 3)
        chk.append(c)
        # ---- B. descriptors from the archived geometry
        row = dict(job=job, model=rec["model"], series=rec.get("series"), pname=rec.get("pname"), pattern=rec.get("pattern"), conf=rec.get("conf"),
                   mode=rec.get("mode"), charge=rec.get("charge"), stage=0, status=status)
        if status == "converged":
            row.update(heat=ARCH.loc[job, "heat"], gnorm=ARCH.loc[job, "gnorm"])      # checked against the .out file in check_outputs.csv
        row["ncycles"] = ARCH.loc[job, "ncycles"]
        sel, sxyz = start_geometry(rec)
        row.update(descriptors(rec, el, xyz, sel, sxyz))
        if "Eu" in el and sel is not None and sel == list(el):
            d = integrity(rec, el, xyz, sel, sxyz, g.get("charges") if status == "converged" else None)
            row.update({k: d[k] for k in ("intact", "n_broken_bonds", "n_new_bonds", "n_CO2_lost", "q_Eu")})
            integ.append(dict(job=job, model=rec["model"], series=rec.get("series"), pname=rec.get("pname"), conf=rec.get("conf"), mode=rec.get("mode"), status=status, **d))
        rows.append(row)
    df = pd.DataFrame(rows)
    df["sys"] = [j.split("_", 1)[1] for j in df.job]; df["superseded"] = df["mode"].eq("ester") & df["series"].eq("main")      # defective ester starts of M6, replaced by ester2
    # ---- vibrational analysis and displacement test
    vib = pd.read_csv(os.path.join(DATA, "vibrational_summary_v2.csv")).set_index("job")
    for col, src in (("vib_n_imag", "n_imag"), ("vib_n_imag_gt20", "n_imag_gt20"), ("vib_n_imag_gt50", "n_imag_gt50"), ("vib_f1", "f1"), ("vib_f2", "f2"),
                     ("vib_gnorm_tangent", "gnorm_tangent"), ("vib_sumrule_max", "sumrule_max"), ("vib_n_mopac_neg", "n_mopac_neg"),
                     ("vib_m1_water_H", "m1_water_H"), ("vib_m1_water_O", "m1_water_O")):
        df[col] = df.job.map(vib[src])
    disp = pd.read_csv(os.path.join(DATA, "displacement_test_results.csv")); dc = disp[disp.status == "converged"]
    P = dc.groupby("parent").agg(disp_n_converged=("job", "size"), disp_dE_min=("dE", "min"), disp_dE_max=("dE", "max"), disp_rmsd_max=("rmsd_heavy_vs_parent", "max"),
                                 disp_max_dEuO=("max_dEuO_first_sphere", "max"), disp_same_first_sphere=("same_first_sphere_atoms", "all"))
    for col in P.columns: df[col] = df.job.map(P[col])
    df["analysis_set"] = (df.status == "converged") & (df.intact == True) & (~df.superseded)      # noqa: E712
    df = df[[c for c in ARCH.reset_index().columns if c in df.columns] + [c for c in df.columns if c not in ARCH.reset_index().columns]]
    df.to_csv(os.path.join(a.out, "results_v2_all_models.csv"), index=False)
    pd.DataFrame(integ).to_csv(os.path.join(a.out, "TableS9_integrity.csv"), index=False)
    C = pd.DataFrame(chk).set_index("job"); C.to_csv(os.path.join(a.out, "check_outputs.csv"))
    # ---- report of link A
    cv = C[C.status == "converged"]; arch = ARCH.loc[cv.index]
    print(f"A. outputs found for {int((C.status != 'missing').sum())} of {len(C)} optimizations; status equal to the archived one in {int((C.status == ARCH.loc[C.index, 'status']).sum())}")
    print(f"   converged runs: {len(cv)}; final geometry in the .out file within {cv.geometry_max_deviation.max():.0e} A of the archived one in {int(cv.geometry_max_deviation.notna().sum())};"
          f" heat of formation within {(cv.heat_out - arch.heat).abs().max():.3f} kcal/mol; gradient norm within {(cv.gnorm_out - arch.gnorm).abs().max():.0e};"
          f" cycles within {int((cv.ncycles_out - arch.ncycles).abs().max())}; charge of Eu equal (0.001) in {int(np.isclose(cv.q_Eu_out, arch.q_Eu, atol=1.5e-3).sum())}")
    print(f"   runs that did not converge (geometry not checkable from the .out file): {int((C.status != 'converged').sum())}")
    print("B. analysis set by model:", df[df.analysis_set].groupby("model").size().to_dict(), "| converged:", df[df.status == "converged"].groupby("model").size().to_dict())
    # ---- C. displacement test: every re-optimization against its .out file
    stages = [os.path.join(SM, "revision", s) for s in ("R7_displace_reopt", "R9_displace_reopt_2")]
    n = 0; nconv = 0; dv = dict(dE=0.0, rmsd=0.0, dEuO=0.0); same = 0; st_ok = 0
    for r in disp.itertuples():
        jd = next((os.path.join(s, r.job) for s in stages if os.path.isdir(os.path.join(s, r.job))), None)
        if jd is None: continue
        stem = [s for s in (r.job, r.job + "_r1", r.job + "_r2") if os.path.exists(os.path.join(jd, s + ".out"))][-1]
        status = analyze.out_status(os.path.join(jd, stem + ".out")); n += 1; st_ok += int(status == r.status)
        if status != "converged": continue
        o = read_out(os.path.join(jd, stem + ".out")); g = G[r.parent]; el = g["el"]; p0 = np.array(g["xyz"], float); x = o["xyz"]; hv = np.array(el) != "H"
        eu = el.index("Eu"); oi = [i for i in range(len(el)) if el[i] == "O"]
        d0 = np.linalg.norm(p0[oi] - p0[eu], axis=1); d1 = np.linalg.norm(x[oi] - x[eu], axis=1); fs = (d0 < 3.0) | (d1 < 3.0)
        dv["dE"] = max(dv["dE"], abs((o["heat"] - ARCH.loc[r.parent, "heat"]) - r.dE)); dv["rmsd"] = max(dv["rmsd"], abs(analyze.kabsch_rmsd(p0[hv], x[hv]) - r.rmsd_heavy_vs_parent))
        dv["dEuO"] = max(dv["dEuO"], abs(float(np.abs(d1 - d0)[fs].max()) - r.max_dEuO_first_sphere)); same += int(bool(((d0 < 3.0) == (d1 < 3.0)).all()) == bool(r.same_first_sphere_atoms)); nconv += 1
    print(f"C. displacement test: {n} of {len(disp)} outputs found; status equal in {st_ok}; converged {nconv}; largest deviations from data/displacement_test_results.csv: "
          f"energy change {dv['dE']:.3f} kcal/mol, RMSD {dv['rmsd']:.0e} A, first-sphere Eu-O change {dv['dEuO']:.0e} A; first-sphere atom set equal in {same}")
    if a.compare:
        x = df.set_index("job"); y = ARCH; bad = {}
        for col in y.columns:
            if col not in x.columns: bad[col] = "missing"; continue
            p, q = x[col], y[col]
            if pd.api.types.is_numeric_dtype(q) and pd.api.types.is_numeric_dtype(p):
                ok = np.isclose(p.astype(float), q.astype(float), rtol=0, atol=1e-9, equal_nan=True)
            else:
                ok = (p.fillna("").astype(str) == q.fillna("").astype(str)).values
            if not ok.all(): bad[col] = int((~ok).sum())
        print(f"compare with data/results_v2_all_models.csv: {len(y.columns)} columns, {len(y)} rows; identical values = {not bad}", bad if bad else "")


if __name__ == "__main__":
    main()
