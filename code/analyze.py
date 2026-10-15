"""Analysis of the Eu/pectin campaign: first coordination sphere descriptors for every job."""
import os, re, json, glob, sys
import numpy as np
import pandas as pd
from mopio import read_aux
from sphere import classify_oxygens

PECTIN_TYPES = ("carboxylate", "carbonyl_ester", "ether_ester", "hydroxyl", "ring_or_glycosidic", "acid_OH", "carbonyl_acid", "alkoxide")
WATER_TYPES = ("water", "hydroxide", "oxonium", "oxide", "water?")


def last_stage(jobdir, job, max_stage=6):
    """Return (name, stage) of the last stage that has an .aux file."""
    best = (job, 0)
    for s in range(1, max_stage + 1):
        if os.path.exists(os.path.join(jobdir, f"{job}_r{s}.aux")):
            best = (f"{job}_r{s}", s)
    return best


def out_status(outfile):
    """Termination status of a MOPAC run from its .out file: converged (gradient criterion met), cycles (cycle budget
    exhausted: "EXCESS NUMBER OF OPTIMIZATION CYCLES"), ended (normal termination otherwise; in this study, SCF
    failure: "UNABLE TO ACHIEVE SELF-CONSISTENCE"), running (no termination line), missing, trust_radius, herbert."""
    if not os.path.exists(outfile):
        return "missing"
    txt = open(outfile, errors="ignore").read()
    if "== MOPAC DONE ==" not in txt and "JOB ENDED NORMALLY" not in txt and "ended normally" not in txt:
        return "running"
    if "CUTOFF" in txt and "IS LESS THAN CUTOFF" in txt:
        return "converged"
    if "TRUST RADIUS NOW LESS THAN" in txt:
        return "trust_radius"
    if "HERBERTS TEST" in txt:
        return "herbert"
    if "EXCESS NUMBER OF OPTIMIZATION CYCLES" in txt or ("CYCLES" in txt and "EXCEEDED" in txt):
        return "cycles"
    return "ended"


def read_xyz(fn):
    lines = open(fn).read().splitlines()
    n = int(lines[0]); el = []; xyz = []
    for l in lines[2:2 + n]:
        t = l.split(); el.append(t[0]); xyz.append([float(t[1]), float(t[2]), float(t[3])])
    return el, np.array(xyz)


def kabsch_rmsd(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return float(np.sqrt(np.mean(np.sum((Pc @ R.T - Qc) ** 2, axis=1))))


def dihedral(p0, p1, p2, p3):
    b0 = p0 - p1; b1 = p2 - p1; b2 = p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n; w = b2 - np.dot(b2, b1n) * b1n
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w))))


def analyze_job(campdir, rec, cutoff=3.0):
    job = rec["job"]; jd = os.path.join(campdir, job)
    name, stage = last_stage(jd, job)
    aux = os.path.join(jd, name + ".aux")
    row = dict(job=job, model=rec["model"], series=rec.get("series"), pname=rec.get("pname"), pattern=rec.get("pattern"),
               conf=rec.get("conf"), mode=rec.get("mode"), charge=rec.get("charge"), stage=stage, status=out_status(os.path.join(jd, name + ".out")))
    if not os.path.exists(aux):
        row["status"] = "missing"
        return row
    a = read_aux(aux)
    if "xyz" not in a or "el" not in a:
        row["status"] = "noxyz"
        return row
    el, xyz = a["el"], a["xyz"]
    row["heat"] = a.get("heat"); row["gnorm"] = a.get("gnorm"); row["ncycles"] = len(a.get("heat_traj", []))
    if "Eu" not in el:
        return row
    eu = el.index("Eu")
    types = classify_oxygens(el, xyz)
    dist = {i: float(np.linalg.norm(xyz[i] - xyz[eu])) for i in types}
    fs = [(dist[i], i, types[i]) for i in types if dist[i] < cutoff]
    fs.sort()
    n_pect = sum(1 for d, i, t in fs if t in PECTIN_TYPES)
    n_wat = sum(1 for d, i, t in fs if t in WATER_TYPES)
    CN = len(fs)
    row.update(CN=CN, n_pectin=n_pect, n_water=n_wat,
               f_pectin=(n_pect / CN if CN else np.nan), f_water=(n_wat / CN if CN else np.nan))
    for t in ("carboxylate", "hydroxyl", "ring_or_glycosidic", "carbonyl_ester", "ether_ester", "acid_OH", "carbonyl_acid", "alkoxide", "hydroxide", "oxonium"):
        row["n_" + t] = sum(1 for d, i, tt in fs if tt == t)
    dp = [d for d, i, t in fs if t in PECTIN_TYPES]; dw = [d for d, i, t in fs if t in WATER_TYPES]
    dc = [d for d, i, t in fs if t == "carboxylate"]
    row.update(d_pectin_mean=(np.mean(dp) if dp else np.nan), d_pectin_min=(min(dp) if dp else np.nan),
               d_water_mean=(np.mean(dw) if dw else np.nan), d_water_min=(min(dw) if dw else np.nan), d_water_max=(max(dw) if dw else np.nan),
               d_carbox_mean=(np.mean(dc) if dc else np.nan), d_all_mean=(np.mean([d for d, i, t in fs]) if fs else np.nan),
               d_all_std=(np.std([d for d, i, t in fs]) if fs else np.nan),
               d_min_any=min(dist.values()) if dist else np.nan)
    # second-sphere / near contacts (3.0-3.6 A)
    row["n_pectin_3.6"] = sum(1 for i in types if types[i] in PECTIN_TYPES and dist[i] < 3.6)
    row["n_water_3.6"] = sum(1 for i in types if types[i] in WATER_TYPES and dist[i] < 3.6)
    # proton-transfer / integrity checks
    n_H = {i: sum(1 for k in range(len(el)) if el[k] == "H" and np.linalg.norm(xyz[k] - xyz[i]) < 1.25) for i in types}
    waters_bad = sum(1 for i in types if types[i] in ("hydroxide", "oxonium", "oxide"))
    row["n_water_defects"] = waters_bad
    row["n_acidOH"] = sum(1 for i in types if types[i] in ("acid_OH",))
    row["collapsed"] = int(any(d < 2.0 for d, i, t in fs))
    # residues: denticity and residue identity (from manifest residue map in MOPAC atom order)
    res = rec.get("residues")
    if res:
        bound_res = []; dent = []
        for k, r in enumerate(res):
            if r.get("kind") not in ("A", "H", "M"):
                continue
            oo = [r.get("O6a"), r.get("O6b")]
            nb = sum(1 for o in oo if o is not None and dist.get(o, 9) < cutoff)
            others = [key for key in ("O2", "O3", "O4", "O5", "O1") if r.get(key) is not None and dist.get(r[key], 9) < cutoff]
            if nb or others:
                bound_res.append(k)
                dent.append({"res": k, "kind": r["kind"], "carboxyl_O": nb, "other_O": others})
        row["n_residues_bound"] = len(bound_res)
        row["residues_bound"] = ",".join(str(k) for k in bound_res)
        row["n_bidentate_carboxylate"] = sum(1 for d in dent if d["carboxyl_O"] == 2)
        row["n_monodentate_carboxylate"] = sum(1 for d in dent if d["carboxyl_O"] == 1)
        row["n_carboxylates_bound"] = sum(1 for d in dent if d["carboxyl_O"] >= 1)
        row["chelate_ring_O"] = sum(1 for d in dent if d["carboxyl_O"] >= 1 and ("O5" in d["other_O"] or "O1" in d["other_O"] or "O4" in d["other_O"]))
        row["hydroxyl_bound"] = sum(1 for d in dent for o in d["other_O"] if o in ("O2", "O3"))
        motif = []
        for d in dent:
            tag = {2: "bi", 1: "mono", 0: ""}[d["carboxyl_O"]]
            if d["other_O"]:
                tag += "+" + "".join(d["other_O"])
            motif.append(f"{d['kind']}{d['res']}:{tag}")
        row["motif"] = " ".join(motif)
        # chain conformational change vs start
        sx = os.path.join(jd, job + "_start.xyz")
        if os.path.exists(sx):
            el0, xyz0 = read_xyz(sx)
            chain_idx = [i for i in range(len(el)) if i < eu and el[i] != "H"]
            if len(el0) == len(el):
                row["chain_rmsd"] = kabsch_rmsd(xyz0[chain_idx], xyz[chain_idx])
                # glycosidic torsions
                tors = []
                for i in range(len(res) - 1):
                    A, B = res[i], res[i + 1]
                    phi = dihedral(xyz[A["O5"]], xyz[A["C1"]], xyz[B["O4"]], xyz[B["C4"]])
                    psi = dihedral(xyz[A["C1"]], xyz[B["O4"]], xyz[B["C4"]], xyz[B["C5"]])
                    phi0 = dihedral(xyz0[A["O5"]], xyz0[A["C1"]], xyz0[B["O4"]], xyz0[B["C4"]])
                    psi0 = dihedral(xyz0[A["C1"]], xyz0[B["O4"]], xyz0[B["C4"]], xyz0[B["C5"]])
                    tors.append((round(phi, 1), round(psi, 1), round(phi0, 1), round(psi0, 1)))
                row["torsions_final"] = ";".join(f"{p:.0f}/{s:.0f}" for p, s, p0, s0 in tors)
                dphi = [abs(((p - p0 + 180) % 360) - 180) for p, s, p0, s0 in tors]
                dpsi = [abs(((s - s0 + 180) % 360) - 180) for p, s, p0, s0 in tors]
                row["max_dtorsion"] = max(dphi + dpsi) if (dphi + dpsi) else np.nan
    # Eu displacement from start
    sx = os.path.join(jd, job + "_start.xyz")
    if os.path.exists(sx):
        el0, xyz0 = read_xyz(sx)
        if len(el0) == len(el):
            row["eu_min_chainO_start"] = float(min(np.linalg.norm(xyz0[i] - xyz0[eu]) for i in types if types[i] in PECTIN_TYPES)) if any(types[i] in PECTIN_TYPES for i in types) else np.nan
    return row


def analyze_campaign(campdir, cutoff=3.0):
    man = json.load(open(os.path.join(campdir, "manifest.json")))
    rows = [analyze_job(campdir, rec, cutoff) for rec in man]
    return pd.DataFrame(rows)


if __name__ == "__main__":
    campdir = sys.argv[1]
    df = analyze_campaign(campdir)
    df.to_csv(os.path.join(campdir, "results.csv"), index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    print(df[["job", "status", "stage", "gnorm", "CN", "n_pectin", "n_water", "f_pectin", "n_carboxylate", "n_hydroxyl", "n_ring_or_glycosidic", "d_carbox_mean", "collapsed", "n_water_defects", "motif"]].to_string())
