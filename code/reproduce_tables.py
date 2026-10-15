"""Recompute Tables 2-6, Tables S2-S5 and the quantitative statements of the manuscript from the data files.

Inputs (same folder as this script, or pass --data <folder>):
    results_final_v2.csv       Table S1: one row per optimization (180 rows), descriptors at cutoff 3.0 A
    scans_semiempirical.json   rigid Eu-O scans (heats of formation) for Table 2 / Table S3 / Figure 3
    final_geoms_final.json     final geometries and Mulliken charges (atom counts, Section 3.10 charges)

Outputs (written to --out, default ./tables_out):
    Table2_calibration.csv, Table3_by_start.csv, Table4_by_pattern.csv, Table5_waters.csv,
    Table6_matched_pairs.csv, TableS2_not_converged.csv, TableS3_EuO_scans.csv, TableS4_controls.csv,
    TableS5_sensitivity.csv, key_numbers.csv

The definitions are those of Sections 2.5-2.6 of the manuscript and S2 of the Supplementary Information:
converged = gradient norm < 1.0 kcal/mol/A (status == "converged"); main series = series == "main";
inner-sphere starts = bid, mono2, ester; exclusion rule: CN < 6 (one structure, PM7_M6_cl_ester).
Requires: Python 3.9+, numpy, pandas, scipy.
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--data", default=HERE, help="folder with results_final_v2.csv, scans_semiempirical.json, final_geoms_final.json")
ap.add_argument("--out", default=os.path.join(os.getcwd(), "tables_out"))
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)

df = pd.read_csv(os.path.join(args.data, "results_final_v2.csv"))
scans = json.load(open(os.path.join(args.data, "scans_semiempirical.json")))
geoms = json.load(open(os.path.join(args.data, "final_geoms_final.json")))

DM = {"A6": 0, "H3blk": 0, "H3alt": 0, "P3": 33, "B3": 50, "L3": 50, "B2": 67, "D2": 67, "S1": 83, "M6": 100}
PATTERNS = ["A6", "H3blk", "H3alt", "P3", "B3", "L3", "B2", "D2", "S1", "M6"]
NAME = {"PM7": "Sparkle/PM7*", "RM1": "RM1"}

conv = df[df.status == "converged"].copy()
main = conv[conv.series == "main"].copy()
inner = main[main["mode"] != "hyd"].copy()
ok = inner[inner.CN >= 6].copy()                       # exclusion rule (Section 2.6)
hyd = main[main["mode"] == "hyd"].copy()
excluded = inner[inner.CN < 6]

KEY = []                                               # (id, value, description)


def key(kid, value, desc):
    KEY.append((kid, value, desc))


# ---------------------------------------------------------------- run statistics (Section 2.5)
key("n_jobs", len(df), "optimizations in the campaign")
key("n_main", int((df.series == "main").sum()), "main-series optimizations")
key("n_sens", int((df.series == "sens_hydration").sum()), "hydration-sensitivity optimizations")
key("n_ctrl", int((df.series == "control").sum()), "control optimizations")
key("n_conv", len(conv), "converged optimizations (gradient norm < 1.0)")
for m in ("PM7", "RM1"):
    key(f"n_conv_{m}", int((conv.model == m).sum()), f"converged, {NAME[m]}")
    key(f"n_jobs_{m}", int((df.model == m).sum()), f"optimizations, {NAME[m]}")
key("n_not_conv", int((df.status != "converged").sum()), "optimizations without convergence (Table S2)")
key("n_cycle_budget", int((df.status == "cycles").sum()), "2.5: optimizations that exhausted the 3000-cycle budget")
key("n_scf_fail", int((df.status == "ended").sum()), "2.5: optimizations that did not reach self-consistence")
nat = {j: len(g["el"]) for j, g in geoms.items()}
nat_main = [nat[j] for j in df[df.series != "control"].job]
nat_ctrl = [nat[j] for j in df[df.series == "control"].job]
key("natoms_min", min(nat_main), "smallest hexasaccharide cluster (atoms)")
key("natoms_max", max(nat_main), "largest hexasaccharide cluster (atoms)")
key("natoms_ctrl_min", min(nat_ctrl), "smallest control cluster (atoms)")
key("natoms_ctrl_max", max(nat_ctrl), "largest control cluster (atoms)")
key("n_excluded_CN_lt_6", len(excluded), "converged inner-sphere structures with CN < 6 (excluded)")

# ---------------------------------------------------------------- Table 2 (scan minima + controls)
r = np.array(scans["r"])


def scan_be(method, lig):
    d = scans[method]
    e = np.array([np.nan if x is None else x for x in d[lig]], float)
    return e - d["Eu3+"] - d[lig + "_ligand"]


def scan_min(method, lig):
    be = scan_be(method, lig)
    k = int(np.nanargmin(be))
    return r[k], be[k]


LIG = [("water", "H2O (calibration target)"), ("methanol", "CH3OH, neutral hydroxyl (control)"),
       ("ether", "(CH3)2O, ring/glycosidic O (control)"), ("ester", "CH3COOCH3, ester C=O (control)"),
       ("acetate", "CH3COO-, bidentate (control)"), ("methoxide", "CH3O-, localised anion (stress test)")]
METHODS = ("Sparkle/PM7", "Sparkle/PM7*", "RM1")
t2 = []
for lig, lab in LIG:
    row = {"Quantity": lab}
    for m in METHODS:
        rm, be = scan_min(m, lig)
        row[m] = f"{rm:.2f} / {be:.0f}"
        key(f"scan_min_{m}_{lig}", f"{rm:.2f} / {be:.0f}", f"Table 2: minimum (A) / depth (kcal/mol), {m}, {lig}")
    t2.append(row)
ctrl = conv[conv.series == "control"].set_index("job")


def cval(job, col):
    return ctrl.loc[job, col] if job in ctrl.index else np.nan


def f3(x):
    return "—" if pd.isna(x) else f"{x:.3f}"


for lab, col, jobs in (("[Eu(H2O)8]3+ optimised: mean Eu-O (A)", "d_water_mean", "Eu_w8"),
                       ("[Eu(H2O)9]3+ optimised: mean Eu-O (A)", "d_water_mean", "Eu_w9"),
                       ("Eu-O(carboxylate), GalA monomer chelate (A)", "d_carbox_mean", "GalA_bid"),
                       ("Eu-O(carboxylate), acetate chelate (A)", "d_carbox_mean", "acetate_bid")):
    row = {"Quantity": lab, "Sparkle/PM7": "—"}
    for m, mm in (("PM7", "Sparkle/PM7*"), ("RM1", "RM1")):
        row[mm] = f3(cval(f"{m}_CTRL_{jobs}", col))
        key(f"ctrl_{m}_{jobs}_{col}", row[mm], f"Table 2: {lab}, {NAME[m]}")
    t2.append(row)
row = {"Quantity": "CN from 10 and 12 waters", "Sparkle/PM7": "—"}
for m, mm in (("PM7", "Sparkle/PM7*"), ("RM1", "RM1")):
    row[mm] = f"{int(cval(f'{m}_CTRL_Eu_w10', 'CN'))}, {int(cval(f'{m}_CTRL_Eu_w12', 'CN'))}"
    key(f"ctrl_{m}_CN_w10_w12", row[mm], f"Table 2: CN from 10 and 12 waters, {NAME[m]}")
t2.append(row)
row = {"Quantity": "ninth-water binding, E(w9) - E(w8) (kcal/mol)", "Sparkle/PM7": "—"}
for m, mm in (("PM7", "Sparkle/PM7*"), ("RM1", "RM1")):
    row[mm] = f"{cval(f'{m}_CTRL_Eu_w9', 'heat') - cval(f'{m}_CTRL_Eu_w8', 'heat'):.1f}"
    key(f"ctrl_{m}_ninth_water", row[mm], f"Table 2: E(w9) - E(w8), {NAME[m]}")
t2.append(row)
pd.DataFrame(t2).to_csv(os.path.join(args.out, "Table2_calibration.csv"), index=False)

# Table S3 (binding energy vs distance) -- same numbers as Figure 3 before re-referencing
s3 = {"r_EuO_A": r}
for m in METHODS:
    for lig, _ in LIG:
        s3[f"{m}|{lig}|BE_kcal"] = np.round(scan_be(m, lig), 3)
pd.DataFrame(s3).to_csv(os.path.join(args.out, "TableS3_EuO_scans.csv"), index=False)

# ---------------------------------------------------------------- Table 3 (by Hamiltonian x start)
t3 = []
for m in ("PM7", "RM1"):
    for mode, lab in (("bid", "bidentate"), ("mono2", "monodentate + bridge"), ("ester", "ester carbonyl"), ("hyd", "hydrated (outer-sphere)")):
        d = (hyd if mode == "hyd" else ok)
        d = d[(d.model == m) & (d["mode"] == mode)]
        if not len(d):
            continue
        other = d.n_pectin - d.n_carboxylate
        row = {"Model": NAME[m], "Start": lab, "n": len(d),
               "CN (range)": f"{d.CN.mean():.2f} ({int(d.CN.min())}–{int(d.CN.max())})",
               "f_pectin": f"{d.f_pectin.mean():.3f} ± {d.f_pectin.std():.3f}",
               "carboxyl O": f"{d.n_carboxylate.mean():.2f}", "other pectin O": f"{other.mean():.2f}",
               "water O": f"{d.n_water.mean():.2f}",
               "Eu-O(carb) (A)": f3(d.d_carbox_mean.mean()) if d.d_carbox_mean.notna().any() else "—",
               "Eu-OH2 (A)": f3(d.d_water_mean.mean())}
        t3.append(row)
        for col in ("n", "CN (range)", "f_pectin", "carboxyl O", "other pectin O", "water O", "Eu-O(carb) (A)", "Eu-OH2 (A)"):
            key(f"T3_{m}_{mode}_{col}", row[col], f"Table 3: {NAME[m]}, {lab}, {col}")
pd.DataFrame(t3).to_csv(os.path.join(args.out, "Table3_by_start.csv"), index=False)

# ---------------------------------------------------------------- Table 4 (by pattern)
t4 = []
for p in PATTERNS:
    for m in ("PM7", "RM1"):
        d = ok[(ok.pname == p) & (ok.model == m)]
        if not len(d):
            continue
        row = {"Pattern": p, "DM (%)": DM[p], "Model": {"PM7": "PM7*", "RM1": "RM1"}[m], "n": len(d),
               "CN": f"{d.CN.mean():.1f}", "mean f_pectin": f"{d.f_pectin.mean():.3f}", "max f_pectin": f"{d.f_pectin.max():.2f}",
               "multi-COO": int((d.n_carboxylates_bound >= 2).sum()), "max residues": int(d.n_residues_bound.max()),
               "waters retained": f"{int(d.n_water.min())}–{int(d.n_water.max())}"}
        t4.append(row)
        for col in ("n", "CN", "mean f_pectin", "max f_pectin", "multi-COO", "max residues", "waters retained"):
            key(f"T4_{p}_{m}_{col}", row[col], f"Table 4: {p}, {NAME[m]}, {col}")
pd.DataFrame(t4).to_csv(os.path.join(args.out, "Table4_by_pattern.csv"), index=False)

# ---------------------------------------------------------------- Table 5 (retained waters by motif class)
ok["motif_class"] = np.where(ok.n_bidentate_carboxylate > 0, "bidentate carboxylate present",
                             np.where(ok.n_monodentate_carboxylate > 0, "monodentate carboxylate only", "ester / backbone O only (M6)"))
t5 = []
for m in ("PM7", "RM1"):
    for cl in ("bidentate carboxylate present", "monodentate carboxylate only", "ester / backbone O only (M6)"):
        d = ok[(ok.model == m) & (ok.motif_class == cl)]
        if not len(d):
            continue
        vc = d.n_water.value_counts().sort_index()
        t5.append({"Model": NAME[m], "Motif class": cl, "n": len(d), "median": int(d.n_water.median()),
                   "range": f"{int(d.n_water.min())}–{int(d.n_water.max())}",
                   "distribution": "; ".join(f"{int(k)} H2O: {int(v)}" for k, v in vc.items())})
    d = hyd[hyd.model == m]
    vc = d.n_water.value_counts().sort_index()
    t5.append({"Model": NAME[m], "Motif class": "hydrated start (no pectin O in first sphere)", "n": len(d), "median": int(d.n_water.median()),
               "range": f"{int(d.n_water.min())}–{int(d.n_water.max())}", "distribution": "; ".join(f"{int(k)} H2O: {int(v)}" for k, v in vc.items())})
for row in t5:
    tag = row["Motif class"].split()[0]
    for col in ("n", "median", "range", "distribution"):
        key(f"T5_{row['Model']}_{tag}_{col}", row[col], f"Table 5: {row['Model']}, {row['Motif class']}, {col}")
pd.DataFrame(t5).to_csv(os.path.join(args.out, "Table5_waters.csv"), index=False)

# ---------------------------------------------------------------- Table 6 (matched pairs at equal DM)
t6 = []
for a, b, lab in (("H3blk", "H3alt", "DM 0 (protonation variants)"), ("B3", "L3", "DM 50"), ("B2", "D2", "DM 67")):
    for m in ("PM7", "RM1"):
        d = ok[(ok.model == m) & (ok["mode"].isin(["bid", "mono2"]))]
        A = d[d.pname == a].set_index(["conf", "mode"])
        B = d[d.pname == b].set_index(["conf", "mode"])
        idx = A.index.intersection(B.index)
        dfp = A.loc[idx].f_pectin - B.loc[idx].f_pectin
        dres = A.loc[idx].n_residues_bound - B.loc[idx].n_residues_bound
        dco = A.loc[idx].n_carboxylates_bound - B.loc[idx].n_carboxylates_bound
        row = {"Pair": f"{a} vs {b} ({lab})", "Model": NAME[m], "pairs": len(idx),
               "dfp >0/=0/<0": f"{int((dfp > 1e-9).sum())} / {int((abs(dfp) <= 1e-9).sum())} / {int((dfp < -1e-9).sum())}",
               "median dfp (range)": f"{dfp.median():+.3f} ({dfp.min():+.3f} to {dfp.max():+.3f})",
               "d residues >0/=0/<0": f"{int((dres > 0).sum())} / {int((dres == 0).sum())} / {int((dres < 0).sum())}",
               "pairs with an extra carboxylate": int((dco > 0).sum())}
        t6.append(row)
        for col in ("pairs", "dfp >0/=0/<0", "median dfp (range)", "d residues >0/=0/<0", "pairs with an extra carboxylate"):
            key(f"T6_{a}{b}_{m}_{col}", row[col], f"Table 6: {a} vs {b}, {NAME[m]}, {col}")
pd.DataFrame(t6).to_csv(os.path.join(args.out, "Table6_matched_pairs.csv"), index=False)

# ---------------------------------------------------------------- statements in Sections 3.1-3.10
bid = ok[ok["mode"] == "bid"]
m2 = ok[ok["mode"] == "mono2"]
kept = m2[(m2.n_monodentate_carboxylate >= 1) & (m2.n_bidentate_carboxylate == 0)]
key("S3.2_n_bid", len(bid), "3.2: converged bidentate starts (CN >= 6)")
key("S3.2_bid_preserved", int((bid.n_bidentate_carboxylate >= 1).sum()), "3.2: bidentate starts that kept the chelate")
key("S3.2_dcarb_bid_all", f"{bid.d_carbox_mean.mean():.3f} ± {bid.d_carbox_mean.std():.3f}", "3.2: Eu-O(carboxylate), bidentate starts, mean ± SD")
key("S3.2_CN_bid", f"{bid.CN.mean():.1f}", "3.2: mean CN, bidentate starts")
key("S3.3_n_mono2", len(m2), "3.3: converged monodentate-plus-bridge starts")
key("S3.3_kept_mono", len(kept), "3.3: retained a single monodentate carboxylate")
key("S3.3_closed_bid", int((m2.n_bidentate_carboxylate >= 1).sum()), "3.3: closed to a bidentate chelate")
key("S3.3_dcarb_mono", f"{kept.d_carbox_mean.mean():.3f} ± {kept.d_carbox_mean.std():.3f}", "3.3: Eu-O(carboxylate), monodentate structures, mean ± SD")
key("S3.3_bid_minus_mono", f"{bid.d_carbox_mean.mean() - kept.d_carbox_mean.mean():.2f}", "3.3: bidentate minus monodentate Eu-O (A)")
key("S3.3_n_backbone", int((m2.n_backbone_O > 0).sum()), "3.3: mono2 structures with a ring/glycosidic O in the first sphere")
key("S3.3_n_hydroxyl", int((m2.n_hydroxyl_O > 0).sum()), "3.3: mono2 structures with a C2/C3 hydroxyl O")
key("S3.3_n_2COO", int((m2.n_carboxylates_bound >= 2).sum()), "3.3: mono2 structures with a second free carboxylate")
b3 = ok[ok.job == "RM1_B3_h2_mono2"].iloc[0]
key("S3.3_B3_example", f"CN {int(b3.CN)}, n_pectin {int(b3.n_pectin)}, n_water {int(b3.n_water)}, f {b3.f_pectin:.2f}, {b3.motif}", "3.3: RM1_B3_h2_mono2")
key("S3.4_n_ok", len(ok), "3.4: inner-sphere structures passing the exclusion criterion")
key("S3.4_n_pectin_distribution", "; ".join(f"{k}: {v}" for k, v in ok.n_pectin.value_counts().sort_index().items()), "3.4: number of pectin O (count of structures)")
key("S3.4_median_n_pectin", int(ok.n_pectin.median()), "3.4: median number of pectin O")
key("S3.4_fpectin_range", f"{ok.f_pectin.min():.2f}–{ok.f_pectin.max():.2f}", "3.4: f_pectin range")
key("S3.4_n_water_distribution", "; ".join(f"{k}: {v}" for k, v in ok.n_water.value_counts().sort_index().items()), "3.4: retained waters (count of structures)")
key("S3.5_n_hyd", f"{len(hyd)} ({int((hyd.model == 'PM7').sum())} PM7*, {int((hyd.model == 'RM1').sum())} RM1)", "3.5: converged hydrated starts")
key("S3.5_hyd_max_fpectin", f"{hyd.f_pectin.max():.3f}", "3.5: largest f_pectin among hydrated starts")
for m in ("PM7", "RM1"):
    h = hyd[hyd.model == m]
    key(f"S3.5_outer_contacts_{m}", int(((h["n_pectin_3.6"] - h.n_pectin) >= 1).sum()), f"3.5: hydrated starts with pectin O at 3.0-3.6 A, {NAME[m]}")
for m in ("PM7", "RM1"):
    d = ok[ok.model == m]
    rho, p = spearmanr(d.pname.map(DM), d.f_pectin)
    key(f"S3.6_spearman_{m}", f"{rho:.2f} (p = {p:.2f}, n = {len(d)})", f"3.6: Spearman DM vs f_pectin, {NAME[m]}")
    pm = d[d.pname != "M6"].groupby("pname").f_pectin.mean()
    key(f"S3.6_pattern_mean_range_{m}", f"{pm.min():.3f}–{pm.max():.3f}", f"3.6: range of pattern means (carboxylate patterns), {NAME[m]}")
m6 = ok[(ok.pname == "M6") & (ok.model == "RM1")]
key("S3.7_M6_RM1", f"CN {int(m6.CN.min())}–{int(m6.CN.max())}, f_pectin {m6.f_pectin.min():.2f}–{m6.f_pectin.max():.2f}, waters {int(m6.n_water.min())}–{int(m6.n_water.max())}", "3.7: M6 with RM1")
key("S3.7_M6_PM7", ", ".join(f"{r.job} CN {int(r.CN)}" for r in excluded.itertuples()), "3.7: excluded M6 structure with Sparkle/PM7*")
allm = pd.concat([ok, hyd])
key("S3.8_water_EuO", f"{allm.d_water_mean.mean():.3f} ± {allm.d_water_mean.std():.3f}", "3.8: Eu-OH2, all converged main-series structures (mean ± SD of per-structure means)")
for m in ("RM1", "PM7"):
    key(f"S3.8_water_EuO_{m}", f"{allm[allm.model == m].d_water_mean.mean():.3f}", f"3.8: Eu-OH2, {NAME[m]}")
key("S3.8_chain_rmsd", f"{ok.chain_rmsd.mean():.2f}", "3.8: mean chain heavy-atom RMSD of the inner-sphere structures (A)")
key("S3.8_max_dtorsion", f"{ok.max_dtorsion.mean():.1f}", "3.8: mean largest glycosidic torsion change of the inner-sphere structures (deg)")
key("S3.8_chain_rmsd_hyd", f"{hyd.chain_rmsd.mean():.2f}", "3.8: mean chain RMSD of hydrated starts (A)")
sens = conv[conv.series == "sens_hydration"]
key("S3.9_sens_converged", ", ".join(sorted(sens.job)), "3.9: converged jobs of the hydration-sensitivity series")
pairs = []
for r4 in sens[sens.job.str.endswith("cw4")].itertuples():
    b0 = sens[sens.job == r4.job.replace("cw4", "cw0")]
    if len(b0):
        b0 = b0.iloc[0]
        pairs.append(f"{r4.job}: motif {'same' if r4.motif == b0.motif else 'DIFFERENT'}, f {r4.f_pectin:.3f} vs {b0.f_pectin:.3f}, dEuO {abs(r4.d_carbox_mean - b0.d_carbox_mean):.3f}")
key("S3.9_sens_pairs", " | ".join(pairs), "3.9: cw4 variant vs cw0 baseline")
q = {j: geoms[j]["charges"][geoms[j]["el"].index("Eu")] for j in ("RM1_B3_h2_mono2", "RM1_CTRL_Eu_w9") if geoms[j].get("charges")}
key("S3.10_Eu_charge_B3", f"{q.get('RM1_B3_h2_mono2', np.nan):+.2f}", "3.10: Mulliken charge on Eu, RM1_B3_h2_mono2 (final SCF of the optimization)")
key("S3.10_Eu_charge_aqua", f"{q.get('RM1_CTRL_Eu_w9', np.nan):+.2f}", "3.10: Mulliken charge on Eu, RM1_CTRL_Eu_w9")

# ---------------------------------------------------------------- Supplementary Tables S2, S4 and S5
# Table S2: optimizations without convergence (descriptors of the last available geometry).
# status "cycles" = the 3000-cycle budget was exhausted; "ended" = MOPAC stopped because self-consistence was not achieved.
nc = df[df.status != "converged"]
REASON = {"cycles": "cycle budget exhausted", "ended": "SCF not converged"}
pd.DataFrame({"Job": nc.job, "Model": nc.model.map(NAME), "Reason": nc.status.map(REASON), "Cycles": nc.ncycles.astype(int), "CN": nc.CN.astype(int),
              "f_pectin": nc.f_pectin.map(lambda x: f"{x:.3f}"), "Motif": nc.motif.fillna("—"), "rerun used": nc.rerun}).to_csv(
    os.path.join(args.out, "TableS2_not_converged.csv"), index=False)
# Table S4: controls
c = df[df.series == "control"]
pd.DataFrame({"Job": c.job, "Model": c.model.map(NAME), "System": c.pname, "Start": c["mode"], "q": c.charge, "Status": c.status,
              "CN": c.CN.astype(int), "Eu-OH2 mean (A)": c.d_water_mean.map(f3), "Eu-O(carb) (A)": c.d_carbox_mean.map(f3),
              "dHf (kcal/mol)": c.heat.map(lambda x: f"{x:.2f}")}).to_csv(os.path.join(args.out, "TableS4_controls.csv"), index=False)
# Table S5: hydration-sensitivity series
h5 = df[df.series == "sens_hydration"]
pd.DataFrame({"Job": h5.job, "Model": h5.model.map(NAME), "Pattern": h5.pname,
              "Variant": h5.job.str.endswith("cw4").map({True: "cw4: 4 H-bonded waters on each free carboxylate not bound to Eu",
                                                         False: "cw0: no water H-bonded to the free carboxylates"}),
              "Status": h5.status.where(h5.status == "converged", "-"), "CN": h5.CN.astype(int), "f_pectin": h5.f_pectin.map(lambda x: f"{x:.3f}"),
              "Eu-O(carb) (A)": h5.d_carbox_mean.map(f3), "Motif": h5.motif}).to_csv(os.path.join(args.out, "TableS5_sensitivity.csv"), index=False)

kdf = pd.DataFrame(KEY, columns=["id", "value", "description"])
kdf.to_csv(os.path.join(args.out, "key_numbers.csv"), index=False)
print(f"tables and {len(kdf)} key numbers written to {args.out}")
