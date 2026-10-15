"""Recompute Tables 2-6, the supplementary tables and the quantitative statements of the revised manuscript (v2).

Inputs (folder --data; default: ../data next to this script):
    results_v2_all_models.csv     one row per optimization (388 rows; models PM7 = Sparkle/PM7+EuO, SRM1 = Sparkle/RM1,
                                  RM1 = RM1 lanthanide model, PM7std = Sparkle/PM7 as distributed); column analysis_set =
                                  converged, heavy-atom connectivity intact (integrity_check.py) and not superseded
    scans_semiempirical_v2.json   rigid Eu-O scans (heats of formation), four models
    final_geoms_v2_all.json       final geometries and net atomic charges
    vibrational_summary_v2.csv    constrained harmonic analysis (constrained_modes.py)
    displacement_test_results.csv displacement test (re-optimizations after a displacement along the imaginary modes)

Outputs (--out, default ./tables_v2): Table2..Table6 CSV, supplementary tables, key_numbers_v2.csv.
Definitions: Sections 2.5-2.6 of the manuscript.  Inner-sphere starts = bid, mono2, ester2 (ester2 = corrected
ester-carbonyl start of the fully esterified chain).  Requires Python 3.9+, numpy, pandas, scipy.
"""
import argparse, json, os
import numpy as np, pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--data", default=os.path.join(HERE, "..", "data"))
ap.add_argument("--out", default=os.path.join(os.getcwd(), "tables_v2"))
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)
rd = lambda f: os.path.join(args.data, f)
df = pd.read_csv(rd("results_v2_all_models.csv"))
scans = json.load(open(rd("scans_semiempirical_v2.json")))
geoms = json.load(open(rd("final_geoms_v2_all.json")))
vib = pd.read_csv(rd("vibrational_summary_v2.csv"))
disp = pd.read_csv(rd("displacement_test_results.csv"))

DM = {"A6": 0, "H3blk": 0, "H3alt": 0, "P3": 33, "B3": 50, "L3": 50, "B2": 67, "D2": 67, "S1": 83, "M6": 100}
PATTERNS = ["A6", "H3blk", "H3alt", "P3", "B3", "L3", "B2", "D2", "S1", "M6"]
NAME = {"PM7": "Sparkle/PM7+EuO", "SRM1": "Sparkle/RM1", "RM1": "RM1", "PM7std": "Sparkle/PM7"}
PROD = ("PM7", "SRM1", "RM1")                      # models used for the structural statistics
SCAN = {"PM7std": "Sparkle/PM7", "PM7": "Sparkle/PM7*", "SRM1": "Sparkle/RM1", "RM1": "RM1"}
INNER = ("bid", "mono2", "ester2")

df = df[~df.superseded].copy()                     # the six defective M6 "ester" starts were replaced by "ester2"
conv = df[df.status == "converged"].copy()
aset = df[df.analysis_set].copy()
main = aset[aset.series == "main"].copy()
ok = main[main["mode"].isin(INNER) & main.model.isin(PROD)].copy()
hyd = main[(main["mode"] == "hyd") & main.model.isin(PROD)].copy()
KEY = []


def key(kid, value, desc):
    KEY.append((kid, value, desc))


def f3(x):
    return "\u2014" if pd.isna(x) else f"{x:.3f}"


def ms(s, nd=3):
    return f"{s.mean():.{nd}f} \u00b1 {s.std():.{nd}f}" if len(s) > 1 else (f"{s.mean():.{nd}f}" if len(s) else "\u2014")


# ---------------------------------------------------------------- run statistics and integrity (Section 2.5)
for m in ("PM7", "SRM1", "RM1", "PM7std"):
    d = df[df.model == m]; c = d[d.status == "converged"]
    key(f"runs_{m}", f"{len(d)} optimizations; {len(c)} converged; {int((d.status == 'cycles').sum())} cycle budget; "
        f"{int((d.status == 'ended').sum())} SCF failure; {int((c.intact == True).sum())} converged and intact; {int(d.analysis_set.sum())} in the analysis set",
        f"2.5: run statistics, {NAME[m]}")
key("runs_total", f"{len(df)} optimizations, {len(conv)} converged, {int(df.analysis_set.sum())} in the analysis set", "2.5: all models (superseded starts not counted)")
for m in ("PM7", "SRM1", "RM1", "PM7std"):
    c = conv[conv.model == m]
    key(f"integrity_{m}", f"{int((c.intact != True).sum())} of {len(c)} converged structures not intact; CO2 lost in {int((c.n_CO2_lost > 0).sum())}", f"2.5: integrity, {NAME[m]}")

# ---------------------------------------------------------------- Table 2 (scans + controls, four models)
r = np.array(scans["r"])


def scan_be(model, lig):
    d = scans[SCAN[model]]
    return np.array([np.nan if x is None else x for x in d[lig]], float) - d["Eu3+"] - d[lig + "_ligand"]


def scan_minima(model, lig):
    """global grid minimum, and inner local minimum (lowest grid point below 2.8 A that is a local minimum)."""
    be = scan_be(model, lig); k = int(np.nanargmin(be))
    loc = [i for i in range(1, len(r) - 1) if be[i] < be[i - 1] and be[i] < be[i + 1] and r[i] < 2.8]
    ki = min(loc, key=lambda i: be[i]) if loc else None
    return (r[k], be[k]), ((r[ki], be[ki]) if ki is not None else None)


LIG = [("water", "H2O"), ("methanol", "CH3OH, neutral hydroxyl"), ("ether", "(CH3)2O, ring/glycosidic O"),
       ("ester", "CH3COOCH3, ester C=O"), ("acetate", "CH3COO-, bidentate"), ("methoxide", "CH3O-, localised anion (stress test)")]
MODELS4 = ("PM7std", "PM7", "SRM1", "RM1")
t2 = []
for lig, lab in LIG:
    row = {"Quantity": "scan: " + lab}
    for m in MODELS4:
        (rg, bg), inner = scan_minima(m, lig)
        txt = f"{inner[0]:.1f} / {inner[1]:.0f}" if inner else "none"
        if inner is None or abs(rg - inner[0]) > 1e-6:
            txt += f" [global: {rg:.1f} / {bg:.0f}]"
        row[NAME[m]] = txt
        key(f"scan_{m}_{lig}", txt, f"Table 2: inner minimum (A) / depth (kcal/mol) [global grid minimum if different], {NAME[m]}, {lig}")
    t2.append(row)
ctrl = aset[aset.series == "control"].set_index("job")


def cval(job, col):
    return ctrl.loc[job, col] if job in ctrl.index else np.nan


for lab, col, jobs in (("[Eu(H2O)8]3+ optimised: mean Eu-O (A)", "d_water_mean", "Eu_w8"), ("[Eu(H2O)9]3+ optimised: mean Eu-O (A)", "d_water_mean", "Eu_w9"),
                       ("Eu-O(carboxylate), GalA monomer chelate (A)", "d_carbox_mean", "GalA_bid"), ("Eu-O(carboxylate), acetate chelate (A)", "d_carbox_mean", "acetate_bid"),
                       ("shortest Eu-O, GalA monomer chelate (A)", "d_min_any", "GalA_bid"), ("shortest Eu-O, acetate chelate (A)", "d_min_any", "acetate_bid")):
    row = {"Quantity": lab}
    for m in MODELS4:
        row[NAME[m]] = f3(cval(f"{m}_CTRL_{jobs}", col)); key(f"ctrl_{m}_{jobs}_{col}", row[NAME[m]], f"Table 2: {lab}, {NAME[m]}")
    t2.append(row)
row = {"Quantity": "CN from 8, 9, 10 and 12 waters"}
for m in MODELS4:
    row[NAME[m]] = ", ".join(str(int(cval(f"{m}_CTRL_Eu_w{n}", "CN"))) for n in (8, 9, 10, 12)); key(f"ctrl_{m}_CN", row[NAME[m]], f"Table 2: CN of the aqua clusters, {NAME[m]}")
t2.append(row)
row = {"Quantity": "ninth-water binding, E(w9) - E(w8) (kcal/mol)"}
for m in MODELS4:
    row[NAME[m]] = f"{cval(f'{m}_CTRL_Eu_w9', 'heat') - cval(f'{m}_CTRL_Eu_w8', 'heat'):.1f}"; key(f"ctrl_{m}_ninth_water", row[NAME[m]], f"Table 2: E(w9) - E(w8), {NAME[m]}")
t2.append(row)
pd.DataFrame(t2).to_csv(os.path.join(args.out, "Table2_calibration.csv"), index=False)
s3 = {"r_EuO_A": r}
for m in MODELS4:
    for lig, _ in LIG:
        s3[f"{NAME[m]}|{lig}|BE_kcal"] = np.round(scan_be(m, lig), 3)
pd.DataFrame(s3).to_csv(os.path.join(args.out, "TableS3_EuO_scans.csv"), index=False)

# ---------------------------------------------------------------- Sparkle/PM7 as distributed (Section 3.2)
cs = conv[conv.model == "PM7std"]
col = cs[cs.d_min_any < 2.0]
key("S3.2_std_converged", f"{len(cs)} of {int((df.model == 'PM7std').sum())}", "3.2: converged, Sparkle/PM7 as distributed")
key("S3.2_std_collapsed", f"{len(col)} of {len(cs)}", "3.2: converged structures with an Eu-O contact below 2.0 A")
key("S3.2_std_collapsed_by_start", "; ".join(f"{k}: {int((g.d_min_any < 2.0).sum())} of {len(g)}" for k, g in cs.groupby("mode")), "3.2: collapsed by start")
key("S3.2_std_dmin", f"{cs.d_min_any.min():.2f}", "3.2: shortest Eu-O contact (A)")
key("S3.2_std_not_intact", f"{int((cs.intact != True).sum())} (all collapsed: {bool((cs[cs.intact != True].d_min_any < 2.0).all())})", "3.2: converged structures with a broken or new heavy-atom bond")
a = aset[(aset.model == "PM7std") & (aset.d_min_any >= 2.0)].set_index("sys"); b = aset[aset.model == "PM7"].set_index("sys"); ii = a.index.intersection(b.index)
key("S3.2_std_vs_EuO_noncollapsed", f"{len(ii)} systems; same CN in {int((a.loc[ii, 'CN'] == b.loc[ii, 'CN']).sum())}; same motif in {int((a.loc[ii, 'motif'].fillna('') == b.loc[ii, 'motif'].fillna('')).sum())}; "
    f"mean Eu-O difference (EuO - std) {(b.loc[ii, 'd_all_mean'] - a.loc[ii, 'd_all_mean']).mean():+.3f} A", "3.2: structures of Sparkle/PM7 without collapse compared with Sparkle/PM7+EuO")
for mode in ("bid", "mono2", "hyd"):
    g = cs[(cs["mode"] == mode) & (cs.series == "main")]
    key(f"S3.2_std_{mode}", f"n {len(g)}; collapsed {int((g.d_min_any < 2.0).sum())}; mean shortest Eu-O {g.d_min_any.mean():.2f}", f"3.2: main series, start {mode}")

# ---------------------------------------------------------------- agreement between the two sparkle models
a = aset[aset.model == "PM7"].set_index("sys"); b = aset[aset.model == "SRM1"].set_index("sys"); ii = a.index.intersection(b.index)
dd = b.loc[ii, "d_all_mean"] - a.loc[ii, "d_all_mean"]
key("pair_PM7_SRM1", f"{len(ii)} systems; same CN {int((a.loc[ii, 'CN'] == b.loc[ii, 'CN']).sum())}; same motif {int((a.loc[ii, 'motif'].fillna('') == b.loc[ii, 'motif'].fillna('')).sum())}; "
    f"same number of pectin O {int((a.loc[ii, 'n_pectin'] == b.loc[ii, 'n_pectin']).sum())}; mean Eu-O difference (Sparkle/RM1 - Sparkle/PM7+EuO) {dd.mean():+.3f} \u00b1 {dd.std():.3f} A", "paired comparison of the two sparkle models")
a = aset[aset.model == "RM1"].set_index("sys"); b = aset[aset.model == "SRM1"].set_index("sys"); ii = a.index.intersection(b.index); im = [i for i in ii if a.loc[i, "series"] == "main"]
key("pair_RM1_SRM1", f"{len(im)} main-series systems; same CN {int((a.loc[im, 'CN'] == b.loc[im, 'CN']).sum())}; same motif {int((a.loc[im, 'motif'].fillna('') == b.loc[im, 'motif'].fillna('')).sum())}; "
    f"RM1 has more pectin O in {int((a.loc[im, 'n_pectin'] > b.loc[im, 'n_pectin']).sum())}, fewer in {int((a.loc[im, 'n_pectin'] < b.loc[im, 'n_pectin']).sum())}", "paired comparison RM1 vs Sparkle/RM1 (intact RM1 structures)")

# ---------------------------------------------------------------- Table 3 (by Hamiltonian x start)
t3 = []
for m in PROD:
    for mode, lab in (("bid", "bidentate"), ("mono2", "monodentate + bridge"), ("ester2", "ester carbonyl (M6)"), ("hyd", "hydrated (outer-sphere)")):
        d = (hyd if mode == "hyd" else ok); d = d[(d.model == m) & (d["mode"] == mode)]
        if not len(d):
            t3.append({"Model": NAME[m], "Start": lab, "n": 0}); continue
        other = d.n_pectin - d.n_carboxylate
        row = {"Model": NAME[m], "Start": lab, "n": len(d), "CN (range)": f"{d.CN.mean():.2f} ({int(d.CN.min())}\u2013{int(d.CN.max())})",
               "f_pectin": ms(d.f_pectin), "carboxyl O": f"{d.n_carboxylate.mean():.2f}", "other pectin O": f"{other.mean():.2f}", "water O": f"{d.n_water.mean():.2f}",
               "Eu-O(carb) (A)": f3(d.d_carbox_mean.mean()) if d.d_carbox_mean.notna().any() else "\u2014", "Eu-OH2 (A)": f3(d.d_water_mean.mean())}
        t3.append(row)
        for c in row:
            if c not in ("Model", "Start"):
                key(f"T3_{m}_{mode}_{c}", row[c], f"Table 3: {NAME[m]}, {lab}, {c}")
pd.DataFrame(t3).to_csv(os.path.join(args.out, "Table3_by_start.csv"), index=False)

# ---------------------------------------------------------------- Table 4 (by pattern)
t4 = []
for p in PATTERNS:
    for m in PROD:
        d = ok[(ok.pname == p) & (ok.model == m)]
        if not len(d):
            t4.append({"Pattern": p, "DM (%)": DM[p], "Model": NAME[m], "n": 0}); continue
        row = {"Pattern": p, "DM (%)": DM[p], "Model": NAME[m], "n": len(d), "CN": f"{d.CN.mean():.1f}", "mean f_pectin": f"{d.f_pectin.mean():.3f}",
               "max f_pectin": f"{d.f_pectin.max():.2f}", "multi-COO": int((d.n_carboxylates_bound >= 2).sum()), "max residues": int(d.n_residues_bound.max()),
               "waters retained": f"{int(d.n_water.min())}\u2013{int(d.n_water.max())}"}
        t4.append(row)
        for c in ("n", "CN", "mean f_pectin", "max f_pectin", "multi-COO", "max residues", "waters retained"):
            key(f"T4_{p}_{m}_{c}", row[c], f"Table 4: {p}, {NAME[m]}, {c}")
pd.DataFrame(t4).to_csv(os.path.join(args.out, "Table4_by_pattern.csv"), index=False)

# ---------------------------------------------------------------- Table 5 (retained waters by motif class)
ok["motif_class"] = np.where(ok.n_bidentate_carboxylate > 0, "bidentate carboxylate present", np.where(ok.n_monodentate_carboxylate > 0, "monodentate carboxylate only", "ester / backbone O only (M6)"))
t5 = []
for m in PROD:
    for cl in ("bidentate carboxylate present", "monodentate carboxylate only", "ester / backbone O only (M6)"):
        d = ok[(ok.model == m) & (ok.motif_class == cl)]
        if not len(d):
            continue
        vc = d.n_water.value_counts().sort_index()
        t5.append({"Model": NAME[m], "Motif class": cl, "n": len(d), "median": int(d.n_water.median()), "range": f"{int(d.n_water.min())}\u2013{int(d.n_water.max())}",
                   "distribution": "; ".join(f"{int(k)} H2O: {int(v)}" for k, v in vc.items())})
    d = hyd[hyd.model == m]
    if len(d):
        vc = d.n_water.value_counts().sort_index()
        t5.append({"Model": NAME[m], "Motif class": "hydrated start (no pectin O in first sphere)", "n": len(d), "median": int(d.n_water.median()),
                   "range": f"{int(d.n_water.min())}\u2013{int(d.n_water.max())}", "distribution": "; ".join(f"{int(k)} H2O: {int(v)}" for k, v in vc.items())})
    else:
        t5.append({"Model": NAME[m], "Motif class": "hydrated start (no pectin O in first sphere)", "n": 0})
for row in t5:
    for c in ("n", "median", "range", "distribution"):
        if c in row:
            key(f"T5_{row['Model']}_{row['Motif class'].split()[0]}_{c}", row[c], f"Table 5: {row['Model']}, {row['Motif class']}, {c}")
pd.DataFrame(t5).to_csv(os.path.join(args.out, "Table5_waters.csv"), index=False)

# ---------------------------------------------------------------- Table 6 (matched pairs at equal DM)
t6 = []
for a_, b_, lab in (("H3blk", "H3alt", "DM 0 (protonation variants)"), ("B3", "L3", "DM 50"), ("B2", "D2", "DM 67")):
    for m in PROD:
        d = ok[(ok.model == m) & (ok["mode"].isin(["bid", "mono2"]))]
        A = d[d.pname == a_].set_index(["conf", "mode"]); B = d[d.pname == b_].set_index(["conf", "mode"]); idx = A.index.intersection(B.index)
        if not len(idx):
            t6.append({"Pair": f"{a_} vs {b_} ({lab})", "Model": NAME[m], "pairs": 0}); key(f"T6_{a_}{b_}_{m}_pairs", 0, f"Table 6: {a_} vs {b_}, {NAME[m]}"); continue
        dfp = A.loc[idx].f_pectin - B.loc[idx].f_pectin; dres = A.loc[idx].n_residues_bound - B.loc[idx].n_residues_bound; dco = A.loc[idx].n_carboxylates_bound - B.loc[idx].n_carboxylates_bound
        row = {"Pair": f"{a_} vs {b_} ({lab})", "Model": NAME[m], "pairs": len(idx),
               "dfp >0/=0/<0": f"{int((dfp > 1e-9).sum())} / {int((abs(dfp) <= 1e-9).sum())} / {int((dfp < -1e-9).sum())}",
               "median dfp (range)": f"{dfp.median():+.3f} ({dfp.min():+.3f} to {dfp.max():+.3f})",
               "d residues >0/=0/<0": f"{int((dres > 0).sum())} / {int((dres == 0).sum())} / {int((dres < 0).sum())}", "pairs with an extra carboxylate": int((dco > 0).sum()),
               "multi-COO in block / alternating": f"{int((A.n_carboxylates_bound >= 2).sum())} of {len(A)} / {int((B.n_carboxylates_bound >= 2).sum())} of {len(B)}"}
        t6.append(row)
        for c in row:
            if c not in ("Pair", "Model"):
                key(f"T6_{a_}{b_}_{m}_{c}", row[c], f"Table 6: {a_} vs {b_}, {NAME[m]}, {c}")
pd.DataFrame(t6).to_csv(os.path.join(args.out, "Table6_matched_pairs.csv"), index=False)

# ---------------------------------------------------------------- statements of Section 3
bid = ok[ok["mode"] == "bid"]; m2 = ok[ok["mode"] == "mono2"]
kept = m2[(m2.n_monodentate_carboxylate >= 1) & (m2.n_bidentate_carboxylate == 0)]
key("bid_n", f"{len(bid)} (" + ", ".join(f"{int((bid.model == m).sum())} {NAME[m]}" for m in PROD) + ")", "bidentate starts in the analysis set")
key("bid_preserved", f"{int((bid.n_bidentate_carboxylate >= 1).sum())} of {len(bid)}", "bidentate starts that kept the chelate")
key("bid_dcarb_all", ms(bid.d_carbox_mean), "Eu-O(carboxylate), bidentate starts, mean \u00b1 SD")
for m in PROD:
    g = bid[bid.model == m]
    key(f"bid_{m}", f"n {len(g)}; Eu-O(carb) {g.d_carbox_mean.mean():.3f}; CN {g.CN.mean():.2f}; waters {int(g.n_water.min())}-{int(g.n_water.max())}; Eu-OH2 {g.d_water_mean.mean():.3f}; shortest Eu-O {g.d_min_any.min():.2f}", f"bidentate starts, {NAME[m]}")
key("bid_CN", f"{bid.CN.mean():.1f}", "mean CN, bidentate starts")
key("mono2_n", f"{len(m2)} (" + ", ".join(f"{int((m2.model == m).sum())} {NAME[m]}" for m in PROD) + ")", "monodentate-plus-bridge starts in the analysis set")
key("mono2_kept", len(kept), "retained a single monodentate carboxylate (no bidentate)"); key("mono2_closed", int((m2.n_bidentate_carboxylate >= 1).sum()), "closed to a bidentate chelate")
key("mono2_dcarb", ms(kept.d_carbox_mean), "Eu-O(carboxylate), monodentate structures"); key("bid_minus_mono", f"{bid.d_carbox_mean.mean() - kept.d_carbox_mean.mean():.2f}", "bidentate minus monodentate Eu-O (A)")
for m in PROD:
    g = m2[m2.model == m]; k = kept[kept.model == m]
    key(f"mono2_{m}", f"n {len(g)}; kept mono {len(k)}; closed bi {int((g.n_bidentate_carboxylate >= 1).sum())}; Eu-O(carb, mono) {k.d_carbox_mean.mean():.3f}; ring/glycosidic O in {int((g.n_backbone_O > 0).sum())}; "
        f"hydroxyl O in {int((g.n_hydroxyl_O > 0).sum())}; second free carboxylate in {int((g.n_carboxylates_bound >= 2).sum())}; mean non-carboxylate pectin O {(g.n_pectin - g.n_carboxylate).mean():.2f}", f"monodentate-plus-bridge starts, {NAME[m]}")
key("mono2_backbone", int((m2.n_backbone_O > 0).sum()), "mono2 structures with a ring/glycosidic O"); key("mono2_hydroxyl", int((m2.n_hydroxyl_O > 0).sum()), "mono2 structures with a hydroxyl O"); key("mono2_2COO", int((m2.n_carboxylates_bound >= 2).sum()), "mono2 structures with a second free carboxylate")
b3 = aset[aset.job == "RM1_B3_h2_mono2"].iloc[0]; b3s = aset[aset.job == "SRM1_B3_h2_mono2"].iloc[0]; b3p = aset[aset.job == "PM7_B3_h2_mono2"].iloc[0]
for lab, x in (("RM1", b3), ("SRM1", b3s), ("PM7", b3p)):
    key(f"B3_h2_mono2_{lab}", f"CN {int(x.CN)}, n_pectin {int(x.n_pectin)}, n_water {int(x.n_water)}, f {x.f_pectin:.2f}, {x.motif}", f"B3_h2_mono2, {NAME[lab]}")
key("ok_n", f"{len(ok)} (" + ", ".join(f"{int((ok.model == m).sum())} {NAME[m]}" for m in PROD) + ")", "inner-sphere structures in the analysis set")
key("ok_CN_lt_6", int((ok.CN < 6).sum()), "inner-sphere structures with CN < 6")
key("ok_n_pectin", "; ".join(f"{k}: {v}" for k, v in ok.n_pectin.value_counts().sort_index().items()), "number of pectin O (count of structures)")
key("ok_n_pectin_median", int(ok.n_pectin.median()), "median number of pectin O"); key("ok_fpectin_range", f"{ok.f_pectin.min():.2f}\u2013{ok.f_pectin.max():.2f}", "f_pectin range")
key("ok_n_water", "; ".join(f"{k}: {v}" for k, v in ok.n_water.value_counts().sort_index().items()), "retained waters (count of structures)")
for m in PROD:
    g = ok[ok.model == m]
    key(f"ok_{m}", f"n {len(g)}; n_pectin " + "; ".join(f"{k}: {v}" for k, v in g.n_pectin.value_counts().sort_index().items()) + f" | f {g.f_pectin.min():.2f}-{g.f_pectin.max():.2f} | waters " + "; ".join(f"{k}: {v}" for k, v in g.n_water.value_counts().sort_index().items()), f"composition, {NAME[m]}")
key("hyd_n", f"{len(hyd)} (" + ", ".join(f"{int((hyd.model == m).sum())} {NAME[m]}" for m in PROD) + ")", "hydrated starts in the analysis set")
key("hyd_max_fpectin", f"{hyd.f_pectin.max():.3f}", "largest f_pectin among hydrated starts"); key("hyd_CN", "; ".join(f"{k}: {v}" for k, v in hyd.CN.value_counts().sort_index().items()), "CN of hydrated starts")
for m in PROD:
    h = hyd[hyd.model == m]
    key(f"hyd_outer_contacts_{m}", f"{int(((h['n_pectin_3.6'] - h.n_pectin) >= 1).sum())} of {len(h)}", f"hydrated starts with pectin O at 3.0-3.6 A, {NAME[m]}")
for m in PROD:
    d = ok[ok.model == m]; rho, p = spearmanr(d.pname.map(DM), d.f_pectin)
    key(f"spearman_{m}", f"{rho:.2f} (p = {p:.2f}, n = {len(d)})", f"Spearman DM vs f_pectin, {NAME[m]}")
    pm = d[d.pname != "M6"].groupby("pname").f_pectin.mean(); key(f"pattern_mean_range_{m}", f"{pm.min():.3f}\u2013{pm.max():.3f}", f"range of pattern means (carboxylate patterns), {NAME[m]}")
    key(f"multiCOO_{m}", f"{int((d.n_carboxylates_bound >= 2).sum())} of {len(d)}", f"structures with two or more free carboxylates bound, {NAME[m]}")
m6 = aset[(aset.pname == "M6") & (aset["mode"] == "ester2")]
for m in ("PM7", "SRM1", "RM1", "PM7std"):
    g = m6[m6.model == m]
    key(f"M6_{m}", f"n {len(g)}; CN {int(g.CN.min())}\u2013{int(g.CN.max())}; f_pectin {g.f_pectin.min():.2f}\u2013{g.f_pectin.max():.2f}; waters {int(g.n_water.min())}\u2013{int(g.n_water.max())}; ester C=O in first sphere {g.n_carbonyl_ester.tolist()}; "
        f"shortest Eu-O {g.d_min_any.min():.2f}\u2013{g.d_min_any.max():.2f}; motifs {g.motif.tolist()}", f"M6 corrected ester start, {NAME[m]}")
allm = pd.concat([ok, hyd])
key("water_EuO", ms(allm.d_water_mean), "Eu-OH2, main-series analysis set (per-structure means)")
for m in PROD:
    g = allm[allm.model == m]; key(f"water_EuO_{m}", f"{g.d_water_mean.mean():.3f}", f"Eu-OH2, {NAME[m]}")
    gb = bid[bid.model == m]; gk = kept[kept.model == m]
    key(f"dcarb_{m}", f"bidentate {gb.d_carbox_mean.mean():.3f} (n {len(gb)}); monodentate {gk.d_carbox_mean.mean():.3f} (n {len(gk)})", f"Eu-O(carboxylate), {NAME[m]}")
    go = ok[ok.model == m]; gh = hyd[hyd.model == m]
    key(f"chain_{m}", f"inner-sphere: rmsd {go.chain_rmsd.mean():.2f}, max torsion change {go.max_dtorsion.mean():.1f}; hydrated: rmsd {gh.chain_rmsd.mean() if len(gh) else float('nan'):.2f}", f"chain deformation, {NAME[m]}")
key("chain_rmsd", f"{ok.chain_rmsd.mean():.2f}", "mean chain RMSD, inner-sphere structures"); key("chain_dtors", f"{ok.max_dtorsion.mean():.1f}", "mean largest torsion change"); key("chain_rmsd_hyd", f"{hyd.chain_rmsd.mean():.2f}", "mean chain RMSD, hydrated starts")
sens = aset[aset.series == "sens_hydration"]
key("sens_in_set", ", ".join(sorted(sens.job)), "hydration-sensitivity jobs in the analysis set")
pairs = []
for r4 in sens[sens.job.str.endswith("cw4")].itertuples():
    b0 = sens[sens.job == r4.job.replace("cw4", "cw0")]; bm = aset[aset.job == r4.job.replace("_cw4", "")]
    if len(b0):
        b0 = b0.iloc[0]; pairs.append(f"{r4.job}: motif {'same' if r4.motif == b0.motif else 'DIFFERENT'}, f {r4.f_pectin:.3f} vs {b0.f_pectin:.3f}, dEuO {abs(r4.d_carbox_mean - b0.d_carbox_mean):.3f}")
key("sens_pairs", " | ".join(pairs), "cw4 variant vs cw0 baseline"); key("sens_n_pairs", len(pairs), "number of cw4/cw0 pairs in the analysis set")
q = {j: geoms[j]["charges"][geoms[j]["el"].index("Eu")] for j in ("RM1_B3_h2_mono2", "RM1_CTRL_Eu_w9") if geoms[j].get("charges")}
key("qEu_B3", f"{q.get('RM1_B3_h2_mono2', np.nan):+.2f}", "net atomic charge on Eu, RM1_B3_h2_mono2"); key("qEu_aqua", f"{q.get('RM1_CTRL_Eu_w9', np.nan):+.2f}", "net atomic charge on Eu, RM1_CTRL_Eu_w9")
rm = conv[conv.model == "RM1"]
key("qEu_RM1_intact", f"{rm[rm.intact == True].q_Eu.min():.2f} to {rm[rm.intact == True].q_Eu.max():.2f}", "net atomic charge of Eu, intact RM1 structures"); key("qEu_RM1_broken", f"{rm[rm.intact != True].q_Eu.min():.2f} to {rm[rm.intact != True].q_Eu.max():.2f}", "net atomic charge of Eu, RM1 structures that fail the integrity criterion")

# ---------------------------------------------------------------- vibrational analysis and displacement test
for m in PROD:
    v = vib[vib.model == m]; good = v[v.sumrule_max <= 0.02]
    key(f"vib_{m}", f"n {len(v)}; no imaginary {int((v.n_imag == 0).sum())}; none below 20i {int((v.n_imag_gt20 == 0).sum())}; none below 50i {int((v.n_imag_gt50 == 0).sum())}; lowest {v.f1.min():.0f}; "
        f"sum rule median {v.sumrule_max.median():.4f}, max {v.sumrule_max.max():.3f}; Hessians within 0.02: {len(good)}, of which none below 20i: {int((good.n_imag_gt20 == 0).sum())}", f"constrained vibrational analysis, {NAME[m]}")
key("vib_unconstrained_neg", f"{int(vib.n_mopac_neg.min())}\u2013{int(vib.n_mopac_neg.max())}", "imaginary frequencies per structure in the standard (unconstrained) analysis")
dc = disp[disp.status == "converged"]
for m in ("PM7", "SRM1", "RM1"):
    g = dc[dc.model == m]; a_ = disp[disp.model == m]
    ret = int(((g.dE.abs() <= 0.1) & (g.rmsd_heavy_vs_parent <= 0.05)).sum())
    key(f"disp_{m}", f"{len(a_)} re-optimizations of {a_.parent.nunique()} structures; converged {len(g)}; returned (|dE| <= 0.1, rmsd <= 0.05) {ret}; lower by > 0.1: {int((g.dE < -0.1).sum())}; lower by > 1: {int((g.dE < -1).sum())}; "
        f"dE min {g.dE.min():.2f}, median {g.dE.median():.2f}; same CN {int((g.CN == g.parent_CN).sum())}; same motif {int((g.motif.fillna('') == g.parent_motif.fillna('')).sum())}; same first-sphere atoms {int(g.same_first_sphere_atoms.sum())}; "
        f"largest first-sphere Eu-O change {g.max_dEuO_first_sphere.max():.3f} A; largest change of mean Eu-O {(g.d_all_mean - g.parent_d_all_mean).abs().max():.3f} A", f"displacement test, {NAME[m]}")
key("disp_all", f"{len(disp)} re-optimizations; {len(dc)} converged; first sphere unchanged in {int(dc.same_first_sphere_atoms.sum())}; returned {int(((dc.dE.abs() <= 0.1) & (dc.rmsd_heavy_vs_parent <= 0.05)).sum())}; lower by > 0.1 kcal/mol {int((dc.dE < -0.1).sum())}", "displacement test, all")

# ---------------------------------------------------------------- supplementary tables
nc = df[(df.status != "converged") | (df.intact != True)]
pd.DataFrame({"Job": nc.job, "Model": nc.model.map(NAME), "Status": nc.status, "Integrity": nc.intact.map({True: "intact", False: "not intact"}), "CN": nc.CN, "shortest Eu-O (A)": nc.d_min_any.map(f3),
              "CO2 lost": nc.n_CO2_lost, "Motif": nc.motif.fillna("\u2014")}).to_csv(os.path.join(args.out, "TableS2_excluded.csv"), index=False)
c = df[df.series == "control"]
pd.DataFrame({"Job": c.job, "Model": c.model.map(NAME), "System": c.pname, "Start": c["mode"], "q": c.charge, "Status": c.status, "Intact": c.intact, "CN": c.CN.astype(int),
              "Eu-OH2 mean (A)": c.d_water_mean.map(f3), "Eu-O(carb) (A)": c.d_carbox_mean.map(f3), "shortest Eu-O (A)": c.d_min_any.map(f3), "dHf (kcal/mol)": c.heat.map(lambda x: f"{x:.2f}")}).to_csv(os.path.join(args.out, "TableS4_controls.csv"), index=False)
h5 = df[df.series == "sens_hydration"]
pd.DataFrame({"Job": h5.job, "Model": h5.model.map(NAME), "Pattern": h5.pname, "Variant": h5.job.str.endswith("cw4").map({True: "cw4", False: "cw0"}), "In analysis set": h5.analysis_set,
              "CN": h5.CN, "f_pectin": h5.f_pectin.map(f3), "Eu-O(carb) (A)": h5.d_carbox_mean.map(f3), "Motif": h5.motif}).to_csv(os.path.join(args.out, "TableS5_sensitivity.csv"), index=False)
kdf = pd.DataFrame(KEY, columns=["id", "value", "description"]); kdf.to_csv(os.path.join(args.out, "key_numbers_v2.csv"), index=False)
print(f"tables and {len(kdf)} key numbers written to {args.out}")
