"""The Eu-O core-core pair term of Sparkle/PM7 and the two-parameter calibration map (Section 2.3 and Figure 4 of the article,
Section S3 and Table S7 of the Supplementary Information).

At fixed geometry a sparkle has no orbitals, so a change of the Eu-O pair parameters XFAC and ALPB changes only the
core-core term   E(r) = Z(Eu) Z(O) (ss|ss) [1 + 2 XFAC exp(-ALPB (r + 0.0003 r^6))].
The scan of any pair is therefore the archived Sparkle/PM7 scan plus an analytic difference.  The script
  1. checks this relation point by point: difference between the archived scans of Sparkle/PM7+EuO and of Sparkle/PM7 as
     distributed (dE_MOPAC) against the analytic difference of the two pair terms (dE_analytic) (check_fallback_vs_scans.csv);
     the reference pair is the one MOPAC 23.2.5 assigns at run time when no Eu-O pair is defined (half of the O-O values);
  2. computes, for a grid of (XFAC, ALPB), the position of the minimum of the rigid Eu-OH2 and Eu-OCH3(-) scans
     (fig05b_calibration_map.csv), which is the map of Figure 4b.

Usage (from code/):  python calibration_map.py [--data ../data] [--out calibration_out] [--compare]
Needs numpy, pandas and scipy.
"""
import argparse, json, os
import numpy as np, pandas as pd
from scipy.interpolate import CubicSpline

A0B, EV, KC = 0.529177210903, 27.2113862, 23.060548
PO9 = {"Eu": 0.23799411, "O": 0.90976817}        # Klopman-Ohno core terms (bohr), MOPAC 23.2.5, PM7
ZZ = 3 * 6                                        # core charges of Eu and O
X_OO, A_OO = 0.356439, 2.457135                   # PM7 O-O pair parameters
X0, A0 = 0.5 * X_OO, 0.5 * A_OO                   # pair assigned at run time to Eu-O (ccrep.F90)
X1, A1 = 2.2075, 2.5000                           # pair of this work (eu_fix.txt)


def gab(r):
    return EV / np.sqrt((r / A0B) ** 2 + (PO9["Eu"] + PO9["O"]) ** 2)


def pair(r, x, a):
    return KC * ZZ * gab(r) * 2 * x * np.exp(-a * (r + 0.0003 * r ** 6))


def main():
    here = os.path.dirname(os.path.abspath(__file__)); ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=os.path.join(here, "..", "data")); ap.add_argument("--out", default=os.path.join(os.getcwd(), "calibration_out"))
    ap.add_argument("--compare", action="store_true"); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    sc = json.load(open(os.path.join(a.data, "scans_semiempirical_v2.json"))); r = np.array(sc["r"])
    binding = lambda model, lig: np.array(sc[model][lig]) - sc[model]["Eu3+"] - sc[model][lig + "_ligand"]
    ligs = ("water", "methanol", "ether", "methoxide"); be0 = {l: binding("Sparkle/PM7", l) for l in ligs}
    rows = []
    for l in ligs:
        d_mopac = binding("Sparkle/PM7*", l) - be0[l]; d_analytic = pair(r, X1, A1) - pair(r, X0, A0)
        rows += [dict(ligand=l, r=ri, dE_MOPAC=m, dE_analytic=an, resid=m - an) for ri, m, an in zip(r, d_mopac, d_analytic)]
    chk = pd.DataFrame(rows).round(6); chk.to_csv(os.path.join(a.out, "check_fallback_vs_scans.csv"), index=False)
    print(f"run-time pair: XFAC = {X0:.7f}, ALPB = {A0:.7f}; {len(chk)} scan points; largest deviation between the difference of the two archived "
          f"scans and the analytic difference: {chk.resid.abs().max():.3f} kcal/mol, over a range of {chk.dE_MOPAC.max() - chk.dE_MOPAC.min():.0f} kcal/mol")
    RF = np.linspace(1.3, 3.5, 551)

    def rmin(lig, x, al):
        y = CubicSpline(r, be0[lig] + pair(r, x, al) - pair(r, X0, A0))(RF)
        return RF[y.argmin()], y.min()
    xs = np.arange(0.0625, 6.0001, 0.0625); als = np.arange(1.0, 4.0001, 0.025)
    g = pd.DataFrame([(x, al, rmin("water", x, al)[0], rmin("methoxide", x, al)[0]) for al in als for x in xs], columns=["XFAC", "ALPB", "rmin_water_A", "rmin_methoxide_A"]).round(4)
    g.to_csv(os.path.join(a.out, "fig05b_calibration_map.csv"), index=False)
    for lab, (x, al) in {"run-time pair": (X0, A0), "this work": (X1, A1)}.items():
        print(f"{lab:14s}", "  ".join(f"{l}: {rmin(l, x, al)[0]:.2f} A / {rmin(l, x, al)[1]:.0f} kcal/mol" for l in ligs))
    if a.compare:
        for name in ("fig05b_calibration_map.csv", "check_fallback_vs_scans.csv"):
            p, q = pd.read_csv(os.path.join(a.out, name)), pd.read_csv(os.path.join(a.data, name))
            num = p.select_dtypes("number").columns
            print(f"{name}: identical values = {p.shape == q.shape and list(p.columns) == list(q.columns) and bool(np.allclose(p[num], q[num], rtol=0, atol=2e-6))}")


if __name__ == "__main__":
    main()
