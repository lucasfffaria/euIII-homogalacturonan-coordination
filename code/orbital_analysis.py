"""Frontier-orbital numbers quoted in Section 3.12 of the article and in the captions of Figures 18 and 19.

Reads the two RM1 single points of mopac_single_points/ (B3mono2_sp: Eu(III) bound to pattern B3, job RM1_B3_h2_mono2;
aqua_sp: [Eu(H2O)9]3+, job RM1_CTRL_Eu_w9; keywords 1SCF GRAPHF AUX VECTORS ALLVEC MULLIK on the optimized geometries)
and prints the HOMO and LUMO energies, the weight of the Eu atomic orbitals in each of them (sum of the squared
coefficients in the orthogonalized basis) and the net atomic charge of Eu.  No MOPAC run and no rendering.

Usage (from code/):  python orbital_analysis.py [--sm ..]        Needs numpy.
"""
import argparse, os, sys
import numpy as np

here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
from orbitals import read_basis, read_mgf

ap = argparse.ArgumentParser(); ap.add_argument("--sm", default=os.path.join(here, "..")); a = ap.parse_args()
for folder, label in (("B3mono2_sp", "Eu(III)-B3 complex (RM1_B3_h2_mono2)"), ("aqua_sp", "[Eu(H2O)9]3+ (RM1_CTRL_Eu_w9)")):
    d = os.path.join(a.sm, "mopac_single_points", folder)
    b = read_basis(os.path.join(d, "sp.aux")); m = read_mgf(os.path.join(d, "sp.mgf"))
    el = list(b["el"]); eu = el.index("Eu"); ao_atom = np.asarray(b["ao_atom"])
    nocc = int(round(m["occ"].sum() / 2)); homo, lumo = nocc - 1, nocc
    w = [float(np.sum(m["C"][k][ao_atom == eu] ** 2)) for k in (homo, lumo)]
    print(f"{label}: HOMO {m['eig'][homo]:.3f} eV, LUMO {m['eig'][lumo]:.3f} eV; Eu weight in HOMO {w[0]:.3f}, in LUMO {w[1]:.3f}; "
          f"net atomic charge of Eu {b['charges'][eu]:+.3f}")
