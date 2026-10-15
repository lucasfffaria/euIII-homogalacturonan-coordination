"""Chemical-integrity criterion of the final structures (Section 2.5 of the article; Tables S2 and S9).

A structure passes when the covalent connectivity of its heavy atoms is that of its starting structure (<job>_start.xyz).
Bond criterion: distance < 1.25 x (sum of covalent radii); bonds to Eu are not considered.

This module holds the bond criterion (bond_matrix) and the XYZ reader.  The criterion is applied to every optimization by
reproduce_descriptors_v2.py, which writes data/TableS9_integrity.csv: for every job, the broken and the new heavy-atom
bonds and, for the residues whose C5-C6 bond is broken (loss of CO2 from a free carboxylate), the geometry and charge of
the detached unit, its distance to Eu and the net atomic charge of Eu.
"""
import numpy as np

COV = {"H": 0.31, "C": 0.76, "O": 0.66, "Eu": 1.98}
FACTOR = 1.25


def read_xyz(path):
    L = open(path).read().splitlines(); n = int(L[0])
    return [l.split()[0] for l in L[2:2 + n]], np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])


def bond_matrix(el, xyz):
    el = np.array(el); R = np.array([COV[e] for e in el])
    D = np.linalg.norm(xyz[:, None] - xyz[None], axis=2)
    M = (D < FACTOR * (R[:, None] + R[None])) & ~np.eye(len(el), dtype=bool)
    M[el == "Eu", :] = False; M[:, el == "Eu"] = False
    return M, D
