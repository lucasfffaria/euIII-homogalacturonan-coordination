"""Rigid ligand geometries and Eu approach vectors for the Eu-O calibration scans (shared by MOPAC and DFT)."""
import numpy as np


def _rdkit_geom(smiles):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    conf = m.GetConformer()
    el = [a.GetSymbol() for a in m.GetAtoms()]
    xyz = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z] for i in range(m.GetNumAtoms())])
    return m, el, xyz


def _unit(v):
    return v / np.linalg.norm(v)


def ligands():
    """Return dict name -> (el, xyz, donor O index or list, charge, approach unit vector(s)).
    Geometry convention: donor O at origin; Eu placed at  r * approach  (bidentate acetate: Eu on the C2 axis, r = Eu-O)."""
    L = {}
    # water: idealized (O-H 0.962 A, 104.4 deg), Eu on the bisector
    L["water"] = dict(el=["O", "H", "H"], xyz=np.array([[0, 0, 0], [0.76, 0.59, 0], [-0.76, 0.59, 0]], float), donor=0, charge=0,
                      approach=np.array([0.0, -1.0, 0.0]), label="H2O (neutral donor)")
    # methanol: same rigid geometry as the earlier scans (Eu-O-C 127 deg, in the H-O-C plane, anti to H)
    el = ["O", "C", "H", "H", "H", "H"]
    xyz = np.array([[0, 0, 0], [-1.42, 0, 0], [0.32, 0.905, 0], [-1.80, -0.50, 0.85], [-1.80, -0.50, -0.85], [-1.80, 1.00, 0]], float)
    L["methanol"] = dict(el=el, xyz=xyz, donor=0, charge=0, approach=_unit(np.array([1.44, -1.92, 0.0])), label="CH3OH (neutral hydroxyl)")
    # methoxide: methanol without the hydroxyl H (localized anionic O; stress test for collapse)
    L["methoxide"] = dict(el=["O", "C", "H", "H", "H"], xyz=xyz[[0, 1, 3, 4, 5]], donor=0, charge=-1, approach=_unit(np.array([1.44, -1.92, 0.0])), label="CH3O- (localized anionic O)")
    # acetate, bidentate: C at origin, O at (0.63, +-1.10); Eu on the C2 axis at x such that Eu-O = r
    el = ["C", "O", "O", "C", "H", "H", "H"]
    xyz = np.array([[0, 0, 0], [0.63, 1.10, 0], [0.63, -1.10, 0], [-1.52, 0, 0], [-1.90, 1.02, 0], [-1.90, -0.51, 0.88], [-1.90, -0.51, -0.88]], float)
    L["acetate"] = dict(el=el, xyz=xyz, donor=[1, 2], charge=-1, approach=None, label="CH3COO- bidentate (delocalized anion)")
    # dimethyl ether (ring / glycosidic O model): MMFF geometry, Eu on the C-O-C bisector, away from the carbons
    m, el, xyz = _rdkit_geom("COC")
    iO = el.index("O"); iC = [i for i in range(len(el)) if el[i] == "C"]
    xyz = xyz - xyz[iO]
    bis = _unit(_unit(xyz[iC[0]]) + _unit(xyz[iC[1]]))
    L["ether"] = dict(el=el, xyz=xyz, donor=iO, charge=0, approach=-bis, label="(CH3)2O (ether: ring / glycosidic O)")
    # methyl acetate (ester carbonyl O): MMFF geometry; Eu in the O=C-O plane, C=O...Eu = 140 deg, on the side away from the ether O
    m, el, xyz = _rdkit_geom("CC(=O)OC")
    # identify carbonyl O (double bond) and ether O
    iOc = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "O" and any(b.GetBondTypeAsDouble() == 2 for b in a.GetBonds())][0]
    iOe = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "O" and a.GetIdx() != iOc][0]
    iCc = [n.GetIdx() for n in m.GetAtomWithIdx(iOc).GetNeighbors()][0]
    xyz = xyz - xyz[iOc]
    v1 = _unit(xyz[iCc])                      # O -> C(carbonyl)
    w = xyz[iOe] - xyz[iCc]                   # C -> ether O
    nrm = _unit(np.cross(v1, w))               # plane normal
    pp = _unit(np.cross(nrm, v1))              # in-plane perpendicular to C=O
    if np.dot(pp, w) > 0:                      # point away from the ether O
        pp = -pp
    th = np.radians(140.0)
    L["ester"] = dict(el=el, xyz=xyz, donor=iOc, charge=0, approach=np.cos(th) * v1 + np.sin(th) * pp, label="CH3COOCH3 (ester carbonyl O)")
    return L


def place_eu(lig, r):
    """Return (el, xyz) of ligand + Eu with Eu-O distance r."""
    el = list(lig["el"]) + ["Eu"]
    if lig["approach"] is None:               # bidentate acetate: Eu on the x axis
        y = lig["xyz"][1][1]; x0 = lig["xyz"][1][0]
        eu = np.array([x0 + np.sqrt(max(r * r - y * y, 0.0)), 0.0, 0.0])
    else:
        eu = r * lig["approach"]
    return el, np.vstack([lig["xyz"], eu])


R_GRID = [1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6]
R_GRID_DFT = [1.5, 1.7, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 3.0, 3.2, 3.6]

if __name__ == "__main__":
    L = ligands()
    for k, v in L.items():
        el, xyz = place_eu(v, 2.4)
        d = [np.linalg.norm(xyz[-1] - xyz[i]) for i in range(len(el) - 1)]
        print(k, "Eu-O(donor) =", round(min(d), 3), "min Eu-X =", round(min(d), 3), "closest non-donor", round(sorted(d)[1 if v["approach"] is not None else 2], 3))
