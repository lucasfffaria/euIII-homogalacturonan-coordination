"""
Geometric assembly of alpha-(1->4)-D-galacturonan oligomers from a 4C1 chair monomer,
with user-defined glycosidic torsions (phi = O5(i)-C1(i)-O4(i+1)-C4(i+1),
psi = C1(i)-O4(i+1)-C4(i+1)-C5(i+1)).
"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms
from rdkit.Geometry import Point3D
from galA import hg_mol, residue_atom_map, ring_axial_check, glycosidic_torsions, TAGS

CHAIR_4C1 = {"O1": (0.75, 1.0), "O2": (0.0, 0.5), "O3": (0.0, 0.5), "O4": (0.75, 1.0), "C6": (0.0, 0.5)}


def is_4C1(row):
    return all(lo <= row[k] <= hi for k, (lo, hi) in CHAIR_4C1.items())


def monomer_4C1(kind="H", seed=0, ntry=300):
    """Return an RDKit mol (with H) of a single GalpA residue in the 4C1 chair, MMFF-optimised."""
    m = hg_mol(kind)
    cids = AllChem.EmbedMultipleConfs(m, numConfs=ntry, randomSeed=seed)
    props = AllChem.MMFFGetMoleculeProperties(m)
    best = None
    for cid in cids:
        ff = AllChem.MMFFGetMoleculeForceField(m, props, confId=cid)
        ff.Minimize(maxIts=2000)
        e = ff.CalcEnergy()
        rm = residue_atom_map(m)
        row = ring_axial_check(m, m.GetConformer(cid), rm)[0]
        if is_4C1(row) and (best is None or e < best[0]):
            best = (e, cid)
    if best is None:
        raise RuntimeError("no 4C1 chair found")
    conf = m.GetConformer(best[1])
    mm = Chem.Mol(m, confId=best[1])
    mm.RemoveAllConformers()
    c = Chem.Conformer(m.GetNumAtoms())
    for i in range(m.GetNumAtoms()):
        c.SetAtomPosition(i, conf.GetAtomPosition(i))
    mm.AddConformer(c, assignId=True)
    return mm


def place_atom(a, b, c, bond, angle_deg, dihedral_deg):
    """Position of atom d such that |cd|=bond, angle(b,c,d)=angle, dihedral(a,b,c,d)=dihedral."""
    ang = np.deg2rad(angle_deg)
    dih = np.deg2rad(dihedral_deg)
    bc = c - b
    bc /= np.linalg.norm(bc)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    d2 = np.array([-bond * np.cos(ang), bond * np.sin(ang) * np.cos(dih), bond * np.sin(ang) * np.sin(dih)])
    return c + d2[0] * bc + d2[1] * m + d2[2] * n


def kabsch(P, Q):
    """Return R, t such that R @ p + t ~ q for rows of P, Q."""
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = Qc - R @ Pc
    return R, t


def heavy_and_H_names(mol, rmap):
    """Map each atom index of a residue to a name: heavy atoms by key; H's as key+'_H'+n."""
    names = {}
    for key, idx in rmap.items():
        if key == "kind" or idx is None:
            continue
        names[idx] = key
        hs = [n.GetIdx() for n in mol.GetAtomWithIdx(idx).GetNeighbors() if n.GetSymbol() == "H"]
        for k, h in enumerate(hs):
            names[h] = f"{key}_H{k}"
    return names


def assemble(pattern, torsions, monomer_cache={}, tau=116.5, bond_C1O4=1.42):
    """Build oligomer mol (with H) with 3D coordinates. torsions: list of (phi, psi) per linkage."""
    n = len(pattern)
    assert len(torsions) == n - 1
    oligo = hg_mol(pattern)
    orm = residue_atom_map(oligo)
    assert [d["kind"] for d in orm] == list(pattern), ([d["kind"] for d in orm], pattern)
    pos = np.zeros((oligo.GetNumAtoms(), 3))
    placed = np.zeros(oligo.GetNumAtoms(), bool)
    prev = None
    for i, kind in enumerate(pattern):
        if kind not in monomer_cache:
            monomer_cache[kind] = monomer_4C1(kind, seed=11)
        mono = monomer_cache[kind]
        mrm = residue_atom_map(mono)[0]
        mnames = heavy_and_H_names(mono, mrm)
        mpos = mono.GetConformer().GetPositions()
        name_to_mono = {v: k for k, v in mnames.items()}
        onames = heavy_and_H_names(oligo, orm[i])
        name_to_oligo = {v: k for k, v in onames.items()}
        if i == 0:
            R, t = np.eye(3), np.zeros(3)
        else:
            phi, psi = torsions[i - 1]
            p_O5 = pos[orm[i - 1]["O5"]]
            p_C1 = pos[orm[i - 1]["C1"]]
            p_O4 = pos[orm[i]["O4"]]  # already placed as O1 of previous residue
            # target C4 and C5 of residue i
            # bond O4-C4 and angle O4-C4-C5 from the monomer
            m_O4, m_C4, m_C5 = mpos[mrm["O4"]], mpos[mrm["C4"]], mpos[mrm["C5"]]
            b_O4C4 = np.linalg.norm(m_C4 - m_O4)
            b_C4C5 = np.linalg.norm(m_C5 - m_C4)
            v1, v2 = m_O4 - m_C4, m_C5 - m_C4
            ang_O4C4C5 = np.degrees(np.arccos(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)))
            p_C4 = place_atom(p_O5, p_C1, p_O4, b_O4C4, tau, phi)
            p_C5 = place_atom(p_C1, p_O4, p_C4, b_C4C5, ang_O4C4C5, psi)
            R, t = kabsch(np.array([m_O4, m_C4, m_C5]), np.array([p_O4, p_C4, p_C5]))
        for name, oidx in name_to_oligo.items():
            if placed[oidx]:
                continue
            if name in name_to_mono:
                pos[oidx] = R @ mpos[name_to_mono[name]] + t
                placed[oidx] = True
            else:
                raise KeyError(name)
        # the O1 of residue i is the O4 of residue i+1 (shared atom) -> set from monomer O1 position
        if i < n - 1:
            o1 = orm[i]["O1"]
            pos[o1] = R @ mpos[mrm["O1"]] + t
            placed[o1] = True
            # ensure the next residue's O4 index is the same atom
            assert orm[i + 1]["O4"] == o1
    assert placed.all(), np.where(~placed)
    conf = Chem.Conformer(oligo.GetNumAtoms())
    for k in range(oligo.GetNumAtoms()):
        conf.SetAtomPosition(k, Point3D(*pos[k]))
    oligo.RemoveAllConformers()
    oligo.AddConformer(conf, assignId=True)
    return oligo


def check_stereo(mol):
    """Compare stereo from 3D with the intended SMILES stereo."""
    m3 = Chem.Mol(mol)
    Chem.AssignStereochemistryFrom3D(m3)
    s3 = Chem.MolToSmiles(Chem.RemoveHs(m3))
    s0 = Chem.MolToSmiles(Chem.RemoveHs(mol))
    return s3 == s0


def mmff_relax(mol, maxIters=5000, dielectric=4.0, distance_dependent=True, fixed_torsions=None):
    props = AllChem.MMFFGetMoleculeProperties(mol)
    props.SetMMFFDielectricModel(2 if distance_dependent else 1)
    props.SetMMFFDielectricConstant(dielectric)
    ff = AllChem.MMFFGetMoleculeForceField(mol, props)
    if fixed_torsions:
        for (a, b, c, d, val) in fixed_torsions:
            ff.MMFFAddTorsionConstraint(a, b, c, d, False, val - 0.5, val + 0.5, 500.0)
    ff.Initialize()
    ff.Minimize(maxIts=maxIters)
    return ff.CalcEnergy()


def helix_params(mol, rmap):
    """Residues per turn and rise per residue from successive residue frames (C1 positions)."""
    pos = mol.GetConformer().GetPositions()
    # use screw transformation between residue frames (C1,C4,O5)
    frames = []
    for d in rmap:
        P = np.array([pos[d["C1"]], pos[d["C4"]], pos[d["O5"]]])
        frames.append(P)
    Rs, ts = [], []
    for i in range(len(frames) - 1):
        R, t = kabsch(frames[i], frames[i + 1])
        Rs.append(R)
        ts.append(t)
    # screw angle from rotation matrices
    angles = [np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))) for R in Rs]
    ang = np.mean(angles)
    # rise: projection of displacement of residue centroid onto axis
    axes = []
    for R in Rs:
        w, v = np.linalg.eig(R)
        ax = np.real(v[:, np.argmin(np.abs(w - 1))])
        axes.append(ax / np.linalg.norm(ax))
    cents = [frames[i].mean(0) for i in range(len(frames))]
    rises = [abs(np.dot(cents[i + 1] - cents[i], axes[i])) for i in range(len(Rs))]
    return 360.0 / ang if ang > 1e-6 else np.inf, float(np.mean(rises))


def write_xyz(mol, fname, comment=""):
    pos = mol.GetConformer().GetPositions()
    with open(fname, "w") as f:
        f.write(f"{mol.GetNumAtoms()}\n{comment}\n")
        for a, p in zip(mol.GetAtoms(), pos):
            f.write(f"{a.GetSymbol():2s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}\n")


if __name__ == "__main__":
    import sys
    tors = [(80.0, -120.0)] * 5
    m = assemble("AAAMMM", tors)
    print("stereo OK:", check_stereo(m))
    rm = residue_atom_map(m)
    print("torsions:", [(round(p), round(q)) for p, q in glycosidic_torsions(m, m.GetConformer(), rm)])
    for row in ring_axial_check(m, m.GetConformer(), rm):
        print({k: round(float(v), 2) for k, v in row.items()}, is_4C1(row))
    e = mmff_relax(m)
    print("MMFF E:", e, "stereo after relax:", check_stereo(m))
    print("torsions after relax:", [(round(p), round(q)) for p, q in glycosidic_torsions(m, m.GetConformer(), rm)])
    print("helix:", helix_params(m, rm))
    write_xyz(m, "/tmp/test_hex.xyz")
