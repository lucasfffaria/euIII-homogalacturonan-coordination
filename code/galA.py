"""
Homogalacturonan (HG) oligomer builder for the Eu(III)/pectin study.

Residue template written from O4 to C1 so that residues concatenate as
"O" + T + "O" + T + ... + "O"   (alpha-(1->4) linkages, reducing end = last).
Stereo tags are determined by enumeration on the monomer so that the
pyranose is alpha-D-galactopyranuronic acid, CIP (2S,3R,4S,5R,6S) in oxane numbering.
"""
import itertools
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdCIPLabeler, rdMolTransforms

# template with placeholders for the five stereocentres (C4, C5, C1, C2, C3 order of appearance)
# ring path: O4-C4(1)-C5(C6OO)-O5-C1( C2(O2)-C3(O3)-ring close )-
TEMPLATE = "[C{a}H]1[C{b}H]({carboxyl})O[C{c}H]([C{d}H](O)[C{e}H]1O)"

CARBOXYL = {"H": "C(=O)O", "A": "C(=O)[O-]", "M": "C(=O)OC"}  # acid, anion, methyl ester

def residue_smiles(tags, kind):
    a, b, c, d, e = tags
    return TEMPLATE.format(a=a, b=b, c=c, d=d, e=e, carboxyl=CARBOXYL[kind])

def oligomer_smiles(pattern, tags):
    """pattern: string like 'AAAMMM' (non-reducing end first)."""
    s = ""
    for k in pattern:
        s += "O" + residue_smiles(tags, k)
    return s + "O"

def cip_of_monomer(tags):
    mol = Chem.MolFromSmiles(oligomer_smiles("H", tags))
    rdCIPLabeler.AssignCIPLabels(mol)
    # identify oxane positions: ring O = pos1; C2 = ring C bearing carboxyl; C6 = anomeric C (ring C bonded to two O)
    ring = [r for r in mol.GetRingInfo().AtomRings() if len(r) == 6][0]
    ringO = [i for i in ring if mol.GetAtomWithIdx(i).GetSymbol() == "O"][0]
    labels = {}
    for i in ring:
        at = mol.GetAtomWithIdx(i)
        if at.GetSymbol() != "C":
            continue
        nO = sum(1 for n in at.GetNeighbors() if n.GetSymbol() == "O")
        hasCarboxyl = any(n.GetSymbol() == "C" and not n.IsInRing() for n in at.GetNeighbors())
        cip = at.GetProp("_CIPCode") if at.HasProp("_CIPCode") else "?"
        labels[i] = (nO, hasCarboxyl, cip)
    # walk the ring from ring O: C2 is the carboxyl-bearing carbon
    # build adjacency in ring
    ringset = set(ring)
    def ring_nb(i):
        return [n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() in ringset]
    c2 = [i for i in ring_nb(ringO) if labels[i][1]][0]
    c6 = [i for i in ring_nb(ringO) if i != c2][0]
    order = [ringO, c2]
    while len(order) < 6:
        nb = [n for n in ring_nb(order[-1]) if n not in order]
        order.append(nb[0])
    # order = O1, C2, C3, C4, C5, C6
    return "".join(labels[i][2] for i in order[1:])

def find_alpha_D_GalpA_tags():
    target = "SRSRS"  # (2S,3R,4S,5R,6S)
    for combo in itertools.product(["@", "@@"], repeat=5):
        if cip_of_monomer(combo) == target:
            return combo
    raise RuntimeError("no stereo combination matched")

TAGS = find_alpha_D_GalpA_tags()


def hg_mol(pattern, add_hs=True):
    smi = oligomer_smiles(pattern, TAGS)
    mol = Chem.MolFromSmiles(smi)
    if add_hs:
        mol = Chem.AddHs(mol)
    return mol


def residue_atom_map(mol):
    """Return list of dicts per residue with indices of key atoms:
    C1..C6, O2,O3,O4,O5, O6a (C=O), O6b (O-/OH/OMe), CH3 (if ester)."""
    ri = mol.GetRingInfo()
    rings = [r for r in ri.AtomRings() if len(r) == 6 and any(mol.GetAtomWithIdx(i).GetSymbol() == "O" for i in r)]
    res = []
    for ring in rings:
        ringset = set(ring)
        ringO = [i for i in ring if mol.GetAtomWithIdx(i).GetSymbol() == "O"][0]
        def ring_nb(i):
            return [n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() in ringset]
        def exo(i, sym):
            return [n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() not in ringset and n.GetSymbol() == sym]
        nbs = ring_nb(ringO)
        c5 = [i for i in nbs if exo(i, "C")][0]
        c1 = [i for i in nbs if i != c5][0]
        c2 = [i for i in ring_nb(c1) if i != ringO][0]
        c3 = [i for i in ring_nb(c2) if i != c1][0]
        c4 = [i for i in ring_nb(c3) if i != c2][0]
        c6 = exo(c5, "C")[0]
        o6 = [n.GetIdx() for n in mol.GetAtomWithIdx(c6).GetNeighbors() if n.GetSymbol() == "O"]
        o6a = [o for o in o6 if mol.GetBondBetweenAtoms(c6, o).GetBondTypeAsDouble() == 2.0][0]
        o6b = [o for o in o6 if o != o6a][0]
        ch3 = [n.GetIdx() for n in mol.GetAtomWithIdx(o6b).GetNeighbors() if n.GetSymbol() == "C" and n.GetIdx() != c6]
        d = dict(C1=c1, C2=c2, C3=c3, C4=c4, C5=c5, C6=c6, O5=ringO,
                 O1=exo(c1, "O")[0], O2=exo(c2, "O")[0], O3=exo(c3, "O")[0], O4=exo(c4, "O")[0],
                 O6a=o6a, O6b=o6b, CH3=(ch3[0] if ch3 else None))
        d["kind"] = "M" if ch3 else ("A" if mol.GetAtomWithIdx(o6b).GetFormalCharge() == -1 else "H")
        res.append(d)
    # order residues from non-reducing end: residue whose O4 is bonded to H (free) is first;
    # chain: C1(res i) - O4(res i+1)
    o4_to_res = {d["O4"]: k for k, d in enumerate(res)}
    c1_to_res = {d["C1"]: k for k, d in enumerate(res)}
    # next residue of res k: the residue whose O4 == O1 of res k (glycosidic O)
    nxt = {}
    for k, d in enumerate(res):
        o1 = d["O1"]
        if o1 in o4_to_res:
            nxt[k] = o4_to_res[o1]
    starts = [k for k in range(len(res)) if k not in nxt.values()]
    order = [starts[0]]
    while order[-1] in nxt:
        order.append(nxt[order[-1]])
    return [res[k] for k in order]


def glycosidic_torsions(mol, conf, resmap):
    """phi = O5(i)-C1(i)-O4(i+1)-C4(i+1); psi = C1(i)-O4(i+1)-C4(i+1)-C5(i+1) (per IUPAC-like defs)."""
    out = []
    for i in range(len(resmap) - 1):
        a, b = resmap[i], resmap[i + 1]
        phi = rdMolTransforms.GetDihedralDeg(conf, a["O5"], a["C1"], b["O4"], b["C4"])
        psi = rdMolTransforms.GetDihedralDeg(conf, a["C1"], b["O4"], b["C4"], b["C5"])
        out.append((phi, psi))
    return out


def set_glycosidic_torsions(mol, conf, resmap, phi, psi):
    for i in range(len(resmap) - 1):
        a, b = resmap[i], resmap[i + 1]
        rdMolTransforms.SetDihedralDeg(conf, a["O5"], a["C1"], b["O4"], b["C4"], phi)
        rdMolTransforms.SetDihedralDeg(conf, a["C1"], b["O4"], b["C4"], b["C5"], psi)


def ring_axial_check(mol, conf, resmap):
    """Return per residue the |cos| of angle between substituent bond and ring normal
    for O1, O2, O3, O4, C6: axial ~1, equatorial ~0.3."""
    pos = conf.GetPositions()
    out = []
    for d in resmap:
        ring = [d["C1"], d["C2"], d["C3"], d["C4"], d["C5"], d["O5"]]
        P = pos[ring]
        c = P.mean(axis=0)
        u, s, vt = np.linalg.svd(P - c)
        n = vt[2]
        row = {}
        for sub, base in (("O1", "C1"), ("O2", "C2"), ("O3", "C3"), ("O4", "C4"), ("C6", "C5")):
            v = pos[d[sub]] - pos[d[base]]
            row[sub] = abs(np.dot(v, n) / np.linalg.norm(v))
        out.append(row)
    return out


if __name__ == "__main__":
    print("alpha-D-GalpA tags:", TAGS)
    print("monomer CIP:", cip_of_monomer(TAGS))
    m = hg_mol("AAAMMM")
    print(Chem.MolToSmiles(Chem.RemoveHs(m)))
    AllChem.EmbedMolecule(m, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    rm = residue_atom_map(m)
    print("residues:", [r["kind"] for r in rm])
    for row in ring_axial_check(m, m.GetConformer(), rm):
        print({k: round(v, 2) for k, v in row.items()})
    print("torsions:", [(round(p), round(q)) for p, q in glycosidic_torsions(m, m.GetConformer(), rm)])
