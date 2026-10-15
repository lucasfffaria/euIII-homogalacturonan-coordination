"""SMILES table of every molecular model of the study (Table S6) and check against the archived start geometries.

For each homogalacturonan hexasaccharide pattern the script writes the stereo-explicit SMILES exactly as assembled by
galA.oligomer_smiles (residue template concatenated through alpha-(1->4) linkages, non-reducing end first), its RDKit
canonical form, the InChIKey, formula and charges.  It then re-derives the SMILES of the chain from the Cartesian
coordinates of every archived starting structure (mopac_optimizations/<job>/<job>_start.xyz; Eu and water removed,
bonds perceived with rdDetermineBonds, stereochemistry assigned from 3D) and compares it with the canonical SMILES.

Usage (from code/):  python smiles_table.py [--runs ../mopac_optimizations] [--out ../data]
Needs RDKit (tested with 2026.09.1), numpy, pandas.
"""
import os, re, sys, argparse
import numpy as np, pandas as pd
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds, rdMolDescriptors
from rdkit.Chem.inchi import MolToInchiKey
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import galA

PATTERNS = [("A6", "AAAAAA"), ("H3blk", "AAAHHH"), ("H3alt", "AHAHAH"), ("P3", "AAMMAA"), ("B3", "AAAMMM"),
            ("L3", "AMAMAM"), ("B2", "AAMMMM"), ("D2", "AMMAMM"), ("S1", "MMAMMM"), ("M6", "MMMMMM")]
# monomer controls (same residue template) and the rigid ligands of the Eu-O scans
CONTROLS = [("GalA", "A", "alpha-D-galacturonate monomer (controls CTRL_GalA_*)"),
            ("GalAMe", "M", "methyl alpha-D-galacturonate monomer (controls CTRL_GalAMe_*)")]
SMALL = [("acetate", "CC(=O)[O-]", "control CTRL_acetate_bid; bidentate scan ligand"), ("water", "O", "scan ligand; explicit solvent"),
         ("methanol", "CO", "scan ligand (neutral hydroxyl)"), ("methoxide", "C[O-]", "scan ligand (stress test)"),
         ("dimethyl ether", "COC", "scan ligand (ring / glycosidic O)"), ("methyl acetate", "CC(=O)OC", "scan ligand (ester carbonyl)")]


def read_xyz(path):
    L = open(path).read().splitlines(); n = int(L[0])
    el = [l.split()[0] for l in L[2:2 + n]]
    xyz = np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])
    return el, xyz


def strip_eu_and_water(el, xyz):
    """Remove Eu and every water molecule (O with two H within 1.15 A and no C within 1.7 A)."""
    el = np.array(el); keep = el != "Eu"
    iO, iH, iC = (np.where(el == s)[0] for s in ("O", "H", "C"))
    for o in iO:
        hs = iH[np.linalg.norm(xyz[iH] - xyz[o], axis=1) < 1.15]
        if len(hs) == 2 and (len(iC) == 0 or np.linalg.norm(xyz[iC] - xyz[o], axis=1).min() > 1.7):
            keep[o] = False; keep[hs] = False
    return list(el[keep]), xyz[keep]


def smiles_from_3d(el, xyz, charge):
    m = Chem.RWMol()
    for e in el:
        m.AddAtom(Chem.Atom(e))
    conf = Chem.Conformer(len(el))
    for i, p in enumerate(xyz):
        conf.SetAtomPosition(i, tuple(map(float, p)))
    m.AddConformer(conf); m = m.GetMol()
    rdDetermineBonds.DetermineBonds(m, charge=charge)
    Chem.SanitizeMol(m); Chem.AssignStereochemistryFrom3D(m)
    return Chem.MolToSmiles(Chem.RemoveHs(m))


def describe(name, smi, role):
    m = Chem.MolFromSmiles(smi)
    return dict(code=name, role=role, formula=rdMolDescriptors.CalcMolFormula(m), charge=Chem.GetFormalCharge(m),
                n_atoms=Chem.AddHs(m).GetNumAtoms(), n_stereocentres=len(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)),
                smiles_as_built=smi, smiles_canonical=Chem.MolToSmiles(m), inchikey=MolToInchiKey(m))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=os.path.join("..", "mopac_optimizations"))
    ap.add_argument("--out", default=os.path.join("..", "data"))
    a = ap.parse_args()
    rows, canon = [], {}
    for code, pat in PATTERNS:
        d = describe(code, galA.oligomer_smiles(pat, galA.TAGS), "HG hexasaccharide")
        d.update(sequence=pat, DM_percent=round(100 * pat.count("M") / 6), n_COO=pat.count("A"), n_COOCH3=pat.count("M"),
                 n_COOH=pat.count("H"), cluster_charge=d["charge"] + 3)
        rows.append(d); canon[code] = (d["smiles_canonical"], d["charge"])
    for code, kind, role in CONTROLS:
        d = describe(code, galA.oligomer_smiles(kind, galA.TAGS), role); d.update(sequence=kind, cluster_charge=d["charge"] + 3)
        rows.append(d); canon[code] = (d["smiles_canonical"], d["charge"])
    for code, smi, role in SMALL:
        rows.append(describe(code, smi, role))
    canon["acetate"] = (Chem.MolToSmiles(Chem.MolFromSmiles("CC(=O)[O-]")), -1)
    # check against every archived start geometry
    chk = []
    for job in sorted(os.listdir(a.runs)):
        p = os.path.join(a.runs, job, job + "_start.xyz")
        if not os.path.isfile(p):
            continue
        code = job.split("_")[1] if "_CTRL_" not in job else job.split("_")[2]
        if code not in canon:
            chk.append(dict(job=job, code=code, checked=False, match="", note="no organic ligand (aqua-ion control)")); continue
        el, xyz = strip_eu_and_water(*read_xyz(p))
        try:
            smi3 = smiles_from_3d(el, xyz, canon[code][1]); ok = smi3 == canon[code][0]; note = ""
        except Exception as e:
            ok = False; note = str(e)[:100]
        chk.append(dict(job=job, code=code, checked=True, match=ok, n_ligand_atoms=len(el), note=note))
    T = pd.DataFrame(rows); C = pd.DataFrame(chk)
    n = C[C.checked].groupby("code").agg(n_start_geometries=("job", "size"), n_match=("match", "sum")).reset_index()
    T = T.merge(n, on="code", how="left")
    cols = ["code", "role", "sequence", "DM_percent", "n_COO", "n_COOCH3", "n_COOH", "charge", "cluster_charge", "formula", "n_atoms",
            "n_stereocentres", "n_start_geometries", "n_match", "smiles_as_built", "smiles_canonical", "inchikey"]
    os.makedirs(a.out, exist_ok=True)
    T[cols].to_csv(os.path.join(a.out, "TableS6_smiles.csv"), index=False)
    C.to_csv(os.path.join(a.out, "TableS6_smiles_check_per_job.csv"), index=False)
    print("alpha-D-GalpA stereo tags:", galA.TAGS, " monomer CIP (oxane numbering):", galA.cip_of_monomer(galA.TAGS))
    print(T[["code", "formula", "charge", "n_atoms", "n_start_geometries", "n_match"]].to_string(index=False))
    print("start geometries checked:", int(C.checked.sum()), " matching:", int(C[C.checked].match.sum()), " not applicable:", int((~C.checked).sum()))


if __name__ == "__main__":
    main()
