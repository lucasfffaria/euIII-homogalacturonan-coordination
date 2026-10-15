# -*- coding: utf-8 -*-
"""
Crystal-structure benchmark of the Eu-O distances (revision v2, Section S5 of the Supplementary Information).

  python xray_benchmark.py build   --cif <folder with the COD .cif files> --out <campaign folder> [--eufix eu_fix.txt]
  python xray_benchmark.py analyze --runs <campaign folder with the MOPAC outputs> --out <folder for the tables>

build    extracts the discrete Eu(III) complex from each crystal structure (unit cell rebuilt with gemmi, molecules found by
         distance criteria across the periodic images), assigns the charge of every ligand (RDKit bond-order perception with the
         smallest total charge that gives a closed-shell Lewis structure), sets the O-H lengths to 0.96 A (water) or 0.97 A
         (hydroxyl, acid, phenol) along the crystal O-H direction and writes one MOPAC input per model, with the O-H lengths
         fixed (H atoms of O-H groups in internal coordinates), as in the article. It also writes the reference tables.
analyze  reads the final geometries (.out), compares every Eu-O distance within 3.0 A in the crystal with the optimized
         value and writes the tables of Section S5.

The 17 structures are listed in SELECTION (COD identifiers). They were chosen, before any calculation, among the COD entries
that contain a discrete complex with at most two Eu atoms, only oxygen donor atoms, no disorder in the complex and a
coordination number of at least six (screening table: data/xray_benchmark/R10_COD_screening.csv), so as to cover the donor
classes of the article. The H atoms of the water molecules of COD 2102807 are not in the CIF and are added here.
Requires numpy, scipy, pandas, gemmi and rdkit.
"""
import argparse, os, re, sys
from itertools import product
import numpy as np, pandas as pd
from scipy.spatial.distance import cdist

SELECTION = [('XB01_aqua8a', '7209926'), ('XB02_aqua8b', '7201013'), ('XB03_aqua9', '2102807'), ('XB04_TFA2_w6', '7028020'),
             ('XB05_ODA3', '7715872'), ('XB06_MAA_dimer', '2217178'), ('XB07_HBA_dimer', '7062635'), ('XB08_DCS_dimer', '1543204'),
             ('XB09_oxalate_dimer', '7044674'), ('XB10_TPA3_MeOH4', '2020541'), ('XB11_ArO3_DME2', '4307878'),
             ('XB12_dik3_diglyme', '1547880'), ('XB13_TMDGA3', '7040442'), ('XB14_dik3_AcOH3', '7114347'), ('XB15_dik3_w2', '4329032'),
             ('XB16_FBA_dimer', '7061033'), ('XB17_ArO3_THF3', '7018331')]
# donor classes that the automatic classification cannot tell apart (assigned from the identity of the ligand)
RELABEL = {'XB11_ArO3_DME2': {'carbonyl/enolate': 'aryloxide'}, 'XB13_TMDGA3': {'carbonyl/enolate': 'amide C=O'},
           'XB12_dik3_diglyme': {'carbonyl/enolate': 'diketonate'},
           'XB14_dik3_AcOH3': {'carbonyl/enolate': 'diketonate', 'ester/acid C=O': 'carboxylic acid C=O'},
           'XB15_dik3_w2': {'carbonyl/enolate': 'diketonate'}, 'XB10_TPA3_MeOH4': {'hydroxyl': 'alcohol OH'},
           'XB16_FBA_dimer': {'ester/acid C=O': 'carboxylic acid C=O'}, 'XB17_ArO3_THF3': {'carbonyl/enolate': 'aryloxide'}}
KEYWORDS = {'PM7': 'PM7 SPARKLE CHARGE={q} AUX GNORM=1.0 CYCLES=3000 THREADS=2 GEO-OK T=3D DUMP=30M BONDS EF LET EXTERNAL=eu_fix.txt',
            'PM7std': 'PM7 SPARKLE CHARGE={q} AUX GNORM=1.0 CYCLES=3000 THREADS=2 GEO-OK T=3D DUMP=30M BONDS EF LET',
            'SRM1': 'RM1 SPARKLE CHARGE={q} AUX GNORM=1.0 CYCLES=3000 THREADS=2 GEO-OK T=3D DUMP=30M BONDS EF LET',
            'RM1': 'RM1 CHARGE={q} AUX GNORM=1.0 CYCLES=3000 THREADS=2 GEO-OK T=3D DUMP=30M BONDS EF LET'}
MODEL_LABEL = {'PM7': 'Sparkle/PM7+EuO', 'PM7std': 'Sparkle/PM7 as distributed', 'SRM1': 'Sparkle/RM1', 'RM1': 'RM1 lanthanide model'}
EU_FIX = "XFAC_Eu O 2.207500\nALPB_Eu O 2.500000\n"
RCOV = {'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'S': 1.05, 'Cl': 1.02, 'Br': 1.20, 'I': 1.39, 'P': 1.07, 'B': 0.84,
        'Si': 1.11, 'Se': 1.20}
IONS = {'Na', 'K', 'Li', 'Cs', 'Rb', 'Ca', 'Sr', 'Ba', 'Mg'}


# ----------------------------------------------------------------------------------------------- crystal -> complex
def load_cell(path):
    import gemmi
    st = gemmi.read_small_structure(path)
    sites = st.get_all_unit_cell_sites()
    el = [s.element.name for s in sites]
    fr = np.array([[s.fract.x, s.fract.y, s.fract.z] for s in sites]) % 1.0
    occ = np.array([s.occ for s in sites]); lab = [s.label for s in sites]
    M = np.array(st.cell.orth.mat.tolist())          # fractional -> Cartesian
    return el, fr, occ, lab, M


def bonded(e1, e2, d):
    if e1 == 'Eu' or e2 == 'Eu':
        o = e2 if e1 == 'Eu' else e1
        return o in ('O', 'N', 'Cl', 'S', 'F', 'Br', 'I') and d < 3.10
    if e1 in IONS or e2 in IONS: return False
    if e1 == 'H' and e2 == 'H': return False
    if e1 not in RCOV or e2 not in RCOV: return False
    return 0.5 < d < 1.25 * (RCOV[e1] + RCOV[e2])


def neighbours(el, fr, M):
    n = len(el); xyz = fr @ M.T; shifts = np.array(list(product((-1, 0, 1), repeat=3))); nb = [[] for _ in range(n)]
    for k, s in enumerate(shifts @ M.T):
        D = cdist(xyz, xyz + s)
        for i, j in zip(*np.where(D < 3.2)):
            if i == j and not shifts[k].any(): continue
            if bonded(el[i], el[j], D[i, j]): nb[i].append((j, tuple(shifts[k]), D[i, j]))
    return nb, xyz


def extract(path, maxn=260):
    """Molecule that contains the first Eu site; None if it is polymeric or larger than maxn atoms."""
    el, fr, occ, lab, M = load_cell(path)
    info = dict(n_cell=len(el), disorder_cell=bool((occ < 0.999).any()))
    if 'Eu' not in el: info['reason'] = 'no Eu site'; return info, None
    if len(el) > 3000: info['reason'] = 'cell too large'; return info, None
    nb, xyz = neighbours(el, fr, M)
    start = el.index('Eu'); seen = {(start, (0, 0, 0))}; first = {start: (0, 0, 0)}; q = [(start, (0, 0, 0))]; poly = False
    while q:
        i, im = q.pop()
        for j, s, d in nb[i]:
            jm = (im[0] + s[0], im[1] + s[1], im[2] + s[2])
            if (j, jm) in seen: continue
            if j in first and first[j] != jm: poly = True; break
            first.setdefault(j, jm); seen.add((j, jm)); q.append((j, jm))
            if len(seen) > maxn: break
        if poly or len(seen) > maxn: break
    if poly: info['reason'] = 'polymeric'; return info, None
    if len(seen) > maxn: info['reason'] = 'too large'; return info, None
    nodes = sorted(seen); idx = {nd: k for k, nd in enumerate(nodes)}
    cel = [el[i] for i, _ in nodes]; cxyz = np.array([xyz[i] + np.array(im) @ M.T for i, im in nodes])
    bonds = set()
    for (i, im) in nodes:
        for j, s, d in nb[i]:
            jm = (im[0] + s[0], im[1] + s[1], im[2] + s[2])
            if (j, jm) in idx: bonds.add(tuple(sorted((idx[(i, im)], idx[(j, jm)]))))
    cocc = np.array([occ[i] for i, _ in nodes])
    info.update(n_atoms=len(cel), n_Eu=cel.count('Eu'), disorder_complex=bool((cocc < 0.999).any()))
    return info, dict(el=cel, xyz=cxyz, bonds=sorted(bonds))


def classify(c):
    """(Eu, donor, class, distance) for every atom bonded to Eu."""
    el = c['el']; n = len(el); adj = [[] for _ in range(n)]; xyz = c['xyz']
    for i, j in c['bonds']: adj[i].append(j); adj[j].append(i)
    donors = []
    for e in [i for i in range(n) if el[i] == 'Eu']:
        for j in adj[e]:
            d = float(np.linalg.norm(xyz[j] - xyz[e])); t = el[j]
            if el[j] == 'O':
                nbr = [k for k in adj[j] if el[k] != 'Eu']; nh = sum(el[k] == 'H' for k in nbr); hv = [k for k in nbr if el[k] != 'H']
                if not hv: t = {2: 'water', 1: 'OH?/water-1H', 0: 'O?/water-0H'}[min(nh, 2)]
                elif len(hv) == 1 and el[hv[0]] == 'C':
                    ci = hv[0]; co = [k for k in adj[ci] if el[k] == 'O']; cn = len(adj[ci])
                    if nh == 1: t = 'hydroxyl' if cn == 4 else ('acid OH' if len(co) == 2 else 'enol/phenol OH')
                    elif cn == 4: t = 'alkoxide'
                    elif len(co) >= 2:
                        other = [k for k in co if k != j][0]; onb = [k for k in adj[other] if el[k] != 'Eu' and k != ci]
                        t = 'carboxylate' if not onb else ('ester/acid C=O' if len(co) == 2 else 'carbonate')
                        if len(co) == 3: t = 'carbonate'
                    else: t = 'carbonyl/enolate'
                elif len(hv) == 2 and all(el[k] == 'C' for k in hv): t = 'ether' if nh == 0 else 'protonated ether?'
                elif len(hv) == 1: t = f'O-{el[hv[0]]}'
                else: t = 'O-other'
            donors.append((e, j, t, d))
    return donors


def frag_charge(el, xyz, bonds):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdDetermineBonds
    RDLogger.DisableLog('rdApp.*')
    rw = Chem.RWMol()
    for e in el: rw.AddAtom(Chem.Atom(e))
    for i, j in bonds: rw.AddBond(int(i), int(j), Chem.BondType.SINGLE)
    conf = Chem.Conformer(len(el))
    for i, p in enumerate(xyz): conf.SetAtomPosition(i, [float(v) for v in p])
    rw.AddConformer(conf)
    for q in (0, -1, -2, -3, -4, 1, 2):
        m = Chem.RWMol(rw)
        try:
            rdDetermineBonds.DetermineBondOrders(m, charge=q); Chem.SanitizeMol(m)
            return q, Chem.MolToSmiles(Chem.RemoveHs(m))
        except Exception:
            continue
    return None, None


def ligands(c):
    el = c['el']; n = len(el); keep = [i for i in range(n) if el[i] != 'Eu']; adj = {i: [] for i in keep}
    for i, j in c['bonds']:
        if el[i] != 'Eu' and el[j] != 'Eu': adj[i].append(j); adj[j].append(i)
    seen = set(); out = []
    for s in keep:
        if s in seen: continue
        comp = []; st = [s]; seen.add(s)
        while st:
            a = st.pop(); comp.append(a)
            for b in adj[a]:
                if b not in seen: seen.add(b); st.append(b)
        comp = sorted(comp); mp = {a: k for k, a in enumerate(comp)}
        bb = [(mp[i], mp[j]) for i, j in c['bonds'] if i in mp and j in mp]
        q, smi = frag_charge([el[a] for a in comp], c['xyz'][comp], bb)
        out.append(dict(atoms=comp, q=q, smiles=smi))
    return out


def add_water_H(c):
    """Two H atoms on every O of an aqua ion whose H atoms are not in the CIF (O-H 0.96 A, H-O-H 104.5 deg, bisector along Eu->O)."""
    el = list(c['el']); xyz = [np.array(p) for p in c['xyz']]; bonds = list(c['bonds']); eu = el.index('Eu')
    Os = [i for i, e in enumerate(el) if e == 'O']; placed = []; half = np.radians(104.5 / 2)
    for o in Os:
        u = xyz[o] - xyz[eu]; u /= np.linalg.norm(u); a = np.cross(u, [1, 0, 0]); a = a if np.linalg.norm(a) > 0.3 else np.cross(u, [0, 1, 0])
        a /= np.linalg.norm(a); b = np.cross(u, a); best = None
        for ang in np.radians(np.arange(0, 180, 10)):
            w = np.cos(ang) * a + np.sin(ang) * b
            h1 = xyz[o] + 0.96 * (np.cos(half) * u + np.sin(half) * w); h2 = xyz[o] + 0.96 * (np.cos(half) * u - np.sin(half) * w)
            others = [xyz[k] for k in Os if k != o] + placed
            score = min(min(np.linalg.norm(h1 - p), np.linalg.norm(h2 - p)) for p in others)
            if best is None or score > best[0]: best = (score, h1, h2)
        for h in best[1:]:
            el.append('H'); xyz.append(h); placed.append(h); bonds.append((o, len(el) - 1))
    return dict(el=el, xyz=np.array(xyz), bonds=sorted(bonds))


def build_structure(name, cif_path):
    info, c = extract(cif_path)
    assert c is not None, (name, info)
    if name == 'XB03_aqua9': c = add_water_H(c)
    el = c['el']; n = len(el); xyz = np.array(c['xyz'], float).copy(); adj = [[] for _ in range(n)]
    for i, j in c['bonds']: adj[i].append(j); adj[j].append(i)
    hfix = {}
    for h in range(n):
        if el[h] != 'H': continue
        assert len(adj[h]) == 1, (name, h); p = adj[h][0]
        if el[p] == 'O':
            hv = [k for k in adj[p] if el[k] not in ('H', 'Eu')]; r = 0.96 if not hv else 0.97
            v = xyz[h] - xyz[p]; xyz[h] = xyz[p] + v / np.linalg.norm(v) * r; hfix[h] = (p, r)
    L = ligands(c); assert all(l['q'] is not None for l in L), name
    don = [(e, j, RELABEL.get(name, {}).get(t, t).replace('O-S', 'sulfonate'), d) for e, j, t, d in classify(c)]
    order = [i for i in range(n) if el[i] != 'H'] + [i for i in range(n) if el[i] == 'H']
    return dict(el=el, xyz=xyz, xyz_cif=np.array(c['xyz'], float), adj=adj, hfix=hfix, donors=don, ligands=L, order=order,
                charge=int(3 * el.count('Eu') + sum(l['q'] for l in L)))


# ----------------------------------------------------------------------------------------------- MOPAC geometry
def angle(a, b, c):
    u = a - b; v = c - b
    return np.degrees(np.arccos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1)))


def int2cart(atoms):
    """MOPAC mixed geometry -> Cartesian. atoms: (el, v1, f1, v2, f2, v3, f3, na, nb, nc); na = 0 means Cartesian.
    Checked against the Cartesian coordinates printed by MOPAC 23.2.5 for the inputs of the article (5e-5 A)."""
    X = []
    for el, v1, f1, v2, f2, v3, f3, na, nb, nc in atoms:
        if na == 0: X.append(np.array([v1, v2, v3])); continue
        A = X[na - 1]; B = X[nb - 1]; C = X[nc - 1]; r = v1; th = np.radians(v2); ph = np.radians(v3)
        bc = A - B; bc /= np.linalg.norm(bc); nrm = np.cross(B - C, bc); nrm /= np.linalg.norm(nrm); m = np.cross(nrm, bc)
        X.append(A + (-r * np.cos(th)) * bc + (r * np.sin(th) * np.cos(ph)) * m + (r * np.sin(th) * np.sin(ph)) * nrm)
    return np.array(X)


def read_mop(path):
    L = open(path).read().splitlines(); atoms = []
    for ln in L[3:]:
        p = ln.split()
        if len(p) >= 7: atoms.append((p[0], float(p[1]), int(p[2]), float(p[3]), int(p[4]), float(p[5]), int(p[6]), *(int(x) for x in p[7:10])))
    return L[:3], atoms


def geo_lines(b):
    """Geometry block: heavy atoms first, then H; H of O-H groups in internal coordinates with the bond length fixed."""
    el = b['el']; xyz = b['xyz']; adj = b['adj']; order = b['order']; new = {o: k + 1 for k, o in enumerate(order)}
    heavy = [i for i in range(len(el)) if el[i] != 'H']; lines = []; nint = 0
    for i in order:
        if i not in b['hfix']:
            x, y, z = xyz[i]; lines.append(f"  {el[i]:<2s} {x:14.8f} +1 {y:12.7f} +1 {z:12.7f} +1     0     0     0"); continue
        p, r = b['hfix'][i]; cand = sorted([k for k in adj[p] if el[k] != 'H'], key=lambda k: el[k] == 'Eu'); done = False
        for nb_ in cand + [k for k in heavy if k != p and k not in cand]:
            a1 = angle(xyz[i], xyz[p], xyz[nb_])
            if not (15 < a1 < 165): continue
            ncs = [k for k in adj[nb_] if el[k] != 'H' and k != p] + sorted([k for k in heavy if k not in (p, nb_)],
                                                                         key=lambda k: np.linalg.norm(xyz[k] - xyz[nb_]))
            for nc_ in ncs:
                if not (25 < angle(xyz[p], xyz[nb_], xyz[nc_]) < 155): continue
                b1 = xyz[p] - xyz[i]; b2 = xyz[nb_] - xyz[p]; b3 = xyz[nc_] - xyz[nb_]; n1 = np.cross(b1, b2); n2 = np.cross(b2, b3)
                dih = np.degrees(np.arctan2(np.cross(n1, n2) @ (b2 / np.linalg.norm(b2)), n1 @ n2))
                for s in (dih, -dih):                      # keep the sign that rebuilds the position
                    X = int2cart([('X', *np.ravel([[v, 1] for v in xyz[p]]), 0, 0, 0), ('X', *np.ravel([[v, 1] for v in xyz[nb_]]), 0, 0, 0),
                                  ('X', *np.ravel([[v, 1] for v in xyz[nc_]]), 0, 0, 0), ('H', r, 0, a1, 1, s, 1, 1, 2, 3)])
                    if np.linalg.norm(X[3] - xyz[i]) < 1e-6:
                        lines.append(f"  H  {r:14.8f} +0 {a1:12.7f} +1 {s:13.7f} +1 {new[p]:5d} {new[nb_]:5d} {new[nc_]:5d}")
                        done = True; nint += 1; break
                if done: break
            if done: break
        assert done, (i, p)
    return lines, nint


def cmd_build(a):
    camp = a.out; os.makedirs(os.path.join(camp, '_reference'), exist_ok=True); mrows = []; refrows = []; maxdev = 0.0
    eufix = open(a.eufix).read() if a.eufix else EU_FIX
    for name, cid in SELECTION:
        b = build_structure(name, os.path.join(a.cif, cid + '.cif')); order = b['order']; new = {o: k for k, o in enumerate(order)}
        lines, nint = geo_lines(b); elo = [b['el'][i] for i in order]; xo = b['xyz'][order]
        for m, kw in KEYWORDS.items():
            job = f'{m}_{name}'; d = os.path.join(camp, job); os.makedirs(d, exist_ok=True)
            title = f"{job}: {MODEL_LABEL[m]}; complex taken from the crystal structure COD {cid} (isolated complex, O-H lengths set to 0.96/0.97 A and kept fixed)"
            open(os.path.join(d, job + '.mop'), 'w', newline='\n').write('\n'.join([kw.format(q=b['charge']), title, ''] + lines) + '\n')
            if m == 'PM7': open(os.path.join(d, 'eu_fix.txt'), 'w', newline='\n').write(eufix)
            _, at = read_mop(os.path.join(d, job + '.mop')); maxdev = max(maxdev, float(np.abs(int2cart(at) - xo).max()))
            mrows.append(dict(job=job, model=m, structure=name, cod_id=cid, charge=b['charge'], natoms=len(elo), n_Eu=elo.count('Eu'),
                              n_OH_fixed=nint, keywords=kw.format(q=b['charge'])))
        with open(os.path.join(camp, '_reference', name + '_xray.xyz'), 'w') as fh:
            note = ' ; H atoms of the water molecules added (not in the CIF)' if name == 'XB03_aqua9' else ''
            fh.write(f"{len(elo)}\n{name} COD {cid} isolated complex, crystal coordinates (same atom order as the .mop files){note}\n")
            for e, p in zip(elo, b['xyz_cif'][order]): fh.write(f"{e:2s} {p[0]:12.5f} {p[1]:12.5f} {p[2]:12.5f}\n")
        for e, j, t, dd in b['donors']:
            if dd <= 3.0: refrows.append(dict(structure=name, cod_id=cid, Eu_index=new[e] + 1, O_index=new[j] + 1, donor_type=t, d_xray=round(dd, 4)))
    pd.DataFrame(mrows).to_csv(os.path.join(camp, 'manifest_R10.csv'), index=False)
    pd.DataFrame(refrows).to_csv(os.path.join(camp, '_reference', 'xray_EuO_distances.csv'), index=False)
    print(f"{len(mrows)} inputs written to {camp}; {len(refrows)} reference Eu-O distances; largest deviation of the written geometry {maxdev:.1e} A")


# ----------------------------------------------------------------------------------------------- analysis
def read_out(path):
    """Elements, final Cartesian geometry (last CARTESIAN COORDINATES block) and final heat of formation of a MOPAC .out file."""
    t = open(path, errors='replace').read(); k = t.rfind('CARTESIAN COORDINATES'); el = []; xyz = []
    for l in t[k:].splitlines()[1:]:
        p = l.split()
        if len(p) == 5 and p[0].isdigit(): el.append(p[1]); xyz.append([float(v) for v in p[2:5]])
        elif el: break
    hm = re.findall(r'FINAL HEAT OF FORMATION\s*=\s*([-\d.]+)\s*KCAL/MOL', t)
    return el, np.array(xyz), (float(hm[-1]) if hm else None)


def read_xyz(path):
    L = open(path).read().splitlines(); n = int(L[0])
    return [l.split()[0] for l in L[2:2 + n]], np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])


def kabsch_rmsd(P, Q):
    P = P - P.mean(0); Q = Q - Q.mean(0); U, S, Vt = np.linalg.svd(P.T @ Q); d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1, 1, d]) @ Vt
    return float(np.sqrt(((P @ R - Q) ** 2).sum() / len(P)))


def bond_matrix(el, x):
    n = len(el); D = cdist(x, x); M = np.zeros((n, n), bool)
    for i in range(n):
        for j in range(i + 1, n):
            if 'Eu' in (el[i], el[j]) or 'H' in (el[i], el[j]): continue
            if el[i] in RCOV and el[j] in RCOV and D[i, j] < 1.25 * (RCOV[el[i]] + RCOV[el[j]]): M[i, j] = M[j, i] = True
    return M


def cmd_analyze(a):
    camp = a.runs; man = pd.read_csv(os.path.join(camp, 'manifest_R10.csv')); REF = pd.read_csv(os.path.join(camp, '_reference', 'xray_EuO_distances.csv'))
    os.makedirs(a.out, exist_ok=True); rows = []; pr = []; geo = []
    for m in man.itertuples():
        d = os.path.join(camp, m.job); stems = [m.job] + [f'{m.job}_r{k}' for k in (1, 2) if os.path.exists(os.path.join(d, f'{m.job}_r{k}.out'))]
        stem = stems[-1]; out = open(os.path.join(d, stem + '.out'), errors='replace').read()
        g = re.findall(r'GRADIENT NORM\s+=\s+([\d.]+)', out); gn = float(g[-1]) if g else np.nan
        ok = ('JOB ENDED NORMALLY' in out) or ('ended normally' in out) or ('MOPAC DONE' in out)
        el, x, heat = read_out(os.path.join(d, stem + '.out')); elc, xc = read_xyz(os.path.join(camp, '_reference', m.structure + '_xray.xyz'))
        assert list(el) == list(elc), m.job
        geo.append((m.job, el, x, heat))
        ref = REF[REF.structure == m.structure]; eus = [i for i, e in enumerate(el) if e == 'Eu']; Os = [i for i, e in enumerate(el) if e == 'O']
        dc = np.array([np.linalg.norm(x[r.O_index - 1] - x[r.Eu_index - 1]) for r in ref.itertuples()]); dx = ref.d_xray.values
        refpairs = {(r.Eu_index - 1, r.O_index - 1) for r in ref.itertuples()}
        gained = [(e, o) for e in eus for o in Os if (e, o) not in refpairs and np.linalg.norm(x[o] - x[e]) < 3.0]; lost = int((dc > 3.0).sum())
        Ms = bond_matrix(el, xc); Mf = bond_matrix(el, x); hv = np.array(el) != 'H'
        rows.append(dict(job=m.job, model=m.model, structure=m.structure, cod_id=m.cod_id, converged=bool(ok and gn <= 1.0), gnorm=gn, restarts=len(stems) - 1,
                         heat=heat, n_ref=len(dx), CN_xray=len(dx) / len(eus), CN_calc=(len(dx) - lost + len(gained)) / len(eus), n_lost=lost, n_gained=len(gained),
                         n_broken=int((Ms & ~Mf).sum() // 2), n_formed=int((Mf & ~Ms).sum() // 2), min_EuO=min(np.linalg.norm(x[o] - x[e]) for e in eus for o in Os),
                         mean_xray=dx.mean(), mean_calc=dc.mean(), MSE=(dc - dx).mean(), MUE=np.abs(dc - dx).mean(), max_abs_err=np.abs(dc - dx).max(),
                         rmsd_heavy=kabsch_rmsd(xc[hv], x[hv]),
                         EuEu_xray=(np.linalg.norm(xc[eus[0]] - xc[eus[1]]) if len(eus) == 2 else np.nan),
                         EuEu_calc=(np.linalg.norm(x[eus[0]] - x[eus[1]]) if len(eus) == 2 else np.nan)))
        pr += [dict(job=m.job, model=m.model, structure=m.structure, Eu_index=r.Eu_index, O_index=r.O_index, donor_type=r.donor_type, d_xray=r.d_xray,
                    d_calc=round(float(dd), 4), err=round(float(dd - r.d_xray), 4)) for r, dd in zip(ref.itertuples(), dc)]
    XB = pd.DataFrame(rows); XP = pd.DataFrame(pr)
    XB.to_csv(os.path.join(a.out, 'xray_benchmark_results.csv'), index=False); XP.to_csv(os.path.join(a.out, 'xray_benchmark_distances.csv'), index=False)
    order = ['PM7std', 'PM7', 'SRM1', 'RM1']

    def stats(g):
        e = g.err
        return pd.Series(dict(n=len(e), MSE=e.mean(), MUE=e.abs().mean(), RMSE=np.sqrt((e ** 2).mean()), max_abs=e.abs().max(), n_gt_0p10=int((e.abs() > 0.10).sum())))
    allm = XP.groupby('model').apply(stats).loc[order]
    extra = XB.groupby('model').agg(n_structures=('job', 'size'), converged=('converged', 'sum'), CN_kept=('CN_calc', lambda s: int((s == XB.loc[s.index, 'CN_xray']).sum())),
                                   shortest_EuO=('min_EuO', 'min'), mean_structure_MUE=('MUE', 'mean'), max_structure_MUE=('MUE', 'max')).loc[order]
    pd.concat([allm, extra], axis=1).round(5).to_csv(os.path.join(a.out, 'TableS13_benchmark_by_model.csv'))
    bt = XP.groupby(['donor_type', 'model']).err.agg(n='size', MSE='mean', MUE=lambda s: s.abs().mean()).round(5).reset_index()
    bt.to_csv(os.path.join(a.out, 'TableS14_benchmark_by_donor.csv'), index=False)
    XB.pivot(index='structure', columns='model', values='MUE')[order].round(5).to_csv(os.path.join(a.out, 'TableS15_benchmark_by_structure.csv'))
    with open(os.path.join(a.out, 'xray_benchmark_final_geometries.xyz'), 'w') as fh:
        for job, el, x, heat in geo:
            fh.write(f"{len(el)}\n{job} heat_of_formation_kcal_mol={heat}\n")
            for e, p in zip(el, x): fh.write(f"{e:2s} {p[0]:12.4f} {p[1]:12.4f} {p[2]:12.4f}\n")
    print(pd.concat([allm, extra], axis=1).round(3).to_string())


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter); sub = ap.add_subparsers(dest='cmd', required=True)
    b = sub.add_parser('build'); b.add_argument('--cif', required=True); b.add_argument('--out', required=True); b.add_argument('--eufix', default=None)
    c = sub.add_parser('analyze'); c.add_argument('--runs', required=True); c.add_argument('--out', required=True)
    a = ap.parse_args(); {'build': cmd_build, 'analyze': cmd_analyze}[a.cmd](a)
