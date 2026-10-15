"""Generator of the MOPAC inputs of the production campaign (mopac_optimizations/): HG hexasaccharides x chain
conformers x starting placements, two Hamiltonians, plus the hydration-sensitivity series and the controls.

Provenance: this is make_campaign_v4.py, which wrote the inputs of an earlier campaign (v4) that was stopped; the
production inputs are those same geometries (identical .mop coordinate blocks and _start.xyz files) run with the
keywords below (CYCLES=3000 and EF LET instead of CYCLES=8000). The keyword lines were updated here to match the
production inputs; nothing else was changed. The starting coordinates depend on the RDKit version (conformer
embedding and MMFF relaxation of the charged residues), so a re-run reproduces the construction but not
necessarily the archived coordinates; the archived inputs in mopac_optimizations/ are the reference.
Usage: python make_campaign.py <output_folder>

Models:
  PM7  : Sparkle/PM7 with the short-range Eu-O core-core correction (EXTERNAL=eu_fix.txt)
  RM1  : RM1 with explicit Eu(III) orbitals (Filho et al. 2014) -- only for net charge >= 0
Every O-H bond is frozen (mixed internal/Cartesian input) to keep the protonation state fixed.
Explicit water: 9 around Eu (TTP vertices at 2.48 A + second shell at 4.6 A) + 2 H-bonded waters per free carboxylate
not bound to Eu (main series); sensitivity series: 0 (cw0) or 4 (cw4) H-bonded waters per free carboxylate.
"""
import os, json, sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from galA import residue_atom_map, glycosidic_torsions
from assemble import assemble, mmff_relax, helix_params, check_stereo, write_xyz
from place_eu import build_complex, water_atoms, fib_sphere, find_site, fill_waters
from mopio import write_mop_mixed

OUT = sys.argv[1] if len(sys.argv) > 1 else "mopac_inputs_regenerated"
os.makedirs(OUT, exist_ok=True)

PATTERNS = {
    "A6": "AAAAAA", "B3": "AAAMMM", "L3": "AMAMAM", "P3": "AAMMAA",
    "B2": "AAMMMM", "D2": "AMMAMM", "S1": "MMAMMM", "M6": "MMMMMM",
    "H3alt": "AHAHAH", "H3blk": "AAAHHH",
}
CONFORMERS = {
    "h3": [(75.0, 105.0)] * 5,
    "h2": [(105.0, 135.0)] * 5,
    "cl": [(105.0, 165.0), (75.0, 105.0), (90.0, 165.0), (75.0, 120.0), (105.0, 135.0)],
}
KEYS = {
    "PM7": "PM7 CHARGE={q} AUX GNORM=1.0 CYCLES=3000 THREADS=2 GEO-OK T=3D DUMP=30M BONDS EF LET EXTERNAL=eu_fix.txt",
    "RM1": "RM1 CHARGE={q} AUX GNORM=1.0 CYCLES=3000 THREADS=2 GEO-OK T=3D DUMP=30M BONDS EF LET",
}
EU_FIX = "XFAC_Eu O 2.207500\nALPB_Eu O 2.500000\n"
manifest = []


def relax_chain(pattern, tors):
    m = assemble(pattern, tors)
    rm = residue_atom_map(m)
    fixed = []
    for i, (phi, psi) in enumerate(tors):
        a, b = rm[i], rm[i + 1]
        fixed += [(a["O5"], a["C1"], b["O4"], b["C4"], phi), (a["C1"], b["O4"], b["C4"], b["C5"], psi)]
    e = mmff_relax(m, maxIters=5000, fixed_torsions=fixed)
    assert check_stereo(m)
    return m, rm, e


def emit(job, model, atoms, q, title, extra):
    d = os.path.join(OUT, job)
    os.makedirs(d, exist_ok=True)
    if model == "PM7":
        open(os.path.join(d, "eu_fix.txt"), "w", newline="\n").write(EU_FIX)
    el = [a[0] for a in atoms]; xyz = np.array([a[1] for a in atoms])
    n_int, perm = write_mop_mixed(os.path.join(d, job + ".mop"), el, xyz, KEYS[model].format(q=q), title)
    with open(os.path.join(d, job + "_start.xyz"), "w") as f:
        f.write(f"{len(el)}\n{title}\n")
        for i in perm:
            f.write(f"{el[i]:2s} {xyz[i][0]:12.6f} {xyz[i][1]:12.6f} {xyz[i][2]:12.6f}\n")
    rec = dict(job=job, model=model, charge=q, natoms=len(el), n_internal=n_int, perm=[int(x) for x in perm])
    rec.update(extra)
    manifest.append(rec)


def residue_index_map(rm, perm):
    """Residue atom indices in the *reordered* (MOPAC) numbering (0-based)."""
    inv = {old: new for new, old in enumerate(perm)}
    out = []
    for k, d in enumerate(rm):
        out.append({kk: (inv[vv] if isinstance(vv, (int, np.integer)) else vv) for kk, vv in d.items() if vv is not None})
    return out


# ---------------------------------------------------------------- main series
for pname, pattern in PATTERNS.items():
    q = 3 - pattern.count("A")
    for cname, tors in CONFORMERS.items():
        m, rm, e = relax_chain(pattern, tors)
        n, h = helix_params(m, rm)
        nA = pattern.count("A")
        if nA == 0:
            modes = ["ester", "hyd"]
        elif nA == 1:
            modes = ["bid", "hyd"]
        else:
            modes = ["bid", "mono2", "hyd"]
        for mode in modes:
            res = build_complex(m, mode, rmap=rm, n_water=9, carbox_waters=2)
            if res is None:
                continue
            atoms, info = res
            for model in ("PM7", "RM1"):
                if model == "RM1" and q < 0:
                    continue
                job = f"{model}_{pname}_{cname}_{mode}"
                title = f"{job} pattern={pattern} conf={cname} start={mode} q={q}"
                extra = dict(info)
                extra.update(dict(series="main", pattern=pattern, pname=pname, conf=cname, mode=mode,
                                  chain_mmff_E=float(e), helix_n=float(n), helix_rise=float(h),
                                  torsions=[(round(p, 1), round(s, 1)) for p, s in glycosidic_torsions(m, m.GetConformer(), rm)]))
                emit(job, model, atoms, q, title, extra)
                manifest[-1]["residues"] = residue_index_map(rm, manifest[-1]["perm"])
            print(pname, cname, mode, len(atoms), "atoms; 1st-shell waters", info["n_first_shell_waters"], "targets", info.get("target_dist"))

# ---------------------------------------------------------------- sensitivity: hydration of carboxylates (0 and 4 waters per carboxylate)
for pname in ("A6", "B3", "L3", "P3"):
    pattern = PATTERNS[pname]; q = 3 - pattern.count("A")
    m, rm, e = relax_chain(pattern, CONFORMERS["h2"])
    for cw in (0, 4):
        atoms, info = build_complex(m, "bid", rmap=rm, n_water=9, carbox_waters=cw)
        for model in ("PM7", "RM1"):
            if model == "RM1" and q < 0:
                continue
            job = f"{model}_{pname}_h2_bid_cw{cw}"
            extra = dict(info); extra.update(dict(series="sens_hydration", pattern=pattern, pname=pname, conf="h2", mode="bid", carbox_waters=cw))
            emit(job, model, atoms, q, f"{job} pattern={pattern} conf=h2 start=bid q={q} carbox_waters={cw}", extra)
            manifest[-1]["residues"] = residue_index_map(rm, manifest[-1]["perm"])
        print(pname, "cw", cw, len(atoms))

# ---------------------------------------------------------------- controls
def ttp_cluster(n):
    center = np.zeros(3)
    atoms = [("Eu", center)]
    if n == 9:
        hh = 2.48 * 0.62; rr = np.sqrt(2.48 ** 2 - hh ** 2)
        pts = []
        for z in (+hh, -hh):
            for k in range(3):
                a = 2 * np.pi * k / 3
                pts.append([rr * np.cos(a), rr * np.sin(a), z])
        for k in range(3):
            a = 2 * np.pi * k / 3 + np.pi / 3
            pts.append([2.58 * np.cos(a), 2.58 * np.sin(a), 0.0])
        pts = np.array(pts)
    else:
        pts = 2.48 * fib_sphere(min(n, 9))
    for p in pts:
        atoms += water_atoms(np.array(p), center)
    if n > 9:
        for p in 4.6 * fib_sphere(n - 9) + np.array([0.1, 0.05, 0.02]):
            atoms += water_atoms(np.array(p), center)
    return atoms

for model in ("PM7", "RM1"):
    for n in (8, 9, 10, 12):
        job = f"{model}_CTRL_Eu_w{n}"
        emit(job, model, ttp_cluster(n), 3, f"{job} [Eu(H2O){n}]3+ gas-phase cluster", dict(series="control", pattern="water", pname="CTRL_water", conf="-", mode=f"w{n}", n_water=n))
    for kind, lab in (("A", "GalA"), ("M", "GalAMe")):
        m = assemble(kind, []); mmff_relax(m, maxIters=2000); rm = residue_atom_map(m)
        q = 3 - (1 if kind == "A" else 0)
        for mode in (["bid", "hyd"] if kind == "A" else ["ester", "hyd"]):
            atoms, info = build_complex(m, mode, rmap=rm, n_water=9, carbox_waters=(2 if kind == "A" else 0))
            job = f"{model}_CTRL_{lab}_{mode}"
            extra = dict(info); extra.update(dict(series="control", pattern=kind, pname="CTRL_" + lab, conf="-", mode=mode))
            emit(job, model, atoms, q, f"{job} monomer control q={q}", extra)
            manifest[-1]["residues"] = residue_index_map(rm, manifest[-1]["perm"])
    # acetate
    ac = Chem.AddHs(Chem.MolFromSmiles("CC(=O)[O-]")); AllChem.EmbedMolecule(ac, randomSeed=3); AllChem.MMFFOptimizeMolecule(ac)
    pos = ac.GetConformer().GetPositions(); el = np.array([a.GetSymbol() for a in ac.GetAtoms()])
    Oidx = [i for i, e in enumerate(el) if e == "O"]
    Cc = [i for i, e in enumerate(el) if e == "C" and sum(1 for j in Oidx if np.linalg.norm(pos[i] - pos[j]) < 1.6) == 2][0]
    eu, f = find_site(pos[Oidx], pos, el, Oidx, soft=[Cc])
    waters, nf = fill_waters(eu, pos, el, Oidx, 9)
    atoms = [(el[i], pos[i]) for i in range(len(el))] + [("Eu", eu)]
    for w in waters: atoms += water_atoms(w, eu)
    emit(f"{model}_CTRL_acetate_bid", model, atoms, 2, f"{model}_CTRL_acetate_bid Eu + acetate + 9 H2O q=2", dict(series="control", pattern="acetate", pname="CTRL_acetate", conf="-", mode="bid", n_water=9, n_first_shell_waters=nf))

json.dump(manifest, open(os.path.join(OUT, "manifest.json"), "w"), indent=1)
print("jobs:", len(manifest), "| PM7:", sum(1 for r in manifest if r["model"] == "PM7"), "RM1:", sum(1 for r in manifest if r["model"] == "RM1"))
