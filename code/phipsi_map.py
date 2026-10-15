"""MMFF94 (distance-dependent dielectric, eps=4r) phi/psi map of the alpha-(1->4) GalpA disaccharide."""
import numpy as np, json, sys
from rdkit import Chem
from rdkit.Chem import AllChem
from galA import residue_atom_map, glycosidic_torsions
from assemble import assemble, mmff_relax

def run_map(kind="M", step=15):
    grid = np.arange(-180, 180, step)
    E = np.full((len(grid), len(grid)), np.nan)
    best = {}
    for i, phi in enumerate(grid):
        for j, psi in enumerate(grid):
            m = assemble(kind * 2, [(float(phi), float(psi))])
            rm = residue_atom_map(m)
            a, b = rm[0], rm[1]
            fixed = [(a["O5"], a["C1"], b["O4"], b["C4"], float(phi)), (a["C1"], b["O4"], b["C4"], b["C5"], float(psi))]
            try:
                e = mmff_relax(m, maxIters=2000, fixed_torsions=fixed)
            except Exception as ex:
                e = np.nan
            E[i, j] = e
    return grid, E

if __name__ == "__main__":
    out = {}
    for kind in ["M", "A", "H"]:
        grid, E = run_map(kind, step=15)
        Emin = np.nanmin(E)
        rel = E - Emin
        out[kind] = {"grid": grid.tolist(), "E": rel.tolist()}
        # local minima
        mins = []
        n = len(grid)
        for i in range(n):
            for j in range(n):
                v = rel[i, j]
                nb = [rel[(i + di) % n, (j + dj) % n] for di in (-1, 0, 1) for dj in (-1, 0, 1) if (di, dj) != (0, 0)]
                if v <= min(nb) and v < 6.0:
                    mins.append((float(grid[i]), float(grid[j]), float(v)))
        mins.sort(key=lambda x: x[2])
        out[kind]["minima"] = mins
        print(kind, "minima (phi, psi, dE kcal/mol):", [(p, q, round(e, 2)) for p, q, e in mins])
    json.dump(out, open("phipsi_map.json", "w"))
