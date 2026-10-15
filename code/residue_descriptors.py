"""Recompute residue-level descriptors from final geometries with a strict carboxylate definition.

Only residues of kind A (free carboxylate, -COO-) count towards carboxylate denticity and multi-carboxylate
motifs; ester (M) and protonated acid (H) contacts through their C6 oxygens are counted separately.
Usage: python residue_descriptors.py results.csv final_geoms.json manifest.json out.csv
"""
import sys, json, numpy as np, pandas as pd

CUT = 3.0


def residue_row(el, xyz, residues, cutoff=CUT):
    xyz = np.asarray(xyz, float)
    eu = list(el).index("Eu")
    d = np.linalg.norm(xyz - xyz[eu], axis=1)
    out = dict(n_residues_bound=0, n_carboxylates_bound=0, n_bidentate_carboxylate=0, n_monodentate_carboxylate=0,
               n_esters_bound=0, n_acids_bound=0, n_backbone_O=0, n_hydroxyl_O=0, motif="", residues_bound="")
    tags = []; bound = []
    for k, r in enumerate(residues):
        kind = r.get("kind")
        if kind not in ("A", "H", "M"):
            continue
        n6 = sum(1 for key in ("O6a", "O6b") if r.get(key) is not None and d[r[key]] < cutoff)
        back = [key for key in ("O5", "O1", "O4") if r.get(key) is not None and d[r[key]] < cutoff]
        hyd = [key for key in ("O2", "O3") if r.get(key) is not None and d[r[key]] < cutoff]
        if n6 == 0 and not back and not hyd:
            continue
        bound.append(k)
        out["n_backbone_O"] += len(back); out["n_hydroxyl_O"] += len(hyd)
        if kind == "A":
            if n6 >= 1:
                out["n_carboxylates_bound"] += 1
            if n6 == 2:
                out["n_bidentate_carboxylate"] += 1; tag = "bi"
            elif n6 == 1:
                out["n_monodentate_carboxylate"] += 1; tag = "mono"
            else:
                tag = ""
        elif kind == "M":
            if n6 >= 1:
                out["n_esters_bound"] += 1
            tag = {0: "", 1: "est", 2: "est2"}[n6]
        else:
            if n6 >= 1:
                out["n_acids_bound"] += 1
            tag = {0: "", 1: "acid", 2: "acid2"}[n6]
        extra = "".join(back + hyd)
        if extra:
            tag = (tag + "+" if tag else "+") + extra
        tags.append(f"{kind}{k}:{tag}")
    out["n_residues_bound"] = len(bound)
    out["residues_bound"] = ",".join(str(k) for k in bound)
    out["motif"] = " ".join(tags)
    return out


if __name__ == "__main__":
    csv, geoms, manifest, outcsv = sys.argv[1:5]
    df = pd.read_csv(csv)
    G = json.load(open(geoms))
    man = {m["job"]: m for m in json.load(open(manifest))}
    new = {}
    for job in df.job:
        rec = man.get(job)
        if job in G and rec and rec.get("residues"):
            new[job] = residue_row(G[job]["el"], G[job]["xyz"], rec["residues"])
    nd = pd.DataFrame.from_dict(new, orient="index")
    # replace old residue-level columns
    for col in nd.columns:
        df[col] = df.job.map(nd[col])
    df = df.drop(columns=[c for c in ("chelate_ring_O", "hydroxyl_bound") if c in df.columns])
    df.to_csv(outcsv, index=False)
    print("written", outcsv, df.shape)
