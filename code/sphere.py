"""First-coordination-sphere analysis of Eu/pectin clusters from MOPAC .aux geometries."""
import re, numpy as np


def classify_oxygens(el, xyz):
    """Return dict idx -> type in {'carboxylate','carbonyl_ester','ether_ester','hydroxyl','ring','glycosidic','water','hydroxide','oxonium'}"""
    n = len(el)
    xyz = np.asarray(xyz)
    heavy = [i for i in range(n) if el[i] not in ("H", "Eu")]
    Hs = [i for i in range(n) if el[i] == "H"]
    def nb(i, sym, cut):
        return [j for j in range(n) if j != i and el[j] == sym and np.linalg.norm(xyz[j] - xyz[i]) < cut]
    types = {}
    for i in range(n):
        if el[i] != "O":
            continue
        hs = nb(i, "H", 1.25)
        cs = nb(i, "C", 1.65)
        if not cs:
            types[i] = {0: "oxide", 1: "hydroxide", 2: "water", 3: "oxonium"}.get(len(hs), "water?")
            continue
        if len(cs) == 2:
            # ring O5 (both carbons in ring) or glycosidic / ester ether O
            c_methyl = [c for c in cs if len(nb(c, "H", 1.2)) == 3]
            types[i] = "ether_ester" if c_methyl else "ring_or_glycosidic"
            continue
        c = cs[0]
        # carbon environment
        c_os = nb(c, "O", 1.65)
        c_hs = nb(c, "H", 1.2)
        if len(c_os) == 2 and len(c_hs) == 0:
            # carboxyl carbon: this O is carboxylate O, carbonyl O, or acid OH, or ester carbonyl/ether
            other = [o for o in c_os if o != i][0]
            other_c = nb(other, "C", 1.65)
            if hs:
                types[i] = "acid_OH"
            elif len(other_c) == 2 or len(nb(other, "H", 1.25)) == 1:
                # other O is ester ether (bonded to CH3) or acid OH -> this is a carbonyl O
                types[i] = "carbonyl_ester" if len(other_c) == 2 else "carbonyl_acid"
            else:
                types[i] = "carboxylate"
        else:
            types[i] = "hydroxyl" if hs else "alkoxide"
    return types


def first_sphere(el, xyz, cutoff=3.0):
    el = list(el); xyz = np.asarray(xyz)
    eu = el.index("Eu")
    types = classify_oxygens(el, xyz)
    out = []
    for i, t in types.items():
        d = float(np.linalg.norm(xyz[i] - xyz[eu]))
        if d < cutoff:
            out.append((d, i, t))
    out.sort()
    return out, types


def summarize(el, xyz, cutoff=3.0):
    fs, types = first_sphere(el, xyz, cutoff)
    n_pect = sum(1 for d, i, t in fs if t not in ("water", "hydroxide", "oxonium", "oxide", "water?"))
    n_wat = sum(1 for d, i, t in fs if t in ("water", "hydroxide", "oxonium", "oxide", "water?"))
    comp = {}
    for d, i, t in fs:
        comp[t] = comp.get(t, 0) + 1
    return dict(CN=len(fs), n_pectin=n_pect, n_water=n_wat, composition=comp,
                d_pectin=[round(d, 3) for d, i, t in fs if t not in ("water",)],
                d_water=[round(d, 3) for d, i, t in fs if t == "water"],
                collapsed=[round(d, 2) for d, i, t in fs if d < 2.1])
