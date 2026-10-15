"""Place Eu(III) and explicit water molecules around a homogalacturonan oligomer for MOPAC inputs."""
import numpy as np
from scipy.optimize import minimize
from galA import residue_atom_map

R_EUO = 2.48        # target Eu-O distance for first-sphere placement
R_SHELL2 = 4.6      # second-shell distance
N_WATER = 9


def fib_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.c_[np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)]


def water_atoms(o_pos, center, tilt_vec=None):
    """Water with O at o_pos, H's pointing away from center (lone pairs toward the ion)."""
    v = o_pos - center
    v = v / np.linalg.norm(v)
    a = np.array([1.0, 0, 0]) if abs(v[0]) < 0.9 else np.array([0, 1.0, 0])
    u = np.cross(v, a); u /= np.linalg.norm(u)
    if tilt_vec is not None:
        u = tilt_vec / np.linalg.norm(tilt_vec)
    ang = np.deg2rad(104.5 / 2)
    h1 = o_pos + 0.96 * (np.cos(ang) * v + np.sin(ang) * u)
    h2 = o_pos + 0.96 * (np.cos(ang) * v - np.sin(ang) * u)
    return [("O", o_pos), ("H", h1), ("H", h2)]


def chain_arrays(mol):
    pos = mol.GetConformer().GetPositions()
    el = np.array([a.GetSymbol() for a in mol.GetAtoms()])
    return pos, el


def carboxylate_residues(rmap):
    return [k for k, d in enumerate(rmap) if d["kind"] == "A"]


def site_objective(x, targets, pos, el, r0=R_EUO, soft=(), weights=None):
    """targets: coordinates of oxygen atoms that should sit at r0 from Eu.
    soft: indices (into pos) of atoms allowed to be closer (carboxyl carbon / partner O of a targeted group)."""
    d = np.linalg.norm(targets - x, axis=1)
    w = np.ones(len(d)) * 3.0 if weights is None else np.asarray(weights)
    f = np.sum(w * (d - r0) ** 2)
    dd = np.linalg.norm(pos - x, axis=1)
    heavy = (el != "H")
    thr = np.where(heavy, 2.9, 2.3)
    for k in soft:
        thr[k] = 2.6 if el[k] != "H" else 2.3
    f += np.sum(np.clip(thr - dd, 0, None) ** 2)
    return f


def find_site(targets, pos, el, exclude_idx, r0=R_EUO, soft=(), weights=None):
    mask = np.ones(len(pos), bool)
    mask[exclude_idx] = False
    idx = np.where(mask)[0]
    p, e = pos[mask], el[mask]
    soft_local = [int(np.where(idx == k)[0][0]) for k in soft if k in idx]
    c = targets[0] if weights is not None else targets.mean(0)
    best = None
    for direction in fib_sphere(60):
        x0 = c + 2.0 * direction
        res = minimize(site_objective, x0, args=(targets, p, e, r0, soft_local, weights), method="Nelder-Mead",
                       options={"xatol": 1e-3, "fatol": 1e-5, "maxiter": 3000})
        if best is None or res.fun < best.fun:
            best = res
    return best.x, best.fun


def carboxylate_pairs(pos, rmap):
    A = carboxylate_residues(rmap)
    pairs = []
    for i in range(len(A)):
        for j in range(i + 1, len(A)):
            ra, rb = rmap[A[i]], rmap[A[j]]
            # closest O-O distance between the two carboxylate groups
            dmin = min(np.linalg.norm(pos[oa] - pos[ob]) for oa in (ra["O6a"], ra["O6b"]) for ob in (rb["O6a"], rb["O6b"]))
            pairs.append((dmin, A[i], A[j]))
    pairs.sort()
    return pairs


def choose_targets(mol, rmap, mode):
    """Return list of (residue index, [oxygen atom indices]) used as coordination targets."""
    pos, el = chain_arrays(mol)
    A = carboxylate_residues(rmap)
    if not A:
        return None
    if mode == "bid" or len(A) == 1:
        # bidentate on the most central free carboxylate
        rc = min(A, key=lambda r: (abs(r - (len(rmap) - 1) / 2), r))
        d = rmap[rc]
        return [(rc, [d["O6a"], d["O6b"]])]
    pairs = carboxylate_pairs(pos, rmap)
    if mode in ("chel", "mono2"):
        dist, ri, rj = pairs[0]
    elif mode == "bridge":
        dist, ri, rj = pairs[-1]
    else:
        raise ValueError(mode)
    # order: the more central residue first (primary, monodentate anchor)
    if abs(rj - (len(rmap) - 1) / 2) < abs(ri - (len(rmap) - 1) / 2):
        ri, rj = rj, ri
    out = []
    for ra, rb in ((ri, rj), (rj, ri)):
        da, db = rmap[ra], rmap[rb]
        cands = [da["O6a"], da["O6b"]]
        o = min(cands, key=lambda k: min(np.linalg.norm(pos[k] - pos[q]) for q in (db["O6a"], db["O6b"])))
        out.append((ra, [o]))
    return out


def ttp_directions():
    """Unit vectors of an ideal tricapped trigonal prism (9 vertices)."""
    h = 0.62; rr = np.sqrt(1 - h ** 2)
    pts = []
    for z in (+h, -h):
        for k in range(3):
            a = 2 * np.pi * k / 3
            pts.append([rr * np.cos(a), rr * np.sin(a), z])
    for k in range(3):
        a = 2 * np.pi * k / 3 + np.pi / 3
        pts.append([np.cos(a), np.sin(a), 0.0])
    return np.array(pts)


def random_rotations(n, seed=0):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        q = rng.normal(size=4); q /= np.linalg.norm(q)
        a, b, c, d = q
        out.append(np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                             [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                             [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]]))
    return out


def fill_waters(eu, pos, el, first_sphere_O, n_water=N_WATER, n_rot=1500):
    """Place water oxygens on the free vertices of a tricapped trigonal prism around Eu
    (vertices occupied by chain oxygens or clashing with the chain are skipped);
    remaining waters go to a second shell."""
    occupied = [pos[i] for i in first_sphere_O]
    heavy = el != "H"
    dirs = ttp_directions()
    best = None
    for R in random_rotations(n_rot):
        verts = eu + R_EUO * (dirs @ R.T)
        free = []
        for v in verts:
            if any(np.linalg.norm(v - o) < 2.3 for o in occupied):
                continue  # vertex taken by a coordinating chain oxygen
            if np.min(np.linalg.norm(pos[heavy] - v, axis=1)) < 2.5: continue
            if np.min(np.linalg.norm(pos[~heavy] - v, axis=1)) < 1.95: continue
            if any(np.linalg.norm(v - o) < 2.6 for o in occupied): continue
            free.append(v)
        # score: number of free vertices, then how well occupied O's sit on vertices
        fit = sum(min(np.linalg.norm(verts - o, axis=1)) for o in occupied)
        score = (len(free), -fit)
        if best is None or score > best[0]:
            best = (score, free)
    n_first_target = max(0, 9 - len(first_sphere_O))
    waters = [np.array(v) for v in best[1][:min(n_first_target, n_water)]]
    n_first = len(waters)
    # second (and, if needed, third) shell
    def ok_second(x):
        allpos = np.vstack([pos, np.array(occupied + waters)]) if (occupied or waters) else pos
        if np.min(np.linalg.norm(allpos - x, axis=1)) < 2.75: return False
        return True
    for r_shell in (R_SHELL2, R_SHELL2 + 1.9, R_SHELL2 + 3.8):
        cand2 = eu + r_shell * fib_sphere(600)
        while len(waters) < n_water:
            scores = []
            for x in cand2:
                if not ok_second(x):
                    continue
                # prefer H-bonding distance (~2.9 A) to the nearest already-placed water
                dmin = min(np.linalg.norm(w - x) for w in waters) if waters else 3.0
                scores.append((abs(dmin - 2.9), tuple(x)))
            if not scores:
                break
            scores.sort()
            waters.append(np.array(scores[0][1]))
        if len(waters) >= n_water:
            break
    return waters, n_first


def hbond_water(o_pos, c_pos, occupied, all_pos, r=2.85):
    """Place a water donating an H-bond to carboxylate oxygen o_pos (carboxyl carbon at c_pos).
    Returns [(el,xyz)...] or None."""
    v = o_pos - c_pos; v /= np.linalg.norm(v)
    best = None
    for d in fib_sphere(300):
        if np.dot(d, v) < 0.25:      # keep on the lone-pair side
            continue
        x = o_pos + r * d
        dmin = np.min(np.linalg.norm(all_pos - x, axis=1))
        if dmin < 2.7:
            continue
        docc = min([np.linalg.norm(x - q) for q in occupied], default=9.0)
        if docc < 2.8:
            continue
        score = np.dot(d, v) + 0.1 * min(dmin, 4.0)
        if best is None or score > best[0]:
            best = (score, x)
    if best is None:
        return None
    x = best[1]
    # orientation: one H toward the carboxylate O, the other roughly perpendicular
    u = o_pos - x; u /= np.linalg.norm(u)
    a = np.array([1.0, 0, 0]) if abs(u[0]) < 0.9 else np.array([0, 1.0, 0])
    w = np.cross(u, a); w /= np.linalg.norm(w)
    h1 = x + 0.96 * u
    ang = np.deg2rad(104.5)
    h2 = x + 0.96 * (np.cos(ang) * u + np.sin(ang) * w)
    return [("O", x), ("H", h1), ("H", h2)]


def hydrate_carboxylates(atoms, rmap, pos, el, per_carboxylate=2):
    """Add H-bonded waters to every free carboxylate (both O's if possible)."""
    allpos = np.array([a[1] for a in atoms])
    occupied = [a[1] for a in atoms if a[0] == "O"]
    eu = [a[1] for a in atoms if a[0] == "Eu"][0]
    added = []
    for d in rmap:
        if d["kind"] != "A":
            continue
        targets = ([d["O6a"], d["O6b"]] * ((per_carboxylate + 1) // 2))[:per_carboxylate]
        targets = [o for o in targets if np.linalg.norm(pos[o] - eu) > 3.2]
        for o in targets:
            w = hbond_water(pos[o], pos[d["C6"]], occupied, allpos)
            if w is None:
                continue
            added += w
            occupied.append(w[0][1])
            allpos = np.vstack([allpos, np.array([x[1] for x in w])])
    return added


def build_complex(mol, mode, rmap=None, n_water=N_WATER, hyd_offset=3.6, carbox_waters=0):
    """Return (atoms list [(el, xyz)], info dict)."""
    if rmap is None:
        rmap = residue_atom_map(mol)
    pos, el = chain_arrays(mol)
    info = {"mode": mode}
    if mode in ("chel", "bridge", "bid", "mono2"):
        tg = choose_targets(mol, rmap, mode)
        if tg is None:
            return None
        target_idx = [o for _, oo in tg for o in oo]
        soft = []
        for r, oo in tg:
            soft += [rmap[r]["C6"], rmap[r]["O6a"], rmap[r]["O6b"]]
        soft = [k for k in soft if k not in target_idx]
        weights = None
        if mode == "mono2" and len(target_idx) == 2:
            d12 = np.linalg.norm(pos[target_idx[0]] - pos[target_idx[1]])
            weights = [3.0, 0.15 if d12 < 7.5 else 0.0]   # firm monodentate anchor; weak pull toward the 2nd carboxylate if reachable
        eu, f = find_site(pos[target_idx], pos, el, target_idx, soft=soft, weights=weights)
        info["target_residues"] = [r for r, _ in tg]
        info["site_obj"] = float(f)
        info["target_dist"] = [round(float(np.linalg.norm(pos[o] - eu)), 2) for o in target_idx]
        # which chain O's are within 2.9 A (counted as first sphere for water filling)
        dO = np.linalg.norm(pos - eu, axis=1)
        first = [i for i in range(len(pos)) if el[i] == "O" and dO[i] < 2.9]
        waters, n_first = fill_waters(eu, pos, el, first, n_water)
    elif mode == "hyd":
        tg = choose_targets(mol, rmap, "bid")
        if tg is None:
            # fully esterified: use ester carbonyl O's of the two central residues
            tg = [(2, [rmap[2]["O6a"]]), (3, [rmap[3]["O6a"]])] if len(rmap) > 3 else [(0, [rmap[0]["O6a"]])]
        target_idx = [o for _, oo in tg for o in oo]
        site, f = find_site(pos[target_idx], pos, el, target_idx)
        # push outward: direction from chain centroid (heavy atoms) to site
        cen = pos[el != "H"].mean(0)
        v = site - pos[target_idx].mean(0)
        v /= np.linalg.norm(v)
        eu = site + hyd_offset * v
        # make sure no chain atom within 4.6 A of Eu (so that 9 waters fit); push further if needed
        for _ in range(20):
            if np.min(np.linalg.norm(pos - eu, axis=1)) >= 4.9:
                break
            eu = eu + 0.3 * v
        info["target_residues"] = [r for r, _ in tg]
        waters, n_first = fill_waters(eu, pos, el, [], n_water)
    elif mode == "ester":
        # Eu near ester carbonyl oxygens (control for fully esterified chain)
        tg = [(2, [rmap[2]["O6a"]]), (3, [rmap[3]["O6a"]])] if len(rmap) > 3 else [(0, [rmap[0]["O6a"]])]
        target_idx = [o for _, oo in tg for o in oo]
        eu, f = find_site(pos[target_idx], pos, el, target_idx)
        info["target_residues"] = [r for r, _ in tg]
        dO = np.linalg.norm(pos - eu, axis=1)
        first = [i for i in range(len(pos)) if el[i] == "O" and dO[i] < 2.9]
        waters, n_first = fill_waters(eu, pos, el, first, n_water)
    else:
        raise ValueError(mode)
    info["n_first_shell_waters"] = n_first
    info["eu_min_chain_O"] = float(min(np.linalg.norm(pos[i] - eu) for i in range(len(pos)) if el[i] == "O"))
    atoms = [(str(el[i]), pos[i]) for i in range(len(pos))]
    atoms.append(("Eu", eu))
    for w in waters:
        atoms += water_atoms(w, eu)
    n_eu_waters = len(waters)
    if carbox_waters:
        extra = hydrate_carboxylates(atoms, rmap, pos, el, per_carboxylate=carbox_waters)
        atoms += extra
        info["n_carboxylate_waters"] = len(extra) // 3
    info["n_water"] = n_eu_waters + info.get("n_carboxylate_waters", 0)
    return atoms, info


def write_mop(fname, atoms, keywords, title):
    with open(fname, "w", newline="\n") as f:
        f.write(keywords + "\n" + title + "\n\n")
        for e, p in atoms:
            f.write(f"{e:2s} {p[0]:12.6f} 1 {p[1]:12.6f} 1 {p[2]:12.6f} 1\n")
