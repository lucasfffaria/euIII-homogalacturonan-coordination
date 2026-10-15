"""Harmonic vibrational analysis on a surface with frozen bond lengths (Section 2.5 of the article, S4 of the SI).

A structure optimized with constrained O-H distances d_k is stationary on the constrained surface: its Cartesian gradient
is g = sum_k lambda_k grad(d_k) (+ a small residual).  The harmonic frequencies of the constrained system are the
eigenvalues of the Hessian of the Lagrangian,  W = H - sum_k lambda_k hess(d_k),  projected on the subspace orthogonal
(in mass-weighted coordinates) to the constrained stretches and to the overall translations and rotations.

Units: coordinates in Angstrom, masses in amu, Hessian in mdyn/A, gradient in kcal/mol/A (MOPAC), frequencies in cm^-1.
Functions: constrained_modes (the analysis), read_force_aux / read_force_gradient / read_mop_constraints (MOPAC readers),
analyse_job (one FORCE job).
"""
import re
import numpy as np

KCAL_PER_MDYNA = 143.9325          # 1 mdyn*A = 143.9325 kcal/mol
FREQ_FACTOR = 1302.7914            # cm^-1 from sqrt(mdyn/A/amu)


def _tr_vectors(xyz, masses):
    """Translations and rotations in mass-weighted Cartesian coordinates (rows), about the centre of mass."""
    n = len(masses); sm = np.sqrt(np.repeat(masses, 3)); com = (masses[:, None] * xyz).sum(0) / masses.sum(); x = xyz - com
    V = []
    for a in range(3):
        v = np.zeros((n, 3)); v[:, a] = 1.0; V.append(v.ravel() * sm)
    for a in range(3):
        e = np.zeros(3); e[a] = 1.0; V.append(np.cross(e, x).ravel() * sm)
    return np.array(V)


def _dist_derivs(xyz, i, j):
    """Gradient (3N) and Hessian (3N x 3N) of the distance |r_i - r_j|."""
    n = len(xyz); u = xyz[i] - xyz[j]; d = np.linalg.norm(u); u = u / d
    g = np.zeros(3 * n); g[3 * i:3 * i + 3] = u; g[3 * j:3 * j + 3] = -u
    B = (np.eye(3) - np.outer(u, u)) / d
    Hd = np.zeros((3 * n, 3 * n))
    for a, sa in ((i, 1), (j, -1)):
        for b, sb in ((i, 1), (j, -1)):
            Hd[3 * a:3 * a + 3, 3 * b:3 * b + 3] = sa * sb * B
    return d, g, Hd


def constrained_modes(xyz, masses, H, grad, pairs, tol_null=1e-8):
    """xyz (N,3) A; masses (N,) amu; H (3N,3N) mdyn/A; grad (3N,) kcal/mol/A; pairs = [(i, j), ...] frozen distances.
    Returns a dict: freq (cm^-1, negative = imaginary), n_imag, lambdas (kcal/mol/A), gnorm_full, gnorm_tangent
    (norm of the gradient left after removing the constraint forces), freq_unconstrained (plain analysis, for reference)."""
    xyz = np.asarray(xyz, float); masses = np.asarray(masses, float); n = len(masses)
    H = 0.5 * (H + H.T); g = np.asarray(grad, float) / KCAL_PER_MDYNA
    sm = np.sqrt(np.repeat(masses, 3)); TR = _tr_vectors(xyz, masses)
    J = np.zeros((len(pairs), 3 * n)); Hd = []
    for k, (i, j) in enumerate(pairs):
        _, J[k], hk = _dist_derivs(xyz, i, j); Hd.append(hk)
    lam = np.linalg.lstsq(J.T, g, rcond=None)[0] if len(pairs) else np.zeros(0)
    g_t = g - J.T @ lam if len(pairs) else g
    W = H.copy()
    for k in range(len(pairs)):
        W -= lam[k] * Hd[k]

    def modes(Hm, extra):
        A = np.vstack([TR] + ([extra / sm] if len(extra) else []))          # vectors to remove, mass-weighted
        U, s, Vt = np.linalg.svd(A, full_matrices=True)
        r = int((s > tol_null * s.max()).sum()); Bn = Vt[r:].T               # orthonormal basis of the complement
        w = np.linalg.eigvalsh(Bn.T @ (Hm / np.outer(sm, sm)) @ Bn)
        return np.sign(w) * FREQ_FACTOR * np.sqrt(np.abs(w)), r

    f_con, r_con = modes(W, J)
    f_unc, _ = modes(H, np.zeros((0, 3 * n)))
    return dict(freq=f_con, n_imag=int((f_con < 0).sum()), n_modes=len(f_con), n_removed=r_con,
                lambdas=lam * KCAL_PER_MDYNA, gnorm_full=float(np.linalg.norm(g) * KCAL_PER_MDYNA),
                gnorm_tangent=float(np.linalg.norm(g_t) * KCAL_PER_MDYNA), freq_unconstrained=f_unc)


# ----------------------------------------------------------------------------------------------- MOPAC readers
_FLOAT = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?")


def _aux_block(txt, key):
    m = re.search(re.escape(key) + r"[^\[\n]*\[\s*(\d+)\s*\]\s*=\s*\n", txt)
    if not m:
        return None
    n = int(m.group(1)); vals = []
    for line in txt[m.end():].splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            if s.startswith("#") and not vals:
                continue
            if vals:
                break
            continue
        if re.match(r"^[A-Z_#.]", s) and not re.match(r"^[-+]?\d", s):
            break
        vals += [v.replace("D", "E").replace("d", "E") for v in _FLOAT.findall(s)]   # tolerates numbers printed without a separating blank
        if len(vals) >= n:
            break
    return np.array(vals[:n], float)


def read_force_aux(path):
    """masses (amu), xyz (A, orientation of the FORCE calculation), H (mdyn/A, un-mass-weighted), MOPAC frequencies."""
    txt = open(path, errors="ignore").read()
    masses = _aux_block(txt, " ISOTOPIC_MASSES")
    xyz = _aux_block(txt, " ATOM_X_FORCE:ANGSTROMS")
    if xyz is None:
        xyz = _aux_block(txt, " ORIENTATION_ATOM_X:ANGSTROMS")
    hm = _aux_block(txt, " HESSIAN_MATRIX:")
    freq = _aux_block(txt, " VIB._FREQ:CM(-1)")
    n = len(masses); n3 = 3 * n
    assert len(hm) == n3 * (n3 + 1) // 2, (len(hm), n3)
    Hmw = np.zeros((n3, n3)); Hmw[np.tril_indices(n3)] = hm; Hmw = Hmw + Hmw.T - np.diag(np.diag(Hmw))
    sm = np.sqrt(np.repeat(masses, 3))
    return masses, xyz.reshape(n, 3), Hmw * np.outer(sm, sm), freq, Hmw


def read_force_gradient(path, natoms):
    """Cartesian derivatives (kcal/mol/A) printed by FORCE in the .out file."""
    L = open(path, errors="ignore").read().splitlines()
    i0 = max(i for i, l in enumerate(L) if "CARTESIAN COORDINATE DERIVATIVES" in l)
    g = []
    for l in L[i0 + 1:]:
        s = l.split()
        if len(s) == 5 and s[0].isdigit():
            g.append([float(v) for v in s[2:5]])
            if len(g) == natoms:
                break
    assert len(g) == natoms, (len(g), natoms)
    return np.array(g).ravel()


def read_mop_constraints(path):
    """Frozen distances of a MOPAC input in mixed coordinates: atoms whose first coordinate has flag 0 and a connectivity
    entry (internal coordinate).  Returns [(i, j)] 0-based (atom, atom it is bonded to) and the element list."""
    L = open(path, errors="ignore").read().splitlines()[3:]
    pairs, el = [], []
    for l in L:
        s = l.split()
        if len(s) < 7:
            continue
        el.append(s[0])
        if len(s) >= 10 and int(s[7]) > 0 and int(s[2]) == 0:
            pairs.append((len(el) - 1, int(s[7]) - 1))
    return pairs, el


def analyse_job(mop, out, aux):
    pairs, el = read_mop_constraints(mop)
    masses, xyz, H, f_mopac, Hmw = read_force_aux(aux)
    grad = read_force_gradient(out, len(el))
    res = constrained_modes(xyz, masses, H, grad, pairs)
    # internal checks of the reading: (1) frequencies recomputed from the Hessian as read reproduce those printed by MOPAC;
    # (2) the frozen distances measured on the FORCE coordinates are the constrained values
    w = np.linalg.eigvalsh(Hmw); f_re = np.sort(np.sign(w) * FREQ_FACTOR * np.sqrt(np.abs(w)))
    res["check_freq_maxdev"] = float(np.abs(np.sort(f_re)[6:] - np.sort(f_mopac)[6:]).max()) if f_mopac is not None and len(f_mopac) == len(f_re) else np.nan
    res["frozen_distances"] = np.array([np.linalg.norm(xyz[i] - xyz[j]) for i, j in pairs])
    res["n_constraints"] = len(pairs); res["natoms"] = len(el)
    return res


def sum_rule(H):
    """Translational sum rule of a Cartesian Hessian (mdyn/A): for an exact Hessian sum_j H[i alpha, j beta] = 0.
    Returns the largest and the root-mean-square violation."""
    n = H.shape[0] // 3; Hs = 0.5 * (H + H.T); S = Hs.reshape(n, 3, n, 3).sum(axis=2)
    return float(np.abs(S).max()), float(np.sqrt((S ** 2).mean()))


if __name__ == "__main__":
    # python constrained_modes.py --runs ../revision/R3_force ../revision/R8_force_new_series [--data ../data]
    # analyses every FORCE job of the given folders whose outputs (<name>_force.out and .aux) are present and, with --data,
    # compares the result with data/vibrational_summary_v2.csv (the summary of all 226 structures of the analysis set).
    import argparse, os
    import pandas as pd
    ap = argparse.ArgumentParser(); ap.add_argument("--runs", nargs="+", required=True); ap.add_argument("--data", default=None); a = ap.parse_args()
    rows = []
    for runs in a.runs:
        for job in sorted(os.listdir(runs)):
            st = os.path.join(runs, job, job + "_force")
            if not (os.path.exists(st + ".aux") and os.path.exists(st + ".out")): continue
            r = analyse_job(st + ".mop", st + ".out", st + ".aux"); f = np.sort(r["freq"]); H = read_force_aux(st + ".aux")[2]; smax, srms = sum_rule(H)
            rows.append(dict(job=job, natoms=r["natoms"], n_constraints=r["n_constraints"], gnorm_tangent=r["gnorm_tangent"], n_modes=r["n_modes"], n_imag=int((f < 0).sum()),
                             n_imag_gt20=int((f < -20).sum()), f1=f[0], f2=f[1], sumrule_max=smax, sumrule_rms=srms))
    T = pd.DataFrame(rows); pd.set_option("display.width", 200)
    print(T.round(4).to_string(index=False))
    if a.data:
        ref = pd.read_csv(os.path.join(a.data, "vibrational_summary_v2.csv")).set_index("job").loc[T.job]
        cols = ["n_constraints", "gnorm_tangent", "n_imag", "n_imag_gt20", "f1", "f2", "sumrule_max"]
        dev = {c: float(np.abs(T[c].values.astype(float) - ref[c].values.astype(float)).max()) for c in cols}
        print(f"{len(T)} FORCE jobs analysed; largest deviation from data/vibrational_summary_v2.csv:", {c: round(v, 6) for c, v in dev.items()})
