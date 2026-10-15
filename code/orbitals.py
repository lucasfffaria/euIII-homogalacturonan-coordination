"""Evaluate MOPAC molecular orbitals / electron density on a grid from the .aux (basis) and .mgf (coefficients).

Slater-type orbitals: R_n(r) = N r^(n-1) exp(-zeta r), N = (2 zeta)^(n+1/2) / sqrt((2n)!)
Real angular parts for s, p, d in MOPAC order: S PX PY PZ X2 XZ Z2 YZ XY.
"""
import re, math
import numpy as np
from mopio import _aux_records, _fortran_float

BOHR = 0.529177210903


def read_basis(aux_file):
    recs = _aux_records(open(aux_file, errors="ignore").read())
    d = {k: v for k, h, v in recs}
    el = d["ATOM_EL"]
    xyz = np.array(d["ATOM_X_OPT"] if "ATOM_X_OPT" in d else d["ATOM_X"], float).reshape(-1, 3)
    aoi = np.array(d["AO_ATOMINDEX"], int) - 1
    sym = d["ATOM_SYMTYPE"]
    zeta = np.array(d["AO_ZETA"], float)
    pqn = np.array(d["ATOM_PQN"], int)
    out = dict(el=el, xyz=xyz, ao_atom=aoi, ao_sym=sym, ao_zeta=zeta, ao_pqn=pqn)
    if "EIGENVECTORS" in d:
        nao = len(aoi)
        v = np.array([_fortran_float(x) for x in d["EIGENVECTORS"]])
        out["eigvec_aux"] = v.reshape(-1, nao)
    if "EIGENVALUES" in d:
        out["eigval_aux"] = np.array([_fortran_float(x) for x in d["EIGENVALUES"]])
    if "MOLECULAR_ORBITAL_OCCUPANCIES" in d:
        out["occ_aux"] = np.array(d["MOLECULAR_ORBITAL_OCCUPANCIES"], float)
    if "ATOM_CHARGES" in d:
        out["charges"] = np.array(d["ATOM_CHARGES"], float)
    if "TOTAL_DENSITY_MATRIX" in d:
        nao = len(aoi)
        tri = np.array([_fortran_float(x) for x in d["TOTAL_DENSITY_MATRIX"]])
        P = np.zeros((nao, nao))
        k = 0
        for i in range(nao):
            for j in range(i + 1):
                P[i, j] = P[j, i] = tri[k]; k += 1
        out["P"] = P
    if "OVERLAP_MATRIX" in d:
        nao = len(aoi)
        tri = np.array([_fortran_float(x) for x in d["OVERLAP_MATRIX"]])
        S = np.zeros((nao, nao))
        k = 0
        for i in range(nao):
            for j in range(i + 1):
                S[i, j] = S[j, i] = tri[k]; k += 1
        out["S"] = S
    return out


def read_mgf(mgf_file):
    """All MO coefficients (orthogonal Löwdin basis), eigenvalues, occupations and the inverse matrix."""
    lines = open(mgf_file, errors="ignore").read().splitlines()
    nat = int(lines[0].split()[0])
    orbs = []
    i = 1 + 2 * nat
    coefs = None
    while i < len(lines):
        l = lines[i]
        if l.strip().startswith("ORBITAL"):
            t = l.split()
            occ = int(t[1]); ev = float(t[-1])
            i += 1
            vals = []
            while i < len(lines) and not lines[i].strip().startswith("ORBITAL") and not lines[i].strip().startswith("INVERSE") and not lines[i].strip().startswith("Keywords"):
                s = lines[i]
                # fixed-width Fortran D format, 5 per line, width 15
                for k in range(0, len(s.rstrip()), 15):
                    tok = s[k:k + 15].strip()
                    if tok:
                        vals.append(float(tok.replace("D", "E")))
                i += 1
            orbs.append((occ, ev, np.array(vals)))
        elif l.strip().startswith("INVERSE_MATRIX"):
            n = int(re.search(r"\[(\d+)x", l).group(1))
            i += 1
            vals = []
            while i < len(lines) and not lines[i].strip().startswith("Keywords"):
                s = lines[i]
                for k in range(0, len(s.rstrip()), 15):
                    tok = s[k:k + 15].strip()
                    if tok:
                        vals.append(float(tok.replace("D", "E")))
                i += 1
            M = np.zeros((n, n)); k = 0
            for a in range(n):
                for b in range(a + 1):
                    M[a, b] = M[b, a] = vals[k]; k += 1
            coefs = M
        else:
            i += 1
    C = np.array([o[2] for o in orbs])
    return dict(occ=np.array([o[0] for o in orbs]), eig=np.array([o[1] for o in orbs]), C=C, Minv=coefs)


def _norm(n, zeta):
    return (2 * zeta) ** (n + 0.5) / math.sqrt(math.factorial(2 * n))


def ao_values(basis, pts, chunk=200000):
    """Return array (npts, nao) of STO values at pts (Angstrom)."""
    xyz = basis["xyz"] / BOHR  # bohr
    P = np.asarray(pts) / BOHR
    nao = len(basis["ao_atom"])
    out = np.zeros((len(P), nao), dtype=np.float32)
    for s in range(0, len(P), chunk):
        p = P[s:s + chunk]
        cache = {}
        for mu in range(nao):
            a = basis["ao_atom"][mu]
            if a not in cache:
                dvec = p - xyz[a]
                r = np.linalg.norm(dvec, axis=1) + 1e-12
                cache[a] = (dvec, r)
            dvec, r = cache[a]
            n = int(basis["ao_pqn"][mu]); z = float(basis["ao_zeta"][mu]); sym = basis["ao_sym"][mu]
            rad = _norm(n, z) * r ** (n - 1) * np.exp(-z * r)
            x, y, zz = dvec[:, 0] / r, dvec[:, 1] / r, dvec[:, 2] / r
            if sym == "S":
                ang = np.full_like(r, math.sqrt(1 / (4 * math.pi)))
            elif sym == "PX":
                ang = math.sqrt(3 / (4 * math.pi)) * x
            elif sym == "PY":
                ang = math.sqrt(3 / (4 * math.pi)) * y
            elif sym == "PZ":
                ang = math.sqrt(3 / (4 * math.pi)) * zz
            elif sym == "X2":   # x2-y2
                ang = math.sqrt(15 / (16 * math.pi)) * (x * x - y * y)
            elif sym == "XZ":
                ang = math.sqrt(15 / (4 * math.pi)) * x * zz
            elif sym == "Z2":
                ang = math.sqrt(5 / (16 * math.pi)) * (3 * zz * zz - 1)
            elif sym == "YZ":
                ang = math.sqrt(15 / (4 * math.pi)) * y * zz
            elif sym == "XY":
                ang = math.sqrt(15 / (4 * math.pi)) * x * y
            else:
                raise ValueError(sym)
            out[s:s + chunk, mu] = (rad * ang).astype(np.float32)
    return out


def make_grid(xyz, margin=3.5, spacing=0.25):
    lo = xyz.min(0) - margin; hi = xyz.max(0) + margin
    nx, ny, nz = [int(np.ceil((hi[i] - lo[i]) / spacing)) + 1 for i in range(3)]
    gx = lo[0] + spacing * np.arange(nx); gy = lo[1] + spacing * np.arange(ny); gz = lo[2] + spacing * np.arange(nz)
    X, Y, Z = np.meshgrid(gx, gy, gz, indexing="ij")
    pts = np.c_[X.ravel(), Y.ravel(), Z.ravel()]
    return pts, (nx, ny, nz), lo, spacing


def mo_on_grid(basis, coef, pts, deortho=None):
    c = coef if deortho is None else deortho @ coef
    A = ao_values(basis, pts)
    return A @ c.astype(np.float32)


def density_on_grid(basis, P, pts, deortho=None, chunk=100000):
    """rho = sum_mu,nu P_mu,nu chi_mu chi_nu  (P in the orthogonal basis; transformed with deortho if given)."""
    Pa = P if deortho is None else deortho @ P @ deortho.T
    rho = np.zeros(len(pts), dtype=np.float32)
    for s in range(0, len(pts), chunk):
        A = ao_values(basis, pts[s:s + chunk])
        rho[s:s + chunk] = np.einsum("pi,ij,pj->p", A, Pa.astype(np.float32), A)
    return rho


def coulomb_esp(charges, xyz, pts, chunk=200000):
    """Point-charge electrostatic potential (kcal/mol per unit charge... returned in atomic units e/bohr)."""
    out = np.zeros(len(pts))
    X = xyz / BOHR
    P = np.asarray(pts) / BOHR
    for s in range(0, len(P), chunk):
        p = P[s:s + chunk]
        v = np.zeros(len(p))
        for q, x in zip(charges, X):
            r = np.linalg.norm(p - x, axis=1) + 1e-6
            v += q / r
        out[s:s + chunk] = v
    return out
