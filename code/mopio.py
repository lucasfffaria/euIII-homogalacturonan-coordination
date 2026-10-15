"""MOPAC input writer with mixed Cartesian / internal coordinates:
all O-H bonds (water and hydroxyl) are written as internal coordinates with the
bond length frozen (flag 0) so that spurious proton transfer to the Eu sparkle cannot occur.
Also: .arc / .aux parsers."""
import re
import numpy as np

OH_BOND = 0.97


def bonded_pairs(el, xyz, scale=1.2):
    cov = {"C": 0.76, "O": 0.66, "H": 0.31, "Eu": 1.0, "N": 0.71}
    n = len(el)
    nb = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if el[i] == "Eu" or el[j] == "Eu":
                continue
            d = np.linalg.norm(xyz[i] - xyz[j])
            if d < scale * (cov[el[i]] + cov[el[j]]):
                nb[i].append(j); nb[j].append(i)
    return nb


def dihedral(p0, p1, p2, p3):
    b0 = p0 - p1; b1 = p2 - p1; b2 = p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    x = np.dot(v, w)
    y = np.dot(np.cross(b1n, v), w)
    return np.degrees(np.arctan2(y, x))


def angle(p0, p1, p2):
    a = p0 - p1; b = p2 - p1
    return np.degrees(np.arccos(np.clip(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b), -1, 1)))


def write_mop_mixed(fname, el, xyz, keywords, title, freeze_OH=True):
    """el: list of symbols; xyz: (n,3). Atoms must be ordered so that every H comes after its O
    and after the two reference atoms used for its angle/dihedral."""
    el0 = list(el); xyz0 = np.asarray(xyz, float)
    # reorder: all non-H atoms first (original order), then all H atoms (original order)
    perm = [i for i, e in enumerate(el0) if e != "H"] + [i for i, e in enumerate(el0) if e == "H"]
    el = [el0[i] for i in perm]; xyz = xyz0[perm]; n = len(el)
    nb = bonded_pairs(el, xyz)
    ieu = el.index("Eu") if "Eu" in el else None
    lines = []
    for i in range(n):
        internal = False
        if freeze_OH and el[i] == "H" and nb[i] and el[nb[i][0]] == "O":
            o = nb[i][0]
            # reference b: a heavy neighbour of O other than H, else Eu, else nearest heavy atom
            cands = [k for k in nb[o] if el[k] != "H" and k < i]
            if cands:
                b = cands[0]
            elif ieu is not None and ieu < i:
                b = ieu
            else:
                heavy = [k for k in range(i) if el[k] not in ("H",) and k != o]
                b = min(heavy, key=lambda k: np.linalg.norm(xyz[k] - xyz[o]))
            # reference c: neighbour of b other than o (must precede i) else nearest heavy atom to b
            cc = [k for k in nb[b] if k != o and k < i and el[k] != "H"] if b != ieu else []
            if cc:
                c = cc[0]
            else:
                heavy = [k for k in range(i) if k not in (o, b) and el[k] != "H"]
                c = min(heavy, key=lambda k: np.linalg.norm(xyz[k] - xyz[b]))
            r = np.linalg.norm(xyz[i] - xyz[o])
            a = angle(xyz[i], xyz[o], xyz[b])
            d = dihedral(xyz[i], xyz[o], xyz[b], xyz[c])
            # guard against collinear references
            if a < 2 or a > 178 or angle(xyz[o], xyz[b], xyz[c]) < 2 or angle(xyz[o], xyz[b], xyz[c]) > 178:
                internal = False
            else:
                lines.append(f"{el[i]:2s} {r:10.5f} 0 {a:10.4f} 1 {d:10.4f} 1 {o+1:5d} {b+1:5d} {c+1:5d}")
                internal = True
        if not internal:
            x, y, z = xyz[i]
            lines.append(f"{el[i]:2s} {x:12.6f} 1 {y:12.6f} 1 {z:12.6f} 1     0     0     0")
    with open(fname, "w", newline="\n") as f:
        f.write(keywords + "\n" + title + "\n\n" + "\n".join(lines) + "\n")
    n_internal = sum(1 for l in lines if not l.rstrip().endswith("0     0     0"))
    return n_internal, perm


def read_arc_geometry(fn):
    """Final Cartesian geometry from a MOPAC .arc file (works for internal-coordinate runs too,
    because the arc stores the geometry in the input format... so prefer the .aux parser)."""
    txt = open(fn, errors="ignore").read()
    lines = txt.splitlines()
    i = [k for k, l in enumerate(lines) if "FINAL GEOMETRY OBTAINED" in l][-1]
    atoms = []
    for l in lines[i + 1:]:
        t = l.split()
        if len(t) >= 7 and re.match(r"^[A-Z][a-z]?$", t[0]):
            atoms.append((t[0], float(t[1]), float(t[3]), float(t[5])))
        elif atoms:
            break
    return atoms


def _aux_records(txt):
    """Split a MOPAC .aux file into (key, header_line, data_tokens) records."""
    recs = []
    key = None; data = []
    for line in txt.splitlines():
        m = re.match(r"^\s*([A-Z_][A-Z_0-9]*)(\[|:|=)", line)
        if m:
            if key is not None:
                recs.append((key, hdr, data))
            key = m.group(1); hdr = line.strip(); data = []
            # inline value after '=' on the same line (scalars)
            if "=" in line and not line.strip().endswith("="):
                data = [line.split("=", 1)[1].strip()]
        elif key is not None and not line.lstrip().startswith("#"):
            data.extend(line.split())
    if key is not None:
        recs.append((key, hdr, data))
    return recs


def _fortran_float(s):
    return float(s.replace("D", "E"))


def read_aux(fn):
    """Parse a MOPAC .aux file: dict with el, xyz (final), charges, heat, gnorm, plus trajectories."""
    txt = open(fn, errors="ignore").read()
    recs = _aux_records(txt)
    out = {"xyz_traj": [], "heat_traj": [], "gnorm_traj": []}
    for key, hdr, data in recs:
        if key == "ATOM_EL":
            out["el"] = data
        elif key == "ATOM_X_UPDATED":
            out["xyz_traj"].append(np.array(data, float).reshape(-1, 3))
        elif key == "ATOM_X_OPT":
            out["xyz"] = np.array(data, float).reshape(-1, 3)
        elif key == "ATOM_X" and "xyz0" not in out:
            out["xyz0"] = np.array(data, float).reshape(-1, 3)
        elif key == "HEAT_OF_FORM_UPDATED":
            out["heat_traj"].append(_fortran_float(data[0]))
        elif key == "GRADIENT_NORM_UPDATED":
            out["gnorm_traj"].append(_fortran_float(data[0]))
        elif key == "HEAT_OF_FORMATION":
            out["heat"] = _fortran_float(data[0])
        elif key == "GRADIENT_NORM":
            out["gnorm"] = _fortran_float(data[0])
        elif key == "ATOM_CHARGES":
            out["charges"] = np.array(data, float)
        elif key == "DIPOLE":
            out["dipole"] = data
        elif key == "EIGENVALUES":
            out["eigenvalues"] = np.array([_fortran_float(x) for x in data])
        elif key == "MOLECULAR_ORBITAL_OCCUPANCIES":
            out["occ"] = np.array(data, float)
    if "xyz" not in out and out["xyz_traj"]:
        out["xyz"] = out["xyz_traj"][-1]
    return out


def read_out_final_xyz(fn):
    """Final Cartesian coordinates from the .out file (CARTESIAN COORDINATES block after optimisation)."""
    txt = open(fn, errors="ignore").read()
    blocks = re.findall(r"CARTESIAN COORDINATES\s*\n\s*\n?(?:\s*NO\.\s+ATOM\s+X\s+Y\s+Z\s*\n\s*\n)?(.*?)(?:\n\s*\n)", txt, re.S)
    if not blocks:
        return None
    atoms = []
    for l in blocks[-1].splitlines():
        t = l.split()
        if len(t) >= 5 and t[0].isdigit():
            atoms.append((t[1], float(t[2]), float(t[3]), float(t[4])))
    return atoms
