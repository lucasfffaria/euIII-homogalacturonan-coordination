"""PyVista renderings of Eu/pectin clusters (publication style)."""
import numpy as np
import pyvista as pv
from scipy.spatial import ConvexHull
pv.OFF_SCREEN = True

COL = {"C": "#5b5b5b", "O": "#e0281e", "H": "#f2f2f2", "Eu": "#f5b400"}
COL_WATER_O = "#4a90e2"     # water oxygens in blue to separate them from pectin oxygens
COL_METHYL = "#8e44ad"      # methyl-ester carbons (violet) to make the esterification pattern visible
RAD = {"C": 0.33, "O": 0.32, "H": 0.19, "Eu": 0.72}
COV = {"C": 0.76, "O": 0.66, "H": 0.31, "Eu": 1.98}


def bonds(el, xyz, scale=1.2):
    n = len(el); out = []
    for i in range(n):
        for j in range(i + 1, n):
            if el[i] == "Eu" or el[j] == "Eu":
                continue
            if np.linalg.norm(xyz[i] - xyz[j]) < scale * (COV[el[i]] + COV[el[j]]):
                out.append((i, j))
    return out


def classify_atoms(el, xyz):
    """Return dict: water_O indices, methyl C indices, first-sphere O set."""
    n = len(el)
    nb = [[] for _ in range(n)]
    for i, j in bonds(el, xyz):
        nb[i].append(j); nb[j].append(i)
    water_atoms = set()
    for i in range(n):
        if el[i] == "O" and all(el[k] == "H" for k in nb[i]):
            water_atoms.add(i); water_atoms.update(nb[i])
    methyl = set()
    for i in range(n):
        if el[i] == "C" and sum(1 for k in nb[i] if el[k] == "H") == 3 and any(el[k] == "O" for k in nb[i]):
            methyl.add(i); methyl.update(k for k in nb[i] if el[k] == "H")
    eu = el.index("Eu") if "Eu" in el else None
    first = set()
    if eu is not None:
        for i in range(n):
            if el[i] == "O" and np.linalg.norm(xyz[i] - xyz[eu]) < 3.0:
                first.add(i)
    return dict(water=water_atoms, methyl=methyl, first=first, nb=nb, eu=eu)


def add_molecule(p, el, xyz, cls, water_alpha=0.35, show_H=True, polyhedron=True, eu_bonds=True, chain_alpha=1.0):
    n = len(el)
    eu = cls["eu"]
    first_water = {i for i in cls["first"] if i in cls["water"]}
    # atoms
    for i in range(n):
        e = el[i]
        if e == "H" and not show_H:
            continue
        in_water = i in cls["water"]
        # first-shell waters (O + their H) are opaque, other waters translucent
        opaque_water = in_water and (i in first_water or any(k in first_water for k in cls["nb"][i]))
        alpha = (1.0 if opaque_water else water_alpha) if in_water else chain_alpha
        color = COL[e]
        if e == "O" and in_water:
            color = COL_WATER_O
        if i in cls["methyl"] and e == "C":
            color = COL_METHYL
        if i in cls["methyl"] and e == "H":
            color = "#d7bde2"
        r = RAD[e]
        p.add_mesh(pv.Sphere(radius=r, center=xyz[i], theta_resolution=28, phi_resolution=28), color=color, smooth_shading=True,
                   opacity=alpha, specular=0.35, specular_power=20)
    # bonds
    for i, j in bonds(el, xyz):
        if not show_H and ("H" in (el[i], el[j])):
            continue
        in_water = i in cls["water"] and j in cls["water"]
        opaque_water = in_water and (i in first_water or j in first_water)
        alpha = (1.0 if opaque_water else water_alpha) if in_water else chain_alpha
        mid = (xyz[i] + xyz[j]) / 2
        for a, b, c in ((i, j, el[i]), (j, i, el[j])):
            col = COL[c]
            if c == "O" and a in cls["water"]:
                col = COL_WATER_O
            if a in cls["methyl"] and c == "C":
                col = COL_METHYL
            p.add_mesh(pv.Tube(pointa=xyz[a], pointb=mid, radius=0.085, n_sides=14), color=col, opacity=alpha, smooth_shading=True)
    if eu is not None and eu_bonds:
        for i in cls["first"]:
            p.add_mesh(pv.Tube(pointa=xyz[eu], pointb=xyz[i], radius=0.045, n_sides=10), color="#c69500", opacity=0.9)
    if eu is not None and polyhedron and len(cls["first"]) >= 4:
        pts = np.array([xyz[i] for i in cls["first"]])
        try:
            hull = ConvexHull(pts)
            faces = np.hstack([[3, *s] for s in hull.simplices])
            poly = pv.PolyData(pts, faces)
            p.add_mesh(poly, color="#f5b400", opacity=0.22, show_edges=True, edge_color="#c69500", line_width=1.2, smooth_shading=False)
        except Exception:
            pass


def frame_view(el, xyz, cls, mode="side", margin=2.5):
    """Camera setup: parallel projection framing the whole cluster.
    mode 'side': look perpendicular to the chain axis and to the Eu-chain vector (chain vertical, Eu to one side);
    mode 'eu'  : look from the Eu side toward the chain (coordination sphere in front)."""
    xyz = np.asarray(xyz)
    heavy = [i for i in range(len(el)) if el[i] != "H"]
    P = xyz[heavy]
    center = P.mean(0)
    chain = np.array([xyz[i] for i in heavy if i not in cls["water"] and el[i] != "Eu"])
    c = chain.mean(0)
    u, s, vt = np.linalg.svd(chain - c)
    axis = vt[0]
    if cls["eu"] is not None:
        v = xyz[cls["eu"]] - c
        v = v - np.dot(v, axis) * axis
        if np.linalg.norm(v) < 1e-3:
            v = vt[1]
        v = v / np.linalg.norm(v)
    else:
        v = vt[1]
    if mode == "side":
        d = np.cross(axis, v); d /= np.linalg.norm(d)
        up = axis
    else:
        d = v
        up = axis
    # extent along the two in-plane directions
    rel = P - center
    e_up = np.abs(rel @ up).max(); e_x = np.abs(rel @ np.cross(d, up)).max()
    scale = max(e_up, e_x * 1.0) + margin
    return dict(position=center + d * 60.0, focal=center, up=up, scale=scale)


def render(el, xyz, out, size=(1600, 1300), view=None, zoom=1.0, mode="side", show_H=True, water_alpha=0.35, title=None, bg="white", polyhedron=True, parallel=True):
    el = list(el); xyz = np.asarray(xyz, float)
    cls = classify_atoms(el, xyz)
    p = pv.Plotter(off_screen=True, window_size=size, lighting="three lights")
    p.set_background(bg)
    add_molecule(p, el, xyz, cls, water_alpha=water_alpha, show_H=show_H, polyhedron=polyhedron)
    p.enable_anti_aliasing("ssaa")
    if view is None:
        view = frame_view(el, xyz, cls, mode=mode)
    p.camera.position = tuple(view["position"]); p.camera.focal_point = tuple(view["focal"]); p.camera.up = tuple(view["up"])
    if parallel:
        p.enable_parallel_projection()
        p.camera.parallel_scale = view["scale"] / zoom
    else:
        p.camera.zoom(zoom)
    if title:
        p.add_text(title, position="upper_left", font_size=12, color="#0b0b0b")
    p.screenshot(out, transparent_background=False)
    p.close()
    return cls


def best_view(el, xyz, cls):
    """Camera looking at Eu from the side opposite to the chain (so the coordination sphere is in front)."""
    xyz = np.asarray(xyz)
    eu = xyz[cls["eu"]]
    chain = np.array([xyz[i] for i in range(len(el)) if i not in cls["water"] and el[i] != "Eu"])
    c = chain.mean(0)
    v = eu - c; v = v / (np.linalg.norm(v) + 1e-9)
    # choose an up vector along the chain's principal axis
    u, s, vt = np.linalg.svd(chain - c)
    axis = vt[0]
    axis = axis - np.dot(axis, v) * v; axis /= np.linalg.norm(axis)
    pos = eu + 22.0 * v
    return [tuple(pos), tuple(eu), tuple(axis)]
