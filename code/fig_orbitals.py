"""Electron density / frontier orbitals / ESP figure for an RM1 Eu-pectin structure.
Usage: python fig_orbitals.py <final_geoms_final.json>:<job> <charge> <out_prefix> [model]
   e.g. python fig_orbitals.py ../data/final_geoms_final.json:RM1_B3_h2_mono2 0 B3mono2
        python fig_orbitals.py ../data/final_geoms_final.json:RM1_CTRL_Eu_w9 3 aqua
Runs a MOPAC single point (1SCF GRAPHF AUX VECTORS ALLVEC MULLIK) in <out_prefix>_sp/, evaluates the STOs on a grid,
prints the HOMO/LUMO energies, the Eu weight of HOMO and LUMO and the Eu Mulliken charge, and renders the orbital and
density/ESP panels of Figure 11 with PyVista (<out_prefix>_homo.png, _lumo.png, _density_esp.png).
MOPAC executable: environment variable MOPAC_EXE (default 'mopac')."""
import sys, os, subprocess, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pyvista as pv
pv.OFF_SCREEN = True
from mopio import read_aux
from orbitals import read_basis, read_mgf, make_grid, mo_on_grid, density_on_grid, coulomb_esp, BOHR
from render import classify_atoms, add_molecule, frame_view

M = os.environ.get("MOPAC_EXE", "mopac")
ENV = dict(os.environ)
if os.path.isabs(M):   # Linux release: libmopac.so lives in ../lib next to bin/mopac
    ENV["LD_LIBRARY_PATH"] = os.path.join(os.path.dirname(os.path.dirname(M)), "lib") + os.pathsep + ENV.get("LD_LIBRARY_PATH", "")


def single_point(el, xyz, charge, workdir, model="RM1", extra=""):
    os.makedirs(workdir, exist_ok=True)
    fn = os.path.join(workdir, "sp.mop")
    with open(fn, "w") as f:
        f.write(f"{model} CHARGE={charge} 1SCF GRAPHF AUX VECTORS ALLVEC MULLIK BONDS {extra}\nsingle point for orbitals\n\n")
        for e, p in zip(el, xyz):
            f.write(f"{e} {p[0]:.6f} 0 {p[1]:.6f} 0 {p[2]:.6f} 0\n")
    subprocess.run([M, "sp.mop"], cwd=workdir, env=ENV, capture_output=True)
    return os.path.join(workdir, "sp.aux"), os.path.join(workdir, "sp.mgf")


def grid_to_image(vals, shape, lo, sp):
    img = pv.ImageData(dimensions=shape, spacing=(sp, sp, sp), origin=tuple(lo))
    img.point_data["v"] = np.asarray(vals, dtype=np.float32).ravel(order="F")
    return img


def render_iso(el, xyz, img, isos, colors, out, title=None, opacity=0.55, scalars_map=None, cmap=None, clim=None, cbar_title=None, view=None):
    cls = classify_atoms(el, xyz)
    p = pv.Plotter(off_screen=True, window_size=(1600, 1300), lighting="three lights")
    p.set_background("white")
    add_molecule(p, el, xyz, cls, water_alpha=0.25, show_H=True, polyhedron=False, eu_bonds=True)
    for iso, col in zip(isos, colors):
        c = img.contour([iso], scalars="v")
        if c.n_points == 0:
            continue
        if scalars_map is not None:
            c = c.sample(scalars_map)
            p.add_mesh(c, scalars="esp", cmap=cmap, clim=clim, opacity=opacity, smooth_shading=True, show_scalar_bar=False,
                       scalar_bar_args=dict(title=cbar_title or "", color="#0b0b0b", title_font_size=18, label_font_size=14, vertical=False, position_x=0.3, position_y=0.04, width=0.4, height=0.06))
        else:
            p.add_mesh(c, color=col, opacity=opacity, smooth_shading=True, specular=0.2)
    p.enable_anti_aliasing("ssaa")
    if view is None:
        view = frame_view(el, xyz, cls, mode="side", margin=1.5)
    p.camera.position = tuple(view["position"]); p.camera.focal_point = tuple(view["focal"]); p.camera.up = tuple(view["up"])
    p.enable_parallel_projection(); p.camera.parallel_scale = view["scale"]
    if title:
        p.add_text(title, position="upper_left", font_size=13, color="#0b0b0b")
    p.screenshot(out)
    p.close()


def main(aux_in, charge, prefix, model="RM1", eu_focus=True):
    if ".json:" in aux_in:
        import json as _json
        fn, job = aux_in.rsplit(":", 1)
        g = _json.load(open(fn))[job]
        el, xyz = list(g["el"]), np.asarray(g["xyz"], float)
    else:
        a = read_aux(aux_in)
        el, xyz = a["el"], a["xyz"]
    work = prefix + "_sp"
    aux, mgf = single_point(el, xyz, charge, work, model=model)
    b = read_basis(aux); m = read_mgf(mgf)
    nocc = int(round(m["occ"].sum() / 2))
    homo, lumo = nocc - 1, nocc
    print("HOMO/LUMO eV:", m["eig"][homo], m["eig"][lumo])
    # grid around the Eu region (within R of Eu); only AOs of atoms within R+6 A contribute (STO decay)
    eu = el.index("Eu")
    R = 7.5
    sel = [i for i in range(len(el)) if np.linalg.norm(xyz[i] - xyz[eu]) < R]
    sub = xyz[sel]
    pts, shape, lo, sp = make_grid(sub, margin=2.0, spacing=0.30)
    print("grid", shape, len(pts))
    near_atoms = {i for i in range(len(el)) if np.linalg.norm(xyz[i] - xyz[eu]) < R + 6.0}
    ao_idx = np.array([k for k in range(len(b["ao_atom"])) if b["ao_atom"][k] in near_atoms])
    bsub = dict(b); 
    for key in ("ao_atom", "ao_sym", "ao_zeta", "ao_pqn"):
        bsub[key] = [b[key][k] for k in ao_idx] if key == "ao_sym" else np.asarray(b[key])[ao_idx]
    from orbitals import ao_values
    A = ao_values(bsub, pts)                      # (npts, nao_sub) float32
    Minv = m["Minv"]
    C_all = (Minv @ m["C"].T).T                   # deorthogonalised coefficients (nmo, nao)
    Csub = C_all[:, ao_idx].astype(np.float32)
    psi_h = (A @ Csub[homo]).reshape(shape, order="C")
    psi_l = (A @ Csub[lumo]).reshape(shape, order="C")
    occ_idx = np.where(m["occ"] > 0)[0]
    rho = np.zeros(len(pts), dtype=np.float32)
    for s0 in range(0, len(occ_idx), 40):
        blk = Csub[occ_idx[s0:s0 + 40]]
        psi = A @ blk.T
        rho += (m["occ"][occ_idx[s0:s0 + 40]].astype(np.float32) * psi ** 2).sum(1)
    rho = rho.reshape(shape, order="C")
    esp = coulomb_esp(b["charges"], b["xyz"], pts).reshape(shape, order="C") * 627.5  # kcal/mol per unit charge
    img_h = grid_to_image(psi_h, shape, lo, sp); img_l = grid_to_image(psi_l, shape, lo, sp)
    img_r = grid_to_image(rho, shape, lo, sp)
    img_e = grid_to_image(esp, shape, lo, sp); img_e.point_data["esp"] = img_e.point_data["v"]
    cls = classify_atoms(el, xyz)
    try:
        view = frame_view(el, np.array(sub), classify_atoms([el[i] for i in sel], sub) if eu_focus else cls, mode="eu", margin=1.0)
    except Exception:
        c = xyz[eu]; rad = float(np.max(np.linalg.norm(sub - c, axis=1))) + 1.2
        view = dict(position=(c + np.array([0.0, 0.0, 3.0 * rad])).tolist(), focal=c.tolist(), up=(0.0, 1.0, 0.0), scale=rad)
    render_iso(el, xyz, img_h, [0.03, -0.03], ["#2a78d6", "#e34948"], prefix + "_homo.png", title=None, view=view)
    render_iso(el, xyz, img_l, [0.03, -0.03], ["#2a78d6", "#e34948"], prefix + "_lumo.png", title=None, view=view)
    render_iso(el, xyz, img_r, [0.02], ["#9ec5f4"], prefix + "_density_esp.png", title=None,
               opacity=0.85, scalars_map=img_e, cmap="coolwarm", clim=(-150, 350), cbar_title="ESP (kcal/mol)", view=view)
    # Eu atomic charge and orbital composition on Eu for HOMO/LUMO
    eu_aos = [k for k in range(len(b["ao_atom"])) if b["ao_atom"][k] == eu]
    for lab, idx in (("HOMO", homo), ("LUMO", lumo)):
        c = m["C"][idx]; w = float(np.sum(c[eu_aos] ** 2))
        print(f"{lab}: Eu weight {w:.3f}")
    print("Eu Mulliken charge:", b["charges"][eu])
    return m["eig"][homo], m["eig"][lumo], b["charges"][eu]


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else "RM1")
