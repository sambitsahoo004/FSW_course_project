#!/usr/bin/env python3
"""
Post-processing for the FSW CFD model (ANSYS Fluent 2023 R2, AA6061).

Reads the Fluent HDF5 case/data files (*.cas.h5 / *.dat.h5) and the report-file
monitors (*-rfile.out) directly - no ANSYS licence needed - and regenerates every
figure used in the README and report, plus a short summary of key numbers.

Usage (from the repository root):
    pip install -r postprocessing/requirements.txt
    python postprocessing/fsw_postprocess.py
    python postprocessing/fsw_postprocess.py --case X.cas.h5 --data X.dat.h5 --out images
"""
import argparse, os
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.colors as mc
from matplotlib.patches import Circle, Rectangle, Arc
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 11,
                     "savefig.dpi": 200})
CLEAN = {"axes.spines.top": False, "axes.spines.right": False}

OMEGA = 74.351      # tool rotation, rad/s (710 rpm)
V_WELD = 0.00133    # welding speed, m/s (79.8 mm/min)


# ----------------------------------------------------------------------------- I/O
def load(case_path, data_path):
    """Return a dict with mesh geometry, cell/face fields and zone ranges."""
    f = h5py.File(case_path, "r")
    m = f["meshes/1"]
    X = m["nodes/coords/4"][:]                                   # node coords, m
    FN = m["faces/nodes/1/nodes"][:].reshape(-1, 3).astype(np.int64) - 1  # tri faces
    c0 = m["faces/c0/1"][:].astype(np.int64) - 1
    c1ds = m["faces/c1/1"]
    c1 = c1ds[:].astype(np.int64) - 1
    c1_lo, c1_hi = int(c1ds.attrs["minId"][0]) - 1, int(c1ds.attrs["maxId"][0])
    ncell = int(m.attrs["cellCount"][0])

    # zone ranges (1-based face ids, inclusive)
    zt = m["faces/zoneTopology"]
    names = zt["name"][0].decode().split(";")
    zones = {names[i]: (int(zt["minId"][i]), int(zt["maxId"][i])) for i in range(len(zt["id"]))}

    # cell centroids: for tetrahedra, mean of the 4 face centroids == cell centroid
    fc = X[FN].mean(1)
    C = np.zeros((ncell, 3)); cnt = np.zeros(ncell)
    np.add.at(C, c0, fc); np.add.at(cnt, c0, 1)
    np.add.at(C, c1, fc[c1_lo:c1_hi]); np.add.at(cnt, c1, 1)
    C /= cnt[:, None]
    v1 = X[FN[:, 1]] - X[FN[:, 0]]; v2 = X[FN[:, 2]] - X[FN[:, 0]]
    A = 0.5 * np.linalg.norm(np.cross(v1, v2), axis=1)

    d = h5py.File(data_path, "r")
    pc = d["results/1/phase-1/cells"]; pf = d["results/1/phase-1/faces"]
    cells = {k: pc[f"SV_{k}/1"][:] for k in ["T", "U", "V", "W", "MU_LAM", "P"]}

    def facefield(name):
        out = np.full(len(FN), np.nan)
        for s in pf[name]:
            ds = pf[name][s]
            if ds.ndim == 1:
                out[int(ds.attrs["minId"][0]) - 1:int(ds.attrs["maxId"][0])] = ds[:]
        return out

    faces = {k: facefield(f"SV_{k}") for k in ["T", "U", "V", "W", "HEAT_FLUX"]}
    res = {k: (d[f"results/residuals/phase-1/{k}/iterations"][:],
               d[f"results/residuals/phase-1/{k}/data"][:])
           for k in ["continuity", "x-velocity", "y-velocity", "z-velocity", "energy"]}
    return dict(X=X, FN=FN, fc=fc, C=C, A=A, zones=zones, cells=cells, faces=faces,
                residuals=res, ncell=ncell)


def zslice(M, name):
    a, b = M["zones"][name]
    return np.arange(a - 1, b)


# -------------------------------------------------------------------- figures
def fig_domain(out):
    with plt.rc_context(CLEAN):
        fig, axs = plt.subplots(1, 2, figsize=(11, 4.4), gridspec_kw={"width_ratios": [1.25, 1]})
        ax = axs[0]; R, yc = 120, 96; xa = np.sqrt(R**2 - yc**2)
        th = np.linspace(-np.arcsin(yc / R), np.arcsin(yc / R), 100)
        ax.plot(R*np.cos(th), R*np.sin(th), color="#c0392b", lw=3, label="Pressure outlet")
        ax.plot(-R*np.cos(th), R*np.sin(th), color="#2471a3", lw=3, label="Velocity inlet (1.33 mm/s)")
        ax.plot([-xa, xa], [yc, yc], color="#7f8c8d", lw=3, label="Side walls (h = 30 W/m²K)")
        ax.plot([-xa, xa], [-yc, -yc], color="#7f8c8d", lw=3)
        ax.add_patch(Circle((0, 0), 12, fc="#f5b041", ec="k", lw=1, label="Tool shoulder, Ø24 mm"))
        ax.add_patch(Circle((0, 0), 3, fc="#6e2c00", ec="k", lw=1, label="Pin, Ø6 mm"))
        ax.annotate("", xy=(62, 0), xytext=(22, 0), arrowprops=dict(arrowstyle="-|>", lw=1.8, color="#2471a3"))
        ax.text(42, -6, "material flow\n(tool frame)", ha="center", va="top", fontsize=8.5, color="#2471a3")
        ax.annotate("", xy=(-62, 0), xytext=(-22, 0), arrowprops=dict(arrowstyle="-|>", lw=1.8, color="k"))
        ax.text(-42, -6, "welding\ndirection", ha="center", va="top", fontsize=8.5)
        ax.add_patch(Arc((0, 0), 36, 36, theta1=20, theta2=160, lw=1.5, color="#1e8449"))
        ax.annotate("", xy=(-17.1, 6.2), xytext=(-16.4, 7.9), arrowprops=dict(arrowstyle="-|>", color="#1e8449", lw=1.5))
        ax.text(0, 24, "ω = 710 rpm (CCW)", ha="center", fontsize=8.5, color="#1e8449")
        ax.text(0, 60, "Advancing side (+y)", ha="center", fontsize=9, style="italic")
        ax.text(0, -66, "Retreating side (−y)", ha="center", fontsize=9, style="italic")
        ax.set_aspect("equal"); ax.set_xlim(-130, 130); ax.set_ylim(-110, 110)
        ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)"); ax.set_title("(a) Top view of the fluid domain (z = 0)")
        ax.legend(loc="lower center", fontsize=7.5, ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.42))
        ax = axs[1]
        ax.add_patch(Rectangle((-40, -6), 80, 6, fc="#d6eaf8", ec="k"))
        ax.add_patch(Rectangle((-12, 0), 24, 7, fc="#f5b041", ec="k"))
        ax.add_patch(Rectangle((-3, -5), 6, 5, fc="#6e2c00", ec="k"))
        ax.plot([-40, 40], [-6, -6], color="#922b21", lw=4)
        ax.text(0, -8.3, "Bottom face: backing-plate contact, h = 10 000 W/m²K", ha="center", fontsize=8.5, color="#922b21")
        ax.text(-26, 1.0, "Top face: h = 30 W/m²K", ha="center", fontsize=8.5)
        ax.text(0, 7.8, "Tool (adiabatic, no-slip, rotating)", ha="center", fontsize=8.5)
        ax.text(36.5, -3.3, "6 mm", fontsize=8.5, ha="right")
        ax.annotate("", xy=(38, 0), xytext=(38, -6), arrowprops=dict(arrowstyle="<->", lw=1))
        ax.annotate("", xy=(5, 0), xytext=(5, -5), arrowprops=dict(arrowstyle="<->", lw=1))
        ax.text(6, -2.7, "pin 5 mm", fontsize=8.5)
        ax.text(-24, -3.3, "AA6061 plate", fontsize=8.5, ha="center")
        ax.set_xlim(-42, 42); ax.set_ylim(-10, 10); ax.set_aspect(2.2)
        ax.set_xlabel("x (mm)"); ax.set_ylabel("z (mm)"); ax.set_title("(b) Section y = 0 (vertical scale ×2.2)")
        plt.tight_layout(); plt.savefig(os.path.join(out, "fig01_domain.png"), bbox_inches="tight"); plt.close()


def top_triangles(M):
    """Top-surface triangles = 'topface' zone + flat shoulder faces of the tool."""
    Xmm = M["X"] * 1000
    tool = zslice(M, "fswtool")
    shoulder = tool[np.all(np.abs(Xmm[M["FN"][tool]][:, :, 2]) < 1e-6, axis=1)]
    return np.concatenate([zslice(M, "topface"), shoulder])


def fig_mesh(M, out):
    Xmm = M["X"] * 1000
    tri = mtri.Triangulation(Xmm[:, 0], Xmm[:, 1], M["FN"][top_triangles(M)])
    with plt.rc_context(CLEAN):
        fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
        for ax, lim, t in zip(axs, [125, 20], ["(a) Entire top surface", "(b) Close-up around the tool shoulder"]):
            ax.triplot(tri, lw=0.25 if lim > 100 else 0.35, color="#1a5276")
            ax.add_patch(Circle((0, 0), 12, fill=False, ec="#e67e22", lw=1.4))
            yl = lim * 0.82 if lim > 100 else lim
            ax.set_aspect("equal"); ax.set_xlim(-lim, lim); ax.set_ylim(-yl, yl)
            ax.set_title(t); ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
        plt.tight_layout(); plt.savefig(os.path.join(out, "fig02_mesh.png"), bbox_inches="tight"); plt.close()


def fig_convergence(M, mon_dir, out):
    cols = {"continuity": "#1f77b4", "x-velocity": "#d62728", "y-velocity": "#2ca02c",
            "z-velocity": "#9467bd", "energy": "#ff7f0e"}
    with plt.rc_context(CLEAN):
        fig, axs = plt.subplots(1, 2, figsize=(11, 3.9))
        ax = axs[0]
        for k, c in cols.items():
            it, D = M["residuals"][k]
            ax.semilogy(it, D[:, 0] / D[:, 1], label=k, color=c, lw=1.4, ls="--" if k == "x-velocity" else "-")
        ax.axhline(1e-6, ls="--", color="k", lw=0.8); ax.text(5, 1.4e-6, "energy criterion 10⁻⁶", fontsize=8)
        ax.axhline(1e-3, ls=":", color="k", lw=0.8); ax.text(5, 1.4e-3, "flow criterion 10⁻³", fontsize=8)
        ax.text(250, 6e-6, "x- and y-velocity\noverlap", fontsize=7.5, color="#555")
        ax.set_xlabel("Iteration"); ax.set_ylabel("Scaled residual"); ax.set_title("(a) Scaled residuals")
        ax.legend(fontsize=8, frameon=False, ncol=2); ax.set_xlim(0, 360)
        ax = axs[1]
        for n, l, c in [("weldforce", "Traverse force $F_x$ (N)", "#1f77b4"),
                        ("lateralforce", "Lateral force $F_y$ (N)", "#2ca02c"),
                        ("thrustforce", "Axial force $F_z$ (N)", "#9467bd")]:
            a = np.loadtxt(os.path.join(mon_dir, f"{n}-rfile.out"), skiprows=3)
            ax.plot(a[:, 0], a[:, 1], label=l, color=c, lw=1.4)
        ax.set_ylim(-20, 30); ax.set_xlabel("Iteration"); ax.set_ylabel("Force (N)")
        ax2 = ax.twinx(); a = np.loadtxt(os.path.join(mon_dir, "torque-rfile.out"), skiprows=3)
        ax2.plot(a[:, 0], -a[:, 1], color="#d62728", lw=1.8, label="Torque $|M_z|$ (N·m)")
        ax2.set_ylim(0, 50); ax2.set_ylabel("Torque (N·m)", color="#d62728"); ax2.spines["right"].set_visible(True)
        h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=8, frameon=False, loc="upper right", ncol=2)
        ax.set_title("(b) Tool force and torque monitors"); ax.set_xlim(0, 360)
        plt.tight_layout(); plt.savefig(os.path.join(out, "fig03_convergence.png"), bbox_inches="tight"); plt.close()


def node_average(M, face_idx, field):
    """Area-weighted average of face values onto nodes (for smooth contouring)."""
    tris = M["FN"][face_idx]; vals = field[face_idx]; w = M["A"][face_idx]
    nv = np.zeros(len(M["X"])); nw = np.zeros(len(M["X"]))
    for k in range(3):
        np.add.at(nv, tris[:, k], vals * w); np.add.at(nw, tris[:, k], w)
    return tris, np.where(nw > 0, nv / np.maximum(nw, 1e-30), 300.0)


def fig_surface_temperature(M, out):
    Xmm = M["X"] * 1000; TF = M["faces"]["T"]
    tt, nT = node_average(M, top_triangles(M), TF)
    bt, bT = node_average(M, zslice(M, "buttomface"), TF)
    tri = mtri.Triangulation(Xmm[:, 0], Xmm[:, 1], tt); btri = mtri.Triangulation(Xmm[:, 0], Xmm[:, 1], bt)
    levels = np.linspace(300, 595, 60)
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.3), gridspec_kw={"width_ratios": [1.25, 1, 1]})
    cs = axs[0].tricontourf(tri, nT, levels=levels, cmap="turbo")
    axs[0].set_aspect("equal"); axs[0].set_title("(a) Top surface, full plate")
    axs[0].set_xlabel("x (mm)"); axs[0].set_ylabel("y (mm)")
    axs[0].add_patch(Circle((0, 0), 12, fill=False, ec="w", lw=0.8, ls="--"))
    axs[0].add_patch(Circle((0, 0), 3, fc="#888", ec="none"))
    for ax, t_, v_, t in [(axs[1], tri, nT, "(b) Top surface, close-up"), (axs[2], btri, bT, "(c) Bottom surface, close-up")]:
        ax.tricontourf(t_, v_, levels=levels, cmap="turbo")
        cl = ax.tricontour(t_, v_, levels=[350, 400, 450, 500], colors="k", linewidths=0.5)
        ax.clabel(cl, fmt="%d", fontsize=7, inline=True)
        ax.add_patch(Circle((0, 0), 12, fill=False, ec="w", lw=1, ls="--"))
        ax.set_xlim(-40, 40); ax.set_ylim(-40, 40); ax.set_aspect("equal"); ax.set_title(t); ax.set_xlabel("x (mm)")
        ax.annotate("", xy=(36, -34), xytext=(20, -34), arrowprops=dict(arrowstyle="-|>", color="w", lw=1.4))
        ax.text(28, -31, "flow", color="w", fontsize=7.5, ha="center")
    axs[1].add_patch(Circle((0, 0), 3, fc="#888", ec="k", lw=0.6))
    cb = fig.colorbar(cs, ax=axs, shrink=0.9, pad=0.015, ticks=np.arange(300, 600, 50))
    cb.set_label("Static temperature (K)")
    plt.savefig(os.path.join(out, "fig04_surface_temperature.png"), bbox_inches="tight"); plt.close()


def build_sections(M):
    """Linear 3-D interpolation of cell + boundary-face data onto section planes."""
    Cmm = M["C"] * 1000; c = M["cells"]; F = M["faces"]
    b = np.arange(M["zones"]["velocityinlet"][0] - 1, len(M["FN"]))  # boundary faces
    fcm = M["fc"][b] * 1000
    s = np.hypot(Cmm[:, 0], Cmm[:, 1]) < 45; sb = np.hypot(fcm[:, 0], fcm[:, 1]) < 45
    pts = np.vstack([Cmm[s], fcm[sb]])
    vals = np.column_stack([np.r_[c[k][s], F[k][b][sb]] for k in ["T", "U", "V", "W"]])
    print("  building 3-D Delaunay interpolant (~30-60 s)...")
    interp = LinearNDInterpolator(Delaunay(pts), vals)
    interp_mu = LinearNDInterpolator(Cmm[s], np.log10(c["MU_LAM"][s]))
    res = {}
    H, Z = np.meshgrid(np.linspace(-30, 30, 601), np.linspace(-6, 0, 121))
    for key, P in [("trans", np.column_stack([np.zeros(H.size), H.ravel(), Z.ravel()])),
                   ("long", np.column_stack([H.ravel(), np.zeros(H.size), Z.ravel()]))]:
        res[key] = (H, Z, interp(P).reshape(*H.shape, 4))
    gx, gy = np.meshgrid(np.linspace(-20, 20, 401), np.linspace(-20, 20, 401))
    for zc in [-0.5, -2.5]:
        P = np.column_stack([gx.ravel(), gy.ravel(), np.full(gx.size, zc)])
        res[zc] = (gx, gy, interp(P).reshape(*gx.shape, 4), interp_mu(P).reshape(gx.shape))
    return res


def fig_sections(res, out):
    levels = np.linspace(300, 595, 60)
    fig, axs = plt.subplots(2, 1, figsize=(10, 5.2))
    for ax, key, lab, t in [
        (axs[0], "trans", "y (mm)   [retreating side ←   → advancing side]", "(a) Transverse section x = 0 (perpendicular to the weld line)"),
        (axs[1], "long", "x (mm)   [leading edge ←   → trailing edge]", "(b) Longitudinal section y = 0 (along the weld line)")]:
        H, Z, o = res[key]; T = o[..., 0].copy(); T[(np.abs(H) < 3) & (Z > -5)] = np.nan
        cs = ax.contourf(H, Z, T, levels=levels, cmap="turbo", extend="both")
        cl = ax.contour(H, Z, T, levels=[400, 450, 500, 550], colors="k", linewidths=0.5); ax.clabel(cl, fmt="%d", fontsize=7)
        ax.add_patch(Rectangle((-3, -5), 6, 5, fc="#888", ec="k", lw=0.6))
        ax.add_patch(Rectangle((-12, 0), 24, 0.9, fc="#bbb", ec="k", lw=0.6, clip_on=False))
        ax.set_xlim(-30, 30); ax.set_ylim(-6, 0.9); ax.set_aspect(1.6)
        ax.set_xlabel(lab); ax.set_ylabel("z (mm)"); ax.set_title(t)
    cb = fig.colorbar(cs, ax=axs, shrink=0.9, pad=0.02, ticks=np.arange(300, 600, 50)); cb.set_label("Static temperature (K)")
    plt.savefig(os.path.join(out, "fig05_sections_temperature.png"), bbox_inches="tight"); plt.close()


def fig_lines(res, out):
    cols = {0: "#c0392b", -3: "#e67e22", -6: "#2471a3"}
    labs = {0: "top surface (z = 0)", -3: "mid-thickness (z = −3 mm)", -6: "bottom surface (z = −6 mm)"}
    with plt.rc_context(CLEAN):
        fig, axs = plt.subplots(1, 2, figsize=(11, 3.9), sharey=True)
        for ax, key, xl, t in [(axs[0], "long", "x (mm)", "(a) Along the weld line (y = 0)"),
                               (axs[1], "trans", "y (mm)", "(b) Across the weld line (x = 0)")]:
            H, Z, o = res[key]; h = H[0]
            for zc in [0, -3, -6]:
                j = np.argmin(np.abs(Z[:, 0] - zc)); T = o[j, :, 0].copy()
                if zc > -5: T[np.abs(h) < 3] = np.nan
                ax.plot(h, T, color=cols[zc], lw=1.8, label=labs[zc])
            ax.axvspan(-12, 12, color="#f5b041", alpha=0.15, lw=0); ax.axvspan(-3, 3, color="#6e2c00", alpha=0.18, lw=0)
            ax.text(0, 305, "pin", ha="center", fontsize=8); ax.text(-7.5, 305, "shoulder", ha="center", fontsize=8)
            ax.set_xlabel(xl); ax.set_title(t); ax.set_xlim(-30, 30); ax.grid(alpha=0.3)
        axs[0].text(-29, 575, "leading", fontsize=8.5, style="italic"); axs[0].text(29, 575, "trailing", fontsize=8.5, style="italic", ha="right")
        axs[1].text(-29, 575, "retreating", fontsize=8.5, style="italic"); axs[1].text(29, 575, "advancing", fontsize=8.5, style="italic", ha="right")
        axs[0].set_ylabel("Static temperature (K)"); axs[0].set_ylim(300, 600)
        h_, l_ = axs[0].get_legend_handles_labels()
        fig.legend(h_, l_, fontsize=9, frameon=False, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.06))
        plt.tight_layout(rect=(0, 0.05, 1, 1)); plt.savefig(os.path.join(out, "fig06_temperature_lines.png"), bbox_inches="tight"); plt.close()


def fig_flow(res, out):
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.2), gridspec_kw={"wspace": 0.42})
    gx, gy, o, mu = res[-2.5]; r = np.hypot(gx, gy)
    Vm = np.sqrt(o[..., 1]**2 + o[..., 2]**2 + o[..., 3]**2) * 1000; Vm[r < 3] = np.nan
    cs = axs[0].contourf(gx, gy, np.clip(Vm, 1, 250), levels=np.logspace(0, np.log10(250), 41),
                         cmap="viridis", norm=mc.LogNorm(1, 250))
    fig.colorbar(cs, ax=axs[0], shrink=0.85, label="|V| (mm/s), log scale", ticks=[1, 3, 10, 30, 100, 250], format="%g")
    U = np.nan_to_num(np.where(r < 3, 0, o[..., 1])); V = np.nan_to_num(np.where(r < 3, 0, o[..., 2]))
    axs[0].streamplot(gx[0], gy[:, 0], U, V, density=1.1, color="w", linewidth=0.5, arrowsize=0.6)
    axs[0].set_title("(a) Velocity, z = −2.5 mm")
    for ax, zc, t in [(axs[1], -0.5, "(b) log₁₀ viscosity, z = −0.5 mm"), (axs[2], -2.5, "(c) log₁₀ viscosity, z = −2.5 mm")]:
        g1, g2, _, m = res[zc]; m = m.copy(); m[np.hypot(g1, g2) < 3] = np.nan
        cs2 = ax.contourf(g1, g2, np.clip(m, 3, 6), levels=np.linspace(3, 6.001, 31), cmap="magma_r"); ax.set_title(t)
    fig.colorbar(cs2, ax=axs[1:], shrink=0.85, label="log₁₀ μ  (μ in Pa·s)", ticks=[3, 4, 5, 6])
    for ax in axs:
        ax.add_patch(Circle((0, 0), 3, fc="#888", ec="k", lw=0.6)); ax.add_patch(Circle((0, 0), 12, fill=False, ec="w", ls="--", lw=0.9))
        ax.set_aspect("equal"); ax.set_xlim(-20, 20); ax.set_ylim(-20, 20); ax.set_xlabel("x (mm)")
    axs[0].set_ylabel("y (mm)")
    plt.savefig(os.path.join(out, "fig07_velocity_viscosity.png"), bbox_inches="tight"); plt.close()


# -------------------------------------------------------------------- summary
def summary(M, mon_dir):
    c = M["cells"]; HF = M["faces"]["HEAT_FLUX"]
    final = {n: np.loadtxt(os.path.join(mon_dir, f"{n}-rfile.out"), skiprows=3)[-1, 1]
             for n in ["weldforce", "lateralforce", "thrustforce", "torque"]}
    P_tool = abs(final["torque"]) * OMEGA
    q = {z: np.nansum(HF[zslice(M, z)]) for z in ["buttomface", "topface", "wall-solid"]}
    print("\n=== Summary ===")
    print(f"Mesh: {M['ncell']:,} tetrahedral cells, {len(M['X']):,} nodes")
    print(f"Peak temperature: {c['T'].max():.1f} K ({c['T'].max()-273.15:.0f} °C)")
    print(f"Viscosity range: {c['MU_LAM'].min():.0f} - {c['MU_LAM'].max():.0f} Pa·s")
    for k, v in final.items(): print(f"{k:>13s}: {v:10.3f}")
    print(f"Tool power |M|·ω = {P_tool:.0f} W;  heat input per length = {P_tool/V_WELD/1000:.0f} J/mm")
    for z, v in q.items(): print(f"Heat through {z:>10s}: {v:9.1f} W")


def main():
    here = os.path.dirname(os.path.abspath(__file__)); root = os.path.dirname(here)
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", default=os.path.join(root, "Case_data_manual.cas.h5"))
    ap.add_argument("--data", default=os.path.join(root, "Case_data_manual.dat.h5"))
    ap.add_argument("--monitors", default=os.path.join(root, "FSW_files", "dp0", "FFF", "Fluent"))
    ap.add_argument("--out", default=os.path.join(root, "images"))
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    print("Loading case/data ..."); M = load(a.case, a.data)
    fig_domain(a.out); fig_mesh(M, a.out); fig_convergence(M, a.monitors, a.out); fig_surface_temperature(M, a.out)
    res = build_sections(M); fig_sections(res, a.out); fig_lines(res, a.out); fig_flow(res, a.out)
    summary(M, a.monitors); print(f"\nFigures written to {a.out}")


if __name__ == "__main__":
    main()
