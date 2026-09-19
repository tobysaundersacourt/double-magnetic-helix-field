#!/usr/bin/env python3
import argparse
import json
from math import pi
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.constants import mu_0
from geometry_utils import (
    strand_phases, strand_overlap_status, normal_plane_surface_contour,
    rotate_xy, rotate_points_about_z,
)


def progress(percent, message):
    print(f"PLOT_PROGRESS {int(percent)} {message}", flush=True)


def read_config(path):
    cfg = {}
    for raw in Path(path).read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        k, v = [s.strip() for s in line.split("=", 1)]
        cfg[k] = v
    floats = ["I", "a", "b", "wavevector", "outer_r", "delta"]
    ints = ["N", "index_cap", "amb_thetas"]
    for k in floats:
        if k in cfg:
            cfg[k] = float(cfg[k])
    for k in ints:
        if k in cfg:
            cfg[k] = int(float(cfg[k]))
    return cfg


def read_field_file(filename):
    with open(filename, "r") as f:
        lines = [f.readline() for _ in range(5)]
    if any(not x for x in lines):
        raise ValueError("Field file must contain exactly five nonempty data lines")

    settings = np.fromstring(lines[0].strip(), sep=",")
    coordinates = np.fromstring(lines[1].strip(), sep=",")
    c_integrals_t = np.fromstring(lines[2].strip(), sep=",")
    field = np.fromstring(lines[3].strip(), sep=",")
    thetas_in = np.fromstring(lines[4].strip(), sep=",")

    if settings.size != 5:
        raise ValueError(f"Expected 5 solver settings, got {settings.size}")
    if coordinates.size % 2 or field.size % 2:
        raise ValueError("Coordinate/field data is malformed")

    wavevector, a, b, J, N = settings
    r, t = coordinates.reshape((-1, 2)).T
    b_r, b_tz = field.reshape((-1, 2)).T
    if len(r) != len(c_integrals_t) or len(r) != len(thetas_in):
        raise ValueError("Solver output arrays have inconsistent lengths")

    c_integrals_z = c_integrals_t[-1] - c_integrals_t
    b_r = b_r * (2 * N * mu_0 * J / pi / wavevector)
    b_t = (
        N * mu_0 * J / pi * np.divide(c_integrals_t, r, out=np.zeros_like(r), where=r != 0)
        + 2 * N * mu_0 * J / pi / wavevector**2
        * np.divide(b_tz, r, out=np.zeros_like(r), where=r != 0)
    )
    b_z = N * mu_0 * J * wavevector / pi * c_integrals_z - 2 * N * mu_0 * J / pi / wavevector * b_tz

    ct = np.cos(t)
    st = np.sin(t)
    b_x = b_r * ct - b_t * st
    b_y = b_r * st + b_t * ct
    x = r * ct
    y = r * st
    b_xy = np.hypot(b_x, b_y)
    bmag = np.sqrt(b_x**2 + b_y**2 + b_z**2)

    return {
        "wavevector": wavevector, "a": a, "b": b, "J": J, "N": int(round(N)),
        "r": r, "t": t, "x": x, "y": y,
        "b_x": b_x, "b_y": b_y, "b_z": b_z,
        "b_xy": b_xy, "bmag": bmag, "thetas_in": thetas_in,
    }


class PolarFieldInterpolator:
    """Piecewise-linear interpolation on the solver's native polar rings.

    This avoids the visible nearest-neighbour ripples without building an
    enormous 2D triangulation of the full point cloud.
    """
    def __init__(self, data, delta):
        self.t = data["t"]
        self.fields = np.column_stack((data["b_x"], data["b_y"], data["b_z"]))
        r = data["r"]

        # Consecutive samples with identical r form one angular ring.
        cuts = np.flatnonzero(np.r_[True, r[1:] != r[:-1], True])
        starts = cuts[:-1]
        ends = cuts[1:]

        # Deduplicate radii at region boundaries using the radial grid index;
        # later rings replace earlier copies, which favors the interior/surface
        # expression at a-b when both are present.
        scale = max(abs(delta), 1e-15)
        ring_by_key = {}
        for s, e in zip(starts, ends):
            if e <= s:
                continue
            rad = float(r[s])
            key = int(round(rad / scale))
            ring_by_key[key] = (rad, int(s), int(e))

        rings = sorted(ring_by_key.values(), key=lambda q: q[0])
        self.radii = np.array([q[0] for q in rings])
        self.rings = [(q[1], q[2]) for q in rings]
        if len(self.radii) < 2:
            raise ValueError("Not enough radial rings to interpolate the field")

        # Sort each angular ring once.  The earlier implementation sorted a
        # ring every time it was sampled, which is needlessly expensive when
        # building a Cartesian x-y preview grid.
        self.ring_t = []
        self.ring_v = []
        for s, e in self.rings:
            tt = self.t[s:e]
            vv = self.fields[s:e]
            if len(tt) > 1:
                order = np.argsort(tt)
                tt = tt[order]
                vv = vv[order]
            self.ring_t.append(tt)
            self.ring_v.append(vv)

    def _angular(self, ring_index, theta):
        tt = self.ring_t[ring_index]
        vv = self.ring_v[ring_index]
        if len(tt) == 1:
            return vv[0].copy()
        tq = float(theta % (2 * pi))
        tt_aug = np.r_[tt, tt[0] + 2 * pi]
        out = np.empty(3)
        for j in range(3):
            v_aug = np.r_[vv[:, j], vv[0, j]]
            out[j] = np.interp(tq, tt_aug, v_aug)
        return out

    def sample_polar(self, radius, theta):
        rq = float(np.clip(radius, self.radii[0], self.radii[-1]))
        hi = int(np.searchsorted(self.radii, rq, side="left"))
        if hi <= 0:
            return self._angular(0, theta)
        if hi >= len(self.radii):
            return self._angular(len(self.radii) - 1, theta)
        lo = hi - 1
        r0, r1 = self.radii[lo], self.radii[hi]
        f0 = self._angular(lo, theta)
        if r1 == r0:
            return f0
        f1 = self._angular(hi, theta)
        w = (rq - r0) / (r1 - r0)
        return (1 - w) * f0 + w * f1

    def sample_xy(self, x, y):
        return self.sample_polar(np.hypot(x, y), np.arctan2(y, x))


def make_xy_field_grid(interp, outer_r, resolution=241):
    """Sample the smooth z=0 field on a regular Cartesian grid.

    This is intentionally a user-facing resampled grid rather than a duplicate
    dump of every native polar solver sample.  The original full-resolution
    compact solver output remains in ``field.csv``.
    """
    resolution = max(33, int(resolution))
    axis = np.linspace(-outer_r, outer_r, resolution)
    X, Y = np.meshgrid(axis, axis, indexing="xy")
    R = np.hypot(X, Y)
    mask = R <= outer_r
    Bx = np.full_like(X, np.nan, dtype=float)
    By = np.full_like(X, np.nan, dtype=float)
    Bz = np.full_like(X, np.nan, dtype=float)

    rows, cols = np.nonzero(mask)
    for iy, ix in zip(rows, cols):
        bx, by, bz = interp.sample_xy(float(X[iy, ix]), float(Y[iy, ix]))
        Bx[iy, ix] = bx
        By[iy, ix] = by
        Bz[iy, ix] = bz

    Bxy = np.hypot(Bx, By)
    Bmag = np.sqrt(Bx**2 + By**2 + Bz**2)
    return {
        "axis": axis, "X": X, "Y": Y, "R": R, "mask": mask,
        "Bx": Bx, "By": By, "Bz": Bz, "Bxy": Bxy, "Bmag": Bmag,
    }


def save_xy_field_csv(path, grid):
    rows, cols = np.nonzero(grid["mask"])
    x = grid["X"][rows, cols]
    y = grid["Y"][rows, cols]
    r = grid["R"][rows, cols]
    theta = np.mod(np.arctan2(y, x), 2*pi)
    table = np.column_stack((
        rows, cols, x, y, r, theta,
        grid["Bx"][rows, cols], grid["By"][rows, cols],
        grid["Bz"][rows, cols], grid["Bxy"][rows, cols],
        grid["Bmag"][rows, cols],
    ))
    np.savetxt(
        path, table, delimiter=",",
        header="grid_row,grid_col,x,y,r,theta,b_x,b_y,b_z,b_xy,b_mag",
        comments="",
    )


def plot_xy_scalar(outpath, grid, values, title, cbar_label, dpi, signed=False):
    fig, ax = plt.subplots(figsize=(7, 6))
    kwargs = {}
    finite = np.isfinite(values)
    if signed and np.any(finite):
        lim = float(np.nanmax(np.abs(values)))
        if lim > 0:
            kwargs.update(vmin=-lim, vmax=lim, cmap="coolwarm")
    mesh = ax.pcolormesh(grid["X"], grid["Y"], values, shading="auto", **kwargs)
    fig.colorbar(mesh, ax=ax, label=cbar_label)
    wall = plt.Circle((0, 0), float(np.nanmax(grid["R"][grid["mask"]])),
                      fill=False, linestyle="--", linewidth=1.0)
    ax.add_patch(wall)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def set_axes_equal_3d(ax):
    ranges = [
        abs(ax.get_xlim3d()[1] - ax.get_xlim3d()[0]),
        abs(ax.get_ylim3d()[1] - ax.get_ylim3d()[0]),
        abs(ax.get_zlim3d()[1] - ax.get_zlim3d()[0]),
    ]
    mids = [
        sum(ax.get_xlim3d()) / 2,
        sum(ax.get_ylim3d()) / 2,
        sum(ax.get_zlim3d()) / 2,
    ]
    rad = max(ranges) / 2
    ax.set_xlim3d(mids[0] - rad, mids[0] + rad)
    ax.set_ylim3d(mids[1] - rad, mids[1] + rad)
    ax.set_zlim3d(mids[2] - rad, mids[2] + rad)


def make_plane_contour(data, interp, contour_points=721):
    """Field sampled on the true helix-normal plane/surface intersection.

    Geometry comes from the same helper used by the live GUI preview, so the
    previewed contour and the post-processed contour are identical.
    """
    k = data["wavevector"]
    geom = normal_plane_surface_contour(
        data["a"], data["b"], k, contour_points=contour_points
    )
    xyz = geom["xyz"]

    B = np.empty_like(xyz)
    for i, (xq, yq, zq) in enumerate(xyz):
        # Undo the screw rotation to sample the native z=0 field, then rotate
        # the field vector itself back into the 3D laboratory frame.
        phi = k * zq
        x0, y0 = rotate_xy(xq, yq, -phi)
        bx0, by0, bz0 = interp.sample_xy(x0, y0)
        bxq, byq = rotate_xy(bx0, by0, phi)
        B[i] = (bxq, byq, bz0)

    bmag = np.linalg.norm(B, axis=1)
    bnormal = B @ geom["nhat"]
    binplane = np.sqrt(np.maximum(bmag**2 - bnormal**2, 0))
    ds = np.linalg.norm(np.diff(xyz, axis=0), axis=1)
    s = np.r_[0.0, np.cumsum(ds)]

    return {
        **geom,
        "B": B, "s": s, "bmag": bmag,
        "bnormal": bnormal, "binplane": binplane,
    }


def generate_plots(field_file, config_file, output_dir, contour_points=721, dpi=160, xy_grid_resolution=241):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    cfg = read_config(config_file)
    data = read_field_file(field_file)
    delta = float(cfg.get("delta", 1e-7))
    outer_r = float(cfg.get("outer_r", np.max(data["r"])))
    I = float(cfg.get("I", 100000.0))

    progress(5, "Building smooth polar interpolator")
    interp = PolarFieldInterpolator(data, delta)

    progress(12, "Sampling the general x-y plane field")
    xy_grid = make_xy_field_grid(interp, outer_r, resolution=xy_grid_resolution)
    save_xy_field_csv(out / "field_xy_plane.csv", xy_grid)

    progress(25, "Solving helix-normal plane/surface contour")
    contour = make_plane_contour(data, interp, contour_points=contour_points)
    xyz = contour["xyz"]
    B = contour["B"]
    phases = strand_phases(data["N"])

    # The reference contour is on strand 0.  The corresponding contours for
    # all other strands are exact rigid rotations by their strand phases.
    all_contours_xyz = [rotate_points_about_z(xyz, float(p)) for p in phases]

    np.savetxt(
        out / "b_normal_plane_contour.csv",
        np.column_stack((
            contour["s"], contour["u"], xyz, B,
            contour["bnormal"], contour["binplane"], contour["bmag"]
        )),
        delimiter=",",
        header="s,surface_angle,x,y,z,b_x,b_y,b_z,b_normal,b_in_plane,b_mag",
        comments="",
    )

    # Save all N new contour geometries as one table.  Field values need only
    # be stored for the reference contour because the N-strand geometry is
    # related by azimuthal symmetry.
    contour_rows = []
    for j, (phase, q) in enumerate(zip(phases, all_contours_xyz)):
        contour_rows.append(np.column_stack((
            np.full(len(q), j),
            np.full(len(q), phase),
            contour["u"],
            q,
        )))
    np.savetxt(
        out / "normal_plane_contours_geometry.csv",
        np.vstack(contour_rows),
        delimiter=",",
        header="strand_index,strand_phase,surface_angle,x,y,z",
        comments="",
    )

    progress(42, "Plotting general x-y field magnitude")
    plot_xy_scalar(
        out / "field_xy_magnitude.png", xy_grid, xy_grid["Bmag"],
        "Magnetic-field magnitude in the x-y plane", "|B| [T]", dpi, signed=False
    )

    progress(54, "Plotting general x-y axial field")
    plot_xy_scalar(
        out / "field_xy_bz.png", xy_grid, xy_grid["Bz"],
        "Axial magnetic field in the x-y plane", "$B_z$ [T]", dpi, signed=True
    )

    progress(66, "Plotting general x-y in-plane field")
    fig, ax = plt.subplots(figsize=(7, 6))
    mesh = ax.pcolormesh(
        xy_grid["X"], xy_grid["Y"], xy_grid["Bxy"], shading="auto"
    )
    fig.colorbar(mesh, ax=ax, label=r"$\sqrt{B_x^2+B_y^2}$ [T]")
    step = max(1, xy_grid_resolution // 22)
    Xq = xy_grid["X"][::step, ::step]
    Yq = xy_grid["Y"][::step, ::step]
    Uq = xy_grid["Bx"][::step, ::step]
    Vq = xy_grid["By"][::step, ::step]
    good = np.isfinite(Uq) & np.isfinite(Vq)
    ax.quiver(Xq[good], Yq[good], Uq[good], Vq[good], angles="xy")
    wall = plt.Circle((0, 0), outer_r, fill=False, linestyle="--", linewidth=1.0)
    ax.add_patch(wall)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title("In-plane magnetic field in the x-y plane")
    fig.tight_layout()
    fig.savefig(out / "field_xy_inplane.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    progress(74, "Plotting magnetic field around the reference contour")
    # Plot physically natural components relative to the helix-normal plane.
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(contour["s"], contour["bmag"], label=r"$|\mathbf{B}|$")
    ax.plot(contour["s"], contour["bnormal"], label=r"$B_{\perp}$")
    ax.plot(contour["s"], contour["binplane"], label=r"$|B_{\parallel}|$")
    ax.set_xlabel("Arc length around contour s [m]")
    ax.set_ylabel("Magnetic field [T]")
    ax.set_title("Magnetic field around the helix-normal contour")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "b_normal_plane_contour.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    # Also retain the laboratory Cartesian components so the complete vector
    # evolution around the contour is visible rather than only its magnitude.
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(contour["s"], B[:, 0], label=r"$B_x$")
    ax.plot(contour["s"], B[:, 1], label=r"$B_y$")
    ax.plot(contour["s"], B[:, 2], label=r"$B_z$")
    ax.set_xlabel("Arc length around contour s [m]")
    ax.set_ylabel("Magnetic field [T]")
    ax.set_title("Cartesian field components around the helix-normal contour")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "b_normal_plane_contour_components.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    progress(82, "Plotting all new contours on the 3D helix geometry")
    k, a, b = data["wavevector"], data["a"], data["b"]
    zmins = [q[:, 2].min() for q in all_contours_xyz]
    zmaxs = [q[:, 2].max() for q in all_contours_xyz]
    zmin, zmax = min(zmins), max(zmaxs)
    zspan = max(zmax-zmin, 2*abs(b), 1e-6)
    zgrid = np.linspace(zmin-0.55*zspan, zmax+0.55*zspan, 150)
    ugrid = np.linspace(0, 2*pi, 70)
    U, Z = np.meshgrid(ugrid, zgrid)

    def strand_surface(phase0):
        X0 = a + b*np.cos(U)
        Y0 = b*np.sin(U)
        P = k*Z + phase0
        return np.cos(P)*X0 - np.sin(P)*Y0, np.sin(P)*X0 + np.cos(P)*Y0

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(projection="3d")
    surf_alpha = min(0.22, max(0.05, 0.55 / max(data["N"], 1)))
    for phase0 in phases:
        Xs, Ys = strand_surface(float(phase0))
        ax.plot_surface(Xs, Ys, Z, alpha=surf_alpha, linewidth=0)

    # Outer cylindrical wall.
    wall_t = np.linspace(0, 2*pi, 70)
    WT, WZ = np.meshgrid(wall_t, np.linspace(zgrid.min(), zgrid.max(), 24))
    WX, WY = outer_r*np.cos(WT), outer_r*np.sin(WT)
    ax.plot_surface(WX, WY, WZ, alpha=0.06, linewidth=0)

    # Centerlines for all strands.
    for j, phase0 in enumerate(phases):
        ax.plot(
            a*np.cos(k*zgrid + phase0), a*np.sin(k*zgrid + phase0), zgrid,
            "--", linewidth=1.0, label="Strand centerlines" if j == 0 else None
        )

    # Draw every new normal-plane/surface contour.
    for j, q in enumerate(all_contours_xyz):
        ax.plot(
            q[:, 0], q[:, 1], q[:, 2],
            linewidth=3.0 if j == 0 else 2.0,
            label="Helix-normal plane/surface contours" if j == 0 else None,
        )

    # Show the reference contour's fixed normal plane.
    center, nhat = contour["center"], contour["nhat"]
    e1 = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(e1, nhat)) > 0.95:
        e1 = np.array([0.0, 1.0, 0.0])
    e1 = e1 - np.dot(e1, nhat) * nhat
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(nhat, e1); e2 /= np.linalg.norm(e2)
    ph = 1.35*max(abs(b), 0.5*zspan)
    q1, q2 = np.meshgrid(np.linspace(-ph, ph, 20), np.linspace(-ph, ph, 20))
    patch = center[:, None, None] + e1[:, None, None]*q1 + e2[:, None, None]*q2
    ax.plot_surface(patch[0], patch[1], patch[2], alpha=0.15, linewidth=0)

    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.set_zlabel("z [m]")
    overlap = strand_overlap_status(a, b, k, data["N"])
    ax.set_title(f"Helix-normal plane/surface contours on all {data['N']} strands")
    if overlap["overlap"] or overlap["touching"]:
        ax.text2D(
            0.02, 0.98, overlap["message"], transform=ax.transAxes,
            va="top", fontsize=9
        )
    ax.legend(loc="upper left")
    set_axes_equal_3d(ax)
    fig.tight_layout()
    fig.savefig(out / "normal_plane_contours_on_helix_3d.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    progress(94, "Writing diagnostics")
    B0 = interp.sample_xy(0.0, 0.0)
    solenoid = mu_0 * I * k / (2*pi)
    overlap = strand_overlap_status(data["a"], data["b"], data["wavevector"], data["N"])
    summary = {
        "strand_count_N": int(data["N"]),
        "neighbor_min_centerline_distance_m": float(overlap["min_centerline_distance"]),
        "neighbor_surface_clearance_m": float(overlap["surface_clearance"]),
        "strands_overlap": bool(overlap["overlap"]),
        "strands_touching": bool(overlap["touching"]),
        "B_z_at_origin_T": float(B0[2]),
        "comparable_solenoid_B_z_T": float(solenoid),
        "normal_plane_contour_closure_error_before_fix_m": contour["closure_before"],
        "max_plane_equation_residual": contour["max_plane_residual"],
        "mean_B_normal_T": float(np.mean(contour["bnormal"])),
        "mean_B_magnitude_T": float(np.mean(contour["bmag"])),
        "solver_samples": int(len(data["r"])),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))

    # General x-y field plots are independent of the old z=0 strand-
    # intersection contour.  The contour plots below refer only to the new
    # helix-normal plane/surface intersection.
    plots = [
        "field_xy_magnitude.png",
        "field_xy_bz.png",
        "field_xy_inplane.png",
        "b_normal_plane_contour.png",
        "b_normal_plane_contour_components.png",
        "normal_plane_contours_on_helix_3d.png",
    ]
    csv_exports = [
        "field_xy_plane.csv",
        "b_normal_plane_contour.csv",
        "normal_plane_contours_geometry.csv",
    ]
    (out / "plots_manifest.json").write_text(json.dumps(plots, indent=2))
    (out / "csv_manifest.json").write_text(json.dumps(csv_exports, indent=2))
    progress(100, "Complete")
    return plots


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--contour-points", type=int, default=721)
    ap.add_argument("--dpi", type=int, default=160)
    ap.add_argument("--xy-grid-resolution", type=int, default=241)
    args = ap.parse_args()
    generate_plots(args.input, args.config, args.output_dir, args.contour_points, args.dpi, args.xy_grid_resolution)


if __name__ == "__main__":
    main()
