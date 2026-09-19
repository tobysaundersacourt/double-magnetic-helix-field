"""Geometry helpers shared by the GUI preview and plotting pipeline."""
from math import pi, sin
import numpy as np
from scipy.optimize import newton, brentq


def strand_phases(N):
    """Return N evenly spaced azimuthal phases in [0, 2*pi)."""
    N = int(N)
    if N <= 0:
        raise ValueError("N must be positive")
    return 2.0 * pi * np.arange(N, dtype=float) / N


def rotate_xy(x, y, angle):
    """Rotate Cartesian x/y coordinates or vector components about +z."""
    c = np.cos(angle)
    s = np.sin(angle)
    return c * x - s * y, s * x + c * y


def rotate_points_about_z(points, angle):
    """Return an (M,3) point array rotated rigidly about the z axis."""
    pts = np.asarray(points, dtype=float)
    out = pts.copy()
    out[:, 0], out[:, 1] = rotate_xy(pts[:, 0], pts[:, 1], angle)
    return out


def rotate_vector_about_z(vector, angle):
    """Rotate one 3-vector rigidly about the z axis."""
    v = np.asarray(vector, dtype=float)
    x, y = rotate_xy(v[0], v[1], angle)
    return np.array([x, y, v[2]], dtype=float)


def normal_plane_surface_contour(a, b, wavevector, contour_points=721):
    """Solve the reference-strand helix-surface / normal-plane intersection.

    The reference strand has centerline position (a, 0, 0) at z=0. Its
    centerline tangent there is proportional to (0, a*k, 1), so the fixed
    plane through (a,0,0) with that normal is orthogonal to the helix tangent.

    The strand surface is the screw rotation of the z=0 circle

        (a + b cos u, b sin u, 0)

    by angle k*z.  The returned curve is the closed branch of the intersection
    continuously connected to z=0 at u=0.

    Contours for the other strands are exact rigid z-rotations of this one by
    their phases 2*pi*j/N; no additional root solve is needed.
    """
    a = float(a)
    b = float(b)
    k = float(wavevector)
    contour_points = int(contour_points)
    if contour_points < 16:
        raise ValueError("contour_points must be at least 16")
    if abs(k) < 1e-15:
        raise ValueError("wavevector is zero")

    def surface_xyz(u, z):
        x0 = a + b * np.cos(u)
        y0 = b * np.sin(u)
        x, y = rotate_xy(x0, y0, k * z)
        return np.array([x, y, z], dtype=float)

    center = np.array([a, 0.0, 0.0])
    normal = np.array([0.0, a * k, 1.0])
    nhat = normal / np.linalg.norm(normal)

    def residual(u, z):
        return float(np.dot(normal, surface_xyz(u, z) - center))

    def residual_dz(u, z):
        xq = surface_xyz(u, z)[0]
        return normal[1] * k * xq + normal[2]

    # Conservative interval for fallback bracketing.  The first term follows
    # from |z| = |a k y| on the plane, with |y| <= a+b for the surface.
    z_bound = abs(a * k) * (abs(a) + abs(b)) + max(abs(b), 1.0) * 1e-6
    z_bound = max(z_bound, 1e-8)

    def fallback_roots(u):
        zz = np.linspace(-z_bound, z_bound, 2501)
        ff = np.array([residual(u, z) for z in zz])
        roots = []
        exact = np.where(np.abs(ff) < 1e-12)[0]
        roots.extend(zz[exact].tolist())
        changes = np.where(ff[:-1] * ff[1:] < 0)[0]
        for j in changes:
            roots.append(
                brentq(
                    lambda z: residual(u, z),
                    zz[j], zz[j + 1], xtol=1e-12, rtol=1e-12,
                )
            )
        if not roots:
            return np.array([])
        roots = np.array(sorted(roots))
        return roots[np.r_[True, np.diff(roots) > 1e-8]]

    uu = np.linspace(0.0, 2.0 * pi, contour_points)
    zz = np.empty_like(uu)
    zz[0] = 0.0
    prev = 0.0
    for i in range(1, len(uu)):
        u = uu[i]
        valid = False
        try:
            zq = newton(
                lambda z: residual(u, z),
                prev,
                fprime=lambda z: residual_dz(u, z),
                tol=1e-12,
                maxiter=50,
            )
            valid = (
                np.isfinite(zq)
                and abs(zq) <= 1.05 * z_bound
                and abs(residual(u, zq)) < 1e-9
            )
        except (RuntimeError, OverflowError, ZeroDivisionError):
            pass
        if not valid:
            roots = fallback_roots(u)
            if len(roots) == 0:
                raise RuntimeError(
                    f"Could not find plane/surface intersection at u={u}"
                )
            zq = roots[np.argmin(np.abs(roots - prev))]
        zz[i] = zq
        prev = zq

    xyz = np.array([surface_xyz(u, z) for u, z in zip(uu, zz)])
    closure_before = float(np.linalg.norm(xyz[-1] - xyz[0]))
    xyz[-1] = xyz[0]
    zz[-1] = zz[0]
    residuals = (xyz - center) @ normal

    return {
        "u": uu,
        "z": zz,
        "xyz": xyz,
        "center": center,
        "normal": normal,
        "nhat": nhat,
        "closure_before": closure_before,
        "max_plane_residual": float(np.max(np.abs(residuals))),
    }


def strand_normal_plane_contours(a, b, wavevector, N, contour_points=721):
    """Return the new normal-plane/surface contour for every strand."""
    reference = normal_plane_surface_contour(a, b, wavevector, contour_points)
    contours = []
    for phase in strand_phases(N):
        contours.append({
            "phase": float(phase),
            "xyz": rotate_points_about_z(reference["xyz"], phase),
            "center": rotate_vector_about_z(reference["center"], phase),
            "normal": rotate_vector_about_z(reference["normal"], phase),
            "nhat": rotate_vector_about_z(reference["nhat"], phase),
            "u": reference["u"],
            "z": reference["z"],
            "closure_before": reference["closure_before"],
            "max_plane_residual": reference["max_plane_residual"],
        })
    return contours


def adjacent_strand_centerline_spacing(a, N):
    """Exact neighboring center spacing in every fixed-z cross-section.

    The modeled strand surfaces are generated by screw-rotating the z=0
    circular cross-sections.  Their centerlines therefore lie on a circle of
    radius ``a`` and neighboring strands differ by an azimuthal phase 2*pi/N.
    At any common z their center separation is the chord

        d = 2 a sin(pi/N).

    Since both strand cross-sections have radius ``b``, their modeled surfaces
    touch at d=2b and overlap at d<2b.
    """
    N = int(N)
    if N <= 1:
        return float("inf")
    return 2.0 * float(a) * sin(pi / N)


def strand_overlap_status(a, b, wavevector, N, rel_tol=1e-9):
    """Describe whether the modeled neighboring strand surfaces touch/overlap.

    ``wavevector`` is accepted so callers can pass the complete helix geometry;
    the overlap criterion is independent of it for this screw-generated surface,
    because all intersecting points necessarily have the same z coordinate.
    """
    del wavevector
    N = int(N)
    if N <= 1:
        return {
            "overlap": False,
            "touching": False,
            "min_centerline_distance": float("inf"),
            "surface_clearance": float("inf"),
            "message": "",
        }

    spacing = adjacent_strand_centerline_spacing(a, N)
    diameter = 2.0 * float(b)
    scale = max(abs(spacing), abs(diameter), 1e-15)
    tol = rel_tol * scale
    overlap = spacing < diameter - tol
    touching = (not overlap) and abs(spacing - diameter) <= tol
    clearance = spacing - diameter

    if overlap:
        msg = (
            "WARNING: neighboring helical strand surfaces overlap. "
            f"Their fixed-z center spacing is {spacing:.6g} m, "
            f"while 2b = {diameter:.6g} m."
        )
    elif touching:
        msg = (
            "WARNING: neighboring helical strand surfaces touch. "
            f"Their fixed-z center spacing ≈ 2b = {diameter:.6g} m."
        )
    else:
        msg = ""

    return {
        "overlap": overlap,
        "touching": touching,
        "min_centerline_distance": spacing,
        "surface_clearance": clearance,
        "message": msg,
    }
