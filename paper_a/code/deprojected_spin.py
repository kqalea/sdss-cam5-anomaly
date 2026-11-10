"""Spin classification using catalogue shape parameters to deproject first.

The baseline classifier finds the centre as the brightest smoothed pixel and the
outer radius from the azimu-thally averaged profile. Both are noisy. The JADES
photometric catalogue gives a centroid (X, Y), a Kron radius, a position angle
and an axis ratio for every object, which removes both guesses.

Deprojection matters because a spiral seen at 60 degrees inclination has its
arms foreshortened along one axis; sampling on ellipses of axis ratio q removes
that distortion and leaves the winding phase intact. Handedness is preserved:
the distortion is a linear map with positive determinant (scaling one axis by
q > 0), so it cannot flip chirality. Only the *quality* of the measurement
improves.

The ellipse parameterisation: with major axis at position angle pa and minor to
major ratio q, the image-plane point for deprojected (r, theta) is

    x = cx + r cos(theta) cos(pa) - q r sin(theta) sin(pa)
    y = cy + r cos(theta) sin(pa) + q r sin(theta) cos(pa)
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates

FWHM_TO_SIGMA = 1.0 / np.sqrt(8.0 * np.log(2.0))


def sample_ellipse(img: np.ndarray, cx: float, cy: float, r_max: float,
                   q: float, pa_deg: float, n_r: int = 64,
                   n_t: int = 180, r_min_frac: float = 0.15
                   ) -> tuple[np.ndarray, np.ndarray]:
    """Sample the image on concentric ellipses; returns (radii, values[n_r, n_t])."""
    q = float(np.clip(q, 0.2, 1.0))          # keep the minor axis usable
    pa = np.radians(pa_deg)
    cpa, spa = np.cos(pa), np.sin(pa)
    radii = np.linspace(r_min_frac * r_max, r_max, n_r)
    theta = np.linspace(0.0, 2.0 * np.pi, n_t, endpoint=False)
    R, T = np.meshgrid(radii, theta, indexing="ij")
    x = cx + R * np.cos(T) * cpa - q * R * np.sin(T) * spa
    y = cy + R * np.cos(T) * spa + q * R * np.sin(T) * cpa
    vals = map_coordinates(img, [y.ravel(), x.ravel()], order=1,
                           mode="constant", cval=np.nan)
    return radii, vals.reshape(n_r, n_t)


def phase_slope(radii: np.ndarray, vals: np.ndarray, m: int = 2
                ) -> tuple[float, float, int]:
    """Signed slope of the m-th azimuthal phase against ln r, plus a strength."""
    theta = np.linspace(0.0, 2.0 * np.pi, vals.shape[1], endpoint=False)
    V = vals - np.nanmean(vals, axis=1, keepdims=True)
    c = np.nansum(V * np.exp(-1j * m * theta)[None, :], axis=1)
    amp = np.abs(c)
    phase = np.angle(c)
    mean_level = np.nanmean(vals, axis=1)

    keep = (amp > 0.25 * np.nanmax(amp)) & (mean_level > 0)
    idx = np.where(keep)[0]
    if len(idx) >= 2:
        runs = np.split(idx, np.where(np.diff(idx) != 1)[0] + 1)
        idx = max(runs, key=len)
    if len(idx) < 8:
        return 0.0, 0.0, len(idx)

    dphi = np.angle(np.exp(1j * np.diff(phase[idx])))
    dlr = np.diff(np.log(radii[idx]))
    w = 0.5 * (amp[idx][1:] + amp[idx][:-1])
    slope = float(np.sum(w * dphi / dlr) / np.sum(w))

    c0 = np.abs(mean_level[idx]) * vals.shape[1]
    strength = float(np.clip(np.median(amp[idx] / np.maximum(c0, 1e-12)), 0.0, 1.0))
    return slope, strength, len(idx)


# Defaults tuned against Shamir (2025)'s labels on the GOODS-S mosaic:
# smooth 0.5 px, inner 25% masked, sample out to R_KRON gives 79% agreement on
# 244 of the 257 catalogue-matched galaxies. The plateau is around 79%; pushing
# the inner mask out past 40% raises the recovered amplitude but loses agreement,
# and m = 1 carries no information at all (49%).
DEFAULT_SMOOTH = 0.5
DEFAULT_R_MIN_FRAC = 0.25


def measure_deproj(img: np.ndarray, cx: float, cy: float, r_max: float,
                   q: float, pa_deg: float, smooth_px: float = DEFAULT_SMOOTH,
                   n_r: int = 64, n_t: int = 180, m: int = 2,
                   r_min_frac: float = DEFAULT_R_MIN_FRAC) -> dict:
    """Classify one galaxy given its catalogue centre, radius, axis ratio and PA."""
    img = np.nan_to_num(img, nan=0.0)
    if smooth_px > 0:
        img = gaussian_filter(img, smooth_px)
    radii, vals = sample_ellipse(img, cx, cy, r_max, q, pa_deg, n_r, n_t, r_min_frac)
    slope, strength, nrad = phase_slope(radii, vals, m=m)
    ok = nrad >= 8 and abs(slope) > 1e-6
    return {"sign": int(np.sign(slope)) if ok else 0, "slope": slope,
            "strength": strength, "n_radii": nrad, "ok": bool(ok)}
