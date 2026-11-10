"""The one implementation of the spin measurement.

Two days of analysis produced twenty-one copies of this calculation, in files that
drifted apart, and every bug that mattered was a copy that did not receive a fix:
the second moments taken over the whole frame, the position-angle offset missing on
COSMOS-Web, the pixel-scale divisions applied twice.

**New code must import from here rather than reimplement.** If a change is needed,
it belongs in this file so that every caller receives it.

The two things that bit hardest are centralised as well: the pixel scales per survey
and the world-to-image position angle conversion.
"""

from __future__ import annotations

import numpy as np

from deprojected_spin import sample_ellipse

# arcsec per pixel, per survey. JADES and COSMOS-Web both sample at 0.03, but they
# are named separately so that a future survey cannot silently inherit the wrong one.
PIXSCALE = {"jades": 0.03, "cosmosweb": 0.03}

# The COSMOS-Web cutout WCS puts the image +x axis at sky position angle -70 degrees,
# so a world position angle maps to the operator frame as theta + 70. Omitting this
# mis-orients the deprojection by 70 degrees. JADES mosaics are north-up with no
# rotation, so the offset there is zero.
PA_IMAGE_OFFSET = {"jades": 0.0, "cosmosweb": +70.0}

KEEP_FRACTION = 0.25   # annulus kept when its m=2 amplitude exceeds this times the peak
MIN_ANNULI = 8         # below this the fit is not attempted


def pa_image(survey: str, theta_world: float) -> float:
    """Convert a world position angle to the operator frame for one survey."""
    return float(theta_world) + PA_IMAGE_OFFSET[survey]


def kron_pixels(survey: str, radius_arcsec: float, lo: float = 5.0, hi: float = 120.0) -> float:
    """Kron radius in pixels, clipped. The catalogue units must already be arcsec."""
    return float(np.clip(radius_arcsec / PIXSCALE[survey], lo, hi))


def measure(img, q: float, pa: float, rk: float) -> dict | None:
    """Signed phase-slope, its sign, the m=2 signal-to-noise and the m=2 concentration.

    Parameters
    ----------
    img : 2-D array
        Deprojected-frame image, already background subtracted.
    q, pa, rk : float
        Axis ratio, operator-frame position angle in degrees, Kron radius in pixels.

    Returns
    -------
    dict with ``slope``, ``sign``, ``snr``, ``m2_frac``, ``n_radii``, or None when the
    galaxy has fewer than MIN_ANNULI usable annuli.
    """
    h, w = img.shape
    radii, vals = sample_ellipse(img, w / 2.0, h / 2.0, rk, q, pa)
    n_t = vals.shape[1]
    theta = np.linspace(0.0, 2 * np.pi, n_t, endpoint=False)
    V = vals - np.nanmean(vals, axis=1, keepdims=True)
    c = np.nansum(V * np.exp(-2j * theta)[None, :], axis=1)
    amp2 = np.abs(c)
    phase = np.angle(c)
    mean_level = np.nanmean(vals, axis=1)
    keep = (amp2 > KEEP_FRACTION * np.nanmax(amp2)) & (mean_level > 0)
    idx = np.where(keep)[0]
    if len(idx) < MIN_ANNULI:
        return None

    x = np.log(radii[idx])
    ph = np.unwrap(phase[idx])
    coef = np.polyfit(x, ph, 1)
    resid = ph - np.polyval(coef, x)

    ct, st = np.cos(2 * theta), np.sin(2 * theta)
    snrs = []
    for i in idx:
        v = V[i]
        if not np.all(np.isfinite(v)):
            continue
        A = 2.0 / n_t * np.sum(v * ct)
        B = 2.0 / n_t * np.sum(v * st)
        amp = float(np.hypot(A, B))
        noise = float(np.std(v - A * ct - B * st))
        if noise <= 0:
            continue
        snrs.append(amp / (noise * np.sqrt(2.0 / n_t)))

    total = np.nansum(np.abs(V), axis=1)
    return {
        "slope": float(coef[0]),
        "sign": int(np.sign(coef[0])),
        "snr": float(np.median(snrs)) if snrs else 0.0,
        "m2_frac": float(np.nanmax(amp2) / max(np.nanmax(total), 1e-12)),
        "n_radii": int(len(idx)),
        "phase_rms": float(np.std(resid)),
    }


def asymmetry(sign) -> tuple[float, float, int]:
    """Mean sign and its binomial error: A = (N+ - N-)/N, sigma = sqrt((1-A^2)/N)."""
    s = np.asarray(sign)
    n = len(s)
    a = (2 * int((s > 0).sum()) - n) / n
    return a, float(np.sqrt(max(1 - a * a, 0) / n)), n


def projection_fit(coss, As, errs):
    """Fit A = A_true cos(theta) + C. Returns (A_true, err, C, err, chi2, dof)."""
    c = np.asarray(coss, float)
    a = np.asarray(As, float)
    e = np.asarray(errs, float)
    X = np.column_stack([c, np.ones(len(c))])
    W = np.diag(1 / e ** 2)
    cov = np.linalg.inv(X.T @ W @ X)
    beta = cov @ (X.T @ W @ a)
    chi2 = float(np.sum(((a - X @ beta) / e) ** 2))
    return (float(beta[0]), float(np.sqrt(cov[0, 0])),
            float(beta[1]), float(np.sqrt(cov[1, 1])),
            chi2, len(c) - 2)
