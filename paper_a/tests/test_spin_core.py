"""spin_core.measure on synthetic two-armed spirals; no data files."""

import sys
from pathlib import Path

import numpy as np

ANALYSIS = Path(__file__).resolve().parents[1]
if str(ANALYSIS) not in sys.path:
    sys.path.insert(0, str(ANALYSIS))

import spin_core  # noqa: E402

SIZE = 256
RK = 44.0
# The cut the pipeline applies after measure(): significance.load(snr_cut=5.0) and
# trend.load(snr > 5). measure() itself has no threshold.
SNR_CUT = 5.0


def _spiral(k=4.0):
    """1 + 0.6 cos(2 theta + k ln r) about the image centre.

    For this image the m=2 phase is k ln r, so measure()'s slope is k and its
    sign is sign(k); mirroring theta -> -theta flips both.
    """
    yy, xx = np.mgrid[0:SIZE, 0:SIZE].astype(float)
    cx = cy = SIZE / 2.0
    r = np.hypot(xx - cx, yy - cy)
    theta = np.arctan2(yy - cy, xx - cx)
    return 1.0 + 0.6 * np.cos(2.0 * theta + k * np.log(np.maximum(r, 1e-6)))


def test_winding_direction_matches_sign():
    for k in (4.0, -4.0):
        d = spin_core.measure(_spiral(k), 1.0, 0.0, RK)
        assert d is not None
        assert d["sign"] == int(np.sign(k))
        assert d["slope"] * k > 0.0
        assert d["snr"] > SNR_CUT


def test_mirror_reverses_sign():
    img = _spiral(4.0)
    d = spin_core.measure(img, 1.0, 0.0, RK)
    assert d is not None
    for mirrored in (np.fliplr(img), np.flipud(img)):
        dm = spin_core.measure(mirrored, 1.0, 0.0, RK)
        assert dm is not None
        assert dm["sign"] == -d["sign"]


def test_180_rotation_preserves_sign():
    """A 180-degree rotation cannot flip an m=2 handedness.

    theta -> theta + pi leaves cos(2 theta) unchanged, so a two-armed spiral is
    invariant and the measured sign is preserved. Only a reflection reverses
    chirality (see test_mirror_reverses_sign). The brief expected a change here;
    the actual, mathematically correct behaviour is asserted instead.
    """
    img = _spiral(4.0)
    d = spin_core.measure(img, 1.0, 0.0, RK)
    dr = spin_core.measure(np.rot90(img, 2), 1.0, 0.0, RK)
    assert d is not None and dr is not None
    assert dr["sign"] == d["sign"]
    assert abs(dr["slope"] - d["slope"]) < 0.05 * abs(d["slope"])


def test_flat_and_noise_stay_below_snr_cut():
    assert spin_core.measure(np.zeros((SIZE, SIZE)), 1.0, 0.0, RK) is None
    flat = np.full((SIZE, SIZE), 100.0)
    d = spin_core.measure(flat, 1.0, 0.0, RK)
    # A flat frame at a positive level still yields a dict: rounding noise in the
    # interpolation passes the relative m=2 mask, so the sign is arbitrary. The
    # result must be rejected by the caller's S/N cut.
    assert d is None or d["snr"] < SNR_CUT
    rng = np.random.default_rng(3)
    dn = spin_core.measure(rng.normal(0.0, 1.0, (SIZE, SIZE)), 1.0, 0.0, RK)
    assert dn is None or dn["snr"] < SNR_CUT
