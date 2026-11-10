"""trend.clusters: single linkage at 1.2 arcsec and the cluster sign."""

import sys
from pathlib import Path

ANALYSIS = Path(__file__).resolve().parents[1]
if str(ANALYSIS) not in sys.path:
    sys.path.insert(0, str(ANALYSIS))

import trend  # noqa: E402

ONE_ARCSEC = 1.0 / 3600.0  # degrees


def _row(ra, dec, sign=1.0, mag=16.0, camcol=1):
    return {"ra": ra, "dec": dec, "sign": sign, "mag": mag, "camcol": camcol}


def test_one_arcsec_pair_merges():
    cl = trend.clusters([_row(10.0, 0.0), _row(10.0, ONE_ARCSEC)])
    assert len(cl) == 1
    assert cl[0]["k"] == 2


def test_two_arcsec_pair_stays_separate():
    cl = trend.clusters([_row(10.0, 0.0), _row(10.0, 2.0 * ONE_ARCSEC)])
    assert len(cl) == 2
    assert sorted(c["k"] for c in cl) == [1, 1]


def test_chain_of_three_merges():
    rows = [_row(10.0, i * ONE_ARCSEC) for i in range(3)]
    cl = trend.clusters(rows)
    assert len(cl) == 1
    assert cl[0]["k"] == 3


def test_cluster_sign_is_member_mean():
    rows = [_row(10.0, 0.0, 1.0, 16.0),
            _row(10.0, ONE_ARCSEC, 1.0, 17.0),
            _row(10.0, 2.0 * ONE_ARCSEC, -1.0, 18.0)]
    cl = trend.clusters(rows)
    assert len(cl) == 1
    assert cl[0]["k"] == 3
    assert abs(cl[0]["s"] - 1.0 / 3.0) < 1e-12
    assert cl[0]["s_best"] == 1.0  # brightest member, mag 16
