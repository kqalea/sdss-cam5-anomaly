"""sdss_more: exact objid parsing and sky dedupe; no network, no data files."""

import sys
from decimal import Decimal
from pathlib import Path

ANALYSIS = Path(__file__).resolve().parents[1]
if str(ANALYSIS) not in sys.path:
    sys.path.insert(0, str(ANALYSIS))

import sdss_more  # noqa: E402


def _row(ra, dec, petro=10.0, objid=1):
    return {"ra": ra, "dec": dec, "petro": petro, "objid": objid}


def test_parse_objid_plain_integer():
    assert sdss_more.parse_objid("1237672793422626838") == 1237672793422626838


def test_parse_objid_scientific_is_decimal_exact():
    text = "1.2376727934226268e18"
    got = sdss_more.parse_objid(text)
    assert got == int(Decimal(text))
    assert got == 1237672793422626800
    # Decimal is exact; a float round-trip would give ...816 instead.
    assert got != int(float(text))


def test_parse_objid_empty_is_zero():
    assert sdss_more.parse_objid("") == 0
    assert sdss_more.parse_objid("   ") == 0


def test_dedupe_drops_within_radius_and_keeps_outside():
    existing = [_row(10.0, 20.0)]
    inside = _row(10.0, 20.0 + 1.5 / 3600.0, petro=9.0, objid=2)
    outside = _row(10.0, 20.0 + 3.0 / 3600.0, petro=8.0, objid=3)
    kept = sdss_more.dedupe([inside, outside], existing)
    assert kept == [outside]


def test_dedupe_collapses_identical_rows():
    r = _row(10.0, 20.0, petro=10.0, objid=7)
    kept = sdss_more.dedupe([dict(r), dict(r)], [])
    assert len(kept) == 1
