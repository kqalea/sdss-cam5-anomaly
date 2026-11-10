"""significance.py: jackknife, cell ids, block bootstrap, binomial sigma."""

import math
import sys
from pathlib import Path

import numpy as np

ANALYSIS = Path(__file__).resolve().parents[1]
if str(ANALYSIS) not in sys.path:
    sys.path.insert(0, str(ANALYSIS))

import significance as sg  # noqa: E402


def test_jk_constant_array_gives_zero_sigma():
    sign = np.full(100, 1.0)
    cid = np.repeat(np.arange(4), 25)
    sigma, used = sg.jk_cells(sign, cid)
    assert used == 4
    assert sigma == 0.0


def test_jk_single_cell_does_not_divide_by_zero():
    sigma, used = sg.jk_cells(np.ones(100), np.zeros(100, dtype=int))
    assert used == 0
    assert math.isnan(sigma)


def test_cell_id_partitions_uniform_grid():
    nra, ndec, reps = 4, 5, 3
    ra = np.array([20.0 + 0.5 * a + 0.01 * r
                   for a in range(nra) for _ in range(ndec) for r in range(reps)])
    dec = np.array([-5.0 + 0.5 * d + 0.01 * r
                    for _ in range(nra) for d in range(ndec) for r in range(reps)])
    cid = sg.cell_id(ra, dec, nra, ndec)
    assert cid.shape == ra.shape
    assert np.issubdtype(cid.dtype, np.integer)
    assert cid.min() == 0
    assert cid.max() == nra * ndec - 1
    assert len(np.unique(cid)) == nra * ndec
    assert np.all(np.bincount(cid, minlength=nra * ndec) == reps)


def test_binomial_sigma_matches_formula():
    sign = np.array([1.0] * 60 + [-1.0] * 40)
    a = sg.amp(sign)
    assert abs(a - 0.2) < 1e-12
    expected = math.sqrt(1.0 - 0.2 ** 2) / math.sqrt(100.0)
    assert abs(sg.binom_sigma(sign) - expected) < 1e-15


def test_boot_cells_constant_blocks_give_zero():
    cid = np.repeat(np.arange(4), 25)
    assert sg.boot_cells(np.ones(100), cid) == 0.0
