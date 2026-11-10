#!/usr/bin/env python
"""Significance audit of the SDSS handedness amplitude.

The amplitude has been quoted at 4.7 sigma, then 5.2, then 3.4. Each number was
computed on a different cutout or with a different error model. This script
makes the error model explicit and reproducible:

  * the amplitude and the binomial error;
  * delete-one-cell jackknife over several sky grids (RA x Dec, RA strips,
    Dec strips) at a range of cell counts;
  * delete-one-run and delete-one-camera-column jackknife;
  * a cell block bootstrap (weights by cell size), which does not rely on the
    equal-influence assumption of the jackknife;
  * the cell chi-squared against the binomial;
  * the camera-column decomposition and the within-sky-cell column contrast;
  * the per-run column contrast, with a sign test across runs;
  * (with two files) the paired sampling difference on the objects both files
    contain.

Usage:
  python analysis/significance.py data/sdss_axis_36k.json \
      --pair data/sdss_axis_s02_256.json
"""
import argparse
import json

import numpy as np


def load(path, snr_cut=5.0):
    rows = json.load(open(path))
    rows = [r for r in rows if r["snr"] > snr_cut]
    ra = np.array([r["ra"] for r in rows])
    dec = np.array([r["dec"] for r in rows])
    sign = np.array([r["sign"] for r in rows], dtype=float)
    camcol = np.array([int(r["camcol"]) for r in rows])
    run = np.array([int(r["run"]) for r in rows])
    mag = np.array([r["mag"] for r in rows])
    rk = np.array([r["rk"] for r in rows])
    return dict(ra=ra, dec=dec, sign=sign, camcol=camcol, run=run,
                mag=mag, rk=rk, n=len(sign), key=[(round(a, 6), round(d, 6))
                for a, d in zip(ra, dec)])


def amp(sign):
    return float(np.mean(sign))


def binom_sigma(sign):
    n = len(sign)
    a = amp(sign)
    return float(np.sqrt(max(1.0 - a * a, 0.0) / n))


def quantile_bins(x, nb):
    """Equal-count bin edges, with repeated edges removed."""
    qs = np.quantile(x, np.linspace(0, 1, nb + 1))
    edges = np.unique(qs)
    return np.clip(np.searchsorted(edges, x, side="right") - 1, 0, len(edges) - 2)


def cell_id(ra, dec, nra, ndec):
    ib = quantile_bins(ra, nra)
    jb = quantile_bins(dec, ndec)
    return ib * ndec + jb


def jk_cells(sign, cid):
    """Delete-one-cell jackknife of the mean. Returns (sigma, n_used)."""
    ids, counts = np.unique(cid, return_counts=True)
    n = len(sign)
    total = sign.sum()
    theta = []
    used = 0
    for c in ids:
        m = cid == c
        nc = m.sum()
        if nc == 0 or nc >= n:
            continue
        theta.append((total - sign[m].sum()) / (n - nc))
        used += 1
    theta = np.array(theta)
    k = len(theta)
    if k < 3:
        return float("nan"), used
    # standard delete-one jackknife variance
    var = (k - 1) / k * np.sum((theta - theta.mean()) ** 2)
    return float(np.sqrt(var)), used


def boot_cells(sign, cid, nboot=4000, seed=7):
    """Block bootstrap that resamples sky cells with replacement."""
    ids = np.unique(cid)
    means = []
    sizes = []
    for c in ids:
        m = cid == c
        if m.sum() == 0:
            continue
        means.append(sign[m].mean())
        sizes.append(m.sum())
    means = np.array(means)
    sizes = np.array(sizes, dtype=float)
    rng = np.random.default_rng(seed)
    k = len(means)
    draws = rng.integers(0, k, size=(nboot, k))
    w = sizes[draws]
    vals = (means[draws] * w).sum(axis=1) / w.sum(axis=1)
    return float(vals.std(ddof=1))


def cell_chi2(sign, cid):
    """Chi-squared of cell means against the binomial null."""
    ids = np.unique(cid)
    chi2 = 0.0
    dof = 0
    for c in ids:
        m = cid == c
        nc = m.sum()
        if nc < 2:
            continue
        ac = sign[m].mean()
        chi2 += (ac ** 2) * nc
        dof += 1
    return chi2, max(dof - 1, 1)


def report(name, d):
    s = d["sign"]
    n = d["n"]
    a = amp(s)
    bs = binom_sigma(s)
    print(f"\n## {name}")
    print(f"N = {n:,}   A = {a:+.4f}   binomial sigma = {bs:.4f}   "
          f"binomial significance = {a / bs:.2f}")
    print("\n| grid | cells | used | jackknife sigma | significance | boot sigma | boot significance | chi2/dof |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|")
    grids = [
        ("RA x Dec 4x4", 4, 4), ("RA x Dec 6x6", 6, 6),
        ("RA x Dec 8x6", 8, 6), ("RA x Dec 10x8", 10, 8),
        ("RA x Dec 14x10", 14, 10), ("RA x Dec 20x14", 20, 14),
    ]
    for label, nra, ndec in grids:
        cid = cell_id(d["ra"], d["dec"], nra, ndec)
        sig, used = jk_cells(s, cid)
        bt = boot_cells(s, cid)
        chi2, dof = cell_chi2(s, cid)
        print(f"| {label} | {nra * ndec} | {used} | {sig:.4f} | {a / sig:.2f} | "
              f"{bt:.4f} | {a / bt:.2f} | {chi2:.1f}/{dof} |")
    # one-dimensional strips
    for label, cid, ncell in [
        ("RA strips 12", quantile_bins(d["ra"], 12), 12),
        ("RA strips 24", quantile_bins(d["ra"], 24), 24),
        ("Dec strips 8", quantile_bins(d["dec"], 8), 8),
        ("Dec strips 16", quantile_bins(d["dec"], 16), 16),
    ]:
        sig, used = jk_cells(s, cid)
        bt = boot_cells(s, cid)
        chi2, dof = cell_chi2(s, cid)
        print(f"| {label} | {ncell} | {used} | {sig:.4f} | {a / sig:.2f} | "
              f"{bt:.4f} | {a / bt:.2f} | {chi2:.1f}/{dof} |")
    # run and camcol
    for label, col in [("runs", d["run"]), ("camcols", d["camcol"])]:
        sig, used = jk_cells(s, col)
        bt = boot_cells(s, col)
        chi2, dof = cell_chi2(s, col)
        print(f"| delete-one {label} | {used} | {used} | {sig:.4f} | {a / sig:.2f} | "
              f"{bt:.4f} | {a / bt:.2f} | {chi2:.1f}/{dof} |")
    return a, bs


def camcol_table(d):
    print("\n### Camera columns")
    print("| camcol | N | A | binomial sigma | sigma | delete-one-jackknife sigma |")
    print("|---:|---:|---:|---:|---:|---:|")
    s = d["sign"]
    for c in sorted(np.unique(d["camcol"])):
        m = d["camcol"] == c
        a = amp(s[m])
        bs = binom_sigma(s[m])
        sig, _ = jk_cells(s[m], cell_id(d["ra"][m], d["dec"][m], 6, 6))
        print(f"| {c} | {m.sum():,} | {a:+.4f} | {bs:.4f} | {a / bs:+.2f} | {sig:.4f} |")
    for c in [5]:
        m = d["camcol"] != c
        a = amp(s[m])
        bs = binom_sigma(s[m])
        sig, used = jk_cells(s[m], cell_id(d["ra"][m], d["dec"][m], 6, 6))
        print(f"| without {c} | {m.sum():,} | {a:+.4f} | {bs:.4f} | {a / bs:+.2f} | {sig:.4f} |")


def within_cell_contrast(d, cid, a_col=5):
    """Amplitude difference between camcol a_col and the rest, inside sky cells.

    Pooled over cells with both groups present; the error is the delete-one-cell
    jackknife over those cells.
    """
    s = d["sign"]
    sel = d["camcol"] == a_col
    ids = np.unique(cid)
    diffs = []
    weights = []
    for c in ids:
        m = cid == c
        n1 = (m & sel).sum()
        n0 = (m & ~sel).sum()
        if n1 < 10 or n0 < 10:
            continue
        diffs.append(s[m & sel].mean() - s[m & ~sel].mean())
        weights.append(min(n1, n0))
    diffs = np.array(diffs)
    weights = np.array(weights, dtype=float)
    if len(diffs) < 4:
        return float("nan"), float("nan"), 0
    est = float(np.average(diffs, weights=weights))
    # weighted mean, delete-one jackknife
    theta = []
    for i in range(len(diffs)):
        keep = np.arange(len(diffs)) != i
        theta.append(np.average(diffs[keep], weights=weights[keep]))
    theta = np.array(theta)
    k = len(theta)
    var = (k - 1) / k * np.sum((theta - theta.mean()) ** 2)
    return est, float(np.sqrt(var)), k


def run_contrast(d, a_col=5, min_n=60):
    """Per-run camcol-a_col contrast, with a sign test across runs."""
    s = d["sign"]
    sel = d["camcol"] == a_col
    runs = np.unique(d["run"])
    rows = []
    for r in runs:
        m = d["run"] == r
        n1 = (m & sel).sum()
        n0 = (m & ~sel).sum()
        if n1 < min_n or n0 < min_n:
            continue
        rows.append((int(r), int(n1), int(n0), s[m & sel].mean(), s[m & ~sel].mean()))
    if not rows:
        return [], 0
    wins = sum(1 for x in rows if x[3] > x[4])
    return rows, wins


def paired(d1, d2, label1="A", label2="B", nra=8, ndec=6):
    """Paired comparison on objects present in both files (same ra/dec key)."""
    key1 = {k: v for k, v in zip(d1["key"], d1["sign"])}
    common = [(k, v) for k, v in zip(d2["key"], d2["sign"]) if k in key1]
    print(f"\n## Paired sampling comparison ({label1} vs {label2})")
    print(f"common objects: {len(common):,} of {d2['n']:,} in {label2} and {d1['n']:,} in {label1}")
    if not common:
        return
    idx = {k: i for i, k in enumerate(d2["key"])}
    keep2 = np.array([idx[k] for k, _ in common])
    s1 = np.array([key1[k] for k, _ in common])
    s2 = np.array([v for _, v in common], dtype=float)
    delta = s2 - s1
    same = float(np.mean(s1 == s2))
    print(f"sign agreement = {same * 100:.1f} per cent")
    print(f"A({label1}) = {s1.mean():+.4f}   A({label2}) = {s2.mean():+.4f}   "
          f"difference = {delta.mean():+.4f}")
    # jackknife on the paired difference over sky cells of the common set
    cid = cell_id(d2["ra"][keep2], d2["dec"][keep2], nra, ndec)
    sig, used = jk_cells(delta, cid)
    bs = boot_cells(delta, cid)
    print(f"paired difference jackknife sigma = {sig:.4f} ({used} cells), "
          f"bootstrap {bs:.4f}, significance of a nonzero shift = {delta.mean() / sig:.2f}")
    print(f"  (an independent-samples comparison would give sigma = "
          f"{np.sqrt(d2['sign'].var() / d2['n'] + d1['sign'].var() / d1['n']):.4f})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sample")
    ap.add_argument("--pair", default=None)
    ap.add_argument("--snr", type=float, default=5.0)
    args = ap.parse_args()

    d = load(args.sample, args.snr)
    a, bs = report(args.sample, d)
    camcol_table(d)

    print("\n### Camera-column contrast within sky cells")
    for label, nra, ndec in [("RA x Dec 6x6", 6, 6), ("RA x Dec 10x8", 10, 8),
                             ("RA x Dec 14x10", 14, 10)]:
        cid = cell_id(d["ra"], d["dec"], nra, ndec)
        est, sig, k = within_cell_contrast(d, cid, 5)
        print(f"| {label} | {k} cells with both | difference (camcol 5 - rest) = "
              f"{est:+.4f} +/- {sig:.4f} | {est / sig:.2f} sigma |")

    rows, wins = run_contrast(d, 5)
    if rows:
        print(f"\n### Per-run camcol-5 contrast (runs with >= 60 objects in both groups): "
              f"{len(rows)} runs, camcol 5 higher in {wins}")
        print("| run | N5 | Nother | A5 | Aother | diff |")
        print("|---:|---:|---:|---:|---:|---:|")
        for r, n1, n0, a1, a0 in rows[:20]:
            print(f"| {r} | {n1} | {n0} | {a1:+.3f} | {a0:+.3f} | {a1 - a0:+.3f} |")

    if args.pair:
        d2 = load(args.pair, args.snr)
        paired(d, d2, label1="native 0.396", label2="0.2/256")


if __name__ == "__main__":
    main()
