"""The column calibration and the corrected amplitude.

Reproduces §5.5 of the paper from the frozen data:
  - the per-column offset relative to the five unaffected columns, per flag and per size bin;
  - the corrected amplitude, with a 300-draw bootstrap over galaxies that includes the
    calibration uncertainty.

Run:  python code/calibrate.py
      python code/calibrate.py --check   (asserts against data/SDSS_RESULT.json)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
BIT = 1 << 48


def load():
    more = [r for r in json.load(open(DATA / "sdss_axis_more.json")) if r["snr"] > 5]
    flags = {int(k): v for k, v in json.load(open(DATA / "flags_all.json")).items()}
    clean = {(round(r["ra"], 5), round(r["dec"], 5)): r["sign"] for r in json.load(open(DATA / "sdss_axis_clean.json"))}
    frozen = json.load(open(DATA / "SDSS_RESULT.json"))
    s, cam, tf, rk = [], [], [], []
    for r in more:
        f = flags.get(r["objid"])
        k = (round(r["ra"], 5), round(r["dec"], 5))
        if f is None or k not in clean:
            continue
        s.append(clean[k]); cam.append(r["camcol"]); tf.append(bool(f[0] & BIT)); rk.append(r["rk"])
    return np.array(s), np.array(cam), np.array(tf), np.array(rk), frozen


def main(argv):
    s, cam, tf, rk, frozen = load()
    q = np.quantile(rk, [1 / 3, 2 / 3])
    edges = [0.0, float(q[0]), float(q[1]), 200.0]

    def bin_of(x):
        return 0 if x < edges[1] else (1 if x < edges[2] else 2)

    bidx = np.array([bin_of(x) for x in rk])

    def offsets(a, c, f, b):
        off = {}
        for fl in (False, True):
            for bb in range(3):
                m = (f == fl) & (b == bb)
                refs = [a[m & (c == k)].mean() for k in [1, 2, 3, 4, 6] if (m & (c == k)).sum() > 20]
                if not refs:
                    continue
                ref = np.mean(refs)
                for k in range(1, 7):
                    mk = m & (c == k)
                    if mk.sum() >= 20:
                        off[(k, fl, bb)] = a[mk].mean() - ref
        return off

    def corrected(idx):
        a, c, f, b = s[idx], cam[idx], tf[idx], bidx[idx]
        off = offsets(a, c, f, b)
        corr = np.array([a[i] - off.get((c[i], f[i], b[i]), 0.0) for i in range(len(a))])
        return corr.mean()

    point = corrected(np.arange(len(s)))
    rng = np.random.default_rng(7)
    boot = np.array([corrected(rng.integers(0, len(s), len(s))) for _ in range(300)])
    err = boot.std()
    lo, hi = np.percentile(boot, [2.5, 97.5])

    # the per-cell table
    print("size-bin edges (rk, px):", np.round(edges, 2))
    print("\n| cam | flag | bin | A | bias | bias_err | N |")
    print("|---|---|---|---:|---:|---:|---:|")
    off = offsets(s, cam, tf, bidx)
    for (k, fl, bb), bias in sorted(off.items()):
        m = (cam == k) & (tf == fl) & (bidx == bb)
        a = s[m].mean(); n = int(m.sum())
        print("| %d | %d | %d | %+.4f | %+.4f | %.4f | %d |" % (k, int(fl), bb, a, bias, np.sqrt((1 - a * a) / n), n))

    print("\ncorrected amplitude %+.4f +/- %.4f (95%% %.4f to %.4f)" % (point, err, lo, hi))

    if "--check" in argv:
        ok = True
        cal = frozen.get("calibration", {})
        good = abs(point - cal.get("corrected_amplitude", 9)) < 5e-4 and abs(err - 0.0089) < 0.002
        ok &= good
        print("\ncheck: corrected %+.4f vs frozen %+.4f, error %.4f vs frozen %.4f -> %s" %
              (point, cal.get("corrected_amplitude", float("nan")), err, cal.get("corrected_error_bootstrap", float("nan")),
               "ok" if good else "MISMATCH"))
        sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main(sys.argv[1:])
