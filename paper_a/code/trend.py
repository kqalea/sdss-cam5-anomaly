"""The trend: amplitude and the camera columns against sample size.

Every sample is collapsed to unique galaxies (single linkage at 1.2 arcsec), because
the rows repeat: 15,952 rows were 9,674 galaxies. The amplitude is reported both as
the brightest member's sign and as the cluster mean, with the delete-one-cell
jackknife over two grids.

Usage: python analysis/trend.py data/sdss_axis_24k.json data/sdss_axis_36k.json data/sdss_axis_more.json
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict

import numpy as np

import spin_core  # noqa: F401  (kept so the import path matches the pipeline)


def load(path):
    rows = [r for r in json.load(open(path)) if r["snr"] > 5]
    return rows


def clusters(rows, link_arcsec=1.2):
    n = len(rows)
    ra = np.array([r["ra"] for r in rows])
    dec = np.array([r["dec"] for r in rows])
    u = np.stack([np.cos(np.radians(dec)) * np.cos(np.radians(ra)),
                  np.cos(np.radians(dec)) * np.sin(np.radians(ra)),
                  np.sin(np.radians(dec))], 1)
    cell = np.radians(0.05)
    g = np.floor(u / cell).astype(np.int64)
    b = defaultdict(list)
    for i in range(n):
        b[(g[i, 0], g[i, 1], g[i, 2])].append(i)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for key, ids in b.items():
        neigh = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    neigh += b.get((key[0] + dx, key[1] + dy, key[2] + dz), [])
        for ii in ids:
            for jj in neigh:
                if jj <= ii:
                    continue
                if float(u[ii] @ u[jj]) > np.cos(np.radians(link_arcsec / 3600.0)):
                    union(ii, jj)
    groups = defaultdict(list)
    for i in range(n):
        groups[find(i)].append(i)
    out = []
    for ix in groups.values():
        s = np.mean([rows[i]["sign"] for i in ix])
        j = min(ix, key=lambda i: rows[i]["mag"])
        out.append(dict(s=float(s), s_best=float(rows[j]["sign"]), ra=rows[j]["ra"], dec=rows[j]["dec"],
                        cc=rows[j]["camcol"], k=len(ix)))
    return out


def jk(vals, ra, dec, nra, ndec):
    qa = np.unique(np.quantile(ra, np.linspace(0, 1, nra + 1)))
    qd = np.unique(np.quantile(dec, np.linspace(0, 1, ndec + 1)))
    ia = np.clip(np.searchsorted(qa, ra, side="right") - 1, 0, len(qa) - 2)
    idd = np.clip(np.searchsorted(qd, dec, side="right") - 1, 0, len(qd) - 2)
    cid = ia * 100 + idd
    n = len(vals)
    tot = vals.sum()
    th = []
    for c in np.unique(cid):
        m = cid == c
        nc = m.sum()
        if nc >= n:
            continue
        th.append((tot - vals[m].sum()) / (n - nc))
    th = np.array(th)
    k = len(th)
    return float(np.sqrt((k - 1) / k * np.sum((th - th.mean()) ** 2)))


def report(path):
    rows = load(path)
    cl = clusters(rows)
    s = np.array([c["s"] for c in cl])
    sb = np.array([c["s_best"] for c in cl])
    ra = np.array([c["ra"] for c in cl])
    dec = np.array([c["dec"] for c in cl])
    cc = np.array([c["cc"] for c in cl])
    sig = jk(s, ra, dec, 8, 6)
    sigb = jk(sb, ra, dec, 8, 6)
    print("\n== %s" % path)
    print("rows %d -> galaxies %d (%.2f rows per galaxy)" % (len(rows), len(cl), len(rows) / len(cl)))
    print("A cluster-mean %+.4f +/- %.4f (%.2f)   A brightest %+.4f +/- %.4f (%.2f)" %
          (s.mean(), sig, s.mean() / sig, sb.mean(), sigb, sb.mean() / sigb))
    line = "camcol:"
    for c in range(1, 7):
        m = cc == c
        line += "  %d N=%d A=%+.3f" % (c, m.sum(), s[m].mean())
    print(line)
    m = cc != 5
    print("without camcol 5: N=%d A=%+.4f" % (m.sum(), s[m].mean()))
    return dict(path=path, n_rows=len(rows), n_gal=len(cl), A=float(s.mean()), sig=float(sig),
                A_best=float(sb.mean()), sig_best=float(sigb), cc={c: float(s[cc == c].mean()) for c in range(1, 7)},
                n_cc={c: int((cc == c).sum()) for c in range(1, 7)})


if __name__ == "__main__":
    out = [report(p) for p in sys.argv[1:]]
    json.dump(out, open("data/trend.json", "w"), indent=1)
    print("\nwrote data/trend.json")
