"""Recompute the paper's canonical numbers from the frozen data.

Run:  python code/analyse.py          (prints the tables)
       python code/analyse.py --check  (also asserts against data/SDSS_RESULT.json)

The data are the two per-galaxy measurement files and the PhotoObj flag bits. Nothing here
queries the network: the inputs are frozen so that the paper's numbers are reproducible by
anyone with this directory.
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
    clean = json.load(open(DATA / "sdss_axis_clean.json"))
    more = [r for r in json.load(open(DATA / "sdss_axis_more.json")) if r["snr"] > 5]
    flags = {int(k): v for k, v in json.load(open(DATA / "flags_all.json")).items()}
    frozen = json.load(open(DATA / "SDSS_RESULT.json"))
    return clean, more, flags, frozen


def main(argv):
    clean, more, flags, frozen = load()
    s = np.array([r["sign"] for r in clean], float)
    cam = np.array([r["camcol"] for r in clean])
    print("N = %d galaxies, A = %+.4f" % (len(s), s.mean()))

    print("\nPer-column amplitude (Table 1)")
    print("| camcol | N | A |")
    print("|---|---:|---:|")
    for c in range(1, 7):
        m = cam == c
        print("| %d | %d | %+.4f |" % (c, m.sum(), s[m].mean()))
    print("| not 5 | %d | %+.4f |" % ((cam != 5).sum(), s[cam != 5].mean()))

    clean_map = {(round(r["ra"], 5), round(r["dec"], 5)): r["sign"] for r in clean}
    rows = []
    for r in more:
        f = flags.get(r["objid"])
        k = (round(r["ra"], 5), round(r["dec"], 5))
        if f is None or k not in clean_map:
            continue
        rows.append((clean_map[k], r["camcol"], bool(f[0] & BIT)))
    ss = np.array([x[0] for x in rows]); cc = np.array([x[1] for x in rows]); tf = np.array([x[2] for x in rows])
    print("\nDetection selection (Table 2, independent-sky sample N = %d)" % len(rows))
    print("| subset | N | A cam5 | A other |")
    print("|---|---:|---:|---:|")
    for name, m in [("flag = 0", ~tf), ("flag = 1", tf)]:
        m5 = m & (cc == 5); mo = m & (cc != 5)
        print("| %s | %d | %+.4f | %+.4f |" % (name, m.sum(), ss[m5].mean(), ss[mo].mean()))
    print("flag = 1 fraction: cam5 %.3f, other %.3f" % (tf[cc == 5].mean(), tf[cc != 5].mean()))

    if "--check" in argv:
        ok = True
        def close(name, got, want, tol=5e-4):
            nonlocal ok
            good = abs(got - want) <= tol
            ok &= good
            print("%-28s got %+.4f  want %+.4f  %s" % (name, got, want, "ok" if good else "MISMATCH"))
        close("overall", s.mean(), frozen["samples"]["combined"]["A"])
        for c in range(1, 7):
            close("camcol %d" % c, s[cam == c].mean(), frozen["columns_combined"][str(c)]["A"])
        close("not 5", s[cam != 5].mean(), frozen["columns_combined"]["not_5"]["A"])
        sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main(sys.argv[1:])
