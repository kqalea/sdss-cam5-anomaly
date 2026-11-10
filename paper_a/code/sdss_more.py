"""More sky: the same fast strip query, offset by half a strip, larger TOP.

The 36k sample used TOP 500 per 5-degree RA strip. This fetches TOP K per strip
on a strip grid shifted by half a strip, so it reaches sky the first run never
touched, and it skips any cutout whose position is already in an existing sample.

Usage:
  python analysis/sdss_more.py --top 700 --workers 40 --out sdss_axis_more.json
"""
from __future__ import annotations

import argparse
import io
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal
from pathlib import Path

import numpy as np
import requests
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent))
import spin_core

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
SQL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"
CUT = "https://skyserver.sdss.org/dr18/SkyServerWS/ImgCutout/getjpeg"
UA = {"User-Agent": "Mozilla/5.0 (research)"}


def parse_objid(text: str) -> int:
    text = text.strip()
    if not text:
        return 0
    if "e" in text.lower() or "." in text:
        return int(Decimal(text))
    return int(text)


def query(top: int, strips: int = 72, shift: float = 2.5) -> list:
    rows = []
    for i in range(strips):
        lo = (360.0 * i / strips + shift) % 360.0
        hi = lo + 360.0 / strips
        if hi <= 360.0:
            cond = "ra BETWEEN %.4f AND %.4f" % (lo, hi)
        else:
            cond = "(ra BETWEEN %.4f AND 360 OR ra BETWEEN 0 AND %.4f)" % (lo, hi - 360.0)
        cmd = (
            "SELECT TOP %d ra, dec, petroRad_r, modelMag_r, objid, mjd, airmass_r, "
            "run, camcol, field FROM PhotoObj "
            "WHERE type = 3 AND %s AND modelMag_r < 16.5 AND petroRad_r > 4" % (top, cond)
        )
        try:
            r = requests.get(SQL, params={"cmd": cmd, "format": "csv"},
                             headers=UA, timeout=120)
        except Exception:
            continue
        if r.status_code != 200:
            continue
        for line in r.text.strip().splitlines()[1:]:
            p = line.split(",")
            if len(p) < 10:
                continue
            try:
                rows.append({"ra": float(p[0]), "dec": float(p[1]),
                             "petro": float(p[2]), "mag": float(p[3]),
                             "objid": parse_objid(p[4]),
                             "mjd": int(float(p[5])) if p[5] else 0,
                             "airmass": float(p[6]) if p[6] else float("nan"),
                             "run": int(float(p[7])) if p[7] else 0,
                             "camcol": int(float(p[8])) if p[8] else 0,
                             "fieldno": int(float(p[9])) if p[9] else 0})
            except ValueError:
                continue
        if (i + 1) % 12 == 0:
            print("  queried %d/%d strips, %d rows" % (i + 1, strips, len(rows)), flush=True)
    return rows


def unit(ra, dec):
    return (np.cos(np.radians(dec)) * np.cos(np.radians(ra)),
            np.cos(np.radians(dec)) * np.sin(np.radians(ra)),
            np.sin(np.radians(dec)))


def key_of(ra, dec):
    x, y, z = unit(ra, dec)
    cell = 0.5 / 3600.0
    return (int(np.floor(x / cell)), int(np.floor(y / cell)), int(np.floor(z / cell))), (x, y, z)


def dedupe(rows, existing, radius_arcsec=2.0):
    cell = radius_arcsec / 3600.0
    grid = {}
    kept = []

    def insert(r):
        k, (x, y, z) = key_of(r["ra"], r["dec"])
        grid.setdefault(k, []).append((x, y, z))

    for r in existing:
        insert(r)
    for r in sorted(rows, key=lambda r: (-r["petro"], r["objid"])):
        k, (x, y, z) = key_of(r["ra"], r["dec"])
        hit = False
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for ux, uy, uz in grid.get((k[0] + dx, k[1] + dy, k[2] + dz), ()):
                        if ux * x + uy * y + uz * z > np.cos(np.radians(radius_arcsec / 3600.0)):
                            hit = True
                            break
                    if hit:
                        break
                if hit:
                    break
            if hit:
                break
        if not hit:
            grid.setdefault(k, []).append((x, y, z))
            kept.append(r)
    return kept


def one(t):
    for attempt in range(3):
        try:
            pars = {"ra": "%.6f" % t["ra"], "dec": "%.6f" % t["dec"],
                    "scale": "0.396", "width": "128", "height": "128"}
            r = requests.get(CUT, params=pars, headers=UA, timeout=45)
            if r.status_code != 200 or len(r.content) < 500:
                return None
            im = Image.open(io.BytesIO(r.content)).convert("L")
            a = np.nan_to_num(np.asarray(im, dtype=float), nan=0.0)
            if a.shape[0] < 64:
                return None
            rk = float(np.clip(t["petro"] / 0.396, 6.0, 60.0))
            d = spin_core.measure(a, 1.0, 0.0, rk)
            if d is None:
                return None
            return {"ra": t["ra"], "dec": t["dec"], "mag": t["mag"],
                    "sign": int(d["sign"]), "snr": float(d["snr"]), "rk": rk,
                    "objid": t["objid"], "mjd": t["mjd"], "airmass": t["airmass"],
                    "run": t["run"], "camcol": t["camcol"], "fieldno": t["fieldno"]}
        except Exception:
            time.sleep(0.5 + attempt)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=700)
    ap.add_argument("--strips", type=int, default=72)
    ap.add_argument("--workers", type=int, default=40)
    ap.add_argument("--out", default="sdss_axis_more.json")
    args = ap.parse_args()

    existing = []
    for name in ("sdss_axis_36k.json", "sdss_axis_more.json"):
        p = DATA / name
        if p.exists():
            existing += json.load(open(p))
    print("existing rows: %d" % len(existing), flush=True)

    t0 = time.time()
    rows = query(args.top, args.strips)
    print("queried %d rows in %.0f s" % (len(rows), time.time() - t0), flush=True)
    uniq = dedupe(rows, existing)
    print("new galaxies after dedupe: %d (from %d rows)" % (len(uniq), len(rows)), flush=True)

    out = []
    t1 = time.time()
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for i, d in enumerate(pool.map(one, uniq), 1):
            if d is not None:
                out.append(d)
            if i % 2000 == 0:
                el = time.time() - t1
                print("  %d/%d measured, %d usable, %.1f/s, eta %.0f min" %
                      (i, len(uniq), len(out), i / max(el, 1), (len(uniq) - i) / max(i / max(el, 1), 0.01) / 60),
                      flush=True)
    path = DATA / args.out
    json.dump(out, open(path, "w"))
    print("wrote %s: %d rows, %.0f s total" % (path, len(out), time.time() - t0))


if __name__ == "__main__":
    main()
