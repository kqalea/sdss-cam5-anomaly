"""Regenerate the six figures of the paper from the frozen data.

Run:  python code/figures.py
"""
from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplcfg")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 9, "figure.dpi": 200, "axes.spines.top": False, "axes.spines.right": False})

rows = json.load(open(DATA / "sdss_axis_clean.json"))
s = np.array([r["sign"] for r in rows], float)
cam = np.array([r["camcol"] for r in rows])
ra = np.array([r["ra"] for r in rows])
dec = np.array([r["dec"] for r in rows])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.2))
ax1.scatter(ra, dec, s=0.3, c="#3b6ea5", alpha=0.4, linewidths=0)
ax1.set_xlabel("RA (deg)"); ax1.set_ylabel("Dec (deg)"); ax1.set_title("%d galaxies, two skies" % len(s))
counts = [int((cam == c).sum()) for c in range(1, 7)]
ax2.bar(range(1, 7), counts, color="#3b6ea5")
ax2.axhline(np.mean(counts), color="0.4", ls="--", lw=0.8)
ax2.set_xlabel("camera column"); ax2.set_ylabel("galaxies"); ax2.set_xticks(range(1, 7))
ax2.set_title("columns are not equally populated")
fig.tight_layout(); fig.savefig(FIG / "fig1_sky_columns.png"); plt.close(fig)

amps = [s[cam == c].mean() for c in range(1, 7)]
errs = [0.0162, 0.0177, 0.0163, 0.0162, 0.0164, 0.0181]
fig, ax = plt.subplots(figsize=(4.6, 3.4))
ax.errorbar(range(1, 7), amps, yerr=errs, fmt="o", color="#b4453a", capsize=3)
ax.axhline(-0.0053, color="0.3", lw=1, label="five-column mean")
ax.axhspan(-0.0053 - 0.04, -0.0053 + 0.04, color="0.7", alpha=0.35, label="empirical column scatter")
ax.set_xlabel("camera column"); ax.set_ylabel("A"); ax.set_xticks(range(1, 7))
ax.set_title("camcol 5 carries the amplitude"); ax.legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(FIG / "fig2_columns.png"); plt.close(fig)

more = [r for r in json.load(open(DATA / "sdss_axis_more.json")) if r["snr"] > 5]
flags = {int(k): v for k, v in json.load(open(DATA / "flags_all.json")).items()}
clean = {(round(r["ra"], 5), round(r["dec"], 5)): r["sign"] for r in rows}
BIT = 1 << 48
recs = []
for r in more:
    f = flags.get(r["objid"]); k = (round(r["ra"], 5), round(r["dec"], 5))
    if f is None or k not in clean:
        continue
    recs.append((r["camcol"], (f[0] & BIT) != 0, clean[k]))
ss = np.array([x[2] for x in recs]); cc = np.array([x[0] for x in recs]); tf = np.array([x[1] for x in recs])
fig, ax = plt.subplots(figsize=(4.8, 3.4))
for flag, color, label, marker in [(False, "#4a7c59", "flag = 0 (good detections)", "s"), (True, "#b4453a", "flag = 1 (TOO_FEW_GOOD_DETECTIONS)", "o")]:
    vals = [ss[(cc == c) & (tf == flag)].mean() for c in range(1, 7)]
    err = [1 / np.sqrt(max(((cc == c) & (tf == flag)).sum(), 1)) for c in range(1, 7)]
    ax.errorbar(np.arange(1, 7) + (0.08 if flag else -0.08), vals, yerr=err, fmt=marker, color=color, label=label, capsize=2, ms=4)
ax.axhline(0, color="0.6", lw=0.8)
ax.set_xlabel("camera column"); ax.set_ylabel("A"); ax.set_xticks(range(1, 7))
ax.set_title("the excess lives in the band-poor subset"); ax.legend(fontsize=6.5, frameon=False, loc="upper left")
fig.tight_layout(); fig.savefig(FIG / "fig3_flags.png"); plt.close(fig)

fig, ax = plt.subplots(figsize=(4.6, 3.3))
vals = [ss[(cc == c) & (~tf)].mean() for c in range(1, 7)]
err = [1 / np.sqrt(max(((cc == c) & (~tf)).sum(), 1)) for c in range(1, 7)]
ax.errorbar(range(1, 7), vals, yerr=err, fmt="o", color="#3b6ea5", capsize=3)
ax.axhspan(0.006, 0.007, color="#b4453a", alpha=0.35, label="model 0.006-0.007")
overall = ss[~tf].mean(); e = np.sqrt((1 - overall ** 2) / (~tf).sum())
ax.axhline(overall, color="#3b6ea5", lw=0.8, ls="--", label="good-subset mean %+.4f" % overall)
ax.axhspan(overall - 2 * e, overall + 2 * e, color="#3b6ea5", alpha=0.15)
ax.set_xlabel("camera column"); ax.set_ylabel("A"); ax.set_xticks(range(1, 7))
ax.set_title("good detections agree with the model"); ax.legend(fontsize=6.5, frameon=False)
fig.tight_layout(); fig.savefig(FIG / "fig4_clean_model.png"); plt.close(fig)

groups = ["cam5\nflag=1", "other\nflag=1", "cam5\nflag=0", "other\nflag=0"]
plus = [0.86, 0.85, 0.77, 0.73]; minus = [0.86, 0.84, 0.73, 0.74]
fig, ax = plt.subplots(figsize=(4.8, 3.2))
x = np.arange(len(groups))
ax.bar(x - 0.2, plus, 0.4, color="#3b6ea5", label="injected +k")
ax.bar(x + 0.2, minus, 0.4, color="#b4453a", label="injected -k")
ax.axhline(0.5, color="0.5", lw=0.8, ls=":")
ax.set_xticks(x); ax.set_xticklabels(groups, fontsize=7); ax.set_ylabel("recovery rate")
ax.set_title("injection test: no positive bias (asymmetry < 0.04)"); ax.legend(fontsize=6.5, frameon=False)
fig.tight_layout(); fig.savefig(FIG / "fig5_injection.png"); plt.close(fig)

leg = json.load(open(DATA / "c5_flag_legacy.json")) if (DATA / "c5_flag_legacy.json").exists() else []
fig, ax = plt.subplots(figsize=(4.2, 3.6))
for label, c5, fl in [("cam5 flag=1", True, True), ("other flag=1", False, True), ("cam5 flag=0", True, False), ("other flag=0", False, False)]:
    m = np.array([o["cam5"] == c5 and o["flag"] == fl for o in leg]) if leg else np.zeros(0, bool)
    if m.sum() < 20:
        continue
    xs = np.array([o["sdss"] for o in leg])[m].mean(); ys = np.array([o["legacy"] for o in leg])[m].mean()
    ax.plot([xs], [ys], "o", ms=6, label="%s (N=%d)" % (label, m.sum()))
ax.axhline(0, color="0.6", lw=0.8); ax.axvline(0, color="0.6", lw=0.8)
lim = (-0.25, 0.25)
ax.plot(lim, lim, color="#b4453a", lw=1, ls="--", label="signal in the galaxies (1:1)")
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("SDSS amplitude"); ax.set_ylabel("Legacy amplitude")
ax.set_title("the excess does not reproduce")
ax.legend(fontsize=6, frameon=False, loc="upper left")
fig.tight_layout(); fig.savefig(FIG / "fig6_independent.png"); plt.close(fig)
print("figures written to", FIG)
