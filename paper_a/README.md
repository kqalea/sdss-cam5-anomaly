# Paper A: code and data

The measurement code and the frozen data for the paper *"A camera-column dependence in SDSS
handedness measurements: localization, calibration and a conditional null"*.

The paper measures the handedness of 26,727 unique SDSS galaxies with the circular m=2 phase
gradient. This directory recomputes the headline numbers and regenerates the five figures without network
access; the raw-cutout fetch itself is in code/sdss_more.py.

## What is here

```
code/
  spin_core.py        the one implementation of the measurement (m=2 phase gradient)
  deprojected_spin.py elliptical sampling used by spin_core
  significance.py     jackknife and error helpers
  trend.py            sample construction helpers
  sdss_more.py        objid parsing and de-duplication
  analyse.py          recomputes Table 1, Table 2 and the canonical numbers
  calibrate.py        recomputes the per-column calibration and the corrected amplitude
  figures.py          regenerates the five figures
tests/
  run_tests.py        18 unit tests (rotation invariance, reflection reversal, de-duplication, ...)
data/
  sdss_axis_clean.json   26,727 galaxies: ra, dec, sign, snr, rk, camcol, sign_grid
  sdss_axis_more.json    the independent-sky sample (15,456 galaxies)
  flags_all.json         PhotoObj flags per objid; bit 48 is TOO_FEW_GOOD_DETECTIONS
  SDSS_RESULT.json       the frozen canonical numbers, including the calibration block
figures/                 the five published figures
```

## Reproduce

```
make check      # recompute the per-column and flag tables, assert against SDSS_RESULT.json
make calibrate  # recompute the column calibration and the corrected amplitude (-0.0116 +/- 0.017)
make test       # run the 18 unit tests
make figures    # regenerate the five figures
```

`make check` prints the per-column table, the detection-selection table, and an ok/MISMATCH
line for every frozen number. It exits non-zero on any mismatch.

## Data provenance

- SDSS DR18 `PhotoObj`: `type = 3`, `modelMag_r < 16.5`, `petroRad_r > 4` arcsec.
- Cutouts: SkyServer `ImgCutout`, 128 px at 0.396 arcsec/px, **without** `opt=G` (which draws a
  coordinate grid, not grayscale).
- The r-band frame pilot uses the SAS `frame-r-*.fits.bz2` products; the channel decomposition
  splits the same cutouts into R, G and B.

## Note on the flags

`TOO_FEW_GOOD_DETECTIONS` is bit 48 of the PhotoObj `flags` mask. The flag split in the paper
is evaluated on the independent-sky sample (15,456 galaxies), which is the only sample for which
the flag bits were fetched.

## Licence

Code: MIT. Data tables: CC BY 4.0 (derived from public SDSS data).
