# Kerr-Cartan — Paper A replication package

**A camera-column dependence in SDSS handedness measurements: localization, calibration and a conditional null**
Yen-Fu Chen (Jecision Digital Co., Ltd., New Taipei City, Taiwan)

## What this is

The measurement code, the frozen data and the figures for the paper. The paper's chain: the SDSS
handedness amplitude is a camera-column dependence (camcol 5 carries +0.1209 against -0.0063 in the
other five), concentrated in TOO_FEW_GOOD_DETECTIONS objects (+0.1787 against -0.0301), localized to
the r-band detector of camcol 5 (the colour decomposition gives 6.6 sigma in the r channel; raw
r-band frames give 2.5 sigma on a 159-galaxy pilot), calibrated per column, flag and size bin; after
the correction the amplitude is -0.0116 +/- 0.017, a conditional null.

## Reproduce

```sh
pip install -r paper_a/requirements.txt
cd paper_a
make check      # per-column and flag tables, asserted against data/SDSS_RESULT.json
make calibrate  # the 36-cell calibration and the corrected amplitude
make test       # 18 unit tests
make figures    # the five figures
```

Pinned: Python 3.15.2, numpy 2.5.4, scipy 1.18.1, matplotlib 3.11.2, astropy 8.0.2, pillow 13.0.0,
requests 2.34.2.

## Contents

```
paper_a/            code, data, figures, tests, Makefile
paper_a_draft.md    the paper source (Markdown)
paper_a.tex         the typeset source
paper_a.pdf         the compiled paper
paper_a_zh.md       the Chinese reading version
CITATION.cff        citation metadata
LICENSE             code MIT, data CC BY 4.0
```

## Cite

If you use this code or the measurements, cite the paper and this repository. The DOI is pending.

## Licence

Code: MIT. Data tables: CC BY 4.0 (derived from public SDSS data).
