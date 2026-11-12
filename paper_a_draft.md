# Paper A (final source)

**Title: A camera-column dependence in SDSS handedness measurements: localization, calibration and a conditional null**

## Abstract

The claimed handedness asymmetry of SDSS spiral galaxies has been debated for fifteen years. We
measure the handedness of 26,727 unique galaxies with one operator, in two disjoint sky samples and
without the coordinate-grid overlay. The amplitude is carried by one camera column: camcol 5 gives
\(+0.1209\) over 5,815 galaxies while the other five give \(-0.0063\). The contrast against the
five-column mean is 7.0 sigma (camcol 5 alone is 7.4 sigma from zero); rescaled by the observed
column overdispersion the contrast is 3.3 sigma. The excess is confined to
galaxies flagged TOO_FEW_GOOD_DETECTIONS: in the independent-sky sample of 15,456 galaxies (58 per
cent of the catalogue) camcol 5 gives \(+0.1787\) against \(-0.0301\) in the other five columns,
while the well-detected objects of that sample (8,182 galaxies) show no column difference and an
amplitude of \(+0.0071 \pm 0.0112\). We treat the column as an instrument with its own zero
point: the other five columns provide the reference, the per-column offset is measured as a
function of the detection flag and the sampling radius (camcol 5 about \(+0.20\) in the
band-poor subset), and after the correction the amplitude of the 15,456-galaxy sample is
\(-0.0116 \pm 0.017\), consistent with zero. No significant net handedness is detected at this
sensitivity;
the camera-column systematic that would otherwise dominate such a measurement is measured,
calibrated and published. Channel decomposition and raw r-band frame measurements localize the anomaly to extended-source
imaging in the r-band detector of camera column 5, although the sub-pixel transfer mechanism
within that detector remains open.

## 1. Introduction

A net handedness in the spins of galaxies is one of the few observables that would connect the
large-scale geometry of the Universe to the physics of structure formation. If galaxies inherit a
preferred sense of rotation from a primordial vorticity field, the excess of one winding over the
other is a direct fossil of that field. The claim is simple to state and hard to establish: the
effect is expected to be small, the measurement is qualitative, and every step between the
primordial field and a pixel in a survey image can introduce a preference of its own.

The observational history is a study in that difficulty. An excess of anticlockwise over clockwise
spirals in the Sloan Digital Sky Survey was reported by Longo (2011) and developed into a series of
claims of a preferred axis and a redshift evolution (Shamir 2012, 2017, 2024).
Countervailing analyses found that the signal responds to the mirror symmetry of the
classification, to the sample selection, and to the imaging pipeline (Land et al. 2008; Hayes et
al. 2017; Iye et al. 2020). The reported amplitudes span a 33-fold range, from about 0.006
to 0.20, and the disagreement is not about the sign of the effect but about which samples and
which pipelines are believed. A recent catalogue-based analysis of 8.5 million DESI Legacy galaxies
found no dipole and a monopole at the \(10^{-3}\) level whose origin the authors attribute to
imaging systematics or sky structure (the DESI chirality release, 2026).

This paper measures where the SDSS amplitude lives and what happens to it when the camera column
is treated as an instrument. We use one operator, 26,727 unique galaxies, two disjoint sky
samples, a fixed pixel scale and no coordinate-grid overlay. The operator is the circular m=2 phase gradient: the RGB cutout is reduced to its luminance
(ITU-R 601), the image is averaged in annuli, the azimuthal mean is removed, the phase of the m=2
Fourier mode is followed with radius after unwrapping the \\(\\pi\\) ambiguity, and the sign of its
least-squares gradient against \\(\\ln r\\) is the handedness. It uses only the pixels, never the
shape model.

Three findings follow. First, the amplitude is carried by one camera column: camcol 5 gives
\(+0.1209\) over 5,815 galaxies while the other five give \(-0.0063\), with a contrast of 7.0 sigma
against the five-column mean (camcol 5 alone is 7.4 sigma from zero; rescaled by the column
overdispersion the contrast is 3.3 sigma). Second, the excess is
confined to galaxies flagged TOO_FEW_GOOD_DETECTIONS: in the independent-sky sample camcol 5 gives
\(+0.1787\) against \(-0.0301\) in the band-poor subsample, and the well-detected objects show no
column difference and \(+0.0071 \pm 0.0112\). Third, the five non-camcol-5 columns provide a
reference: the camcol-5 band-poor offset is measured at about \(+0.20\) as a function of detection
quality and size, and after the correction the sample amplitude is \(-0.0116 \pm 0.0089\),
consistent with zero.

The result has a corollary. The selection that creates the excess, poor detection, is the same
selection that keeps those objects out of the bright, resolved samples on which volunteer
classifications are built; a defect confined to poorly detected objects is therefore invisible to
those samples by construction. This paper concerns the low-quality end of the SDSS catalogue; it
does not adjudicate the bright-sample analyses, which have no access to these objects.

The paper is organised as follows. Section 2 describes the data. Section 3 defines the operator.
Section 4 presents the null tests. Section 5 gives the results and the calibration. Section 6
collects the systematic tests, Section 7 discusses the implications, and Section 8 concludes.

## 2. Data

The galaxies are drawn from the SDSS DR18 PhotoObj table with `type = 3`, `modelMag_r < 16.5`
and `petroRad_r > 4` arcsec. Two samples are used. The first consists of 72 RA strips of 5 degrees
with TOP 500 per strip. The second is an independent sky: 72 strips shifted by half a strip width,
TOP 700 per strip, with every position within 2 arcsec of the first sample dropped before the
cutouts are fetched.

Rows are not galaxies. The discovery sample contains 35,971 rows for 11,271 galaxies and the
independent sample 35,387 rows for 15,456 galaxies; the same sky is imaged in overlapping runs, so
the same galaxy enters the catalogue several times. Duplicates are merged by single linkage at
1.2 arcsec. The combined catalogue has 26,727 unique galaxies. The discovery sample gives
\(A = +0.0294\) (3.1 sigma), the independent sample \(+0.0155\) (1.9 sigma), and the combined
catalogue \(+0.0214 \pm 0.0061\) (binomial; the sky-cell jackknife gives 0.0075, 2.8 sigma).

Cutouts are 128 by 128 pixels at 0.396 arcsec per pixel, a 50.7 arcsec field, from the SkyServer
ImgCutout service. The sampling radius is \(r_k = \mathrm{clip}(\mathrm{petroRad}_r/0.396, 6, 60)\)
pixels. **The SkyServer `opt=G` option draws a coordinate grid**, not grayscale; the discovery
sample was first fetched with it, and the clean refetch of both samples changes 10.8 per cent of the
signs and lowers the camcol-5 amplitude from 0.1391 to 0.1209.

The six columns are not equally populated (camcol 5 has 5,815 galaxies, the others 3,891 to
4,476), an artefact of the column-major ordering of `objid` against the TOP cut. Equalising the
columns to 3,891 galaxies each leaves the camcol-5 amplitude unchanged (+0.117, +0.128, +0.128 in
three draws).

## 3. The operator

The handedness is measured by the circular m=2 phase gradient. For each galaxy the cutout is
sampled in annuli from 1 to \(r_k\) pixels; the azimuthal mean is removed at every radius; the m=2
Fourier coefficients are computed in each annulus; the phase of the mode is followed with
\(\ln r\); the sign of the phase gradient is the handedness. An amplitude signal-to-noise of 5 is
required. The operator never uses the axis ratio, the position angle or the shape model. The operator
and the full analysis are released with this paper (Section 9).

Three checks fix it. Rotation preserves the sign and reflection reverses it (unit-tested). The
same galaxies measured at 0.2 arcsec per pixel agree with the 0.396 arcsec measurement at 91.2 per
cent. On a sample of 254 galaxies with visual classifications the attenuation of the operator
relative to the visual sample is 0.31 to 0.37.

## 4. Null tests

**Repeat observations.** SDSS images the same sky in adjacent runs. Matching positions within
1.2 arcsec yields 21,464 pairs in the full measurement; the two measurements agree 0.845 of the
time for camcol 5 and 0.855 outside it. In an auxiliary bright test catalogue of 36,000 raw cutouts (17,347 repeat pairs before
deduplication) the agreement reaches 0.931; this is the sample used in the earlier version of the
analysis.

**Volunteer classifications.** A Galaxy Zoo sample in the same sky (9,986 volunteer labels; 8,891 with an operator measurement)
gives the volunteer-label amplitude \(-0.027 \) (2.6 sigma); our operator on those images gives \(-0.004\). The operator
reproduces the volunteer labels on the bright sample at 0.664 to 0.72 depending on the sampling
radius.

**Reflection and rotation.** A reflected image returns the opposite sign and a rotated image the
same sign; the test is run on every batch.

## 5. Results

### 5.1 The amplitude against sample size
Discovery sky \(+0.0294\) over 11,271 galaxies (3.1 sigma); independent sky \(+0.0155\) over
15,456 (1.9 sigma); combined \(+0.0214 \pm 0.0061\) (3.5 sigma). The two samples are disjoint in
objects but drawn from the same runs and nights, so they test the stability of the measurement
across objects, not the independence of the camera.

### 5.2 The camera column
| camcol | N | A | jackknife sigma |
|---|---:|---:|---:|
| 1 | 4,476 | +0.0018 | 0.0162 |
| 2 | 4,375 | -0.0487 | 0.0177 |
| 3 | 4,128 | -0.0334 | 0.0163 |
| 4 | 4,042 | +0.0079 | 0.0162 |
| **5** | **5,815** | **+0.1209** | 0.0164 |
| 6 | 3,891 | +0.0460 | 0.0181 |

The five other columns scatter with \(\chi^2 = 17.9\) for four degrees of freedom; their empirical
scatter is 0.037, larger than the jackknife error. Against the galaxy-weighted mean of the five,
\(-0.0063 \pm 0.0080\), camcol 5 is 7.0 sigma; rescaled for the overdispersion it is about
3.3 sigma. Both numbers are quoted.

### 5.3 The detection selection
In the independent-sky sample (15,456 galaxies, 58 per cent of the catalogue) camcol 5 contains
57.1 per cent of objects flagged TOO_FEW_GOOD_DETECTIONS against 44.3 per cent elsewhere, and the
excess is concentrated in that subset (the unflagged camcol-5 residual is +0.0327 +/- 0.026).

| subset | N | A cam5 | A other five |
|---|---:|---:|---:|
| flag = 0 | 8,182 | +0.0327 | +0.0016 |
| flag = 1 | 7,274 | +0.1787 | -0.0301 |

The flag is not a uniform amplifier: per column in flag = 1 the amplitudes are
\(-0.029, -0.077, -0.033, -0.062, +0.179, +0.058\). The setting is an interaction between the
column and the detection quality. A within-run comparison (camcol 5 against the other columns of
the same run) gives \(+0.144 \pm 0.032\); this localizes the effect to the column rather than to
the night or the weather, and is not an independent detection of handedness.

### 5.4 The well-detected subset
In the independent-sky sample the six columns are \(+0.026, -0.020, -0.026, +0.011, +0.033, +0.022\), consistent within
\(\pm 0.03\), and the overall amplitude is
\[ A = +0.0071 \pm 0.0112 \qquad (0.63\ \sigma). \]
The two-sided 95 per cent interval is \([-0.015, +0.029]\); amplitudes above +0.03 are excluded, but
a small positive signal of 0.006 to 0.007 is not. The bound refers to the pixel measurement; after
the 0.31 to 0.37 operator dilution of Section 3 it corresponds to an intrinsic sky asymmetry below
about 0.08 to 0.10. It is consistent with
a small intrinsic asymmetry of 0.006 to 0.007 and with zero; it is a bound, not a detection.

### 5.5 The column calibration
The five non-camcol-5 columns provide the reference. The offset is measured per column, per detection
flag, and in three sampling-radius bins (edges at the terciles of \(r_k\)), with the per-cell
binomial error.

| cam | flag = 1 bias (small / middle / large) |
|---|---|
| 1 | +0.013 / +0.017 / -0.032 |
| 2 | -0.038 / +0.022 / -0.114 |
| 3 | -0.018 / +0.141 / -0.093 |
| 4 | -0.024 / -0.126 / +0.037 |
| **5** | **+0.266 +/- 0.034 / +0.229 +/- 0.048 / +0.114 +/- 0.038** |
| 6 | +0.067 / -0.053 / +0.200 |

The flag = 0 cells used in the same correction are not uniformly positive for camcol 5:
(-0.063 +/- 0.059, +0.003 +/- 0.039, +0.118 +/- 0.045) in the three bins; the offset is therefore
a band-poor phenomenon.

The camcol-5 band-poor cells are consistently positive; the fifteen non-camcol-5 cells scatter at
1.3 to 2.1 times their nominal errors (chi2 = 32.7 for 14 degrees of freedom, p = 0.003), so the
calibration reference is itself not clean. Applying the offset to each galaxy leaves a sample amplitude of
**-0.0116 +/- 0.017** on the 15,456-galaxy sample. The statistical error (0.0089, a 300-draw
bootstrap) is combined in quadrature with the intrinsic between-column SEM (0.0149, from the observed
scatter 0.0372 minus the measurement noise),
because the five reference columns are themselves inconsistent at chi2 = 17.9 for four degrees of
freedom; the 95 per cent interval is [-0.046, +0.023]. The in-sample correction is not independent: the offsets are estimated
from the same galaxies that are then corrected, so the corrected mean is close to the reference by
construction. We therefore validate it out of sample: the independent-sky sample is split in two
halves by RA, the offsets are estimated on one half and applied to the other, giving
-0.0034 +/- 0.0116 and -0.0215 +/- 0.0117, both consistent with zero. The corrected amplitude is the sample-weighted mean of the five-column reference: the camcol-5
galaxies contribute weight but not sign information, so the effective sign-information sample is
12,098 galaxies. The null is conditional on the instrumental interpretation; a signal with the
camcol-5 x flag x size signature is removed by construction. The numbers are frozen in
data/SDSS_RESULT.json. The correction is a zero-point adjustment; it does not change the per-galaxy
signs.

## 6. Systematics
Each test below excludes one explanation of the column offset; the measurements are recorded in the
project notes. The coordinate grid changes 10.8 per cent of signs (camcol 5 falls 0.1391 to
0.1209). Run and night are excluded by the within-run contrast between camcol 5 and the other columns of the
same run (+0.1453 +/- 0.0285, 5.1 sigma; removing run means leaves +0.1189). Inclination is matched by
`expAB_r` (medians 0.569 against 0.566), although the effect is stronger in face-on galaxies and
this test is therefore weak. The recorded selection is matched on magnitude, \(r_k\), S/N, Dec and
RA, giving a contrast of \(+0.115 \pm 0.028\) (4.1 sigma, grid signs) and \(+0.165 \pm 0.046\)
(3.6 sigma) in the four-quantile version. Frame attribution through the SAS WCS reproduces the
cutout sign at 0.90. The sample over-representation is removed by equalising the columns. The
injection test (300 injections per group; Fig. 5) recovers the injected winding at 0.86 in the camcol-5
band-poor group, 0.85 elsewhere in flag = 1, 0.77 in the camcol-5 well-detected group and 0.73
elsewhere, with no significant difference between the +k and -k injections (largest group asymmetry
0.04); a positive sign-recovery bias above about 0.06 is excluded at 95 per cent for the camcol-5
band-poor group.
The photometric scale (a global gain) and the JPEG quantisation (recompression from q90 to q40)
leave the measured amplitude unchanged.

**Band and detector localization.** To distinguish cross-band registration in the three-band
rendering from a single detector, the clean cutouts are decomposed into colour channels over 1,000
matched flag=1 pairs. The contrast is absent in the B channel (g band, \(-0.005 \pm 0.044\)),
marginal in R (i band, \(+0.088 \pm 0.045\), 2.0 sigma), and carried by the G channel (r band:
camcol 5 +0.238 against -0.055 elsewhere, contrast \(+0.293 \pm 0.044\), 6.6 sigma). The luminance
combination dilutes it to \(+0.138 \pm 0.045\) (3.1 sigma) by averaging in the unperturbed g
channel. A pilot measurement of the raw r-band calibrated frames (159 galaxies, 79 in camcol 5)
confirms the signal outside the JPEG: after correcting the parity inversion between the FITS array
convention and the north-up cutouts (the per-galaxy sign agreement is 0.70 after the flip), the
unstretched r-band images give \(+0.291\) in camcol 5 against \(-0.100\) elsewhere, a contrast of
\(+0.391 \pm 0.156\) (2.5 sigma). The 669 field stars show no significant asymmetry (-0.09 +/- 0.13 in
camcol 5 against +0.12 +/- 0.17 elsewhere; difference -0.21 +/- 0.21), a weak null that does not
exclude a PSF effect at radii larger than the 1-8 pixels used. The anomaly is therefore
localized to extended-source imaging in the r-band detector of camera column 5 (r5), not to the
PSF, to the JPEG rendering or to cross-band registration.



## 7. Discussion

**The published amplitude spread.** A sample that includes band-poor detections from a camera
column with a detection-related offset reports an inflated positive amplitude; a bright, resolved
selection (such as the volunteer catalogues) excludes that regime by construction and recovers a
near-zero monopole. The 33-fold spread of the published values are consistent
with differences in selection, before any question about the sky is asked.

**Low-quality and bright samples are different populations.** The excess lives in objects that the
bright, resolved samples do not contain; this paper therefore localizes a systematic in the
low-quality end of the catalogue and does not settle the bright-sample analyses.

**The open question.** Channel decomposition and raw FITS measurements localize the excess to
extended-source imaging in a single detector (r5), while excluding the stellar PSF, the JPEG
compression and cross-band registration. Why the TOO_FEW_GOOD_DETECTIONS regime amplifies the phase
gradient specifically within the r5 detector, whether through a threshold-entry selection bias or a
sub-pixel charge-transfer or readout asymmetry at low surface brightness, is not resolved by the
survey products alone.

## 8. Conclusions

1. One operator, two disjoint sky samples, 26,727 galaxies, no coordinate grid.
2. The handedness amplitude is carried by one camera column (contrast 7.0 sigma; 3.3 sigma rescaled for the
   column overdispersion), and the excess is confined to galaxies flagged
   TOO_FEW_GOOD_DETECTIONS.
3. The well-detected objects show no column difference and +0.0071 +/- 0.0112: a bound of
   |A| < 0.030 at two sigma, not a detection.
4. The column offsets are calibrated per column, per flag and per size bin; the camcol-5
   band-poor offset is about +0.20.
5. After the correction the amplitude is -0.0116 +/- 0.017: no significant net handedness is
   detected at this sensitivity.
6. The mechanism behind the column offset is not identified; the calibration and the null are
   published so that the measurement can be repeated on data that carry neither.
## 9. Data and code availability

This research received no external funding.

The analysis code and the derived measurements are released as a self-contained package: the
operator, the sample construction and the single-linkage deduplication, the jackknife and bootstrap
errors, the column calibration, the figure scripts and eighteen unit tests, together with the frozen
per-galaxy measurements and the calibration table. Four Makefile targets reproduce the paper:
\`make check\` recomputes the per-column and flag tables, \`make calibrate\` recomputes the column
calibration and the corrected amplitude, \`make test\` runs the unit tests and \`make figures\`
regenerates the figures. The package is available from the public repository
(https://github.com/kqalea/sdss-cam5-anomaly) and is archived on Zenodo at https://doi.org/10.5281/zenodo.24150030 (doi:10.5281/zenodo.24150030).

