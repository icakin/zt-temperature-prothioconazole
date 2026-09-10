# Panel A — un-anchored thermal optimum: methods note

**Claim.** The growth thermal optimum and magnitude of *Z. tritici* are recovered from sequence
and external reference data, with nothing fitted to the *Z. tritici* growth curve.

**Peak position (optimum).** Per-enzyme catalytic optima were not anchored to the measured
optimum. Instead each enzyme's Topt was set from its sequence-predicted melting temperature
(TemStaPro) minus a single enzyme Topt–Tm gap **transferred from an independent organism**: the
GECKO *E. coli* model (eciML1515), whose enzyme optima (Li–Engqvist predictor) and melting
temperatures (Leuenberger meltome) give a median gap of **12.4 °C** (bootstrap 95% CI 12.0–12.7,
n = 1368). Applied to *Z. tritici*'s sequence Tm, this predicts the organism optimum at **24 °C**,
matching the measured value.

*Why a constant gap, not a regression.* The *E. coli* Topt–Tm regression is weak (r = 0.64) and
*Z. tritici*'s Tm lies below the *E. coli* training range, so the fitted slope extrapolates
unreliably; the constant median gap is the parsimonious, robust transfer statistic and is
consistent with the physical expectation that the catalysis-vs-denaturation gap is approximately
Tm-independent.

*Tm-scale consistency.* *E. coli* (OGT 37 °C, meltome Tm 55.6 °C) and *Z. tritici* (OGT 24 °C,
TemStaPro Tm 36.7 °C) imply dTm/dOGT = 1.45 °C/°C, within the literature range (~1–2 °C per °C
OGT), indicating the TemStaPro pseudo-Tm sits on a scale consistent with the measured meltome.
(A direct calibration — TemStaPro run on the meltome proteins — remains future work.)

**Magnitude.** Absolute turnover was set by rescaling the DLTKcat predictions (median 3.3 s⁻¹) to
the empirical enzyme-kcat distribution (Bar-Even et al. 2011; central-metabolism median ~10–79
s⁻¹), an external correction of the predictor's documented low bias — not a fit to growth. This
recovers the growth magnitude to ~1.5×; the residual reflects the class-dependent width of the
reference distribution and remaining predictor bias.

**Shape and collapse.** The curve's breadth and its high-temperature collapse emerge from the
sequence-predicted Tm distribution; the in-vivo inactivation offset is the one thermal parameter
still calibrated to the control curve.

**Uncertainty.** The transferred gap is tightly determined (±0.4 °C → ±0.4 °C on the peak,
formally), but the honest outer uncertainty on the predicted optimum is the **±~6 °C spread across
independent sequence-based predictors** (DLTKcat catalytic optima 18 °C; Zeldovich IVYWREL 16 °C;
Tome OGT 28.5 °C; *E. coli*-transferred gap 24 °C). The measured 24 °C falls within this envelope;
the *E. coli*-consistent transfer is the point estimate.

**Honest limitations.** (1) The constant-gap transfer and the TemStaPro/meltome scale equivalence
are assumptions that materially affect the point estimate; the direct Tm-scale calibration is not
yet done. (2) The magnitude carries a ~1.5× residual. (3) The optimum is predicted within the
method envelope, not to <2 °C. The anchored model remains available as the more precise (but
non-predictive) alternative.

**Sources.** Bar-Even et al. 2011, *Biochemistry* 50:4402; Zeldovich et al. 2007, *PLoS Comput
Biol* 3:e5; Li, Engqvist et al. 2019 (Tome), *ACS Synth Biol* 8:1411; Leuenberger et al. 2017
meltome; enzyme Topt–Tm meta-analysis, *Bioscience Reports* 2021, 41:BSR20210336.
