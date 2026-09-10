# v4 verdict — transcriptome-driven growth prediction (both approaches FAIL)

Preregistration: `docs/prereg_v4_transcriptome_models.md` (frozen before runs).
Scripts: 84 (Approach B), 85 (Approach A). All outputs in `tables/revision/v4/`.
Disposition per prereg: **archived here; neither enters the manuscript.**

## Approach B — control-trained ribosome forecast: FAIL (B1)

Fit log(growth) = a + b·s on the three control conditions (s = mean log₂
expression of the 117 KEGG map03010 ribosome genes), predict the three drug
conditions blind.

| condition | predicted | observed [95% CrI] | inside | log-resid |
|---|---|---|---|---|
| 15 °C + drug | 0.0344 | 0.0155 [0.0108, 0.0187] | NO | +0.80 |
| 21 °C + drug | 0.0361 | 0.0266 [0.0255, 0.0272] | NO | +0.31 |
| 27 °C + drug | 0.0366 | 0.0331 [0.0235, 0.0428] | yes | +0.10 |

- **B1 FAIL** (P(B1) = 0.000 over 500 joint RNA × growth bootstraps).
  B2 (attenuation direction) passes (P = 0.69), but the gate required both.
- Identical verdict under the broader 153-gene annotation-OR module
  (first run, archived in 84_log history) — the module unification
  (§7 amendment 1) did not affect the outcome.
- Leave-one-temperature-pair-out: total squared log-residual 0.63; systematic
  structure (overpredicts at 15 °C, underpredicts at 27 °C).

**Interpretation.** The strong ribosome–growth correlation across all six
conditions (r = 0.93) is NOT a law identifiable from the controls alone: control
growth is hump-shaped (0.042/0.052/0.041) while the ribosome score declines
nearly monotonically, so the control-trained slope (0.23 per log₂ unit) is far
too shallow to capture the drug effect at 15 °C. The drug removes substantially
more growth than its ribosomal repression 'buys' at cold temperature — the
correlation is carried by the drug conditions, not by a shared transcript–growth
mapping. This bounds the translation-capacity hypothesis: it remains a
consistent qualitative account of the ATTENUATION (B2 held), but it does not
quantitatively predict growth from ribosome expression.

## Approach A — E-Flux on plain ztGEM: FAIL (A1, A2, A4)

Frozen spec ran exactly as preregistered (V = feasibility-anchored at
control-15 °C; AND=min/OR=sum; one global normalization).

- **A1 FAIL:** 7 of 11 resolvable pairs discordant; Kendall τ = −0.067.
- **A2 FAIL:** predicted drug effect POSITIVE at 21 °C (δ₂₁ = +0.17);
  negative at 15 and 27 °C.
- A3 pass (δ₂₇ = −0.37 > δ₁₅ = −0.50) — vacuous given A1/A2.
- **A4 FAIL on all three prongs:** LOTPO SSR 3.14 vs ribosome-only baseline
  0.63 (5× worse); observed τ below the 95th percentile of BOTH the
  condition-label permutation null (0.60) and the degree-preserving
  gene-to-GPR scramble null (0.60) — and in fact below their medians:
  the network structure contributes no growth information beyond chance.

**Interpretation.** The post-treatment transcriptome, mapped to reaction
capacities through this GEM, does not contain a recoverable growth signal.
Combined with v3 (drug-from-stoichiometry gate failure, sign audit), the
modelling book on predicting the drug response from this data is now closed
from both ends: neither mechanism-driven (enzyme kinetics) nor data-driven
(expression-constrained) GEM variants predict the growth phenotype. The
manuscript's existing framing — condition-constrained flux analysis as an
interpretive layer with explicit caveats — is the correct and final use of the
model for this dataset.

## Manuscript consequences
1. NONE required by these failures (both analyses stay out, per prereg).
2. Module-definition correction (independent of v4, found during it): the
   sample-level interaction statistics are now computed on the unified
   117-gene KEGG map03010 module (matching Fig. 5C): F(2,39) = 12.3,
   P = 7.3 × 10⁻⁵, contrast +0.87 ± 0.38 log₂, LOO max P = 2.7 × 10⁻⁴
   (previously reported from a 153-gene annotation-OR set mislabelled
   "117"). Same conclusion; numbers corrected in manuscript.qmd.
3. OPTIONAL honesty sentence for Discussion (author's call, drafted per
   prereg): "A ribosome-expression model calibrated on the untreated
   conditions alone did not quantitatively predict growth under the drug
   (overpredicting it at 15 °C), so translational capacity indexes, but does
   not by itself quantitatively determine, the growth response."
