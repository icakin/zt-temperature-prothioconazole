# Modelling stress tests (post-review) — results and consequences

## 1. Leave-one-temperature-out (LOTO) for the CYP51 safety-margin model
Fit to six temperatures, predict the seventh, seven times (scripts/41).

| T (°C) | EC50 obs | obs 95% CI | LOTO pred | within CI |
|---|---|---|---|---|
| 15 | 1.12 | 0.69–1.70 | 3.81 | no |
| 18 | 3.42 | 1.79–7.12 | 1.28 | no |
| 21 | 2.04 | 1.28–3.61 | 1.51 | yes |
| 24 | 1.40 | 0.99–2.02 | 1.69 | yes |
| 26 | 1.85 | 1.13–3.36 | 2.38 | yes |
| 27 | 5.00 | 3.20–11.17 | 2.77 | no |
| 28 | 2.66 | 1.64–4.59 | 3.82 | yes |

Out-of-sample log-EC50 correlation r = -0.14; 4/7 within the observed CI.
Held-out 27 °C: predicted 2.77 vs observed 5.00 (predicted 2.5-fold rise from
15 °C against 4.5-fold observed). Fitting without 15 °C collapses Ea to 0.
CONCLUSION: the three-parameter margin model INTERPOLATES the EC50-temperature
trend; it has no demonstrable out-of-sample predictive skill at a held-out
temperature. Contributing factors: only seven EC50 estimates, one clear outlier
(18 °C), and the 27 °C estimate is itself an extrapolation beyond the tested
dose range with a very wide CI (3.2-11.2). The in-sample statements (AICc
superiority over binding-only and temperature-invariant nulls; capture of the
trend) remain valid; "prediction" language is not.

## 2. AOX structural identifiability
Closing r_AOX under the published constraint set (scripts/42).

- All six conditions remain FEASIBLE with AOX closed.
- The network instead routes oxygen through cytochrome c oxidase (r_0438),
  over-produces ATP (ATP-synthase flux rises, e.g. 1.09 -> 3.42 at 27 °C + drug)
  and dissipates it; total pFBA flux rises by only 0-4.4%.
- The previously reported FVA min = max was computed with substrate uptake
  already fixed at its parsimonious optimum, which makes AOX appear uniquely
  determined; it is unique GIVEN that objective hierarchy, not topologically.
CONCLUSION: AOX is the PARSIMONIOUS route by which the network reconciles the
measured physiology, not the only feasible one. What survives, and is arguably
the more robust statement: the measured physiology (high oxygen consumption with
low growth) REQUIRES substantial energy dissipation; the model selects AOX as
the least-cost implementation, and AOX transcript abundance independently tracks
the inferred flux (rho = 0.99). The discrimination between dissipation routes
comes from the transcriptomics, not from the model alone.

## Consequences for the manuscript
Claim changes required:
- "required substantial respiratory rerouting through AOX" -> "the parsimonious
  flux solution routes substantial oxygen consumption through AOX"
- Methods FVA/uniqueness sentence -> state the closure test and its outcome
- Discussion "AOX flux is required to reconcile" -> "parsimonious reconciliation"
- Add the LOTO result (supplement) and avoid any out-of-sample framing for the
  safety-margin model
