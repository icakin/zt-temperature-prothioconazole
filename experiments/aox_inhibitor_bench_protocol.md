# AOX inhibitor experiment — bench protocol (potency-restoration design)

Narrative companion to `AOX_inhibitor_plate_layout.xlsx` (the xlsx is the
authoritative well-by-well sheet). Power figures from
`scripts/87_power_aox_experiment.py` (vial-level SD of log growth = 0.171, ~17% CV).

## The headline question
Does blocking the alternative oxidase (AOX) **restore prothioconazole potency at
supra-optimal temperature**? Finding the efficacy loss (already in the paper) and
then *reversing* it pharmacologically is the high-impact result — mechanism plus a
translational lever (azole + AOX-inhibitor co-application for warm conditions).

## Reagents and roles
- **n-propyl gallate (nPG)** — AOX inhibitor; PRIMARY probe (direct
  *Zymoseptoria/Mycosphaerella* precedent).
- **N,2-dihydroxybenzamide (SHAM)** — AOX inhibitor, second chemotype; CONFIRMATION.
- **sodium azide** — cytochrome-oxidase inhibitor; ACUTE respiration/validation
  ONLY, never in a growth vial.
nPG and SHAM are two chemotypes, not independent target-specific probes; their
convergence is pharmacological support, not gene-level proof.

## PRIMARY read-out — EC50 rescue
Prothioconazole full dose–response (Hill curve) at **15 °C and 27 °C, each ± nPG**.
Signature that makes it a headline: **EC50 falls sharply with nPG at 27 °C but
barely moves at 15 °C** — a temperature-specific rescue, i.e. an
inhibitor × temperature interaction on log EC50. Confirm the 27 °C rescue with SHAM.

## SECONDARY read-outs (same plates, free)
- **Growth difference-in-differences** per temperature:
  I_T = [log r(drug+nPG) − log r(drug)] − [log r(nPG) − log r(0)] — subtracts the
  inhibitor's own effect. Powered for a material (~15–20%) effect, not equivalence.
- **Respiration:** is the drug-induced *extra* respiration nPG/SHAM-sensitive?
  (the respiratory-paradox mechanism.)

## Design / vial counts (see xlsx for the map)
- **Pilot** (1 plate, 27 °C + 2 mg/L drug): nPG titration 0/0.05/0.1/0.25/0.5/1 mM
  + matched cell-free blanks; SHAM 0.5/1/2 mM + blanks; acute azide validation
  (vehicle / azide / azide+nPG / azide+SHAM) + blanks = 24 wells. Fixes the nPG
  working dose (lowest dose that blocks with no control-growth cost and clean
  blanks) and proves blockade (azide-resistant respiration collapses under both
  nPG and SHAM).
- **Main run** = 2 plates (15 °C, 27 °C). Per plate: 8 prothioconazole doses
  (0–4 mg/L) × ±nPG = 16 culture vials + 8 cell-free blanks = 24. One run = one
  independent culture-start.
- **SHAM confirmation:** one plate at 27 °C, same layout, SHAM for nPG.

## Golden rule (every vial)
Same total solvent: 0.4% DMSO = 20 µL/5 mL, split 10 µL drug channel + 10 µL
inhibitor channel; back-fill vehicle so no-drug and no-inhibitor vials still carry
their 10 µL. All stocks are 500× (10 µL → 1× in 5 mL); recipes in the xlsx.

## Replication (decisive)
Split-culture blocking: each independent inoculum feeds ALL arms; culture/day =
biological replicate, vial = technical. ~8–12 independent culture-starts per arm
for a clean EC50 shift; ~24/arm to detect a 15% growth change, ~14/arm for 20%.
Confirming the model's <5% near-neutrality would need ~67/arm — not attempted; the
growth arm is framed as **falsifying** a material effect, never confirming
equivalence. Start 3–4 runs, read the effect size, extend.

## Controls & autoxidation
Cell-free blanks on every plate (medium ± drug ± inhibitor, all combinations used;
ideally spent-medium+inhibitor). nPG/SHAM burn O₂ chemically → do NOT subtract one
mean blank slope; fit a matched background term dO/dt = −q(O₂)·B(t) − a_j(O,T) per
dose, and reject any dose whose abiotic O₂ loss rivals the biological effect.
Matched DMSO in every vial; temperature-specific sensor calibration; dark handling
(nPG photo-oxidises); solubility check at 15 °C; an OD/dry-weight biomass check
independent of DO if at all possible.

## Claim ceiling
Permitted: "AOX-inhibitor-sensitive potency/respiration, strongest at supra-optimal
temperature; co-applied AOX inhibition restores prothioconazole potency at 27 °C."
NOT permitted from this design: exact AOX flux shares; "azide-resistant = AOX by
definition"; AOX-*gene* causality (needs WT vs ΔAOX vs complemented — the decisive
future design). Watch the known gallate–azole interaction: read the rescue as
AOX-attributable only if SHAM reproduces it.
