# Model robustness package (post-review analyses)

Three analyses strengthening the GEM section, run after manuscript completion.
Scripts: scripts/37-40. Outputs: tables/gem/surface_*, medium_sensitivity.csv,
anchor_tpc_shift*.csv; figures/FigureS12.{png,pdf}.

## 1. Dose x temperature surface (37)
The network + Arrhenius CYP51 capacity + competitive inhibition (Ea, Ki from
the EC50 fit; thermal baseline from the control TPC only) predicts the full
relative-growth surface: r = 0.83 vs the Hill-characterised surface over the
49 treated conditions. Because Ki and Ea were fitted to the seven EC50s, this
is a test of surface SHAPE, not EC50 placement. Residual structure is
informative: the model reproduces cool synergy (5/6 confident Bliss signs at
15 C) and high-dose antagonism, but predicts no interaction at supra-optimal
low doses (model delta ~ 0 at 27-28 C, 0.06-0.5 mg/L) where the data show
antagonism. That low-dose supra-optimal antagonism is precisely the component
attributable to the transcriptional convergence mechanism (Figs 3-5) rather
than CYP51 kinetics: enzyme-level constraints capture part of the interaction,
and their failure localises where the regulatory mechanism operates.

## 2. Anchoring sensitivity (39)
Enzyme catalytic optima re-anchored at 20/22/24/26 C: the emergent curve's
peak tracks the anchor 1:1 (as expected; the optimum is anchored), while
CTmax (36.7-36.8 C), rmax (0.216/h) and the collapse shape are invariant.
The high-temperature failure and magnitude are genuinely sequence-derived.
Rising-limb Ea varies modestly (0.41-0.57 eV across +/-4 K).

## 3. Medium sensitivity (38)
Condition-constrained AOX shares recomputed with 38 amino-acid/nucleobase
uptakes opened (0.5 mmol/gDW/h each) as a rich-YMS proxy. The qualitative
pattern survives and strengthens: cool controls remain 0% AOX; drug and
supra-optimal conditions increase (15_2: 14->35%, 21_2: 25->58%, 27_0: 4->43%,
27_2: 42->72%). The sucrose-minimal medium used in the paper is therefore the
CONSERVATIVE choice for the AOX conclusion. P/O pattern unchanged. Note: the
published FVA table (ztGEM_v04_required_AOX_fitted_maintenance.csv) shows
min = max for AOX in all six conditions under the published constraint set,
confirming the Methods uniqueness statement.
