# Claude Code prompt — native sector/calibrated figures + run calibrated & sectors (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`).

NOTE TO THE USER: to run this unattended, launch Claude Code in an auto-approving
mode (accept edits + allow commands, e.g. the "accept all" permission mode or
`--dangerously-skip-permissions`). The prompt below instructs the agent not to pause,
but tool-permission prompts are controlled by the launch mode, not the prompt text.

---

```
Work AUTONOMOUSLY from start to finish. Do not pause for confirmation, do not ask
for approval between steps, make reasonable decisions yourself, and commit your own
changes. Only stop if a step hard-fails in a way you cannot resolve, and then report
exactly what blocked you and what you tried. Print a concise summary at the end.
Do not change core scientific/numerical code (mmrt, enzyme_cost cost math,
providers, tpc, sensitivity, decomposition, control); the figure work is additive.

Read src/etcgem/plotting.py, sensitivity.py (SensitivityResult: it holds `samples`
= input DataFrame and `descriptors` = descriptor DataFrame), sectors.py, cli.py and
strains/eciML1515/strain.yaml + experiments/{calibrated,sectors}.yaml first.

STEP 1 - add two first-class figures to src/etcgem/plotting.py
(a) plot_descriptor_intervals(result, out_dir, fname="descriptor_intervals.png"):
    the calibrated-uncertainty summary. For the key descriptors present in
    result.descriptors (Topt_C, rmax, CTmax_C, niche_width_C, Ea_eV), draw a
    horizontal point-interval plot: ensemble median with the 2.5-97.5 percentile
    interval per descriptor (one row each; normalise/scale per descriptor or use
    separate panels so they are readable together). This is the headline view of
    calibrated uncertainty and works for any sweep.
(b) plot_sector_tradeoff(result, out_dir, fname="sector_tradeoff.png"): only when
    "f_metab" is a column of result.samples (else return None and skip). Join
    result.samples[["f_metab","f_maint"]] with result.descriptors by row and
    scatter rmax vs f_metab coloured by f_maint (add a second panel for CTmax_C vs
    f_metab). This shows the interior growth optimum from metabolic<->biosynthesis
    co-limitation.
Wire both into plot_all(result, out_dir): always call descriptor_intervals; call
sector_tradeoff only when f_metab is present (guard with a hasattr/column check).
These are NEW files only; existing figures (tpc_ensemble, sensitivity_heatmap,
descriptor_distributions) must be byte-for-byte unchanged. Quick-test on the toy
strain: `etcgem sweep --strain _toy --experiment quick` then confirm
descriptor_intervals.png is produced (sector_tradeoff.png absent, since no f_metab).
Commit: "add descriptor-interval and sector-tradeoff figures".

STEP 2 - run the calibrated experiment WITH figures (M1.2)
Run (do NOT pass --no-plots): `etcgem sweep --strain eciML1515 --experiment calibrated`
Confirm it writes strains/eciML1515/outputs/sweep_calibrated/ containing
tpc_ensemble.png, descriptor_distributions.png, descriptor_intervals.png,
descriptors.csv, sensitivity_spearman.csv, summary.json, and resolved_config.yaml
showing envelope_sampling mode (correlated or posterior). If the experiment is set
to mode: posterior, confirm it read strains/eciML1515/dltkcat/fits.csv (per-enzyme
Topt_sd/dCp_sd). This run may take a few minutes; let it finish.

STEP 3 - run the sectors experiment WITH figures, WITHOUT changing the baseline
Sectors must be active only for THIS run, so the standard-pool baseline (default,
calibrated, decompose, control, dltkcat) is unchanged. In order of preference:
(a) Add a proteome_sectors override to experiments/sectors.yaml:
      proteome_sectors: {enabled: true, atpm_reaction: null}
    Then run `etcgem sweep --strain eciML1515 --experiment sectors` and CHECK the
    run's resolved_config.yaml shows proteome_sectors.enabled: true AND that the run
    actually activated sectors (e.g. a reported auto-calibrated translation_coeff, or
    f_metab/f_maint appearing as swept inputs with a non-trivial sensitivity). If so,
    good.
(b) If the experiment-level override does NOT activate sectors (the CLI only reads
    strain-level proteome_sectors), then: set proteome_sectors.enabled: true in
    strains/eciML1515/strain.yaml, run the sectors sweep, THEN REVERT
    strain.yaml back to enabled: false (commit the revert). This leaves the baseline
    strain with sectors disabled.
Either way the sectors run must produce strains/eciML1515/outputs/sweep_sectors/
with tpc_ensemble.png, sensitivity_heatmap.png, sector_tradeoff.png, descriptors.csv,
samples.csv (with f_metab/f_maint columns), summary.json. If sectors cannot be made
to activate at all, report why and stop (do not fabricate a run).
Commit: "run calibrated and sectors experiments (with figures)".

STEP 4 - verify + summarise
- List the contents of sweep_calibrated/ and sweep_sectors/ and confirm the expected
  PNGs + CSVs exist.
- Confirm the baseline strain still has proteome_sectors disabled (grep strain.yaml).
- Print a short summary: what figures were added, the two run output dirs, the
  calibrated envelope_sampling mode used, and whether the sector run showed an
  interior optimum (peak rmax at intermediate f_metab) from sector_tradeoff data.

CONSTRAINTS
- Autonomous: no approval requests between steps; commit your own work.
- Additive figures only; do not alter existing figures or core cost/science code.
- Never pass --no-plots for the calibrated/sectors runs (the report needs the PNGs).
- Leave strains/eciML1515/strain.yaml with proteome_sectors disabled at the end.
```
