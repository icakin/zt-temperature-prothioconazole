# Claude Code prompt — config-model migration: defaults / strain / experiment + two-tier CLI

Self-contained. Run from the project root (`.../MICROADAPT/etcGEMs`). The earlier
src-layout reorg is ALREADY done (package is `src/etcgem`, strains live under
`strains/`). This prompt only changes the CONFIG MODEL and CLI surface so a strain
is runnable with no experiment, and method settings are shared, not duplicated.
Do NOT change any scientific/numerical code.

---

```
The project has already been reorganised to a src-layout (`src/etcgem`, console
command `etcgem`, per-strain folders under `strains/`). Read the current code
before changing anything: src/etcgem/cli.py, src/etcgem/config.py,
strains/eciML1515/config.yaml, configs/*.yaml. This task changes only the config
model and the CLI; reuse all scientific modules (mmrt, enzyme_cost, providers,
tpc, sensitivity, plotting, dltkcat) unchanged. Commit a snapshot first.

CURRENT STATE (verify by reading)
- strains/<name>/config.yaml is self-contained: it mixes organism biophysics
  (provider {type, model_file, prot_prefix, pool_id, biomass_rxn,
  default_Topt_offset, default_dCp, target_fraction, pool_scale}, T0_C,
  temperature_grid) with METHOD settings (sensitivity block, solver_timeout,
  crit_frac).
- configs/ holds toy.yaml and example_gecko.yaml (non-strain).
- cli.py exposes `sweep` (--strain NAME | --config PATH) and `dltkcat`; config.py
  injects provider.model_path and output_dir from the strain folder.
- strains/eciML1515/dltkcat/ already contains input/output/fits csvs and
  outputs/ already has results — PRESERVE all of these.

TARGET CONFIG MODEL
- defaults.yaml (project root): universal METHOD defaults only —
    solver_timeout, crit_frac, and a fallback temperature_grid.
- strains/<name>/strain.yaml: ORGANISM ONLY and self-sufficient to build + run
  the model — the provider block (as today, incl. thermal defaults + pool_scale),
  T0_C, and a strain-appropriate temperature_grid. NO sensitivity / method keys.
- experiments/<name>.yaml: OPTIONAL method overlays. Create:
    experiments/default.yaml       -> {kind: sweep, sensitivity: {...moved from the
                                       strain config...}}
    experiments/quick.yaml         -> a small-n sweep for smoke tests
    experiments/decomposition.yaml -> {kind: decompose, decomposition: {...}}
- Merge precedence when composing a run: defaults <- strain <- experiment.
  Organism keys come only from the strain; method keys from defaults, overridden
  by the experiment. config.resolve(strain, experiment=None) returns the merged
  dict; keep injecting provider.model_path and output_dir; add
  dump_resolved(cfg, out_dir) writing resolved_config.yaml into every run's output
  folder.

MIGRATION STEPS
1. Create defaults.yaml (solver_timeout, crit_frac, fallback temperature_grid).
2. For each strains/<name>/config.yaml: split into strain.yaml (provider + T0_C +
   temperature_grid) and move its `sensitivity` block into experiments/default.yaml;
   drop solver_timeout/crit_frac from the strain (now in defaults). Delete the old
   config.yaml. Keep provider internals exactly as they are so build_provider is
   unchanged.
3. Convert configs/toy.yaml -> strains/_toy/strain.yaml (provider type toy, T0_C,
   temperature_grid). Put example_gecko.yaml into strains/_template/ as a template
   strain.yaml with placeholders. Delete the configs/ directory.
4. Rework the CLI into two tiers (argparse subcommands):
   Strain-only (NO experiment):
     etcgem build --strain NAME              # build provider; print+save model summary
     etcgem tpc   --strain NAME [--fits [PATH]]   # nominal TPC + descriptors + plot
     etcgem fba   --strain NAME --temp C [--fits [PATH]]   # single solve at C
     etcgem dltkcat prep|parse --strain NAME ...  # as today (strain only)
   Strain + experiment:
     etcgem sweep     --strain NAME --experiment EXP [--fits [PATH]] [--resume] [--seconds N] [--no-plots]
     etcgem decompose --strain NAME --experiment EXP [--no-plots]
   Behaviour:
     * output dirs: build->outputs/build, tpc->outputs/tpc, fba->outputs/fba,
       sweep->outputs/sweep_EXP, decompose->outputs/decompose_EXP (all under
       strains/NAME/outputs/). Dump resolved_config.yaml into each.
     * --fits default path = strains/NAME/dltkcat/fits.csv; apply via
       dltkcat.apply_fits_to_provider before running.
     * decompose dispatches to a decomposition module IF present; if absent, print
       "decomposition module not installed yet" and exit non-zero.
     * you may keep a `--config PATH` escape hatch on sweep for ad-hoc configs.
5. Update config.py: load_defaults(), load_strain(name), load_experiment(name),
   resolve(name, experiment=None), dump_resolved(). Keep build_provider identical.

VERIFY (run all; report)
1. Strain-only with NO experiment works:
   etcgem build --strain _toy ; etcgem tpc --strain _toy
   etcgem build --strain eciML1515 ; etcgem tpc --strain eciML1515 ;
   etcgem fba --strain eciML1515 --temp 37
2. Experiment runs: etcgem sweep --strain _toy --experiment quick ;
   etcgem sweep --strain eciML1515 --experiment default (reduce n for a smoke run).
   Each writes outputs + resolved_config.yaml.
3. etcgem dltkcat prep --strain eciML1515 still works and the existing
   strains/eciML1515/dltkcat/ + outputs/ are untouched.
4. python -c "import etcgem; from etcgem import providers, dltkcat, sensitivity".
5. No references to the old self-contained config.yaml remain; README + docs
   updated to defaults/strain/experiment and the two-tier CLI.
Commit as "config model: defaults/strain/experiment + two-tier CLI".

ALSO
- Delete the leftover empty `tpc_pipeline/` wrapper directory at the project root
  if it still exists (it is vestigial from before the reorg).

CONSTRAINTS
- No change to numerical behaviour; provider construction stays byte-for-byte.
- Strains must run with no experiment; experiments are optional overlays.
- Preserve existing DLTKcat artifacts and outputs. Keep it reviewable; snapshot first.
```
