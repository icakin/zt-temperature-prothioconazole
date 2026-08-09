# Running etcgem on your Mac — step by step

All commands are for the macOS Terminal, run from the project root
(`.../MICROADAPT/etcGEMs`). Copy-paste as-is (the path is quoted because the
folder name has spaces).

```bash
cd "/Users/g.yvon-durocher/Library/CloudStorage/OneDrive-UniversityofExeter/Documents/work/MICROADAPT/etcGEMs"
```

---

## Part A — Set up and run the pipeline

### 1. Create a virtual environment and install

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -e .          # installs the `etcgem` CLI + deps
```

This installs cobrapy, numpy, scipy, pandas, matplotlib, pyyaml. The GLPK LP
solver ships with cobrapy — nothing else to install. (Re-activate later with
`source .venv/bin/activate` whenever you open a new Terminal.)

### 2. Smoke-test offline (no model file needed)

```bash
etcgem tpc   --strain _toy                       # nominal curve (no experiment)
etcgem sweep --strain _toy --experiment quick    # a small sweep
```

`_toy` is a synthetic strain (cobrapy's `e_coli_core`, no model file). The sweep
writes `strains/_toy/outputs/sweep_quick/` with `tpc_ensemble.png`,
`sensitivity_heatmap.png`, `descriptors.csv`, etc. If those appear, the install
works.

### 3. Run the real E. coli sweep (eciML1515 strain)

Config is layered as **defaults ← strain ← experiment**: method settings live in
`defaults.yaml` + `experiments/`, organism biophysics in
`strains/eciML1515/strain.yaml`. The strain ships with its model at
`strains/eciML1515/model/eciML1515_batch.xml`. Inspect or solve it first:

```bash
etcgem build --strain eciML1515            # model summary -> outputs/build
etcgem tpc   --strain eciML1515            # nominal TPC   -> outputs/tpc
etcgem fba   --strain eciML1515 --temp 37  # single solve  -> outputs/fba
```

Then the full sweep:

```bash
etcgem sweep --strain eciML1515 --experiment default
```

Takes a few minutes (120 samples x 53 temperatures). Outputs land in
`strains/eciML1515/outputs/sweep_default/`: the TPC ensemble, sensitivity heatmap,
descriptor distributions, `summary.json`, and the exact merged `resolved_config.yaml`.
Edit `strains/eciML1515/strain.yaml` for organism knobs (temperature grid,
`pool_scale`, `default_dCp`) or `experiments/default.yaml` for method knobs
(sample count, parameter ranges).

For long runs on a time-limited machine, checkpoint instead:

```bash
etcgem sweep --strain eciML1515 --experiment default --resume --seconds 60   # call repeatedly until "ALL DONE"
```

That is the whole forward pipeline. Part B adds real per-enzyme thermal data.

**Adding another strain:** `cp -r strains/_template strains/NAME`, drop the
enzyme-constrained model into `strains/NAME/model/`, set `provider.model_file`
in `strains/NAME/strain.yaml`, then `etcgem build --strain NAME` and
`etcgem sweep --strain NAME --experiment default`.

---

## Part B — Add DLTKcat thermal parameters (optional, advanced)

DLTKcat is a separate deep-learning repo with heavy dependencies (PyTorch,
RDKit). Keep it in its **own environment** so it doesn't clash with etcgem.

### 1. Generate DLTKcat's input from the model

With `.venv` active, from the project root:

```bash
etcgem dltkcat prep --strain eciML1515 --tmin 5 --tmax 55 --n 11
```

Produces `strains/eciML1515/dltkcat/input.csv` with columns `rxn_id, enz, sub,
bigg, mnx, Temp_C, Temp_K` (~2340 reactions x 11 temps). Keep the `rxn_id`
column through the next steps so predictions map back.

### 2. Set up DLTKcat and predict

```bash
git clone https://github.com/SizheQiu/DLTKcat.git   # outside this repo
cd DLTKcat
# create its environment, then install: pytorch, scikit-learn, rdkit, pandas
```

Then follow DLTKcat's own workflow — its **`/code/GEMs.ipynb` notebook is the
worked example**. In short:

1. `convert_input("input.csv", enz_col='enz', sub_col='sub')` — adds `smiles`
   and `seq` columns (resolved from the UniProt id + substrate name).
2. Normalise temperature to `Temp_K_norm` and `Inv_Temp_norm` (their
   `code/feature_functions.py` / notebook show the exact training normalisation).
3. Predict → a table of log10(kcat), one row per (reaction, temperature). Make
   sure `rxn_id` and `Temp_C` are carried into the output. Save it as
   `strains/eciML1515/dltkcat/output.csv`.

### 3. Fit MMRT and run the sweep with real thermal params

Back in the project root with `.venv` active:

```bash
etcgem dltkcat parse --strain eciML1515 --temp-col Temp_C --kcat-col pred_log10kcat
etcgem sweep --strain eciML1515 --experiment default --fits
```

`parse` reads `strains/eciML1515/dltkcat/output.csv`, fits Topt/dCp per reaction
(log10 input by default; add `--no-log10` if your output is raw kcat) and writes
`strains/eciML1515/dltkcat/fits.csv`. `sweep … --fits` (default path
`strains/eciML1515/dltkcat/fits.csv`) applies those parameters — keeping the
model's calibrated enzyme costs — and runs the sweep into
`strains/eciML1515/outputs/sweep_default/`. You can also apply fits to a single
nominal curve without a sweep: `etcgem tpc --strain eciML1515 --fits`.

Check the `ok`/`r2` columns in `fits.csv`: DLTKcat is noisy (R²≈0.6), so
reactions with a poor or non-peaked fit are flagged and fall back to the default
thermal knobs instead of injecting garbage.

---

## Troubleshooting

- **`command not found: etcgem`** — activate the venv (`source .venv/bin/activate`)
  and `pip install -e .` from the project root. `python -m etcgem ...` works too.
- **Model file not found** — confirm `strains/eciML1515/model/eciML1515_batch.xml`
  is fully downloaded from OneDrive (not a cloud-only placeholder) and that
  `provider.model_file` in the strain config matches its filename.
- **Sweep feels slow** — lower `n_samples` or the temperature-grid `n` in the
  strain config; `solver_timeout` (seconds) caps any single hard LP. Or use
  `--resume --seconds N`.
- **Want faster solves** — optional: install the `gurobipy` or CPLEX solver;
  cobrapy will pick it up automatically. Not required.
