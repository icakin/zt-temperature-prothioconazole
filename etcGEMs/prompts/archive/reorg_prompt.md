# Claude Code prompt — reorganise the etcGEM TPC project

Paste the block below into Claude Code, run from the project root
(`.../MICROADAPT/etcGEMs`). It restructures the project; it must NOT change any
scientific logic.

---

```
You are reorganising an existing Python research project in place. Refactor the
folder structure and packaging ONLY — do not change any scientific/algorithmic
code behaviour. Preserve every function's logic; only move files, rename the
package, consolidate the runnable scripts, update imports, paths, configs, and
docs. Verify with tests at the end.

PROJECT ROOT
/Users/g.yvon-durocher/Library/CloudStorage/OneDrive-UniversityofExeter/Documents/work/MICROADAPT/etcGEMs

SAFETY FIRST
- If this is not already a git repo, run `git init` and commit the current state
  ("pre-reorg snapshot") so the change is revertible. Otherwise commit current
  state first. Use `git mv` for moves where possible.

CURRENT STATE (what exists now)
- etcGEMs/tpc_pipeline/               <- outer wrapper folder
  - tpc_pipeline/                     <- the Python package (double-nesting to fix)
    - __init__.py, mmrt.py, enzyme_cost.py, providers.py, tpc.py,
      sensitivity.py, plotting.py, config.py, dltkcat.py, resume.py, __main__.py
  - configs/ (example_toy.yaml, example_gecko.yaml, eciML1515.yaml, _eci_sandbox.yaml)
  - outputs/ (toy_run/, eciML1515_run/, eciML1515_dltkcat_input.csv, maybe eciML1515_run_dltkcat/)
  - prompts/ (dltkcat_sweep_prompt.md, reorg_prompt.md)
  - run_with_fits.py, requirements.txt, README.md, RUNBOOK.md
- etcGEMs/eciML1515_batch.xml         <- model file, loose at root
- etcGEMs/key_refs/                   <- papers (PDFs, docx)

TARGET STRUCTURE (create under etcGEMs/, the new project root)
- pyproject.toml            # package name = etcgem; deps from requirements.txt;
                            #   console_scripts entry point "etcgem = etcgem.cli:main"
- README.md                 # move+update from tpc_pipeline/README.md
- src/etcgem/               # the package, RENAMED from tpc_pipeline to etcgem
    __init__.py, mmrt.py, enzyme_cost.py, providers.py, tpc.py,
    sensitivity.py, plotting.py, config.py, dltkcat.py, cli.py
- strains/
    eciML1515/
      model/eciML1515_batch.xml          # move etcGEMs/eciML1515_batch.xml here
      config.yaml                        # from configs/eciML1515.yaml, de-absolutised (see below)
      dltkcat/                           # move outputs/eciML1515_dltkcat_input.csv -> dltkcat/input.csv
      outputs/
        default/                         # move old outputs/eciML1515_run/* here
        dltkcat/                         # move old outputs/eciML1515_run_dltkcat/* here if present
    _template/                           # empty skeleton: model/ config.yaml dltkcat/ outputs/
- configs/toy.yaml          # from configs/example_toy.yaml (non-strain example)
- prompts/                  # move etcGEMs/tpc_pipeline/prompts/ up to etcGEMs/prompts/
- refs/                     # rename etcGEMs/key_refs/ -> etcGEMs/refs/
- docs/RUNBOOK.md           # move from tpc_pipeline/RUNBOOK.md
After moving, delete the now-empty etcGEMs/tpc_pipeline/ wrapper.

RENAME THE PACKAGE (tpc_pipeline -> etcgem)
- Move the package to src/etcgem/. Internal relative imports (from .mmrt import ...)
  keep working. Update every remaining reference to the old name: `python -m
  tpc_pipeline...`, `import tpc_pipeline`, `from tpc_pipeline import ...`, and any
  mentions in README, RUNBOOK, docstrings, configs, and prompts/*.md.
  Grep the whole tree for "tpc_pipeline" and update all hits.

CONSOLIDATE SCRIPTS INTO ONE CLI (src/etcgem/cli.py)
Fold the behaviour of the old __main__.py, resume.py, and run_with_fits.py into a
single argparse CLI. Keep the underlying functions; just unify the entry points.
Subcommands and flags:
  etcgem sweep  --strain NAME  [--fits [PATH]] [--resume] [--seconds N] [--no-plots]
  etcgem sweep  --config PATH  [...]                # for the non-strain toy config
  etcgem dltkcat prep   --strain NAME [--tmin --tmax --n]
  etcgem dltkcat parse  --strain NAME [--pred PATH] [--key --temp-col --kcat-col --no-log10]
  etcgem dltkcat fit|csv|targets ...                # keep existing dltkcat subcommands
Behaviour:
  - `--strain NAME` resolves: root = strains/NAME; config = strains/NAME/config.yaml;
    model = strains/NAME/model/<file from config>; and it sets the run output dir to
    strains/NAME/outputs/default (or .../outputs/dltkcat when --fits is given).
    The CLI should inject these resolved absolute paths into the loaded config dict
    (override provider.model_path and output_dir) so config.py stays unchanged.
  - `--fits` with no path defaults to strains/NAME/dltkcat/fits.csv; applies via
    dltkcat.apply_fits_to_provider before the sweep (the old run_with_fits logic).
  - `--resume` uses the old resume.py checkpointing loop with --seconds budget.
  - `dltkcat prep` writes strains/NAME/dltkcat/input.csv; `parse` reads
    strains/NAME/dltkcat/output.csv by default and writes fits.csv there.
Also expose `python -m etcgem` (add a __main__.py that calls cli.main).

STRAIN CONFIG (strains/eciML1515/config.yaml)
- Start from the old configs/eciML1515.yaml but remove the absolute model_path and
  output_dir (the CLI injects them). Keep: provider type gecko, prot_prefix,
  pool_id, default_Topt_offset, default_dCp, pool_scale, T0_C, temperature_grid,
  sensitivity block, solver_timeout, crit_frac. Add a `model_file:
  eciML1515_batch.xml` key the CLI joins to strains/NAME/model/.
- Create strains/_template/config.yaml as the same with placeholder values and a
  comment explaining how to add a new strain (copy _template, drop model in model/).
- Delete the sandbox-only configs/_eci_sandbox.yaml if present.

PACKAGING
- pyproject.toml with [build-system] setuptools, [project] name="etcgem", the deps
  currently in requirements.txt, and [project.scripts] etcgem = "etcgem.cli:main".
  Use setuptools package discovery under src/ (package-dir src). Keep a thin
  requirements.txt in sync for non-editable users.

VERIFY (do all of these; report results)
1. `pip install -e .` succeeds in a fresh venv.
2. `etcgem sweep --config configs/toy.yaml` (or `--strain` if you make toy a
   strain) runs and writes its figures. 
3. `etcgem sweep --strain eciML1515 --seconds ...` or a reduced-sample smoke run
   completes and writes strains/eciML1515/outputs/default/ figures + summary.json.
4. `etcgem dltkcat prep --strain eciML1515` writes strains/eciML1515/dltkcat/input.csv
   with the expected columns and ~2300 reactions.
5. `python -c "import etcgem; from etcgem import providers, dltkcat, mmrt"` works.
6. Grep confirms no remaining "tpc_pipeline" references anywhere.
Then update README.md, docs/RUNBOOK.md, and prompts/dltkcat_sweep_prompt.md to the
new package name, paths, and `etcgem ... --strain` commands. Commit as
"reorganise project: src-layout, etcgem package, per-strain folders, unified CLI".

CONSTRAINTS
- Do NOT alter numerical/scientific behaviour of mmrt, enzyme_cost, providers,
  tpc, sensitivity, dltkcat. This is a structure/packaging refactor only.
- If any move hits a permission error (OneDrive), report it rather than forcing.
- Keep changes reviewable; commit the pre-reorg snapshot first.
```
