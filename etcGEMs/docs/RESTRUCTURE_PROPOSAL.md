# Repository audit and restructure proposal

_Prepared for review. Nothing here has been executed. It proposes a target
structure, a naming convention, and a phased, low-risk cleanup._

## 1. What has grown messy (findings)

**A. Heavy legacy / external directories.**

| Dir | Size | Tracked in git | Status |
|-----|------|----------------|--------|
| `DLTKcat/` | 867 MB | no (gitignored) | external vendored clone of the kcat(T) predictor; fits already applied and saved |
| `tpc_pipeline/` | 305 MB | **0 files tracked** | the original scaffold, fully superseded by `src/etcgem/`; only its (ignored) venv makes it big |
| `refs/` | 37 MB | no (gitignored) | reference PDFs; fine as-is |

**B. Root-level clutter.** Stray files sitting in the project root: `.DS_Store`,
`verify_control.log` / `verify_decomp.log` / `verify_eci.log` / `verify_p2.log`,
`*.bak*`, `collaborator_email_draft.md`, and a stray `outputs/toy_run/` (10 tracked
files, separate from the real per-strain outputs). `push_updates.sh` and
`defaults.yaml` also live loose at the root.

**C. Configuration sprawl (4 locations).** Config is scattered across the root
(`defaults.yaml`), the strain (`strains/eciML1515/config.yaml` **and**
`strain.yaml`), `configs/` (`example_gecko.yaml`, `toy.yaml`), and `experiments/`
(13 YAMLs). Notably `strains/eciML1515/config.yaml` looks **superseded by**
`strain.yaml` (strain.yaml says "ORGANISM ONLY", defers method settings to
`defaults.yaml` + `experiments/`, and marks `pool_scale` deprecated; config.yaml
still carries the old `pool_scale: 0.9`).

**D. Output-directory sprawl.** `strains/eciML1515/outputs/` holds ~40 run folders
(260 tracked files) mixing four different things:

- **Canonical (the report needs these):** `sweep_default`, `sweep_dltkcat_ext`,
  `decompose_decomposition_recast`, `control_tuned`, `elasticity_tuned`,
  `decompose_tuned`, `sweep_calibrated`, `sweep_sectors`, `proteome_sectors`,
  `validation_trusted`, `anatomy`, `calibration_vanderlinden_v3`, `ablation_*`.
- **Stale / superseded (labelled so in `assemble.py` itself):**
  `decompose_decomposition_sectors` (SUPERSEDED), `control_control` (superseded by
  `control_tuned`), `validation` (legacy), `percurve_validation` (RETIRED),
  `elasticity_elasticity` (STALE, pre-reconciliation), `calibration_vanderlinden`
  + `_v2` (superseded by `_v3`), `calibration_noll_minimal`, `calibration_phase1`.
- **Throwaway quick runs:** `control_control_quick`, `decompose_decomposition_quick`,
  `sweep_calibrated_quick`, `sweep_quick`.
- **Intermediate diagnostics:** `model_realism_p1b`, `ceiling_diagnostic`,
  `pool_reconciliation`, plus misc `default`, `build`, `tpc`, `fba`.

**E. Ugly / non-general names.** Doubled tags (`elasticity_elasticity`,
`decompose_decomposition_*`) come from the CLI writing `<command>_<experiment>` while
the experiment YAMLs repeat the command name. Version suffixes
(`calibration_vanderlinden` / `_v2` / `_v3`) and an inconsistent emergent-vs-tuned
scheme (`control_control` vs `control_tuned`) compound it.

**F. `assemble.py` carries dead wiring.** It still maps `elasticity` → the STALE
dir, `control` → `control_control`, `decompose` → the SUPERSEDED sectors dir, etc.,
even though the report now renders from the `_tuned` and `_recast` dirs.

**G. `prompts/` (~35 files, 384 KB).** Valuable provenance, but a flat pile of
executed one-offs mixed with the few still-active prompts.

## 2. Proposed target structure

```
etcGEMs/
├── README.md  pyproject.toml  requirements.txt
├── configs/                      # ALL config consolidated here
│   ├── defaults.yaml             # (moved from root)
│   ├── examples/                 # example_gecko.yaml, toy.yaml
│   └── experiments/              # the 13 experiment overlays (moved from experiments/)
├── scripts/                      # push_updates.sh, any helper scripts
├── src/etcgem/                   # the package (unchanged; it is clean)
├── strains/
│   └── eciML1515/
│       ├── strain.yaml           # single canonical strain descriptor (config.yaml removed)
│       ├── model/  media/  thermal/  proteomics/  dltkcat/   # RAW DATA — untouched
│       └── outputs/
│           ├── <curated canonical runs only, renamed scheme>
│           └── _archive/         # superseded/diagnostic runs, gitignored
├── reports/etcgem/               # report.qmd, supplementary.qmd, assets/, assemble.py (pruned)
├── prompts/
│   ├── README.md                 # index (updated)
│   ├── <active/pending prompts>
│   └── archive/                  # executed historical prompts
└── docs/                         # RUNBOOK.md, this proposal
```

External, kept outside the tracked tree and documented in the README:
`DLTKcat/` and `refs/` (already gitignored). `tpc_pipeline/` removed.

## 3. Naming convention (outputs)

`<analysis>_<modelstate>`, with modelstate in {`emergent`, `tuned`} where the
distinction applies:

| Now | Proposed |
|-----|----------|
| `elasticity_elasticity` (stale) | `elasticity_emergent` (regenerate) |
| `elasticity_tuned` | `elasticity_tuned` (keep) |
| `decompose_decomposition_recast` | `decompose_emergent` |
| `decompose_tuned` | `decompose_tuned` (keep) |
| `control_control` | `control_emergent` (regenerate) or archive |
| `control_tuned` | `control_tuned` (keep) |
| `calibration_vanderlinden_v3` | `calibration_vanderlinden` (only current kept) |
| `sweep_default` / `sweep_dltkcat_ext` / `sweep_calibrated` / `sweep_sectors` | keep (clear enough) |
| `validation_trusted` | `validation` (after retiring the legacy one) |
| `anatomy`, `ablation_*`, `proteome_sectors` | keep |

Also fix the CLI so `<command>_<experiment>` no longer doubles (rename experiment
overlays so `decompose` + `sectors` → `decompose_sectors`, not
`decompose_decomposition_sectors`).

## 4. Action table

| Item | Action |
|------|--------|
| `.DS_Store`, `verify_*.log`, `*.bak*` | delete (already gitignored) |
| `collaborator_email_draft.md` | move to `docs/correspondence/` or delete |
| root `outputs/toy_run/` | delete (stray) |
| `push_updates.sh` | move to `scripts/` |
| `defaults.yaml` | move to `configs/` |
| `configs/`, `experiments/` | consolidate under `configs/` (examples + experiments) |
| `strains/eciML1515/config.yaml` | remove (superseded by `strain.yaml`) — after confirming |
| `tpc_pipeline/` | remove from working tree (legacy, 0 tracked) |
| stale/superseded output dirs (§1D) | `git mv` to `outputs/_archive/` (or delete) |
| quick + diagnostic runs | delete or archive (regenerable) |
| canonical output dirs | `git mv` to the §3 naming; update `assemble.py` in lockstep |
| `assemble.py` | prune dead entries; point only at canonical/renamed dirs |
| `prompts/*` (executed) | `git mv` to `prompts/archive/`; keep active at top; update README |

## 5. Phased plan (safe order)

1. **Safety.** Do it on a branch; confirm `quarto render` currently succeeds so we
   have a known-good baseline to diff against.
2. **Root + junk.** Delete stray files; create `scripts/`; move `push_updates.sh`.
3. **Config consolidation.** Move `defaults.yaml` + examples + experiments under
   `configs/`; update the CLI's config search paths; remove the superseded
   `config.yaml`; re-run one CLI command to confirm resolution still works.
4. **Legacy removal.** Remove `tpc_pipeline/` and root `outputs/toy_run/`.
5. **Outputs curation.** `git mv` stale/quick/diagnostic runs to `_archive/`
   (gitignored); rename canonical runs to the new scheme.
6. **Report rewire.** Update `assemble.py` paths + any `report.qmd`/`supplementary.qmd`
   asset references; re-run `assemble.py`; `quarto render` BOTH; diff the PDFs
   against the baseline to confirm nothing broke.
7. **Prompts + docs.** Archive executed prompts; update `prompts/README.md` and the
   top-level `README.md` structure section.

## 6. Risks and safeguards

- **The report build is the main dependency.** `assemble.py` and the `.qmd` files
  reference exact output-dir names and asset paths. Renames must update
  `assemble.py` in the same commit, and we re-render + diff the PDF to verify.
- **Preserve history:** use `git mv`, not delete-and-readd.
- **Never touch raw data:** `model/`, `media/`, `thermal/`, `proteomics/`,
  `dltkcat/` inputs are precious and stay put.
- **Keep expensive results:** `calibration_vanderlinden_v3` cost ~4.5 h — keep it
  (renamed), never delete.
- **Regenerable vs precious:** all analysis outputs can be re-run from prompts, so
  archiving/deleting them is low-risk; the only care is not to break the report.

## 7. Decisions needed before executing

1. Superseded/quick/diagnostic outputs: **archive** (keep in-repo under `_archive/`)
   or **delete** (they are regenerable)?
2. `tpc_pipeline/` (305 MB legacy): **remove** from the working tree, or keep?
3. `DLTKcat/` (867 MB): keep vendored in place, or relocate outside the repo and
   reference it by path in the README?
4. Executed prompts: **archive** into `prompts/archive/`, or leave flat?

Once these are settled I can turn the plan into a single phased restructure prompt
for Claude Code to execute on a branch, with the render-and-diff verification built in.
