# etcGEM report (Quarto → PDF/docx)

A Quarto report that **embeds pre-generated pipeline outputs** (figures + CSV/JSON
tables) and renders to PDF and docx for sharing with collaborators. It does **not**
re-run the pipeline or any scientific code — rendering is fast because every figure
is a static PNG and every table just reads a CSV.

## Layout

```
reports/ecoli_tpc/
  _quarto.yml            project: formats pdf+docx, output-dir _output
  header.tex             LaTeX preamble (line numbers, caption styling, title block)
  references.bib         method references
  nature-communications.csl
  assemble.py            copies run outputs -> assets/ (stable names) + derived figs
  report.qmd             main report (sections 1–9)
  supplementary.qmd      full sensitivity / control / identifiability tables + config
  assets/figures/        populated by assemble.py (git-tracked)
  assets/tables/         populated by assemble.py (git-tracked)
  _output/               rendered report.pdf / report.docx (git-ignored)
```

## Prerequisites

- **Quarto CLI** (`quarto --version`).
- **LaTeX**: either a system TeX (with `xelatex`/`pdflatex`) or
  `quarto install tinytex`.
- In the project venv: `pip install jupyter ipykernel pandas tabulate`
  (`great_tables` is an optional upgrade for prettier tables).

## Workflow

1. Run the pipeline so the outputs exist (build / sweep / decompose / control /
   dltkcat). The report reads these run dirs (see the top of `assemble.py`):

   | asset group | default source run dir | produced by |
   |-------------|------------------------|-------------|
   | sweep       | `strains/eciML1515/outputs/sweep_default` | `sweep --experiment default` |
   | dltkcat     | `strains/eciML1515/outputs/sweep_dltkcat_ext` | `sweep --experiment dltkcat_ext --fits …/fits_ext.csv` |
   | decompose   | `strains/eciML1515/outputs/decompose_decomposition_sectors` | `decompose --experiment decomposition_sectors` |
   | control     | `strains/eciML1515/outputs/control_control` | `control --experiment control` (full) |
   | calibrated  | `strains/eciML1515/outputs/sweep_calibrated` | `sweep --experiment calibrated` |
   | sectors     | `strains/eciML1515/outputs/sweep_sectors` | `sweep --experiment sectors` |

   **Run the FULL experiments WITH plots** so the PNGs exist — the `*_quick`
   experiments run with `--no-plots` and emit no figures. For the M1.2 / sector
   sections:

   ```bash
   # (optional, one-off) calibrate the MMRT curvature so nominal Ea hits a target;
   # this rewrites provider.default_dCp (and may extend the grid) in strain.yaml:
   etcgem calibrate-dcp --strain eciML1515 --target-ea 0.65

   etcgem sweep     --strain eciML1515 --experiment default                 # -> sweep_default
   etcgem sweep     --strain eciML1515 --experiment dltkcat_ext \
                    --fits strains/eciML1515/dltkcat/fits_ext.csv            # -> sweep_dltkcat_ext
   etcgem decompose --strain eciML1515 --experiment decomposition_sectors   # -> decompose_decomposition_sectors
   etcgem control   --strain eciML1515 --experiment control                 # -> control_control (full; slow)
   etcgem sweep     --strain eciML1515 --experiment calibrated              # -> sweep_calibrated
   etcgem sweep     --strain eciML1515 --experiment sectors                 # -> sweep_sectors
   ```

   The `sectors` experiment (and `decomposition_sectors`) activate the proteome
   sectors for that run only (via an experiment-level `proteome_sectors` override that
   merges onto the strain, which stays disabled at baseline). If you skip the sectors
   run, `assemble.py` warns and §8 simply renders with a missing figure — remove §8 if
   you want a clean omission.

2. Assemble the assets (copies figures/tables under stable names and generates the
   two derived figures):

   ```bash
   python reports/ecoli_tpc/assemble.py
   ```

3. Render:

   ```bash
   cd reports/ecoli_tpc
   quarto render                    # report + supplementary, pdf + docx
   # or just the PDF of the main report:
   quarto render report.qmd --to pdf
   ```

   Output lands in `reports/ecoli_tpc/_output/` (`report.pdf`, `report.docx`,
   `supplementary.pdf`, `supplementary.docx`).

## Notes

- **Editing sources**: change the run dirs at the top of `assemble.py` to point at
  whichever runs you produced. Prose in `report.qmd` is intentionally short
  placeholder text for the author to expand; figures and tables are wired to real
  assets so it renders immediately.
- Only `_output/` is git-ignored; the `.qmd`, `header.tex`, `references.bib`, CSL,
  `assemble.py` and the copied `assets/` are tracked.
