#!/usr/bin/env python3
"""Assemble the etcGEM report assets.

Copies the chosen pipeline run's figures and tables into
``reports/ecoli_tpc/assets/`` under STABLE names, so report.qmd / supplementary.qmd
reference fixed filenames regardless of which experiment/run-dir produced them.
Also generates two DERIVED figures (calibrated-vs-default descriptor comparison
and the sector trade-off) that the pipeline itself does not emit.

The document never re-runs the pipeline -- it only embeds these assets. Missing
sources are skipped with a warning (assembly never fails on a missing run), and a
summary of what was copied/generated is printed at the end.

Run from the project root:  ``python reports/ecoli_tpc/assemble.py``
"""
from __future__ import annotations

import json
import os
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --- locations (edit source run dirs here) ---------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))          # project root
STRAIN = os.path.join(ROOT, "strains", "eciML1515", "outputs")

RUNS = {
    # Canonical runs only (superseded/quick/diagnostic runs are archived under
    # outputs/_archive/). The report renders from the TUNED-model + calibration dirs.
    "sweep":      os.path.join(STRAIN, "sweep_default"),          # summary.json + provenance
    "dltkcat":    os.path.join(STRAIN, "sweep_dltkcat_ext"),      # DLTKcat fits applied
    "calibrated": os.path.join(STRAIN, "sweep_calibrated"),       # M1.2; kept for resolved_config provenance
    "proteome":   os.path.join(STRAIN, "proteome_sectors"),       # temperature proteomics
    "ablation":   os.path.join(STRAIN),                           # ablation_* live in outputs/
    "anatomy":    os.path.join(STRAIN, "anatomy"),                # reference-point model anatomy
    "valid_trust": os.path.join(STRAIN, "validation"),            # Van Derlinden (MG1655, BHI) validation (renamed from validation_trusted)
    # P2 v3 calibration + P3 TUNED-model analyses (rich BHI, reconciled pool, growth law)
    "calibration_v3": os.path.join(STRAIN, "calibration_vanderlinden"),  # renamed from calibration_vanderlinden_v3
    "elasticity_tuned": os.path.join(STRAIN, "elasticity_tuned"),
    "decompose_tuned":  os.path.join(STRAIN, "decompose_tuned"),
    "control_tuned":    os.path.join(STRAIN, "control_tuned"),
}

FIG_DIR = os.path.join(HERE, "assets", "figures")
TBL_DIR = os.path.join(HERE, "assets", "tables")

# (run_key, source_filename, stable_dest_name)
FIGURES = [
    ("sweep",      "tpc_ensemble.png",             "tpc_ensemble.png"),
    ("sweep",      "sensitivity_heatmap.png",      "sensitivity_heatmap.png"),
    ("sweep",      "descriptor_distributions.png", "descriptor_distributions.png"),
    # TUNED-model decomposition + control + calibration (P3); replace the stale versions
    ("decompose_tuned", "decomp_iqr_bands.png",    "decomp_iqr_bands.png"),
    ("decompose_tuned", "decomp_variance.png",     "decomp_variance.png"),
    ("control_tuned",   "control_coefficient_bar.png", "control_thermal.png"),
    ("calibration_v3",  "prior_vs_posterior_tpc.png", "prior_vs_posterior_tpc.png"),
    ("calibration_v3",  "corner.png",              "corner_v3.png"),
    ("dltkcat",    "tpc_ensemble.png",             "dltkcat_ensemble.png"),
    ("proteome",   "sector_fractions_vs_T.png",    "proteome_sector_fractions.png"),
    ("proteome",   "usage_pred_vs_meas.png",       "proteome_usage_pred_vs_meas.png"),
    ("proteome",   "sector_pred_vs_meas.png",      "proteome_sector_pred_vs_meas.png"),
    ("ablation",   "ablation_comparison.png",      "ablation_comparison.png"),
    ("valid_trust", "validation_trusted_curves.png",  "validation_trusted_curves.png"),
    ("elasticity_tuned", "elasticity_tornado_rmax.png", "elasticity_tornado_rmax.png"),
    ("elasticity_tuned", "elasticity_tornado_CTmax_C.png", "elasticity_tornado_CTmax_C.png"),
    ("elasticity_tuned", "elasticity_heatmap.png",   "elasticity_heatmap.png"),
    ("anatomy",    "reference_tpc.png",            "reference_tpc.png"),
    ("anatomy",    "enzyme_param_densities.png",   "enzyme_param_densities.png"),
    ("anatomy",    "example_enzyme_kcatT.png",     "example_enzyme_kcatT.png"),
]

TABLES = [
    ("sweep",      "descriptors.csv",           "descriptors.csv"),
    ("sweep",      "sensitivity_spearman.csv",  "sensitivity_spearman.csv"),
    ("sweep",      "summary.json",              "summary.json"),
    # TUNED-model decomposition + control + calibration corrections (P3)
    ("decompose_tuned", "decomposition_variance.csv",      "decomposition_variance.csv"),
    ("decompose_tuned", "decomposition_iqr_magnitude.csv", "decomposition_iqr_magnitude.csv"),
    ("decompose_tuned", "decompose_summary.json",          "decompose_summary.json"),
    ("control_tuned",   "thermal_control.csv",       "thermal_control.csv"),
    ("control_tuned",   "identifiability.csv",       "identifiability.csv"),
    ("control_tuned",   "thermal_control_annotated.csv", "thermal_control_annotated.csv"),
    ("control_tuned",   "identifiability_annotated.csv", "identifiability_annotated.csv"),
    ("control_tuned",   "control_top_enzymes.csv",       "control_top_enzymes.csv"),
    ("calibration_v3",  "demanded_corrections.csv",  "demanded_corrections.csv"),
    ("calibration_v3",  "summary.json",              "calibration_v3_summary.json"),
    ("proteome",   "sector_fractions_vs_T.csv", "proteome_sector_fractions.csv"),
    ("proteome",   "validation_correlations.csv","proteome_validation_correlations.csv"),
    ("ablation",   "ablation_summary.csv",      "ablation_summary.csv"),
    ("valid_trust", "validation_trusted_table.csv",   "validation_trusted_table.csv"),
    ("valid_trust", "validation_trusted_summary.json", "validation_trusted_summary.json"),
    ("elasticity_tuned", "elasticity_table.csv",  "elasticity_table.csv"),
    ("elasticity_tuned", "elasticity_bands.json", "elasticity_bands.json"),
]

# resolved_config.yaml for supplementary provenance: first available wins.
PROVENANCE_ORDER = ["calibrated", "sweep", "proteome", "anatomy"]

KEY_DESCRIPTORS = ["Topt_C", "rmax", "CTmax_C", "niche_width_C"]


def _copy(src_dir, src_name, dest_dir, dest_name, copied, missing, tag):
    src = os.path.join(src_dir, src_name)
    if not os.path.exists(src):
        missing.append(f"{tag}:{src_name}  ({src})")
        return False
    os.makedirs(dest_dir, exist_ok=True)
    shutil.copy2(src, os.path.join(dest_dir, dest_name))
    copied.append(f"{dest_name}  <- {tag}/{src_name}")
    return True


def main():
    copied, missing = [], []
    # local enzyme-identity join (gene / protein name) for the control tables
    try:
        import annotate_enzymes
        annotate_enzymes.main()
    except Exception as e:
        print(f"[assemble] enzyme annotation skipped ({e})")
    for key, src_name, dest in FIGURES:
        _copy(RUNS[key], src_name, FIG_DIR, dest, copied, missing, key)
    for key, src_name, dest in TABLES:
        _copy(RUNS[key], src_name, TBL_DIR, dest, copied, missing, key)

    # resolved_config.yaml for provenance (first available run)
    for key in PROVENANCE_ORDER:
        if _copy(RUNS[key], "resolved_config.yaml", TBL_DIR,
                 "resolved_config.yaml", copied, missing, key):
            break

    print("=" * 70)
    print(f"COPIED / GENERATED ({len(copied)}):")
    for c in copied:
        print("  +", c)
    if missing:
        print(f"\nSKIPPED / MISSING ({len(missing)}):")
        for m in missing:
            print("  - WARNING:", m)
    print("=" * 70)
    print(f"assets -> {os.path.relpath(os.path.join(HERE, 'assets'), ROOT)}")


if __name__ == "__main__":
    main()
