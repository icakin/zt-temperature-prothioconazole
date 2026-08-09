# Source: Noll et al. 2023 (Data in Brief) — E. coli K-12 NCM3722

Raw growth data downloaded from the Data in Brief paper:
"Data on the influence of temperature on the growth of Escherichia coli in a
minimal medium containing glucose as the sole carbon source ... growth yields
and rates at each temperature from 27 to 45 °C."
- DOI / ScienceDirect: https://www.sciencedirect.com/science/article/pii/S2352340923001555
- PubMed: 37006390

Strain: E. coli K-12 NCM3722 (prototrophic K-12; close relative of MG1655 / iML1515).
Medium: defined glucose minimal medium (glucose-replete), aerobic.
Temperature range: 27–45 °C.
Contents: (describe the files you placed here — OD curves and/or computed growth rates/yields.)

Purpose: model-matched (K-12, glucose-minimal, replete) TPC for eciML1515 validation/calibration.
To be extracted into ../../ecoli_tpc_curves.csv and ../../ExpGrowth.csv as a new curve.

## Extracted curve
`noll2023_ncm3722_tpc.csv` — Table 1 of the paper, transcribed: per-temperature
growth rate (umax, 1/h) and yield (C mol biomass / C mol glucose), each as mean +/- SD
over n_wells (28-40 growing wells), for 27-45 C. Growth rate at 27 C not reported (blank).
Strain E. coli K-12 NCM3722 (delta motA). Medium = modified MOPS + 0.5 g/L glucose +
5 ug/ml each of L-Met/His/Arg/Pro/Thr/Trp (trace amino acids). Kept as a standalone,
trusted, primary-source curve; NOT merged into the older Smith-derived compilation
(ExpGrowth.csv / ecoli_tpc_curves.csv), whose strains/media/robustness are uncertain.
