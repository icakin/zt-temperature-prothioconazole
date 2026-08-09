# Claude Code prompt — ingest temperature proteomics to ground temperature-dependent proteome allocation (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). Autonomous.

RUN ORDER: run this AFTER `align_with_mres_unfolding_model_prompt` has finished and
committed (it edits the thermal model). This prompt is complementary: the unfolding
model grounds the ENVELOPE in per-enzyme Tm; this grounds temperature-dependent
ALLOCATION in a measured E. coli proteome. It builds on the proteome-sector machinery
(src/etcgem/sectors.py) and must be BACKWARD-COMPATIBLE (off by default).

NOTE TO USER: launch in an auto-approving mode so it runs unattended.

---

```
Work AUTONOMOUSLY end to end; commit parts; print a summary. Read src/etcgem/
sectors.py, providers.py (from_gecko, enzyme_id/gene handling), enzyme_cost.py
(set_allocation), tpc.py (compute_tpc temperature loop) first. Keep everything
backward-compatible: when the new empirical-allocation option is off, behaviour and
numbers are identical.

DATA: a temperature-resolved E. coli proteome from
github.com/DeyuWang-itp/protein_allocation:
  data/tem_proteomic.csv       -> per-protein abundance across temperature
  data/proteomic_diff_medium.csv -> per-protein abundance across media (secondary)
Download/clone these. tem_proteomic.csv columns (verified):
  Accession (b-number, e.g. b3987), Description (has GN=<gene>, OS=E. coli K12),
  genename, COG (one-letter functional category), Localization,
  abundance columns: LB16_1_norm, LB16_2_norm, LB25_{1,2}_norm, LB30_{1,2},
  LB37_{1,2}, LB43_{1,2}  (LB at 16/25/30/37/43 C, two replicates each);
  Glucose{25,30,37}_{1,2}, Glycerol{25,30,37}_{1,2}; synthesis_rate, dna_sequence,
  protein_sequence, essentiality.
Use the LB series as the temperature axis (16,25,30,37,43 C); average replicates.

PART 1 - empirical temperature-dependent SECTOR fractions (data product)
- Map each protein to a coarse sector via its COG category (make the mapping a small,
  documented, configurable dict; sensible default):
    biosynthesis/ribosome = {J}
    chaperone/stress       = {O}
    metabolic              = {C,G,E,F,H,I,P,Q}
    housekeeping/other     = everything else (K,L,D,M,N,T,U,V,W,Y,Z,S,R,blank)
- Compute MASS-weighted sector fractions at each LB temperature: weight each protein's
  averaged abundance by its molecular weight (from protein_sequence length x ~110 Da,
  or a proper MW calc), sum per sector, normalise to 1. Output
  data/sector_fractions_vs_T.csv (rows: 16/25/30/37/43 C; cols: f_metab, f_bio,
  f_chaperone, f_other) and a figure sector_fractions_vs_T.png. Expect (and report)
  the chaperone sector rising toward 43 C.

PART 2 - map proteins to eciML1515 enzymes
- Build a mapping from the proteomics identifiers (b-number / genename) to our enzyme
  ids (EnzymeEntry.enzyme_id = UniProt). Use the model's gene-protein associations
  (iML1515 genes are b-numbers; GECKO maps genes->UniProt for prot_ ids) and/or a
  UniProt id-mapping / the Description gene name. Report COVERAGE: how many eciML1515
  enzymes have a measured abundance profile across the 5 temperatures.

PART 3 - wire temperature-dependent allocation into the model (opt-in)
- Add a strain/experiment option, e.g. `allocation_from_data: sector_fractions_vs_T.csv`.
  When set (and proteome_sectors enabled), make the sector fractions TEMPERATURE-
  DEPENDENT: in the compute_tpc temperature loop, set the sector allocation at each T
  by interpolating the empirical f_sector(T) (linear in T; clamp outside 16-43 C),
  calling set_allocation. This replaces the temperature-INDEPENDENT sectors with the
  measured temperature dependence (the missing ingredient). When the option is absent,
  behaviour is unchanged (hard backward-compat gate).

PART 4 - VALIDATION (the payoff: predicted vs measured allocation)
- For each LB temperature, solve the model and compute predicted per-enzyme proteome
  mass usage (cost_i(T) * |v_i|, i.e. the enzyme mass drawn from the pool). Compare to
  the MEASURED per-enzyme mass (abundance x MW) for mapped enzymes:
  * per-temperature correlation (Spearman and log-Pearson) predicted vs measured;
  * predicted vs measured SECTOR fractions across temperature.
  Write validation_enzyme_usage.csv and figures: usage_pred_vs_meas.png (scatter per
  T with R^2/rho) and sector_pred_vs_meas.png (predicted vs measured f_sector(T)).
  Report whether the model reproduces the observed reallocation (chaperone ramp,
  metabolic-sector shift). This is a real, data-grounded test of the allocation axis
  (H1.2/H1.3) rather than an in-silico assertion.

PART 5 - report + docs
- Add a short report subsection (and Next-steps note): temperature-dependent proteome
  allocation is now grounded in measured E. coli proteomics [@Wang2026], showing
  <report the measured chaperone ramp and the predicted-vs-measured agreement>. Note
  this is E. coli-specific; per-strain proteomes or sequence predictors are the
  library-scale route. Cite Wang2026 (already added) and note the data source.
- Add an `etcgem` CLI helper (e.g. `etcgem proteome-sectors --strain eciML1515
  --data <path>`) or a small script that produces Part 1 + Part 4 outputs.

VERIFY (report all)
1. Mapping coverage (enzymes with a measured temperature profile).
2. The measured chaperone-sector fraction at 43 C vs 30 C (is it higher?).
3. Predicted-vs-measured correlation per temperature; sector agreement.
4. Backward-compat gate: with allocation_from_data unset, nominal/sweep descriptors
   are byte-identical (1e-9) to current code.

CONSTRAINTS
- Backward compatible: opt-in only; disabled -> identical numbers (item 4 is a gate).
- Reuse sectors.py/set_allocation, compute_tpc, providers; do not change core cost/MMRT
  or the unfolding model.
- Handle identifier mapping and replicate averaging carefully; document the COG->sector
  choice and the MW weighting.
- Autonomous; commit parts separately: "ingest temperature proteomics -> empirical
  sector fractions", "map proteome to eciML1515 + temperature-dependent allocation",
  "validate predicted vs measured allocation + report".
```
