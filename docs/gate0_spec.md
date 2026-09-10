# Gate 0 specification: registered feasibility audit for the v2 architecture

**Status: frozen upon I. Cakin's sign-off. The script (`scripts/57_gate0_feasibility.py`)
executes once, on real data, only after that sign-off; the branch taken is recorded
without modifying the script. All quantities computed here are architecture-design
evidence, never later presented as independent tests (charter §10).**

Companion to `docs/v2_model_design_charter.md`. This is a registered feasibility
audit preceding architecture freeze — not a blinded preregistration; the untreated
TPC is published and known.

---

## 1. Pinned literature values (charter open item 2 — now fixed)

**k_R(T), ribosomal elongation kinetics.** Baseline elongation rate in
*S. cerevisiae* is well established (~5.6–10 aa s⁻¹ at 30 °C; BioNumbers
BNID 104314/107785; Riba et al. 2019 PNAS). A directly measured elongation-rate
temperature curve for fungi is not available, so the temperature dependence is
declared as a **generic biosynthetic Q10 assumption, scenario-fixed and never
fitted**: primary Q10 = 2.0 (equivalent Arrhenius Ea ≈ 0.53 eV at 288–300 K),
band {1.5, 3.0}. Only *ratios* k_R(T)/k_R(21) enter the model, so the absolute
baseline rate cancels; the Q10 band is the entire kinetic assumption.

**P_max^lit, metabolic pool upper bound.** From the GECKO lineage
(Sánchez et al. 2017 Mol Syst Biol; Domenzain et al. 2022 Nat Commun): total
protein ≈ 0.5 g gDW⁻¹ in yeast, enzyme mass fraction f ≈ 0.5, saturation σ ≤ 1.
**P_max^lit = 0.5 × 0.5 × 1.0 = 0.25 g gDW⁻¹.** The profile in the build stage is
[P_min, 0.25]; Gate 0 only checks P_min's plausibility against it.

## 2. Frozen primary configuration (the only one that branches)

Control (0 mg L⁻¹) samples only, at 15/21/27 °C. Ribosomal sector = KEGG
map03010 gene list exactly as built by `scripts/25_kegg_gsea_tables.py` (the
published module, 117 genes). Metabolic sector = genes carried by ecZtGEM draw
reactions. Shares are **mass-weighted**: φ_s,c = Σ_{g∈s} MW_g·m_g,c / Σ_g MW_g·m_g,c
with m = DESeq2 size-factor-normalised counts (linear scale; vsd is not used for
shares) and MW from the Ensembl pep FASTA. Growth values: published Hill r₀
posterior medians and 95% CIs (μ15 = 0.0369 [0.0324, 0.0404], μ21 = 0.0523
[0.0460, 0.0590], μ24 = 0.0629, μ27 = 0.0392 [0.0356, 0.0443] h⁻¹).
Diagnostic-only (non-voting): unweighted shares; alternative sector lists.

**Input preparation (performed once, documented here; outputs committed to the
repo).** The user's R installation lacks DESeq2, so normalised counts were
extracted from `tables/rnaseq/dds.rds` in the analysis container by base-R slot
access (counts ÷ the tximport `normalizationFactors` matrix — exactly
`counts(dds, normalized=TRUE)` for this object). The dds is keyed by newer
Ensembl stable IDs (EFMGRG) while all annotation files (KEGG, pep FASTA) are
Mycgr3-keyed; a 1:1 identifier bridge was reconstructed by joining on the
per-gene `baseMean` stored in both the dds and the published DE tables
(9,209/9,209 matched, zero collisions at 8-decimal precision; 9,207 carry
Mycgr3 names, the remaining 2 retain EFMGRG and fall outside all sectors).
Committed artefacts: `tables/rnaseq/normalized_counts.csv` (Mycgr3-keyed,
9,209 × 45) and `tables/rnaseq/id_bridge_efmgrg_mycgr3.csv`. Input files and
their SHA-256 hashes are recorded in the output JSON.

## 3. Computed quantities

1. φ_R and φ_met at 15/21/27 °C (per-sample, then condition means with
   replicate SD), both weightings; **r = φ_R(27)/φ_R(21)** and its bootstrap CI
   (resampling the 7–8 control replicates per temperature).
2. **Cold-limb bracket:** is μ15/μ21 (obs) inside [Q10=3 ratio, Q10=1.5 ratio] =
   [0.518, 0.784]? Pre-computed: observed 0.706 → expected PASS; the audit
   verifies with CI overlap.
3. **27 °C ceiling test (measured temperature):** ceiling ratio = r·k(27)/k(21)
   for each Q10 in {1.5, 2, 3}, with r's CI propagated, against μ27/μ21 obs
   [CI]. Classify: BINDS (obs within ceiling CI), SLACK (ceiling lower bound >
   obs upper), VIOLATION (obs lower > ceiling upper across the whole band).
4. **24 °C compatibility inequalities** (analytic, interpolated temperature —
   design note, not a mechanism falsification): report r_min(peak) = 2·1.203/
   (k24/k21) − 1 and r_max(bind27) = 0.749/(k27/k21) per Q10, and whether the
   measured r lies in neither, one, or both intervals.
5. **P_min:** stage ecZtGEM_full.xml to the execution environment, FBA at its
   built prot_pool_exchange bound (uniform-baseline kcats — documented
   crudeness), scale linearly: P_min = bound × 0.0523/μ_max. Order-of-magnitude
   check only.
6. **DLTKcat coverage:** fraction of metabolic-sector draw enzymes with a
   DLTKcat prediction.

## 4. Decision tree (verbatim from charter §10, thresholds now numeric)

| # | condition (primary config) | branch |
|---|---|---|
| 1 | cold-limb bracket FAIL (obs ratio CI wholly outside [0.518, 0.784]) | translation is not the cold-limb mechanism → **global-Ea fallback architecture** (Ea ∈ {0.3, 0.55, 0.8} eV, fixed) |
| 2 | bracket PASS and no measured-temperature violation | **retain translation ceiling** as primary |
| 3 | 27 °C SLACK at primary Q10 = 2 | falling limb not attributed to translation; build proceeds with sector-limitation test; residual decline pre-labelled "unexplained" |
| 4 | 27 °C VIOLATION across the whole Q10 band | **reject the transcript-to-active-ribosome proxy** (charter Outcome C) |
| 5 | measured r < r_min(peak) at all Q10 (24 °C incompatibility realised) | reject proxy + linear interpolation **jointly** (Outcome B); no nonlinear interpolation fitted as rescue |
| 6 | P_min > 0.25 g gDW⁻¹, or P_min < 0.001 | metabolic-pool module flagged implausible; carried as a caveat into the prereg |
| 7 | DLTKcat coverage < 0.60 of metabolic-sector enzyme mass | claims restricted to the covered subnetwork; missing enzymes get the pre-written fallback (median kcat(T) curve) |

Branches 3, 5, 6, 7 are annotations that shape the prereg; 1, 2 and 4 select the
architecture. Multiple branches can co-fire; the output JSON records all that do.
No branch may be added, and no threshold moved, after the script runs on real data.

## 5. Execution, outputs, and record

Execution environment: the cloud container (the Mac VM lacks scipy/cobra), with
inputs staged from the repo and hashes recorded. Outputs →
`tables/revision/v2gate0/`: `gate0_results.json` (all quantities, hashes,
branches fired), `gate0_sector_shares.csv`, `gate0_log.txt`. The script has been
smoke-tested on synthetic data only; at sign-off time it has never seen the real
expression matrix. One run; the JSON is the immutable record.
