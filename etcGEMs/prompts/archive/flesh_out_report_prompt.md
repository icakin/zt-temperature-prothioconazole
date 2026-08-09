# Claude Code prompt — flesh out the etcGEM report (remove placeholders, full write-up)

Run from the project root (`.../MICROADAPT/etcGEMs`). Rewrites
reports/etcgem/report.qmd into a complete, collaborator-ready report and
re-renders it. Autonomous; commit at the end.

NOTE TO USER: launch in an auto-approving mode so it runs unattended.

---

```
Work AUTONOMOUSLY. Turn reports/etcgem/report.qmd into a complete, well-written
scientific report: remove EVERY placeholder (the "*(Placeholder ...)*" notes and the
placeholder abstract) and flesh out all prose. Keep the existing figure embeds and
the executable python table chunks working (do not break any @fig-/@tbl- crossref or
change asset paths); you may expand captions. Commit at the end and re-render.

HARD RULES ON ACCURACY (this is a grant report — do not invent anything)
- Every quantitative claim in the prose must come from the ACTUAL result files.
  Before writing each results section, READ the relevant files under
  reports/etcgem/assets/tables/ (summary.json, sensitivity_spearman.csv,
  decomposition_table.csv, thermal_control.csv, identifiability.csv,
  calibrated_summary.json, calibrated_descriptors.csv, sectors_samples.csv,
  sectors_descriptors.csv, descriptors.csv) and state the real numbers you find.
  If a number is not in the files, do not state it.
- Only cite references that exist in references.bib; add any missing ones with
  correct BibTeX/DOIs (candidates below). Do not cite a paper you have not added.
- Preserve scientific honesty: keep and expand the caveats listed below.

STRUCTURE (reorganise the .qmd into these sections; keep all current figures/tables,
re-homed into the results sections)

1. Background and objectives
   - Motivation: the broader aim is to predict microbial thermal-performance-curve
     (TPC) parameters from genomes (for scaling to community respiration / carbon-
     cycle feedbacks). A cell's TPC emerges from its metabolism, so an enzyme- and
     temperature-constrained genome-scale model (etc-GEM) is the natural predictor.
     Keep this self-contained (a collaborator without the grant should follow it).
   - Strain(s) and why: we work on *E. coli* via the enzyme-constrained model
     eciML1515 (a GECKO formulation of the manually curated iML1515 reconstruction).
     Explain the choice: E. coli/iML1515 is the best-curated bacterial GEM with the
     best-characterised kcats, so it establishes and de-risks the method before
     scaling to a large, phylogenetically diverse bacterial isolate library.
   - Objectives of the in-silico work — state them as numbered objectives the rest
     of the report circles back to, e.g.:
       O1  reproduce a realistic E. coli growth TPC and define its components;
       O2  determine which model components set the thermal ENVELOPE (Topt, CTmax,
           breadth, Ea) and whether a few enzymes dominate (sequence-predictability);
       O3  determine what sets the rate MAGNITUDE / baseline rate B0 (allocation);
       O4  quantify the SEPARABILITY of a genome-set envelope vs an allocation-set
           magnitude (the key hypothesis);
       O5  assess IDENTIFIABILITY: which parameters are inferable from growth alone
           vs require omics;
       O6  move from a nominal parameter scan to CALIBRATED uncertainty using
           sequence/temperature-aware kcat predictions.

2. The model (etc-GEM) — comprehensive
   - Describe every component and how they link:
       * metabolic network + GPRs from iML1515 [Monk2017];
       * enzyme constraints / shared protein pool (GECKO) [Sanchez2017; Domenzain2022],
         solved with cobrapy [Ebrahim2013]; models from the SysBioChalmers ecModels
         container (github.com/SysBioChalmers/ecModels);
       * temperature layer: per-enzyme kcat(T) by macromolecular rate theory (MMRT),
         negative dCp giving a peaked Arrhenius response [Hobbs2013]; the peak-
         normalisation convention (reference kcat = each enzyme's maximum, so warming
         only ever raises cost) — explain why (avoids super-efficiency, keeps the pool
         binding, gives a single-peaked TPC);
       * proteome allocation: a single shared pool, refined into coarse-grained
         metabolic/biosynthesis/maintenance sectors after the growth laws
         [Basan2015; Scott2010]. Note these sectors are temperature-INDEPENDENT here,
         whereas coarse-grained models of temperature-dependent growth make allocation
         itself temperature-dependent (reallocation to chaperone/stress proteins at
         extremes) [Mairet2021; Wang2026] -- a limitation revisited in the conclusion;
       * kcat(T) data source: DLTKcat temperature-dependent turnover predictions
         [Qiu2024] (github.com/SizheQiu/DLTKcat), fitted to per-enzyme MMRT
         parameters; the in-vitro/in-vivo kcat gap motivating this [Davidi2016;
         Heckmann2020];
       * the temperature-constrained GEM lineage this builds on: the Bayesian etc-GEM
         programme [Li2021; Pettersen2023].
   - Explain how the LINKED structure maps onto the experimental design: because the
     TPC emerges from (envelope x allocation) under a binding proteome budget, the
     model exposes two natural, separable axes to perturb (thermal envelope
     parameters vs proteome allocation) — which is exactly what the sensitivity /
     decomposition / control experiments interrogate to meet O2-O5.
   - Cite provenance explicitly (papers + the two GitHub repos above).

3. Experimental design (in-silico) and key assumptions
   - Describe each analysis and what objective it serves:
       * global sensitivity: Latin-hypercube sweep over envelope knobs (dTopt,
         topt spread, dCp scaling) and allocation (pool budget / sectors), with
         Spearman indices and descriptor distributions (O2, O3);
       * allocation-vs-envelope variance decomposition: crossed two-group design with
         a Shapley split of each descriptor's variance (O4);
       * per-enzyme thermal control + identifiability: a cheap proteome-wide screen
         plus targeted finite-difference control coefficients; identifiability as a
         first-order control-magnitude proxy (O2, O5);
       * DLTKcat-calibrated envelope and calibrated (correlated / posterior) sampling
         with a one-factor shared fraction rho (O6);
       * proteome-sector allocation trade-off (O3).
   - State the KEY ASSUMPTIONS plainly: peak-normalised reference kcats (kcat_ref =
     enzyme maximum); pool operating in the binding regime; MMRT curvature governs
     breadth/Ea and its default is a modelling choice (calibratable to a target Ea);
     the one-factor rho is a simple stand-in for the true thermostability covariance;
     variance fractions are defined relative to the swept input ranges; identifiability
     is first-order (not full Fisher information); DLTKcat skill is limited (log10
     R^2 ~ 0.6); growth-only (respiration/CUE not modelled here); all results are
     structural / in-silico predictions pending phenotype/omics calibration.

4+. Results (one subsection per result; keep the current figures and table chunks).
   For EACH result: (i) describe the figure/table, (ii) report the actual numbers and
   any stats read from the files, (iii) interpret mechanistically, (iv) explicitly
   circle back to the objective/hypothesis it addresses. Cover, in order:
     - Nominal TPC (fig tpc_ensemble; tbl-nominal from summary.json) -> O1. Report the
       nominal Topt, rmax, CTmin/CTmax, breadth, Ea, and note whether Topt/rmax are
       realistic for E. coli.
     - Global sensitivity (figs sensitivity_heatmap, descriptor_distributions; tbl-sens
       from sensitivity_spearman.csv) -> O2/O3. Report the actual Spearman values
       (which inputs drive which descriptors) and draw the envelope-vs-allocation split.
     - Allocation vs envelope decomposition (figs decomp_*; tbl-decomp from
       decomposition_table.csv) -> O4. Report phi_A/phi_E per descriptor and the
       interaction term; state whether the partition is separable/additive, and note
       which allocation axis was used (budget vs sectors) and that a scalar-budget axis
       makes the separation partly by construction.
     - Per-enzyme control + identifiability (figs control_*; tbl-control, tbl-ident)
       -> O2/O5. Name the top thermal-determinant enzyme(s) from thermal_control.csv;
       report the identifiable fraction from identifiability.csv AND be explicit about
       whether it is proteome-wide or top-K-conditioned (read the file / summary to
       tell). Tie to sequence-predictability (few determinants) and the omics need.
     - DLTKcat-calibrated envelope (fig dltkcat_ensemble) -> O6. Compare to the default
       run.
     - Calibrated uncertainty (figs calibrated_ensemble, calibrated_vs_default;
       tbl-calibrated) -> O6. Report calibrated descriptor medians and IQRs and how
       they compare to the default hand-set-range spread (read both CSVs).
     - Proteome-sector trade-off (figs sector_tradeoff, sectors_sensitivity) -> O3.
       Report the f_metab at which rmax peaks (interior optimum) from the sectors data.

Final section - Conclusion and next steps
   - CONCLUSION (write carefully; this is the key interpretive message of the report).
     The decomposition shows a clear division of labour -- allocation sets the rate
     magnitude, the thermal envelope sets curve position and shape -- but this is
     PARTLY A CONSEQUENCE OF MODEL STRUCTURE, not a purely emergent biological finding,
     and must be stated as such:
       (i) peak-normalisation fixes each enzyme's reference kcat to its maximum, so the
           peak growth rate is by construction nearly invariant to the thermal envelope
           and scales with the proteome budget; and
       (ii) the allocation axis is temperature-independent, so it cannot by itself
           reshape the temperature response.
     Consequently, fitting the model to an empirical GROWTH TPC would calibrate the
     emergent curve (and could set the MMRT curvature to the observed conserved Ea) but
     would NOT resolve the envelope-versus-allocation partition, which is
     under-determined from growth data alone [Pettersen2023]. Resolving the partition
     requires data that observe the internal state -- temperature-resolved proteomics
     (allocation) and 13C fluxomics (in vivo kcat) across multiple temperatures --
     which is exactly what motivates those omics measurements.
     Moreover, real proteome allocation is itself temperature-dependent: coarse-grained
     proteome-allocation models of temperature-dependent growth show cells reallocate
     toward chaperone / stress-response proteins at temperature extremes
     [Mairet2021; Wang2026]. Adding a temperature-dependent chaperone/stress allocation
     sector, and calibrating the sectors to quantitative proteomics as those studies do,
     would let allocation contribute to curve SHAPE (especially the high-temperature
     limit), blurring the clean separation and giving a more realistic, identifiable
     partition. The in-silico decomposition should therefore be read as a demonstration
     that the etc-GEM operationalises the genome-set / allocation-set MEASUREMENT -- the
     informative quantities being the interaction terms and the genome-to-genome
     variation in the split -- not as evidence that the separation is a fixed law.
   - NEXT STEPS (concrete, prioritised): calibrate MMRT curvature to the measured
     conserved Ea; add a temperature-dependent chaperone/stress allocation sector and
     calibrate the sectors to quantitative proteomics [Mairet2021; Wang2026]; run the
     H1.3 decomposition on the mechanistic (temperature-dependent) allocation axis;
     extend identifiability proteome-wide if not already; add respiration and CUE TPCs
     (the community-relevant currency); scale the analyses across the bacterial isolate
     library to estimate the genome-set vs allocation-set ceiling across taxa; and
     replace the in-silico priors with phenotype/proteome/flux-calibrated,
     hierarchical-Bayesian posteriors. Note which are already scaffolded vs future work.

REFERENCES to ensure exist in references.bib (add with correct DOIs if missing):
  Monk2017 (iML1515, Nat Biotechnol), Sanchez2017, Domenzain2022 (GECKO 2.0, Nat
  Commun), Ebrahim2013 (cobrapy, BMC Syst Biol), Hobbs2013, Li2021, Pettersen2023,
  Qiu2024, Basan2015, Scott2010, Davidi2016, Heckmann2020, Machado2021,
  Mairet2021 (Optimal proteome allocation and the temperature dependence of microbial
  growth laws, npj Syst Biol Appl, DOI 10.1038/s41540-021-00172-y; code
  github.com/fmairet/Temperature_Allocation), Wang2026 (A proteome optimal allocation
  model for elucidating effects of temperature on bacterial growth, npj Syst Biol Appl,
  DOI 10.1038/s41540-026-00693-4; code github.com/DeyuWang-itp/protein_allocation).

WRITING STYLE
- Clear, precise scientific prose; define terms on first use; no hype; hedge in-silico
  claims appropriately. Use the crossrefs and citations. Keep the abstract to one
  tight paragraph reflecting the real findings.
- The abstract MUST telegraph the conclusion's key caveat: state (briefly) that the
  apparent division of labour is partly a consequence of model structure (peak-
  normalised kcats and a temperature-independent allocation budget) and is not
  resolvable from growth data alone -- distinguishing genome-set from allocation-set
  contributions will require temperature-resolved proteome and flux calibration and a
  temperature-dependent allocation model.

FINISH
- `python reports/etcgem/assemble.py` (refresh assets), then
  `cd reports/etcgem && quarto render`; confirm _output/report.pdf builds with no
  unresolved crossrefs/citations. Print a short summary of what changed.
- Commit: "flesh out etcGEM report: background, model, design, results, next steps".

CONSTRAINTS
- Do not fabricate numbers or citations; state only what the result files and
  references.bib support.
- Do not change pipeline/scientific code or the analysis outputs; this is writing only.
- Keep every existing figure and table chunk functional.
```
