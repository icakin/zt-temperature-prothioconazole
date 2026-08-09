# Claude Code prompt — add the "structural separation / growth-fit / temperature-dependent allocation" caveat to the report (targeted)

Run from the project root (`.../MICROADAPT/etcGEMs`). The report at
reports/etcgem/report.qmd is ALREADY fully written (flesh-out has run). This is a
TARGETED augmentation: add two references, insert one interpretive point into the
existing "Interpretation and caveats" section, add one item to "Next steps", and
telegraph it in the abstract. Do NOT rewrite the report or restate what is already
there. Autonomous; re-render and commit.

---

```
Work autonomously; commit at the end. The report reports/etcgem/report.qmd is already
complete (sections: # Background and objectives, # The model..., # In-silico
experimental design..., # Results, # Interpretation and caveats, # Next steps). Make
only the following targeted edits and do not duplicate existing content. Read the
current abstract, the "# Interpretation and caveats" section and the "# Next steps"
section first.

1. references.bib — add two entries (fetch correct BibTeX / DOIs):
   - Mairet2021: Mairet F, Gouze J-L, de Jong H. "Optimal proteome allocation and the
     temperature dependence of microbial growth laws." npj Systems Biology and
     Applications (2021), DOI 10.1038/s41540-021-00172-y.
     (code: github.com/fmairet/Temperature_Allocation)
   - Wang2026: "A proteome optimal allocation model for elucidating effects of
     temperature on bacterial growth." npj Systems Biology and Applications (2026),
     DOI 10.1038/s41540-026-00693-4. (code: github.com/DeyuWang-itp/protein_allocation)

2. "# Interpretation and caveats" — ADD one new paragraph (place it near the top of
   the section, right after the existing summary paragraph, before the bulleted
   caveats). It must make explicit — clearly and without hype — that the clean division
   of labour is PARTLY STRUCTURAL, not a purely emergent biological result:
   - rmax ~100% allocation is largely BY CONSTRUCTION: peak-normalisation fixes each
     enzyme's reference kcat to its maximum, so at the organismal optimum enzyme cost
     is (to first order) independent of the thermal envelope and simply scales with the
     proteome budget; and the allocation axis is temperature-INDEPENDENT, so it cannot
     by itself reshape the temperature response. Hence the "envelope sets shape,
     allocation sets magnitude" split is partly a property of the modelling framework.
   - Fitting the model to an empirical GROWTH TPC would calibrate the emergent curve
     but would NOT resolve the envelope-vs-allocation partition, which is
     under-determined from growth data alone [@Pettersen2023]; resolving it needs data
     that observe the internal state — temperature-resolved proteomics (allocation) and
     13C fluxomics (in vivo kcat) across multiple temperatures.
   - Real proteome allocation is itself temperature-dependent: coarse-grained
     proteome-allocation models of temperature-dependent growth show reallocation to
     chaperone / stress-response proteins at temperature extremes [@Mairet2021;
     @Wang2026]; a temperature-dependent chaperone/stress sector, calibrated to
     proteomics, would let allocation contribute to curve SHAPE and make the partition
     identifiable. So the decomposition should be read as operationalising the
     genome-set / allocation-set MEASUREMENT (the informative quantities being the
     interaction terms and the across-genome variation in the split), not as evidence
     that the separation is a fixed biological law.

3. Abstract — ADD one concise clause telegraphing the above (do not lengthen it much):
   that the apparent division of labour is partly a consequence of model structure
   (peak-normalised kcats and a temperature-independent allocation budget) and is not
   resolvable from growth data alone — distinguishing genome-set envelope from
   allocation-set magnitude will require temperature-resolved proteome and flux
   calibration and a temperature-dependent allocation model.

4. "# Next steps" — ADD one item (keep the existing numbered list; append/renumber):
   "Constrain the proteome-sector allocation with measured sector fractions."
   Coarse-grained proteome-allocation models of temperature-dependent growth
   [@Mairet2021; @Wang2026] provide E. coli sector fractions (metabolic / ribosome /
   chaperone) and their temperature dependence, which can calibrate the sector budgets
   (P_total, f_metab, f_bio) and add a temperature-dependent chaperone/stress sector —
   partially lifting the temperature-independent-allocation limitation. Note that
   per-enzyme allocation within the metabolic sector (hence the full envelope/allocation
   split) still requires temperature-resolved PER-PROTEIN proteomics (the WP1 proteome
   tier), because coarse sector fractions cannot pin ~2500 individual enzymes.

FINISH
- Re-render: `cd reports/etcgem && quarto render`; confirm _output/report.pdf builds
  with the new citations resolving (Mairet2021, Wang2026) and no broken refs.
- Commit: "report: add structural-separation / growth-fit / temperature-dependent
  allocation caveat and refs (Mairet2021, Wang2026)".

CONSTRAINTS
- Targeted edits only; do not rewrite existing prose or restate existing caveats.
- Do not fabricate numbers; cite only refs now present in references.bib.
- Do not touch pipeline/scientific code or analysis outputs (writing only).
```
