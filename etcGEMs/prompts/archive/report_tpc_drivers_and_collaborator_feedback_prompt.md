# Claude Code prompt — report improvements: (1) generalise "what carries the magnitude" into a full emergent-vs-tuned "what drives each TPC feature" discussion; (2) take on collaborator (Nick) feedback — accessibility, per-enzyme provenance, enzyme identities, abundance handling, compensation (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). REPORT-ONLY: prose + one small enzyme-
annotation join to regenerate the control/identifiability tables. Do NOT re-run the Bayesian,
sensitivity, decomposition, or control analyses, and do NOT change the model. Growth only.

NOTE TO USER: launch in an auto-approving mode.

CONTEXT: P3 has run; the report is on the tuned model at rich BHI. Two improvements:
(1) Deepen the interpretation of what the emergent-vs-tuned DIFFERENCE tells us about the drivers of
    EVERY TPC feature (not just magnitude).
(2) Address a collaborator's (non-modeller) questions so the report is clearer and the enzyme tables
    are interpretable.

---

```
Work AUTONOMOUSLY; commit in parts; print a summary. Read first: reports/etcgem/report.qmd (esp.
"Calibration: what the data demand" incl. the "What carries the magnitude" subsection ~L577; "The
complete model" / equations sections; "Which enzymes matter ..." control/identifiability section
~L789 with tbl-control/tbl-ident), reports/etcgem/supplementary.qmd,
strains/eciML1515/outputs/calibration_vanderlinden_v3/ (emergent-vs-tuned numbers + demanded_
corrections.csv), reports/etcgem/assets/tables/{thermal_control.csv,identifiability.csv}, and the
code that defines the per-enzyme inputs: src/etcgem/{providers.py, enzyme_cost.py, kcat_temp.py or
the kcat(T)/MMRT module, unfolding/Tm module, config.py} plus strain.yaml and any GECKO enzyme->gene
map, so the data-provenance description and enzyme identities are ACCURATE (describe what is actually
implemented, do not guess sources).

PART A - generalise "What carries the magnitude" -> "What the tuning reveals about each TPC feature"
- Broaden the "What carries the magnitude" subsection (in the Calibration section) from a magnitude-
  only discussion into a full treatment of what the EMERGENT -> TUNED difference (the demanded
  corrections) tells us about the driver of EACH curve feature. This is a SECOND, independent line of
  evidence that corroborates the tuned-model decomposition (which is structural): here we read the
  drivers off what the data actually had to CORRECT.
- For each feature, state the emergent value, the tuned value, the lever that moved (from
  demanded_corrections.csv), and the inference:
  * MAGNITUDE (r_max 1.04 -> 2.14 vs 2.40): the data demanded catalytic CAPACITY (sigma railed to its
    ceiling + modest kcat_scale), nothing in the envelope — so height is a capacity quantity. The
    a-priori shortfall was purely capacity, not shape.
  * UPPER LIMIT (CT_max 51.8 -> 45.9): the data demanded dTm ~ -5.6 K and essentially nothing else —
    so the hot collapse is set by protein STABILITY (T_m); the emergent melting proteome sat a few K
    too high and the tuning corrected only that axis.
  * OPTIMUM (T_opt 39.9 -> 38.4 vs 40): dTopt/topt_scale barely moved — the sequence/kinetic-grounded
    optimum was ALREADY right a priori; T_opt needed no real correction.
  * RISING LIMB / E_a: dCp was weakly constrained and the discrepancy term absorbed the possibly-
    flattened digitized E_a rather than the fit chasing it — so the cold-side shape is envelope/
    kinetic and this single curve does not pin it.
- THE HEADLINE INSIGHT to make explicit: the grounded envelope got the SHAPE and POSITION of the TPC
  essentially right out of the box (T_opt, cold-side E_a, rounded peak); the data only had to correct
  the CAPACITY (magnitude) and fine-tune the UPPER LIMIT (T_m). That is a strong statement that the
  a-priori, sequence/structure-grounded envelope is trustworthy and that growth data mainly inform
  the catalytic budget. Note this AGREES with the tuned-model decomposition (magnitude owns height;
  stability owns T_opt/CT_max; kinetics owns the cold side) and say so — two independent arguments,
  same conclusion.
- Keep it honest re the sigma ceiling / growth-law f_metab squeeze caveat already in the report; do
  not overclaim. Rename the subsection to reflect the broader scope (e.g. "What the tuning reveals
  about each feature of the curve").

PART B - accessibility: a plain-language "what goes in, and where each input comes from" (Nick Q1)
- Add a short, non-modeller-friendly passage (in "The complete model" or as a clearly-signposted
  "Inputs and data provenance" subsection) that plainly states, per ENZYME, what the model uses and
  its SOURCE, exactly as implemented (read the code to be accurate):
  * per-enzyme FLUX: model-PREDICTED by the FBA/optimisation (not measured);
  * per-enzyme turnover kcat(T): the temperature-dependent turnover (state the actual method used —
    e.g. DLTKcat-based kcat(T) / MMRT curvature) and its source;
  * per-enzyme thermal OPTIMUM / curvature: state the actual source (sequence/ML-predicted optima;
    Li et al. approach; DLTKcat) as implemented;
  * per-enzyme thermal UNFOLDING (native fraction / T_m): state the actual source (the melting-
    proteome / Leuenberger et al. data via the MRes/Li two-state approach) as implemented.
- Answer Nick's exact question in the text: YES, these are applied PER ENZYME individually (each
  enzyme carries its own flux, kcat(T), optimum/curvature and T_m), and cite the provenance of each.
  Be precise and do not attribute a source we do not actually use.

PART C - enzyme identities in the control/identifiability tables (Nick Q2)
- The control/identifiability tables currently show only rxn_id + enzyme_id (UniProt). Build a small
  annotation join so a reader can tell WHAT each enzyme is: for each enzyme_id/rxn_id add GENE NAME
  (e.g. glmS), EC NUMBER, and a short PROTEIN/ENZYME NAME, keeping the UniProt accession. Source the
  annotations LOCALLY from the model (cobra iML1515 gene .name / gene_reaction_rule, reaction
  annotation EC-code) and the GECKO enzyme->gene map; do NOT fetch from the web. Write the augmented
  table to assets/tables/ and update tbl-control (and, if present, the identifiable-enzyme list) to
  display the identities. Add one sentence naming the top few enzymes in plain English (e.g. GF6PTA /
  P17169 = glucosamine--fructose-6-phosphate aminotransferase, gene glmS, EC 2.6.1.16). If any
  annotation is missing, leave a blank rather than inventing it.

PART D - enzyme abundance handling + the compensatory-expression question (Nick Q3)
- Add a short, plain explanation of how enzyme ABUNDANCE is handled: abundances are NOT individually
  measured; each enzyme's proteome mass usage = flux / kcat(T) x MW, drawn from a SHARED proteome
  pool that the optimiser allocates; at the SECTOR level we supply measured (medium-matched) proteome
  allocation, but per-enzyme abundance is an optimisation OUTCOME, not an input.
- Address Nick's compensation intuition honestly: as temperature rises and per-enzyme kcat falls /
  enzymes unfold, MORE enzyme mass is needed to carry the same flux, so the optimiser DOES draw more
  of the pool toward those steps (an implicit, mass-balance compensation) UNTIL the pool binds — but
  there is NO regulatory SENSING ("detect low product -> upregulate"); it is the least-cost optimum
  under the pool constraint, not active feedback regulation. State this distinction clearly (it is a
  genuine limitation worth naming).

PART E - build + verify
- Re-run assemble.py; quarto render BOTH report and supplementary; confirm no unresolved crossrefs
  and the augmented enzyme table renders. Keep the single-supplement structure. Keep emergent vs
  tuned distinct; do not undo the de-foregrounding (no new results in the setup sections beyond the
  plain provenance description, which is method not results).

VERIFY (report all)
1. The "what drives each feature" subsection now covers magnitude, CT_max, T_opt and E_a via the
   emergent->tuned corrections, states the headline (shape/position right a priori; data correct
   capacity + upper limit), and notes agreement with the decomposition.
2. A plain-language per-enzyme inputs/provenance passage exists and is ACCURATE to the code (flux
   predicted; kcat(T), optima, T_m sources named as actually implemented); Nick's "per enzyme
   individually?" question is answered explicitly.
3. tbl-control (and the identifiable list) now show gene name + EC + protein name + UniProt; the top
   enzymes are named in prose; annotations sourced locally (no web).
4. Abundance-handling + compensation-vs-sensing explanation added and honest.
5. Both PDFs build; one supplement; no dangling refs; setup sections still de-foregrounded.

CONSTRAINTS
- Report-only + the local enzyme-annotation join; no model change, no re-run of Bayesian/sensitivity/
  decomposition/control. Annotations sourced locally (no web fetch).
- Describe ONLY the provenance actually implemented; do not attribute unused sources.
- Autonomous; commit in parts: "report: generalise magnitude subsection into emergent-vs-tuned drivers of every TPC feature",
  "report: plain-language per-enzyme inputs & data provenance (collaborator feedback)",
  "report: enzyme identities (gene/EC/name/UniProt) in control & identifiability tables",
  "report: enzyme-abundance handling + compensation-vs-sensing clarification".
```
