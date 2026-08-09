# Claude Code prompt — P1: rich (BHI) reference operating point + emergent validation against Van Derlinden (MG1655) (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). This is **P1** of the MG1655 rebuild:
switch the reference operating point to a rich medium (BHI) and validate the EMERGENT model
(nothing fit to growth) against the one trusted, exact-strain curve. Growth only.
**No Bayesian tuning and no re-run of the sensitivity/decomposition/identifiability here**
(those are rebuilt on the tuned model in P3).

NOTE TO USER: launch in an auto-approving mode so it runs unattended.

WHY / DATA DECISION: we now validate on the exact GEM strain only. **Van Derlinden 2012**
(E. coli K-12 MG1655, CGSC #6300; BHI broth; 7-46 C; peak ~2.4 h-1) is the sole validation curve.
Noll (NCM3722, wrong strain) and Erdos (unpublished) are dropped from the main analysis (Erdos
may return later as a rich cross-check). Van Derlinden's 7-46 C span covers the cold rising limb
and the hot collapse. It is figure-digitized (no per-point SD).
File: strains/eciML1515/thermal/sources/vanderlinden2012_intfoodmicro/vanderlinden2012_mg1655_bhi_tpc.csv

MEDIUM DECISION (evidence-based): BHI is not encoded in silico in the literature — even the
Rothia GEM (Van Impe/ASM, spectrum.04006-23) used LB/TSB as the defined rich stand-ins and BHI
only as an in-vitro baseline. So we treat BHI as AVAILABILITY via a curated rich component list
(strains/eciML1515/media/BHI_media.csv, already built: our LB list + tissue extras, glucose
included, TSB plant sugars excluded, asn/gln omitted by hydrolysate convention). Do NOT pin
uptake rates (that is a plain-FBA workaround; our enzyme constraints cap growth).

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first:
src/etcgem/{providers.py (set_medium: glucose_minimal + LB branches, lb_media_csv path),
validation.py (CurveSpec, CURVES, _set_medium, load_curve),sectors.py (medium-matched
allocation),config.py, cli.py}, strains/eciML1515/strain.yaml (default_medium_proteome),
strains/eciML1515/media/BHI_media.csv, the Van Derlinden CSV above, and
reports/etcgem/report.qmd (model/anatomy/validation/provenance sections). Do NOT run or modify
the sensitivity/elasticity, decomposition, control/identifiability, or calibration code.

PART A - wire "BHI" as a selectable rich medium
- Extend set_medium with a medium="BHI" option that opens the BHI_media.csv component uptakes as
  AVAILABILITY (the same code path as "LB" reading LB_media.csv; read BHI_media.csv's `Name`
  column of EX_<met>_e ids). Close other carbon sources not in the BHI set; do NOT pin uptake
  rates (uptake_ub stays large; enzyme constraints limit growth). Report components opened vs
  missing from the model. Keep glucose_minimal and LB intact.
- Sector allocation for BHI: the DeyuWang proteomics has LB/Glucose/Glycerol, not BHI, so BHI
  uses the LB (rich) medium-matched sector fractions. Wire BHI -> LB sector curve.

PART B - set the RICH (BHI) reference operating point + regenerate anatomy
- Make the reference/nominal operating point the rich medium (BHI + LB rich sector allocation),
  replacing glucose-minimal as the reference (set default_medium_proteome / the anatomy medium
  accordingly). Keep glucose-minimal available as a non-default option.
- Regenerate the model-anatomy section AT the rich reference: the reference TPC (BHI, absolute
  rate, descriptors marked), the per-enzyme parameter distributions (these are medium-independent
  - unchanged), and the example kcat(T) panel. Note the operating point is now rich (BHI).

PART C - emergent validation against Van Derlinden (raw absolute rates)
- Make Van Derlinden the primary (and only main) validation curve: add a CurveSpec for it
  (medium "BHI", strain MG1655, source CSV above) and REMOVE noll_minimal from the main
  validation; keep erdos_LB only as an optional, clearly-labelled secondary cross-check (or omit).
- Predict the EMERGENT model TPC under BHI across Van Derlinden's temperatures (7-46 C) and score
  on RAW ABSOLUTE rate (1/h): abs R^2 / RMSE, predicted vs observed Topt, rmax, CTmax, and the
  rising-limb Ea (vs the ~0.85 eV benchmark). Expect magnitude UNDER-prediction (emergent ~1.0
  vs observed ~2.4) and read it as the a-priori gap, not a tuning failure (no tuning here).
- Figure: emergent model TPC vs Van Derlinden data points, raw absolute rate, 7-46 C.

PART D - report (front half only; leave downstream analyses for P3)
- Update the model/anatomy, validation and data-provenance sections: reference operating point is
  now rich (BHI); validation is against Van Derlinden (MG1655, exact strain, 7-46 C, digitized).
  Document the medium handling: BHI as availability via a curated rich list, with the evidence
  that even the Rothia GEM did not encode BHI in silico (LB/TSB proxies) - cite it. State the
  caveats plainly: rich medium is a LESS stringent test than minimal (relaxed enzyme constraints);
  Van Derlinden is figure-digitized (no SD); BHI is approximated by a curated LB-plus-tissue list;
  and an MG1655 defined-minimal TPC remains the gold-standard future test.
- Do NOT rewrite the sensitivity/decomposition/identifiability sections in P1 - add a one-line note
  that they will be rebuilt on the tuned model (P3). Re-run assemble.py; quarto render; confirm the
  PDF builds.

VERIFY (report all)
1. medium="BHI" works: components opened vs missing; the emergent rmax under BHI vs under
   glucose_minimal (rich should be higher).
2. Reference operating point is now rich (BHI); anatomy regenerated at it.
3. Emergent validation vs Van Derlinden: abs R^2/RMSE, pred vs obs Topt/rmax/CTmax/Ea; confirm
   the a-priori magnitude gap (~1.0 vs ~2.4).
4. Van Derlinden is the sole main validation curve (Noll removed); PDF builds.

CONSTRAINTS
- Emergent model only: nothing fit to growth in P1. Medium is availability; do NOT pin uptakes.
- Absolute rates are primary. Keep glucose_minimal available (do not delete) but rich is the
  reference now.
- Do NOT touch the sensitivity/decomposition/control/calibration code or their report sections.
- Be honest about the Van Derlinden digitization (no SD) and the rich-medium-is-easier caveat.
- Autonomous; commit in parts: "add BHI rich medium (availability, curated list)",
  "rich (BHI) reference operating point + anatomy regen",
  "emergent validation vs Van Derlinden (MG1655, absolute rates)",
  "report: rich reference + Van Derlinden validation + medium provenance".
```
