# Claude Code prompt — comprehensive rewrite of the report's "Model equations, parameters and methods" section (Li-2021 style, full equations, structure + departures). REPORT-ONLY (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). PURE report writing on report.qmd (and
supplementary.qmd only if a long derivation is better placed there). No code, no analyses, no
output regeneration. SAFE TO RUN CONCURRENTLY with the P2 Bayesian (report files only; P2 owns
code + calibration outputs). Edit report files only.

GOAL: rewrite the "Model equations, parameters and methods" section into a comprehensive,
publication-grade description that (1) writes EVERY governing equation out in full with every term
defined, (2) clearly articulates the model's hierarchical structure and how the components link,
(3) states the design choices, (4) explicitly and fairly marks where our model FOLLOWS vs DEPARTS
FROM / UPGRADES the published Li et al. 2021 yeast etcGEM (Nat Commun 12:190), and (5) reflects the
current in-silico experimental design. Match Li et al.'s STYLE: a numbered constraint system,
then each temperature-dependence equation in turn with all symbols defined, then parameterisation,
then the FBA/inference. Read the Li 2021 methods (the PDF the user provided) to mirror that style.

---

```
Work AUTONOMOUSLY; commit report changes only; print a summary. Read first, and make the written
equations MATCH THE ACTUAL CODE (do not invent): src/etcgem/{enzyme_cost.py (the constraints,
_costs/_costs_unfolding, base_cost, kcat_scale/kappa_scale/dTm/tm_scale, NGAM path, Perturbation),
unfolding.py (native fraction, dGu, convergence temps, kcat(T)/rel_kcat, ngam_T),
sectors.py (sector partition, biosynthesis cap translation_coeff*v_bio<=f_bio*P_total, metabolic
pool, growth-law coupling), tpc.py (descriptors), calibration.py (likelihood/emcee)},
strains/eciML1515/strain.yaml, and reports/etcgem/report.qmd (current "Model equations, parameters
and methods" section + the parameter-provenance table). Do NOT modify code or run analyses.

WRITE THE SECTION with these parts, all as rendered LaTeX with numbered, cross-referenceable
equations (use Quarto equation labels), every symbol defined in prose:

1. THE ENZYME-CONSTRAINED MODEL (the constraint system, Li Eq-1 style). State it as a system:
   steady state S v = 0; per-reaction capacity 0 <= v_i <= kcat_i(T)*[E]_N,i(T) (flux limited by
   the NATIVE enzyme); the proteome-pool constraint (OURS is sector-partitioned — see part 4); and
   the temperature-dependent maintenance v_ATPM >= NGAM(T). Define v (flux vector), S (stoichiometry),
   [E]_N,i (native enzyme conc.), kcat_i(T). Note it is built on the GECKO enzyme-constrained
   reconstruction of iML1515 (eciML1515) for E. coli K-12 MG1655.

2. TEMPERATURE-DEPENDENT TURNOVER kcat(T) (macromolecular rate theory / transition-state; Li Eq 2-3).
   kcat_i(T) = (k_B T / h) * exp(-dG‡_i(T)/(R T));
   dG‡_i(T) = dH‡_i + dCp‡_i (T - T0) - T (dS‡_i + dCp‡_i ln(T/T0)).
   Define k_B, h, R, T0, and the transition-state enthalpy/entropy/heat-capacity dH‡, dS‡, dCp‡.
   State our anchoring choice: turnover is expressed RELATIVE to its value at each enzyme's optimum
   T_opt,i (rel_kcat = kcat(T)/kcat(T_opt)=1 at T_opt), with the absolute level carried in base_cost
   (a DESIGN CHOICE; note it).

3. TWO-STATE UNFOLDING / NATIVE FRACTION (Li Eq 4-6). N <-> U equilibrium:
   f_N,i(T) = 1 / (1 + exp(-dG_u,i(T)/(R T)));
   dG_u,i(T) = dH*_i + dCp_u,i (T - T_H*) - T dS*_i - T dCp_u,i ln(T/T_S*), with convergence
   temperatures T_H* = 373.5 K, T_S* = 385 K. Define dH*, dS*, dCp_u. State that denatured enzyme
   still consumes proteome but does not catalyse (the cost is inflated by 1/f_N), which sets the
   high-temperature falling limb.

4. EFFECTIVE PER-FLUX ENZYME COST + THE PROTEOME-SECTOR PARTITION (our structure). Give
   c_i(T) = base_cost_i / (rel_kcat_i(T) * f_N,i(T)) with base_cost_i = MW_i / kcat_i(T_opt).
   Then the SECTOR partition (Basan/Scott growth laws) — OUR UPGRADE over Li's single pool: total
   protein P_total split into fractions {f_metab, f_bio, f_maint, f_chaperone} summing to 1, with
   the metabolic pool budget P_metab = sigma * f_metab * P_total. Constraints:
     metabolic pool:   sum_i c_i(T) |v_i| <= P_metab
     biosynthesis cap: kappa * v_bio <= f_bio * P_total   (kappa = translation_coeff)
     maintenance:      v_ATPM >= NGAM(T) * (f_maint / f_maint_nom)   [T-dependent, sector-scaled]
   Define sigma (in-vivo saturation), P_total, kappa, v_bio (biomass flux), and each sector fraction;
   note the fractions are MEASURED and medium- and temperature-matched (f_sector(m,T) from proteomics).

5. THE GROWTH-LAW COUPLING (our upgrade; a config-selectable mode). Proteome-conserving Scott-2010
   allocation: f_bio(mu) = f_bio_0 + s*mu and f_metab(mu) = f_metab_0 - s*mu (same slope s;
   mu = v_bio). Show the LINEARISED constraints (kappa - s*P_total) v_bio <= f_bio_0 P_total and
   sum_i c_i|v_i| + s*(sigma P_total) v_bio <= f_metab_0 (sigma P_total). Explain that the coupling
   (proteome shifts metab->ribosomes as growth rises) makes the ribosome<->metabolism trade-off,
   NOT a static cap, set the maximal rate.

6. TEMPERATURE-DEPENDENT MAINTENANCE NGAM(T) (E. coli form; departs from Li's yeast Eq 16).
   NGAM(T) = a (1 - b exp[(-E_m / k_B)(1/T_ref - 1/T)]), floored at T_ref, with a≈8.5, b≈0.62,
   E_m≈0.5 eV, T_ref=298 K (ported from the MRes E. coli model). Define terms.

7. OBJECTIVE + DESCRIPTORS. FBA maximises mu = v_bio at each temperature, giving the organismal TPC
   mu(T). Define the descriptors used throughout: T_opt = argmax mu; r_max = max mu; CT_min/CT_max
   (mu = crit_frac * r_max on the cold/hot side); niche width = CT_max - CT_min; and the rising-limb
   activation energy E_a from a Boltzmann-Arrhenius fit (ln mu vs 1/(k_B T)).

8. THE PERTURBATION / CORRECTION PARAMETERS (used by the in-silico experiments and the Bayesian).
   Explain these are GLOBAL scalars that transform the WHOLE grounded per-enzyme distribution, not
   per-enzyme free parameters: dTopt (shift optima), topt_scale (spread), dCp_scale (curvature),
   dTm (shift melting temps), tm_scale (spread), kcat_scale (turnover level), kappa_scale (translation
   efficiency), f_metab/f_maint (allocation), ngam_scale/ngam_steepness (maintenance amplitude/T-
   dependence). Give each its defining transformation.

9. PARAMETER PROVENANCE. Keep/tidy the provenance table (every value grounded in independent data or
   a stated literature value, none fit to growth in the emergent model): T_opt (Li-Engqvist/Tome),
   T_m (Leuenberger/Meltome), dCp (DLTKcat + MMRT prior), sector fractions (proteomics, per medium/T),
   P_total, sigma, growth-law slope, NGAM constants (source + coverage).

STRUCTURE / HIERARCHY paragraph: state the model as a hierarchy — grounded per-enzyme thermal
parameters -> per-enzyme kcat(T) and native fraction -> effective per-flux cost -> sector-partitioned
proteome constraints + growth-law coupling + NGAM(T) -> FBA growth mu(T) -> TPC descriptors -> (for
analysis) the global correction knobs. Make the linkage between components explicit.

RELATIONSHIP TO / DEPARTURES FROM Li et al. 2021 (a clear, fair subsection). We FOLLOW Li for the
thermal core (transition-state kcat(T), two-state native fraction, T-dependent maintenance, the
enzyme-constrained pool idea — ours is a direct descendant). We DEPART / UPGRADE by: (i) E. coli
iML1515 / MG1655 vs their yeast ecYeast7.6; (ii) grounding the per-enzyme thermal parameters in
independent data so the TPC is an a-priori prediction, rather than fitting all per-enzyme parameters
to the growth curve; (iii) replacing the single undivided protein pool with a MEASURED, medium- and
temperature-matched proteome-SECTOR partition and an explicit biosynthesis/translation cap; (iv) the
growth-rate-dependent (proteome-conserving) allocation / growth law; (v) a likelihood-based (emcee)
calibration of a low-dimensional set of GLOBAL correction knobs, kept as a separate inverse analysis,
vs their SMC-ABC over per-enzyme parameters. State each departure and WHY (the design rationale).

IN-SILICO EXPERIMENTAL DESIGN (update to the current design). Describe, at methods level, each
analysis and what it tests: the emergent a-priori validation against the exact-strain curve; the
equal-perturbation ELASTICITY (standardised step, magnitude-based); the allocation/kinetics/
stability/maintenance DECOMPOSITION (grouped variance shares + reachable-curve bands); per-enzyme
CONTROL + IDENTIFIABILITY; and the BAYESIAN calibration (likelihood, priors, emcee) with the
emergent-vs-tuned separation. Keep it comprehensive and current.

FINALISE: keep it self-contained and readable (define acronyms), match the single-supplement
structure, ensure every equation is numbered/cross-referenced and every symbol defined. quarto
render report (and supplementary if touched); confirm the PDF builds with no unresolved crossrefs.

VERIFY (report all)
1. Every governing equation (parts 1-8) is written in full with all symbols defined; equations
   numbered/cross-referenced; the written forms match the code.
2. The hierarchy/linkage paragraph and the "departures from Li et al." subsection are present and
   accurate (single-pool->sector partition; growth law; grounded params; likelihood vs ABC; E. coli).
3. The in-silico design is comprehensive and reflects the current analyses.
4. report.pdf builds with no unresolved crossrefs; single-supplement structure preserved.

CONSTRAINTS
- Report-only: edit report.qmd (and supplementary.qmd only if needed). No code, analyses, or output
  regeneration. Safe to run alongside P2 (do not touch code/calibration outputs).
- The equations MUST match the actual implementation (read the modules); do not invent forms.
- Match Li et al.'s equation-first, term-defined style; be fair and specific about departures.
- Autonomous; commit in parts: "report: full model equations (ec constraints, kcat(T), unfolding,
  sector partition, growth law, NGAM)", "report: hierarchy + departures-from-Li subsection",
  "report: updated in-silico experimental design".
```
