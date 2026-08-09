# zt_ipo323 — status (2026-08-08)

Mirrors an etcgem strain folder (cf. strains/eciML1515 and cauris_etcgem).

REAL (grounded in project 12129 data):
- model/ecZtGEM_lite.xml — enzyme-constrained conversion of ztGEM v0.3
  (Yeast9 template, diamond-BBH orthology to IPO323, Zt-native GPRs, AOX
  added as r_AOX, 4-reaction gap-fill). NOVELTY (checked 2026-08-08):
  automated CarveFungi drafts of Z. tritici exist (6 strains incl. IPO323;
  Zenodo 7413265; extracted to ../../gem/carvefungi_zt/) but have 3.4x
  fewer genes, absurd growth rates, and NO CYP51 or AOX — claim is "first
  curated, enzyme- and temperature-constrained, experimentally validated
  Zt GEM" (see gem/ztGEM_v04_report.md, CarveFungi section).
  Enzyme MWs computed from the IPO323 proteome (10,931 proteins).
- thermal/ExpGrowth.csv — measured control TPC, 7 temperatures (15-28 C),
  posterior-median biomass-specific growth rates from the Bayesian
  physiology pipeline (zt_tpc_measured.csv adds 95% CrIs).
- media/YMS_media.csv — current exchange set (sucrose C source, YMS proxy).
- Fitted maintenance available: ATP_maint = 79.2 * mu (see
  ../gem/ztGEM_v04_maintenance_fit.json); ngam_reaction r_4046.

LAYER-1 GROUNDING — COMPLETE (2026-08-08):
- thermal/BestParamsTopt.csv: DONE. TemStaPro (ProtT5, CPU, 7.5 h on the
  IPO323 MacBook) over all 742 model enzymes; pseudo-Tm 35.1-64.3 C
  (mean 37.6). strain.yaml now has thermal_model: unfolding +
  enzyme_params uncommented. Result: emergent Topt moved 28.0 -> 24.0 C
  (= the measured optimum, nothing fitted) and CTmax came off the grid
  edge to 36.8 C (real collapse from the enzymes' own melting points).
- dltkcat/: DONE. output.csv = 25,916 DLTKcat predictions (2,356 of 3,330
  arm reactions, 71%; the rest are lipids with no resolvable SMILES and
  keep defaults). Run entirely offline via run_dltkcat_zt.py: sequences
  from the IPO323 proteome, SMILES from ModelSEED (BiGG/MNX/KEGG/name
  matching; input_enriched.csv + substrate_smiles_map.csv). fits.csv:
  77 MMRT fits passed etcgem's ok filter, of which 17 were UNPHYSICAL
  (fitted Topt >= enzyme Tm - these poison unfolding mode: growth pins at
  ~1e-4 at all T) and are marked ok=False. Rule applied: Topt in [2,45] C
  and Topt <= Tm - 2. With the 60 clean fits: TPC unchanged at
  Topt 24.0 / rmax 0.216 / CTmax 36.8 (DLTKcat reshapes temperature
  response only; kcat magnitudes untouched by design).
  >> etcGEMs BUG for Gab: apply_fits_to_provider should reject fits with
  Topt >= Tm before overriding unfolding-mode thermal params.
- Enzyme assignment: full GECKO expansion in model/ecZtGEM_full.xml
  (3,504 isozyme/complex arms, real MWs); strain.yaml points at it.

REMAINING PLACEHOLDER:
- kcat magnitudes: uniform 25 s^-1 at T0 = 21 C (rmax ~3.4x high). This is
  deliberate - absolute rates are the `etcgem calibrate` stage (Gab).
- p_total / sigma: E. coli/yeast literature values; replace with fungal
  measurements if available.
- Falling limb 24-28 C still too flat vs measurement (model loses 4%,
  data lose ~37% by 27 C) - calibration target, plus possible Tm-offset
  tuning (temstapro_to_bestparams_zt.py --center/--slope).

VALIDATION TARGET: the measured TPC in thermal/ExpGrowth.csv; the emergent
etc-GEM TPC should reproduce its shape (Topt ~24 C, CTmax ~29-30 C) BEFORE
any calibration, per the etcgem two-stage logic.

SCIENCE HOOK: temperature x prothioconazole. CYP51 = r_0317
(Mycgr3G110231 + CPR Mycgr3G75805); AOX = r_AOX (Mycgr3G72918). The
single-enzyme safety-margin analysis (../gem/cyp51_margin_model.py,
Ea 0.47 eV, GEM-recovered EC50(T)) is the hand-built special case this
etc-GEM should reproduce and refine per-enzyme.

FIRST RUN THROUGH THE etcgem PIPELINE (2026-08-07, Claude cloud session):
- `etcgem build --strain zt_ipo323` and `etcgem tpc --strain zt_ipo323` both
  run end-to-end on this folder (etcGEMs @ main).
- Emergent nominal TPC (nothing fit to the growth curve): Topt 28.0 C,
  CTmax 35 C (grid edge; mmrt fallback — expect a real collapse once the
  TemStaPro Tm table + unfolding mode are in), rising-limb Ea 0.59 eV
  (measured control growth E = 0.36 [0.13, 0.58] eV — consistent),
  rmax 0.006 /h vs measured 0.052 — magnitude low with uniform placeholder
  kcats, which is exactly what the etcgem calibration stage corrects.
- Two model fixes were required and are baked into ecZtGEM_lite.xml:
  (1) SBML group names de-whitespaced (optlang constraint naming);
  (2) r_4046 (NGAM) bounds freed to (0, 1000) — it was fixed at 0.7, which
  the per-group enzyme budgets cannot satisfy (use growth-coupled
  maintenance, ATP_maint = 79.2 * mu, from the v0.4 fit);
  (3) four pure-reverse costed transports flipped to forward convention.
