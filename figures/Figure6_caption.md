**Figure 6. An enzyme- and temperature-constrained genome-scale metabolic model
accounts for the temperature × prothioconazole interaction in *Zymoseptoria tritici*.**
The model (ecZtGEM) couples the curated metabolic network to a proteome budget in which
each enzyme carries a sequence-predicted turnover number (DLTKcat) and a sequence-predicted
thermal stability (TemStaPro); five global parameters are calibrated to the physiology data
(a turnover scale and an in-vivo melting-temperature offset for the thermal curve; a reference
potency, its temperature dependence and a Hill coefficient for the drug response), while all
per-enzyme kinetic and thermal parameters, the network stoichiometry, and the transcript data
are independent of the growth measurements.
**(A)** Growth thermal performance curve. The dashed grey curve is the *emergent* prediction
with nothing fitted to growth — its optimum position (24 °C) and high-temperature collapse
arise from the enzyme parameters alone; the green curve additionally calibrates two global
knobs (turnover scale, in-vivo Tm offset) to recover the absolute magnitude. Points, measured
control growth (mean ± 95% CI). The collapse near 29 °C reflects in-vivo thermal failure ~7–10 °C
below sequence-predicted unfolding, consistent with in-cell thermal proteome profiling.
**(B)** Condition-constrained flux analysis (measured growth and respiration imposed). The
model-inferred alternative-oxidase (AOX) share of oxygen consumption rises with combined drug
and heat and tracks the independently measured AOX transcript abundance (Spearman ρ = 0.99);
the transcript is not used to constrain the flux.
**(C)** Temperature-dependent potency. Prothioconazole enters as graded competitive inhibition
of CYP51 (Hill coefficient n = 0.75); lines are the model, points the measured relative-growth
surface. This replaces a hard-bound (step) formulation and reproduces the smooth measured dose
response (in-sample r = 0.89).
**(D)** Model-derived effective ATP yield per oxygen atom collapses from ~1.3 in the cool
controls to ~0.07 under combined drug and heat, a mechanistic basis for the measured
carbon-use-efficiency decline.
Panels A (magnitude) and C are calibrated fits and B/D are condition-constrained inferences;
the model is presented as an explanatory account of the observed interaction rather than a
forecast of unseen temperatures (out-of-sample EC₅₀ prediction is limited by measurement noise;
see Supplementary).
