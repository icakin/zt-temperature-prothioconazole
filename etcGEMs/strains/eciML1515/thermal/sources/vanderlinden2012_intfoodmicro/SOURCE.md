# Source: Van Derlinden & Van Impe 2012, Int J Food Microbiol 158:73-78

"Modeling growth rates as a function of temperature: Model performance evaluation with
focus on the suboptimal temperature range." PDF: 1-s2.0-S0168160512002632-main.pdf.
Strain: E. coli K-12 MG1655 (CGSC #6300, E. coli Genetic Stock Center, Yale). EXACT GEM strain.
Medium: Brain Heart Infusion (BHI) broth, aerobic (rich, undefined composition).
Data: 157 mu_max(T) estimates, 7-46 C, compiled from static + dynamic experiments.

## Extracted curve
`vanderlinden2012_mg1655_bhi_tpc.csv` - digitized from Fig 1 (row 1, left panel). The y-axis is
sqrt(mu_max); I read the CENTRE of the red data cloud at each temperature and back-transformed to
mu_max (1/h) by squaring. Peak ~2.4 1/h near 40 C; growth from ~7 to ~46 C. Figure-digitized by
eye (approx +/- 0.05 sqrt-units); no raw table or SD. Both the read sqrt value and the squared
rate are stored. EXACT MG1655 strain, broad range, rich (BHI) medium.
