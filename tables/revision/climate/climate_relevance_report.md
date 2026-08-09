# Climate relevance of the measured thermal bands

ERA5 hourly 2 m air temperature, 10 NW-European winter-wheat sites,
spray window 04-01 to 07-15, 1961-1990 vs 1996-2025.

## Mean spray-window hours per season, per site

            band  early_mean_hours  recent_mean_hours  change_hours  pct_change  sites_increasing  n_sites  paired_t     paired_p
   below_assayed            1767.5             1582.3        -185.2       -10.5                 0       10    -8.984 8.663518e-06
   synergy_15_20             570.5              675.7         105.2        18.4                10       10    12.210 6.639739e-07
   neutral_20_24             157.4              207.4          50.0        31.7                10       10     5.665 3.074920e-04
antagonism_24_28              41.6               63.5          21.9        52.5                 9       10     4.216 2.251246e-03
   above_assayed               7.0               15.2           8.2       118.0                 9       10     2.679 2.524100e-02
          mean_C              12.3               13.4           1.1         8.9                10       10    11.497 1.108204e-06


Supra-optimal (>24 C) spray-window hours: 49 -> 79 per season (+62%).
Pooled trend: +10.4 h per decade (p = 0.00031, R2 = 0.20).

## Per-site trend in supra-optimal hours

                site  slope_hours_per_decade        p    r2
       Cambridge, UK                    7.03 0.012694 0.102
           Ghent, BE                   12.56 0.002803 0.144
           Halle, DE                   21.31 0.000001 0.339
        Oak Park, IE                    0.74 0.551167 0.006
Orleans (Beauce), FR                   23.37 0.000201 0.214
          Rennes, FR                   14.53 0.001252 0.166
        Roskilde, DK                    4.27 0.009658 0.110
      Rothamsted, UK                    5.98 0.015147 0.097
         Uppsala, SE                    5.00 0.079180 0.052
      Wageningen, NL                    8.86 0.010646 0.107
    ALL SITES (mean)                   10.36 0.000309 0.202

## Caveats

2 m air temperature is not canopy temperature. Hours are unweighted by
infection risk or actual spray timing. Bands outside 15-28 C are reported
but not interpreted, as they were not assayed.
