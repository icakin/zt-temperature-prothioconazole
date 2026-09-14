#!/usr/bin/env Rscript
# =============================================================================
# 43_azide_rates.R -- initial O2-consumption rates for the azide / SHAM
# respirometry (Figure 6c-d, Figure S17).
#
# Three independent SensorDish plates (biological replicates) at 27 C, each with
# 10 cell-free wells, 7 wells of 0.5 mM sodium azide and 7 wells of azide +
# 1 mM SHAM (layout in data/oxygen/azide_layout.csv). For every well the rate
# is the negative slope of a straight line fitted to dissolved O2 over the
# initial window: from the post-equilibration maximum to +10 h. The plateau that
# follows in azide-only wells is excluded by construction.
#
# Inputs : data/oxygen/azide_rep{1,2,3}_Oxygen.xlsx   PreSens exports
#          data/oxygen/azide_layout.csv               well -> condition
# Outputs: tables/aox/azide_well_rates.csv    one row per well
#          tables/aox/azide_plate_means.csv   plate x condition means (the unit
#                                             of replication for Figure 6d)
# Run from the repository root:  Rscript scripts/43_azide_rates.R
# =============================================================================
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
source(file.path(ROOT, "scripts/aox_common.R"))

tr <- load_azide_plates(ROOT)
rates <- do.call(rbind, lapply(split(tr, list(tr$plate, tr$well), drop = TRUE), function(g) {
  f <- initial_rate(g$t_h, g$o2)
  data.frame(plate = g$plate[1], well = g$well[1], condition = as.character(g$condition[1]),
             rate_mgL_h = f$rate, r2 = f$r2, window_start_h = f$t_start, window_end_h = f$t_end,
             stringsAsFactors = FALSE)
}))
rates <- rates[order(rates$plate, match(rates$condition, AZ_ORDER), rates$well), ]
rownames(rates) <- NULL

means <- aggregate(rate_mgL_h ~ plate + condition, rates, mean)
means <- reshape(means, idvar = "plate", timevar = "condition", direction = "wide")
names(means) <- sub("rate_mgL_h\\.", "", names(means))
means <- means[, c("plate", AZ_ORDER)]

dir.create(file.path(ROOT, "tables/aox"), showWarnings = FALSE, recursive = TRUE)
write.csv(rates, file.path(ROOT, "tables/aox/azide_well_rates.csv"), row.names = FALSE)
write.csv(means, file.path(ROOT, "tables/aox/azide_plate_means.csv"), row.names = FALSE)

tt <- t.test(means$azide, means$`azide+SHAM`, paired = TRUE)
cat(sprintf("azide respirometry: %d wells, %d plates, mean R2 = %.3f\n",
            nrow(rates), nrow(means), mean(rates$r2)))
print(round(means, 4), row.names = FALSE)
cat(sprintf("azide vs azide+SHAM, paired t on plate means: diff %+.4f mg/L/h, P = %.4f\n",
            mean(means$azide - means$`azide+SHAM`), tt$p.value))
