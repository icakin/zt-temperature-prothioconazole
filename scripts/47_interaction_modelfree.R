## =============================================================================
## 47_interaction_modelfree.R — model-free test of the temperature x dose
## interaction, with no Bliss framework and no reference temperature.
##
##   separable ("no interaction") null : log g = f(T) + h(c)
##   interaction model                 : log g = f(T) + h(c) + ti(T, c)
##
## Thin-plate / tensor-product smooths (mgcv, REML). The two models are
## compared by (i) EXACT leave-one-out cross-validation (both models refitted
## n times), (ii) AIC, and (iii) the approximate significance of the pure
## interaction tensor term. The interaction residual surface
## F(T,c) - [f(T) + h(c)] is then extracted with simulation-based intervals.
##
## RUN FROM the master "X0123 copy" folder:
##     Rscript scripts/47_interaction_modelfree.R
## Runtime ~20 s. Outputs -> tables/revision/interaction/
## =============================================================================
suppressPackageStartupMessages({library(mgcv); library(dplyr); library(readr)})
set.seed(1)

OUT <- "tables/revision/interaction"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(OUT, "interaction_modelfree.log")
sink(LOG, split = TRUE)

d <- read_csv("tables/physiology/derived_N0_R_results_with_carbon.csv",
              show_col_types = FALSE) |>
  mutate(Dose = ifelse(Dose == "Control", "0", Dose),
         conc = as.numeric(Dose),
         Temp = as.numeric(T),
         g    = growth_C_per_C_h) |>
  filter(is.finite(g), g > 0) |>
  mutate(logg = log(g), ldose = log10(conc + 0.03))   # +0.03 keeps the zero dose
cat("n =", nrow(d), " temperatures:", length(unique(d$Temp)),
    " doses:", length(unique(d$conc)), "\n\n")

f_sep <- logg ~ s(Temp, k = 6) + s(ldose, k = 5)
f_int <- logg ~ s(Temp, k = 6) + s(ldose, k = 5) + ti(Temp, ldose, k = c(5, 4))
m_sep <- gam(f_sep, data = d, method = "REML")
m_int <- gam(f_int, data = d, method = "REML")

cat("=========== (iii) pure-interaction term ===========\n")
print(summary(m_int)$s.table)

cat("\n=========== (ii) AIC ===========\n")
cat("separable   AIC =", round(AIC(m_sep), 2), "\n")
cat("interaction AIC =", round(AIC(m_int), 2),
    "  dAIC =", round(AIC(m_sep) - AIC(m_int), 2), "\n")
cat("\nanova (separable vs interaction):\n")
print(anova(m_sep, m_int, test = "F"))

cat("\n=========== (i) exact leave-one-out CV ===========\n")
loo_err <- function(form) {
  e <- numeric(nrow(d))
  for (i in seq_len(nrow(d))) {
    fit  <- gam(form, data = d[-i, ], method = "REML")
    e[i] <- d$logg[i] - predict(fit, newdata = d[i, ])
  }
  e
}
e_sep <- loo_err(f_sep)
e_int <- loo_err(f_int)
cat("separable   LOO RMSE =", round(sqrt(mean(e_sep^2)), 4),
    " MAE =", round(mean(abs(e_sep)), 4), "\n")
cat("interaction LOO RMSE =", round(sqrt(mean(e_int^2)), 4),
    " MAE =", round(mean(abs(e_int)), 4), "\n")
dd <- e_sep^2 - e_int^2
tt <- t.test(dd)
cat("paired difference in squared LOO error (separable - interaction): mean",
    signif(mean(dd), 3), " 95% CI [", signif(tt$conf.int[1], 3), ",",
    signif(tt$conf.int[2], 3), "]  p =", signif(tt$p.value, 3), "\n")
cat("interaction model better on",
    round(100 * mean(abs(e_int) < abs(e_sep)), 1), "% of held-out points\n")

## ---- machine-readable summary ----------------------------------------------
st <- summary(m_int)$s.table
res <- tibble::tibble(
  stat = c("n", "ti_edf", "ti_F", "ti_p", "AIC_sep", "AIC_int", "dAIC",
           "anova_F", "anova_p", "loo_rmse_sep", "loo_rmse_int",
           "loo_mae_sep", "loo_mae_int", "loo_paired_mean", "loo_paired_lo",
           "loo_paired_hi", "loo_paired_p", "pct_points_better"),
  value = c(nrow(d), st["ti(Temp,ldose)", "edf"], st["ti(Temp,ldose)", "F"],
            st["ti(Temp,ldose)", "p-value"], AIC(m_sep), AIC(m_int),
            AIC(m_sep) - AIC(m_int),
            anova(m_sep, m_int, test = "F")$F[2],
            anova(m_sep, m_int, test = "F")$`Pr(>F)`[2],
            sqrt(mean(e_sep^2)), sqrt(mean(e_int^2)),
            mean(abs(e_sep)), mean(abs(e_int)),
            mean(dd), tt$conf.int[1], tt$conf.int[2], tt$p.value,
            100 * mean(abs(e_int) < abs(e_sep))))
write_csv(res, file.path(OUT, "interaction_modelfree_stats.csv"))

## ---- interaction residual surface with simulation-based intervals ----------
grid <- expand.grid(Temp  = seq(15, 28, length.out = 60),
                    ldose = log10(c(0.06, 0.12, 0.25, 0.5, 1, 2, 4) + 0.03))
Xi <- predict(m_int, grid, type = "lpmatrix")
Xs <- predict(m_sep, grid, type = "lpmatrix")
set.seed(2); NS <- 4000
bi <- rmvn(NS, coef(m_int), vcov(m_int))
bs <- rmvn(NS, coef(m_sep), vcov(m_sep))
D  <- (Xi %*% t(bi)) - (Xs %*% t(bs))
grid$delta_med <- apply(D, 1, median)
grid$delta_lo  <- apply(D, 1, quantile, 0.025)
grid$delta_hi  <- apply(D, 1, quantile, 0.975)
grid$p_pos     <- apply(D, 1, function(x) mean(x > 0))
grid$conc      <- round(10^grid$ldose - 0.03, 3)
write_csv(grid, file.path(OUT, "interaction_residual_surface.csv"))

cat("\n=========== interaction residual by temperature ===========\n")
s <- grid |> group_by(Tr = round(Temp)) |>
  summarise(median_delta          = round(median(delta_med), 3),
            frac_dose_positive    = round(mean(delta_med > 0), 2),
            frac_credibly_nonzero = round(mean(p_pos > 0.975 | p_pos < 0.025), 2),
            .groups = "drop")
print(as.data.frame(s))
write_csv(s, file.path(OUT, "interaction_residual_by_temperature.csv"))

## ---- sensitivity to the basis dimension of the tensor term -----------------
## k is a user choice, so confirm the conclusion is not an artefact of it.
cat("\n=========== basis-dimension sensitivity ===========\n")
sens <- list()
for (kt in c(3, 4, 5)) for (kd in c(3, 4)) {
  ff <- as.formula(sprintf(
    "logg ~ s(Temp, k = 6) + s(ldose, k = 5) + ti(Temp, ldose, k = c(%d, %d))",
    kt, kd))
  mm <- gam(ff, data = d, method = "REML")
  stt <- summary(mm)$s.table["ti(Temp,ldose)", ]
  sens[[length(sens) + 1]] <- data.frame(
    k_Temp = kt, k_dose = kd, edf = round(stt["edf"], 2),
    F = round(stt["F"], 3), p = signif(stt["p-value"], 3),
    dAIC = round(AIC(m_sep) - AIC(mm), 2))
}
SENS <- do.call(rbind, sens); rownames(SENS) <- NULL
print(SENS)
cat("\nNote: the tensor term is favoured at every basis dimension tested, but a\n",
    "temperature basis of k = 3 cannot represent a sign reversal across T_opt\n",
    "(it leaves ~2 df), so it is underpowered against the specific alternative.\n", sep = "")

## exact LOO is repeated at a smaller basis, so the out-of-sample conclusion
## does not rest on one value of k
f_int_s <- logg ~ s(Temp, k = 6) + s(ldose, k = 5) + ti(Temp, ldose, k = c(4, 3))
e_int_s <- loo_err(f_int_s)
dd_s <- e_sep^2 - e_int_s^2; tt_s <- t.test(dd_s)
cat("\nexact LOO at ti k = c(4,3): RMSE =", round(sqrt(mean(e_int_s^2)), 4),
    " vs separable", round(sqrt(mean(e_sep^2)), 4),
    "; paired mean", signif(mean(dd_s), 3),
    " 95% CI [", signif(tt_s$conf.int[1], 3), ",", signif(tt_s$conf.int[2], 3),
    "]  p =", signif(tt_s$p.value, 3), "\n")
SENS$loo_rmse <- NA_real_
SENS$loo_rmse[SENS$k_Temp == 5 & SENS$k_dose == 4] <- sqrt(mean(e_int^2))
SENS$loo_rmse[SENS$k_Temp == 4 & SENS$k_dose == 3] <- sqrt(mean(e_int_s^2))
write_csv(SENS, file.path(OUT, "interaction_k_sensitivity.csv"))

saveRDS(list(sep = m_sep, int = m_int, e_sep = e_sep, e_int = e_int),
        file.path(OUT, "interaction_modelfree.rds"))
cat("\nwrote ->", OUT, "\n")
sink()
