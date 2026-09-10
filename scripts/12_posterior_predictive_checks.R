## F_ppc_local.R — posterior predictive checks (run from the master "X0123 copy" folder)
## (review §11: "Document posterior predictive checks")
suppressMessages({library(brms); library(bayesplot)})
set.seed(42)
MD <- "models/physiology"
OUT <- "tables/revision/ppc"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

models <- list(
  hill_growth        = file.path(MD, "growth/fit_hill.rds"),
  loglinear_resp     = file.path(MD, "respiration/fit_loglinear_respiration.rds"),
  tpc_growth         = file.path(MD, "brms_sharpe_schoolfield_tpc_growth_C_per_C_h_by_dose.rds"),
  tpc_respiration    = file.path(MD, "brms_sharpe_schoolfield_tpc_respiration_C_per_C_h_by_dose.rds"),
  cue_arrhenius      = file.path(MD, "brms_arrhenius_log_CUE_common_slope.rds")
)

pbayes <- function(yrep_stat, y_stat) mean(yrep_stat >= y_stat)
summ <- list()
for (nm in names(models)) {
  fit <- readRDS(models[[nm]])
  y   <- brms::get_y(fit)
  yrep <- posterior_predict(fit, ndraws = 500)
  ## Bayesian p-values for key statistics
  ps <- c(mean = pbayes(rowMeans(yrep), mean(y)),
          sd   = pbayes(apply(yrep,1,sd), sd(y)),
          min  = pbayes(apply(yrep,1,min), min(y)),
          max  = pbayes(apply(yrep,1,max), max(y)))
  ## coverage of 95% predictive intervals
  qlo <- apply(yrep, 2, quantile, 0.025); qhi <- apply(yrep, 2, quantile, 0.975)
  cov95 <- mean(y >= qlo & y <= qhi)
  summ[[nm]] <- data.frame(model = nm, n_obs = length(y),
                           p_mean = ps["mean"], p_sd = ps["sd"],
                           p_min = ps["min"], p_max = ps["max"],
                           coverage95 = cov95,
                           max_rhat = max(brms::rhat(fit), na.rm=TRUE))
  ## figure: density overlay + stat hist + observed-vs-predicted
  png(file.path(OUT, paste0("ppc_", nm, ".png")), width = 2200, height = 700, res = 200)
  p1 <- pp_check(fit, ndraws = 60) + ggplot2::ggtitle(paste0(nm, ": density overlay"))
  p2 <- pp_check(fit, type = "stat_2d", stat = c("mean","sd"), ndraws = 500) +
        ggplot2::ggtitle("mean vs sd")
  ep <- posterior_epred(fit, ndraws = 200)
  fitm <- colMeans(ep)
  d <- data.frame(obs = y, fit = fitm)
  p3 <- ggplot2::ggplot(d, ggplot2::aes(fit, obs)) + ggplot2::geom_point(size=0.8, alpha=0.6) +
        ggplot2::geom_abline(slope=1, intercept=0, colour="red") +
        ggplot2::theme_classic() + ggplot2::ggtitle("observed vs posterior mean")
  print(patchwork::wrap_plots(p1, p2, p3, nrow = 1))
  dev.off()
  cat(nm, "done. p(mean)=", round(ps["mean"],3), " p(sd)=", round(ps["sd"],3),
      " cov95=", round(cov95,3), " maxRhat=", round(max(brms::rhat(fit),na.rm=TRUE),4), "\n")
}
res <- do.call(rbind, summ)
write.csv(res, file.path(OUT, "ppc_summary.csv"), row.names = FALSE)
print(res)
