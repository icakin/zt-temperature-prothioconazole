## G_sensitivity_part2.R — finish the LOO comparison
## RUN FROM the physiology/ folder:  Rscript ../revision/G_sensitivity_part2.R
## Refits the Hill model on the SAME 156-row dataset as the alternatives
## (the saved fit_hill.rds was trained on 157 rows, so loo_compare refused),
## then compares all three forms with LOO (reloo for high-pareto-k points).
suppressMessages({library(brms); library(readr); library(dplyr)})
set.seed(42)
OUT_T <- "tables/physiology/sensitivity"; OUT_M <- "models/physiology/sensitivity"

dat <- read_csv("tables/physiology/fig1_extra_dose_response_by_temp.csv", show_col_types=FALSE) %>%
  mutate(conc = Prothioconazole_mg_L, rate = Rate_raw,
         temperature = factor(Temperature_C),
         rep_id = factor(paste(Temperature_C, Replicate, sep="_"))) %>%
  filter(grepl("Growth", Trait))

f_hill <- bf(rate ~ r0 / (1 + (conc/ec50)^n),
             r0 ~ temperature + (1|rep_id), ec50 ~ temperature, n ~ temperature, nl=TRUE)
p_hill <- c(prior(normal(0.04,0.02), nlpar="r0", lb=0),
            prior(student_t(3,0,0.01), class="sd", nlpar="r0"),
            prior(lognormal(0,1), nlpar="ec50", lb=0),
            prior(normal(1.5,1), nlpar="n", lb=0.1))
fit_hill <- brm(f_hill, data=dat, prior=p_hill, family=gaussian(), chains=4, iter=4000,
                warmup=2000, cores=4, seed=42,
                control=list(adapt_delta=0.95, max_treedepth=12),
                file=file.path(OUT_M,"fit_hill_refit"))
fit_w1 <- readRDS(file.path(OUT_M, "fit_weibull1.rds"))
fit_ex <- readRDS(file.path(OUT_M, "fit_exponential.rds"))

lh <- loo(fit_hill, reloo=TRUE); lw <- loo(fit_w1, reloo=TRUE); le <- loo(fit_ex, reloo=TRUE)
cmp <- loo_compare(list(Hill=lh, Weibull1=lw, Exponential=le))
capture.output(print(cmp, simplify=FALSE), file=file.path(OUT_T,"loo_comparison.txt"))
print(cmp)
cat("\nEC50(27C) medians on identical data — Hill refit:",
    round(unname(quantile(as_draws_df(fit_hill)$b_ec50_Intercept +
                          as_draws_df(fit_hill)$b_ec50_temperature27, 0.5)), 3), "\n")
cat("Max Rhat Hill refit:", round(max(rhat(fit_hill), na.rm=TRUE), 4), "\n")
