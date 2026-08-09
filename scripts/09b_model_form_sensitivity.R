## G_sensitivity_local.R — dose–response model-form sensitivity (boss review §11)
## RUN FROM the physiology/ folder:   Rscript G_sensitivity_local.R
## (or source it from RStudio with the working directory set to physiology/)
##
## Refits the growth dose–response with two alternative functional forms
## (Weibull type 1, exponential decay), compares EC50-equivalents per
## temperature and the 27 vs 15 C fold-change under each form, and runs a LOO
## comparison against the existing Hill fit (models/growth/fit_hill.rds).
## Outputs -> tables/sensitivity/ and models/sensitivity/
## Expected runtime: ~10-20 min total (two nonlinear brms fits, 4 chains each).

suppressMessages({library(brms); library(readr); library(dplyr)})
set.seed(42)
TB <- "tables/physiology"; MD <- "models/physiology/growth"
OUT_T <- "tables/physiology/sensitivity"; OUT_M <- "models/physiology/sensitivity"
for (d in c(OUT_T, OUT_M)) dir.create(d, recursive=TRUE, showWarnings=FALSE)

dat <- read_csv(file.path(TB, "fig1_extra_dose_response_by_temp.csv"), show_col_types=FALSE) %>%
  mutate(conc = Prothioconazole_mg_L, rate = Rate_raw,
         temperature = factor(Temperature_C),
         rep_id = factor(paste(Temperature_C, Replicate, sep="_"))) %>%
  filter(grepl("Growth", Trait))
cat("n =", nrow(dat), "\n")

## --- Weibull type 1: rate = r0 * exp(-(conc/lam)^k) --------------------------
f_w1 <- bf(rate ~ r0 * exp(-(conc/lam)^k),
           r0 ~ temperature + (1|rep_id), lam ~ temperature, k ~ temperature, nl=TRUE)
p_w1 <- c(prior(normal(0.04,0.02), nlpar="r0", lb=0),
          prior(student_t(3,0,0.01), class="sd", nlpar="r0"),
          prior(lognormal(0,1), nlpar="lam", lb=0),
          prior(normal(1,1), nlpar="k", lb=0.1))
fit_w1 <- brm(f_w1, data=dat, prior=p_w1, family=gaussian(), chains=4, iter=4000,
              warmup=2000, cores=4, seed=42,
              control=list(adapt_delta=0.95, max_treedepth=12),
              file=file.path(OUT_M,"fit_weibull1"))

## --- Exponential decay: rate = r0 * exp(-conc/lam) ---------------------------
f_ex <- bf(rate ~ r0 * exp(-conc/lam),
           r0 ~ temperature + (1|rep_id), lam ~ temperature, nl=TRUE)
p_ex <- c(prior(normal(0.04,0.02), nlpar="r0", lb=0),
          prior(student_t(3,0,0.01), class="sd", nlpar="r0"),
          prior(lognormal(0,1), nlpar="lam", lb=0))
fit_ex <- brm(f_ex, data=dat, prior=p_ex, family=gaussian(), chains=4, iter=4000,
              warmup=2000, cores=4, seed=42,
              control=list(adapt_delta=0.95, max_treedepth=12),
              file=file.path(OUT_M,"fit_exponential"))

fit_hill <- readRDS(file.path(MD, "fit_hill.rds"))

## --- EC50-equivalents per temperature per posterior draw ---------------------
temps <- levels(dat$temperature)
draw_param <- function(fit, par) {
  dr <- as_draws_df(fit)
  int <- dr[[paste0("b_", par, "_Intercept")]]
  sapply(temps, function(t) {
    cn <- paste0("b_", par, "_temperature", t)
    if (cn %in% names(dr)) int + dr[[cn]] else int
  })
}
q <- function(x) c(median=median(x), lo=unname(quantile(x,0.025)), hi=unname(quantile(x,0.975)))
ec50_tab <- list()
## Hill: ec50 direct
e <- draw_param(fit_hill, "ec50")
for (i in seq_along(temps)) ec50_tab[[length(ec50_tab)+1]] <-
  data.frame(model="Hill (log-logistic)", temperature=temps[i], t(q(e[,i])))
## Weibull-1: EC50 = lam * (ln 2)^(1/k)
lam <- draw_param(fit_w1, "lam"); k <- draw_param(fit_w1, "k")
for (i in seq_along(temps)) ec50_tab[[length(ec50_tab)+1]] <-
  data.frame(model="Weibull-1", temperature=temps[i], t(q(lam[,i]*log(2)^(1/k[,i]))))
## Exponential: EC50 = lam * ln 2
lam <- draw_param(fit_ex, "lam")
for (i in seq_along(temps)) ec50_tab[[length(ec50_tab)+1]] <-
  data.frame(model="Exponential", temperature=temps[i], t(q(lam[,i]*log(2))))
EC <- do.call(rbind, ec50_tab)
names(EC) <- c("model","temperature","EC50_median","EC50_lo","EC50_hi")
write.csv(EC, file.path(OUT_T, "ec50_by_model_form.csv"), row.names=FALSE)
print(EC)

## --- fold-change 27 vs 15 per model form (posterior) -------------------------
fold <- function(kind) {
  if (kind=="hill") { e <- draw_param(fit_hill,"ec50")
    x27 <- e[,temps=="27"]; x15 <- e[,temps=="15"] }
  else if (kind=="w1") { l<-draw_param(fit_w1,"lam"); kk<-draw_param(fit_w1,"k")
    x27 <- l[,temps=="27"]*log(2)^(1/kk[,temps=="27"])
    x15 <- l[,temps=="15"]*log(2)^(1/kk[,temps=="15"]) }
  else { l<-draw_param(fit_ex,"lam"); x27 <- l[,temps=="27"]*log(2); x15 <- l[,temps=="15"]*log(2) }
  r <- x27/x15; c(q(r), p_gt1 = mean(r>1))
}
FC <- rbind(data.frame(model="Hill (log-logistic)", t(fold("hill"))),
            data.frame(model="Weibull-1", t(fold("w1"))),
            data.frame(model="Exponential", t(fold("ex"))))
names(FC) <- c("model","fold_median","fold_lo","fold_hi","p_fold_gt1")
write.csv(FC, file.path(OUT_T, "ec50_foldchange_27v15_by_model_form.csv"), row.names=FALSE)
print(FC)

## --- LOO comparison ----------------------------------------------------------
lh <- loo(fit_hill); lw <- loo(fit_w1); le <- loo(fit_ex)
cmp <- loo_compare(list(Hill=lh, Weibull1=lw, Exponential=le))
capture.output(print(cmp), file=file.path(OUT_T,"loo_comparison.txt"))
print(cmp)
cat("Max Rhat: Hill", round(max(rhat(fit_hill),na.rm=TRUE),4),
    " W1", round(max(rhat(fit_w1),na.rm=TRUE),4),
    " Exp", round(max(rhat(fit_ex),na.rm=TRUE),4), "\n")
cat("\nDone. Key question: does EC50(27C) > EC50(15C) hold under every model form",
    "\n(FC table, p_fold_gt1 ~ 1), and does LOO prefer the Hill form?\n")
