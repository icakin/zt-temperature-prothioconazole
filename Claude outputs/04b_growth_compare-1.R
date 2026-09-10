# =============================================================================
# 04b_growth_compare.R — growth rates from YOUR hand-set windows, straight to
# the dose x temperature comparison. No auto-trimming, no stable-r, no guessing.
#
# Reads : tables/physiology/Oxygen_All_Long.csv
#         tables/physiology/manual_fit_windows.csv   (written by 03)
# Fits  : O(t) = O2_0 + (K/r)*(1 - exp(r*t))   [same model/bounds as 04]
# Writes: tables/physiology/growth_by_window.csv        (one row per curve)
#         tables/physiology/growth_comparison.csv       (dose x temp summary)
#
# Also reports RESPIRATION. In the model O(t)=O2_0+(K/r)(1-exp(rt)), K is the O2
# consumption rate at the window start and equals q*N0 (biomass-specific
# respiration x initial biomass). Every vial here was inoculated at the same
# OD (0.0005), so N0 is constant across wells and K is proportional to q --
# i.e. K is a valid RELATIVE measure of biomass-specific respiration WITHIN
# this experiment. It is not an absolute respiration rate.
# K/r is reported as the respiration-to-growth ratio (the CUE-style cost index).
# =============================================================================
suppressPackageStartupMessages({library(minpack.lm)})

.this <- tryCatch({
  if (requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path))
    dirname(rstudioapi::getActiveDocumentContext()$path)
  else dirname(sys.frame(1)$ofile)
}, error=function(e) getwd())
base_dir   <- dirname(normalizePath(.this, mustWork=FALSE))
tables_dir <- file.path(base_dir,"tables","physiology")
cat("base_dir:", base_dir, "\n")

long <- read.csv(file.path(tables_dir,"Oxygen_All_Long.csv"))
wpath <- file.path(tables_dir,"manual_fit_windows.csv")
if (!file.exists(wpath)) stop("No manual_fit_windows.csv — set windows in 03 first.")
win <- read.csv(wpath)
cat("curves with a window set:", nrow(win), "\n")

resp_model <- function(r,K,t,O2_0) O2_0 + (K/r)*(1-exp(r*t))

fit_window <- function(tt, yy) {
  t0 <- tt - min(tt)
  sl <- suppressWarnings(median(diff(yy)/diff(t0), na.rm=TRUE))
  K0 <- if (is.finite(sl)) max(abs(sl),1e-6) else 1e-3
  ft <- try(nlsLM(Oxygen ~ resp_model(r,K,Time0,O2_0),
        data=data.frame(Time0=t0,Oxygen=yy),
        start=list(r=1e-3,K=K0,O2_0=yy[1]),
        lower=c(r=1e-6,K=1e-10,O2_0=min(yy)-1),
        upper=c(r=0.15, K=1,    O2_0=max(yy)+1),
        control=nls.lm.control(maxiter=500)), silent=TRUE)
  if (inherits(ft,"try-error")) return(c(NA,NA,NA,NA))
  co <- coef(ft); p <- resp_model(co[["r"]],co[["K"]],t0,co[["O2_0"]])
  sst <- sum((yy-mean(yy))^2)
  c(co[["r"]], co[["K"]],
    if (sst>1e-12) 1-sum((yy-p)^2)/sst else NA,
    sqrt(mean((yy-p)^2)))
}

res <- do.call(rbind, lapply(seq_len(nrow(win)), function(i){
  w <- win[i,]
  s <- long[long$T==w$T & as.character(long$Dose)==as.character(w$Dose) &
            toupper(long$Replicate)==toupper(w$Replicate),]
  s <- s[is.finite(s$Time) & is.finite(s$Oxygen),]
  s <- s[s$Time>=w$fit_start & s$Time<=w$fit_end,]
  s <- s[order(s$Time),]
  if (nrow(s)<6) return(NULL)
  f <- fit_window(s$Time, s$Oxygen)
  data.frame(T=w$T, Dose=as.character(w$Dose), Replicate=toupper(w$Replicate),
             fit_start=w$fit_start, fit_end=w$fit_end, n=nrow(s),
             r_per_min=f[1], r_per_h=f[1]*60, K=f[2], r2=f[3], rmse=f[4],
             rt=f[1]*(w$fit_end-w$fit_start))
}))
write.csv(res, file.path(tables_dir,"growth_by_window.csv"), row.names=FALSE)

lev <- c("Control","0.0125","0.025","0.05","0.1","0.2","0.35")
lev <- c(lev[lev %in% res$Dose], setdiff(unique(res$Dose), lev))
res$Dose <- factor(res$Dose, levels=lev)

summ <- list()
for (tt in sort(unique(res$T))) {
  d <- res[res$T==tt,]
  cat("\n", strrep("=",76), "\n", tt, " C   growth rate r (per hour)\n", strrep("=",76), "\n", sep="")
  m <- tapply(d$r_per_h, list(d$Dose,d$Replicate), function(x) x[1])
  ctrl <- if ("Control" %in% rownames(m)) mean(m["Control",],na.rm=TRUE) else NA
  tab <- data.frame(dose=rownames(m), m, check.names=FALSE)
  tab$mean <- round(rowMeans(m,na.rm=TRUE),4)
  tab$sd   <- round(apply(m,1,sd,na.rm=TRUE),4)
  tab$pct_ctrl <- if (is.finite(ctrl)) round(100*tab$mean/ctrl) else NA
  tab$p_vs_ctrl <- NA
  if ("Control" %in% rownames(m)) {
    for (dz in setdiff(rownames(m),"Control")) {
      a <- m[dz,]; b <- m["Control",]; ok <- is.finite(a)&is.finite(b)
      if (sum(ok)>=2) tab$p_vs_ctrl[tab$dose==dz] <- round(t.test(a[ok],b[ok],paired=TRUE)$p.value,3)
    }
  }
  print(tab, row.names=FALSE)

  # ---- respiration (K) and respiration:growth (K/r) ------------------------
  for (metric in c("K","K_over_r")) {
    d$val <- if (metric == "K") d$K else d$K / d$r_per_h
    mm <- tapply(d$val, list(d$Dose, d$Replicate), function(x) x[1])
    cc <- if ("Control" %in% rownames(mm)) mean(mm["Control",], na.rm=TRUE) else NA
    tb <- data.frame(dose = rownames(mm),
                     mean = signif(rowMeans(mm, na.rm=TRUE), 4),
                     sd   = signif(apply(mm, 1, sd, na.rm=TRUE), 3),
                     pct_ctrl = if (is.finite(cc)) round(100*rowMeans(mm,na.rm=TRUE)/cc) else NA)
    tb$p_vs_ctrl <- NA
    if ("Control" %in% rownames(mm))
      for (dz in setdiff(rownames(mm),"Control")) {
        a <- mm[dz,]; b <- mm["Control",]; ok <- is.finite(a) & is.finite(b)
        if (sum(ok) >= 2) tb$p_vs_ctrl[tb$dose==dz] <- round(t.test(a[ok],b[ok],paired=TRUE)$p.value,3)
      }
    cat(sprintf("\n  --- %s ---  %s\n", metric,
        if (metric=="K") "(respiration proxy: O2 consumption rate at window start)"
        else "(respiration : growth ratio)"))
    print(tb, row.names=FALSE)
    dd2 <- d[is.finite(d$val) & d$Dose != "Control",]
    if (nrow(dd2) > 3) {
      ld <- log10(as.numeric(as.character(dd2$Dose)))
      lm2 <- lm(dd2$val ~ ld)
      cat(sprintf("    trend vs log-dose: slope = %+.5g, p = %.4f\n",
          coef(lm2)[2], summary(lm2)$coefficients[2,4]))
    }
  }

  # omnibus: dose effect, blocked by culture
  dd <- d[is.finite(d$r_per_h),]
  if (nlevels(droplevels(dd$Dose))>1) {
    a <- anova(lm(r_per_h ~ factor(Dose) + factor(Replicate), data=dd))
    cat(sprintf("\n  ANOVA dose (blocked by culture): F(%d,%d) = %.2f,  p = %.4f\n",
                a[1,"Df"], a[nrow(a),"Df"], a[1,"F value"], a[1,"Pr(>F)"]))
    cat(sprintf("  -> %s\n", if (a[1,"Pr(>F)"]<0.05) "DOSE EFFECT present" else "no detectable dose effect"))
  }
  tab$T <- tt; summ[[as.character(tt)]] <- tab
}
write.csv(do.call(rbind,summ), file.path(tables_dir,"growth_comparison.csv"), row.names=FALSE)
cat("\nwritten: growth_by_window.csv  and  growth_comparison.csv\n")
