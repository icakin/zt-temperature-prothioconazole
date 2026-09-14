#!/usr/bin/env Rscript
# =============================================================================
# check_npg.R -- look at the nPG growth rates, curve by curve.
#
#   Rscript scripts/check_npg.R              # nPG (default)
#   Rscript scripts/check_npg.R sham         # same report for SHAM
#
# Run from the repository root. Reads:
#   data/aox/<inh>_{15,27}_Oxygen.csv
#   tables/aox/<inh>_fit_windows.csv     <- copy your re-trimmed windows here
#
# Prints, per curve: growth rate, fit quality, oxygen actually consumed inside
# the window, and the window itself. Then the per-culture dose slopes and the
# 15 vs 27 C comparison. Finally two diagnostics that tell you whether the
# trimming is driving the answer.
# =============================================================================
suppressPackageStartupMessages(library(minpack.lm))
INH <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "npg"
resp <- function(t,r,K,O0) O0 + (K/r)*(1-exp(r*t))

fit1 <- function(tt,yy){
  t0 <- tt-min(tt); sl <- median(diff(yy)/diff(t0))
  f <- try(nlsLM(yy ~ resp(t0,r,K,O0),
        start=list(r=1e-3,K=min(max(abs(sl),1e-6),1),O0=yy[1]),
        lower=c(r=1e-6,K=1e-10,O0=min(yy)-1), upper=c(r=.15,K=1,O0=max(yy)+1),
        control=nls.lm.control(maxiter=200,ftol=1e-12,ptol=1e-12)), silent=TRUE)
  if (inherits(f,"try-error")) return(NULL)
  p <- coef(f); pr <- resp(t0,p[["r"]],p[["K"]],p[["O0"]])
  list(r=unname(p[["r"]])*60, R2=1-sum((yy-pr)^2)/sum((yy-mean(yy))^2))
}

win <- read.csv(sprintf("tables/aox/%s_fit_windows.csv",INH), stringsAsFactors=FALSE)
rows <- list()
for (TP in c(15,27)) {
  d <- read.csv(sprintf("data/aox/%s_%d_Oxygen.csv",INH,TP), check.names=FALSE)
  for (cv in setdiff(names(d), c("Time","T"))) {
    dose <- sub("_R[0-9]+$","",cv); rep <- as.integer(sub("^.*_R","",cv))
    w <- win[win$T==TP & win$Dose==dose & win$Replicate==paste0("R",rep),]
    if (!nrow(w)) { message("no window: ",TP," ",cv); next }
    y <- suppressWarnings(as.numeric(d[[cv]]))
    m <- is.finite(y) & d$Time>=w$fit_start[1] & d$Time<=w$fit_end[1]
    if (sum(m) < 30) { message("too few points: ",TP," ",cv); next }
    f <- fit1(d$Time[m], y[m]); if (is.null(f)) { message("fit failed: ",TP," ",cv); next }
    rows[[paste0(TP,cv)]] <- data.frame(temp=TP,
      dose=ifelse(dose=="Control",0,as.numeric(dose)), culture=rep,
      r=round(f$r,4), R2=round(f$R2,4),
      O2used=round(max(y[m])-min(y[m]),2),
      start=w$fit_start[1], end=w$fit_end[1], mins=round(w$fit_end[1]-w$fit_start[1]))
  }
}
D <- do.call(rbind, rows); D <- D[order(D$temp,D$dose,D$culture),]

cat(sprintf("\n=============== %s : %d curves ===============\n", toupper(INH), nrow(D)))
for (tp in c(15,27)) {
  cat(sprintf("\n---------- %d C ----------\n", tp))
  x <- D[D$temp==tp,]
  cat(sprintf("%8s %4s %8s %7s %8s %7s %7s   %s\n",
      "dose","cul","r (1/h)","R2","O2 used","start","end","note"))
  for (i in seq_len(nrow(x))) cat(sprintf("%8.4g %4d %8.4f %7.4f %8.2f %7.0f %7.0f   %s\n",
      x$dose[i], x$culture[i], x$r[i], x$R2[i], x$O2used[i], x$start[i], x$end[i],
      ifelse(x$O2used[i] < 6, "LOW O2 USE - rate poorly identified", "")))
  mr <- tapply(x$r, x$dose, mean)
  cat(sprintf("  mean r by dose:  %s\n",
      paste(sprintf("%g mM = %.4f", as.numeric(names(mr)), mr), collapse="   ")))
}

cat("\n=============== per-culture dose slopes ===============\n")
sl <- list()
for (tp in c(15,27)) {
  x <- D[D$temp==tp & D$dose>0,]
  v <- sapply(split(x,x$culture), function(g) coef(lm(r ~ log10(dose), g))[2])
  sl[[as.character(tp)]] <- v
  cat(sprintf("%d C: %s | mean %+.5f  sd %.5f  P(vs 0) = %.3f\n",
    tp, paste(sprintf("%+.5f",v),collapse="  "), mean(v), sd(v), t.test(v)$p.value))
}
tt <- t.test(sl[["27"]], sl[["15"]], var.equal=FALSE)
cat(sprintf("\n27 C minus 15 C: %+.4f   95%% CI [%+.4f, %+.4f]   Welch P = %.3f\n",
  mean(sl[["27"]])-mean(sl[["15"]]), tt$conf.int[1], tt$conf.int[2], tt$p.value))

cat("\n=============== is the trimming driving the answer? ===============\n")
for (tp in c(15,27)) {
  x <- D[D$temp==tp,]
  c1 <- cor(x$end, x$r); c2 <- cor(x$O2used, x$r)
  cat(sprintf("%d C: r vs window END  %+.2f   |   r vs O2 consumed  %+.2f%s\n", tp, c1, c2,
      ifelse(abs(c1)>0.7|abs(c2)>0.7,"   <-- still trimming-driven","")))
}
cat("\nAim for both correlations well below 0.7. Anything above means the fitted\nrate is following where you stopped rather than how fast the culture grew.\n\n")
