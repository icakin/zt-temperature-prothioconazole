# =============================================================================
# 51_figure1.R — Figure 1 with new panel G (activation energies).
# Panels A-F identical to 24_figure1.R; G per agreed design:
#   paired point-intervals (growth/respiration) per dose, within-dose connector,
#   neutral intervals, dose colour on symbols only, right-hand across-dose
#   contrast + |E_CUE| summary, and an aligned Delta-E strip with zero line.
# Redrawn at 180 mm x 238 mm. Manual panel tags (strip untagged).
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(ggridges); library(patchwork)
})
source("scripts/fig_style.R")
DAT <- "tables/physiology"; OUT <- "tables/revision/fig1G"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
KEV <- 11604.51812; TREF <- 293.15; KB <- 8.617e-5
dose_levels <- c("Control","0.06","0.12","0.25","0.5","1","2","4")
ridge_levels <- rev(dose_levels)
norm_dose <- function(x){x<-as.character(x); x[x %in% c("0","Control","control")]<-"Control"
  sub("^0\\.50$","0.5", sub("^1\\.00$","1", x))}

raw <- read_csv(file.path(DAT,"derived_N0_R_results_with_carbon.csv"), show_col_types=FALSE) %>%
  mutate(Dose = factor(norm_dose(Dose), levels = dose_levels))
dg <- read_csv(file.path(DAT,"posterior_sharpe_schoolfield_tpc_growth_C_per_C_h_by_dose.csv"),
               show_col_types=FALSE) %>% mutate(Dose = factor(norm_dose(Dose), levels=dose_levels))
dr <- read_csv(file.path(DAT,"posterior_sharpe_schoolfield_tpc_respiration_C_per_C_h_by_dose.csv"),
               show_col_types=FALSE) %>% mutate(Dose = factor(norm_dose(Dose), levels=dose_levels))
dc <- read_csv(file.path(DAT,"posterior_arrhenius_log_CUE_common_slope.csv"),
               show_col_types=FALSE) %>% mutate(Dose = factor(norm_dose(Dose), levels=dose_levels))

# ---------- panels A-F (verbatim logic from 24_figure1.R) ----------
Tg <- seq(14.8, 28.6, length.out = 140)
ss_band <- function(d) bind_rows(lapply(split(d, d$Dose), function(dd){
  TK <- Tg + 273.15
  ln <- outer(seq_len(nrow(dd)), seq_along(TK), function(i,j)
    dd$lnB0[i] - dd$E[i]*KEV*(1/TK[j]-1/TREF) -
      log1p(exp(dd$Eh[i]*KEV*(1/dd$Th[i]-1/TK[j]))))
  q <- apply(exp(ln), 2, quantile, c(.025,.5,.975), na.rm=TRUE)
  tibble(Dose=dd$Dose[1], T=Tg, lo=q[1,], med=q[2,], hi=q[3,])}))
set.seed(1)
sub_d <- function(d,n=2000) d %>% group_by(Dose) %>% slice_sample(n=n) %>% ungroup()
bg <- ss_band(sub_d(dg)); br <- ss_band(sub_d(dr))
dose_scale_c <- scale_colour_manual(values=DOSE_RAMP, limits=dose_levels,
                                    name=expression("Prothioconazole (mg L"^-1*")"))
dose_scale_f <- scale_fill_manual(values=DOSE_RAMP, limits=dose_levels,
                                  name=expression("Prothioconazole (mg L"^-1*")"))
tpc_panel <- function(band, pts, ycol, ylab) ggplot() +
  geom_ribbon(data=band, aes(T,ymin=lo,ymax=hi,fill=Dose), alpha=0.11) +
  geom_line(data=band, aes(T,med,colour=Dose), linewidth=0.55) +
  geom_point(data=pts, aes(T,.data[[ycol]],colour=Dose), size=0.75, alpha=0.75, stroke=0) +
  scale_y_log10() + dose_scale_c + dose_scale_f +
  guides(fill="none", colour=guide_legend(nrow=1, override.aes=list(alpha=1,linewidth=1.4,size=0))) +
  labs(x="Temperature (\u00b0C)", y=ylab) + theme_pub()
pA <- tpc_panel(bg, raw, "growth_C_per_C_h", expression("Growth rate (C C"^-1*" h"^-1*")")) + labs(tag="A")
pC <- tpc_panel(br, raw, "respiration_C_per_C_h", expression("Respiration rate (C C"^-1*" h"^-1*")")) + labs(tag="C")
ridge_panel <- function(d, vcol, xlab){
  ctrl_med <- median(d[[vcol]][d$Dose=="Control"])
  ggplot(d, aes(.data[[vcol]], factor(Dose, levels=ridge_levels), fill=Dose, colour=Dose)) +
    geom_density_ridges(scale=1.25, alpha=0.55, linewidth=0.35,
                        quantile_lines=TRUE, quantiles=2, rel_min_height=0.005) +
    geom_vline(xintercept=ctrl_med, linetype="22", colour=INK2, linewidth=0.35) +
    dose_scale_c + dose_scale_f + guides(colour="none", fill="none") +
    labs(x=xlab, y=NULL) + theme_pub() + theme(legend.position="none")}
pB <- ridge_panel(dg, "lnB0", expression("ln(B"[0]*"), growth")) + labs(tag="B")
pD <- ridge_panel(dr, "lnB0", expression("ln(B"[0]*"), respiration")) + labs(tag="D")
xg <- 1/(KB*TREF) - 1/(KB*(Tg+273.15))
cue_band <- bind_rows(lapply(split(dc, dc$Dose), function(dd){
  ln <- outer(dd$alpha, xg, function(a,x) a) + outer(dd$E, xg, function(e,x) e*x)
  q <- apply(exp(ln), 2, quantile, c(.025,.5,.975), na.rm=TRUE)
  tibble(Dose=dd$Dose[1], x=xg, lo=q[1,], med=q[2,], hi=q[3,])}))
raw_cue <- raw %>% mutate(x = 1/(KB*TREF) - 1/(KB*(T+273.15)))
pE <- ggplot() +
  geom_ribbon(data=cue_band, aes(x,ymin=lo,ymax=hi,fill=Dose), alpha=0.11) +
  geom_line(data=cue_band, aes(x,med,colour=Dose), linewidth=0.55) +
  geom_point(data=raw_cue, aes(x,CUE,colour=Dose), size=0.75, alpha=0.75, stroke=0) +
  scale_y_log10() + dose_scale_c + dose_scale_f +
  guides(fill="none", colour=guide_legend(nrow=1, override.aes=list(alpha=1,linewidth=1.4,size=0))) +
  labs(x=expression(frac(1,k*T[ref]) - frac(1,k*T)), y="CUE") + theme_pub() + labs(tag="E")
pF <- ridge_panel(dc, "alpha", expression(alpha*"  (equiv. ln(CUE) at T"[ref]*")")) + labs(tag="F")

# ---------- panel G ----------
qs <- function(x) c(med=median(x), lo95=unname(quantile(x,.025)), hi95=unname(quantile(x,.975)),
                    lo66=unname(quantile(x,.17)), hi66=unname(quantile(x,.83)))
Eg <- dg %>% group_by(Dose) %>% summarise(as_tibble(as.list(qs(E)))) %>% mutate(trait="Growth")
Er <- dr %>% group_by(Dose) %>% summarise(as_tibble(as.list(qs(E)))) %>% mutate(trait="Respiration")
S <- bind_rows(Eg, Er) %>%
  mutate(xi = as.numeric(factor(Dose, levels=dose_levels)),
         x = xi + ifelse(trait=="Growth", -0.17, 0.17))
conn <- S %>% select(Dose, xi, trait, med) %>%
  pivot_wider(names_from=trait, values_from=med)
# draw-wise per-dose and across-dose contrasts (pair by .draw)
dd <- inner_join(dg %>% select(Dose,.draw,Eg=E), dr %>% select(Dose,.draw,Er=E),
                 by=c("Dose",".draw")) %>% mutate(dE = Er - Eg)
per_dose <- dd %>% group_by(Dose) %>%
  summarise(as_tibble(as.list(qs(dE))), p_pos=mean(dE>0)) %>%
  mutate(xi = as.numeric(factor(Dose, levels=dose_levels)))
across <- dd %>% group_by(.draw) %>% summarise(m=mean(dE)) %>% pull(m)
cueEd <- abs(dc$E)
sumG <- tibble(x=c(9.45,10.05), lab=c("Delta*bar(italic(E))[a]", "group('|',italic(E)[CUE],'|')"),
               med=c(median(across), median(cueEd)),
               lo95=c(quantile(across,.025), quantile(cueEd,.025)),
               hi95=c(quantile(across,.975), quantile(cueEd,.975)),
               lo66=c(quantile(across,.17), quantile(cueEd,.17)),
               hi66=c(quantile(across,.83), quantile(cueEd,.83)))
cat(sprintf("across-dose contrast: %.3f [%.3f, %.3f] | |E_CUE| %.3f [%.3f, %.3f]\n",
            sumG$med[1], sumG$lo95[1], sumG$hi95[1], sumG$med[2], sumG$lo95[2], sumG$hi95[2]))
cat("per-dose P(dE>0):", sprintf("%s=%.2f", per_dose$Dose, per_dose$p_pos), "\n")

XLIM <- c(0.45, 10.5); DIVX <- 8.85
pG <- ggplot(S) +
  geom_segment(data=conn, aes(x=xi-0.17, xend=xi+0.17, y=Growth, yend=Respiration),
               colour="#CFCFCF", linewidth=0.3) +
  geom_linerange(aes(x=x, ymin=lo95, ymax=hi95), colour=INK2, linewidth=0.3, alpha=0.75) +
  geom_linerange(aes(x=x, ymin=lo66, ymax=hi66), colour=INK2, linewidth=0.8, alpha=0.9) +
  geom_point(aes(x=x, y=med, shape=trait, fill=Dose), size=1.8, colour=INK, stroke=0.35) +
  geom_vline(xintercept=DIVX, colour="#DDDDDD", linewidth=0.4) +
  geom_linerange(data=sumG, aes(x=x, ymin=lo95, ymax=hi95), colour="#D64B21", linewidth=0.3) +
  geom_linerange(data=sumG, aes(x=x, ymin=lo66, ymax=hi66), colour="#D64B21", linewidth=0.8) +
  geom_point(data=sumG, aes(x=x, y=med), shape=23, size=1.8, fill="#D64B21", colour=INK, stroke=0.35) +
  scale_shape_manual(values=c(Growth=21, Respiration=24), name=NULL) +
  scale_fill_manual(values=DOSE_RAMP, limits=dose_levels, guide="none") +
  scale_x_continuous(breaks=c(1:8, sumG$x), limits=XLIM,
                     labels=c(dose_levels, parse(text=sumG$lab))) +
  labs(x=NULL, y=expression(italic(E)[a]*" (eV)"), tag="G") +
  guides(shape=guide_legend(override.aes=list(fill="grey65"))) +
  theme_pub() +
  theme(legend.position=c(0.115,0.96), legend.direction="horizontal",
        legend.key.size=unit(2.6,"mm"),
        axis.text.x=element_text(size=BASE-0.5))
pGs <- ggplot(per_dose) +
  geom_hline(yintercept=0, linetype="22", colour=INK2, linewidth=0.3) +
  geom_linerange(aes(x=xi, ymin=lo95, ymax=hi95), colour=INK2, linewidth=0.3, alpha=0.75) +
  geom_linerange(aes(x=xi, ymin=lo66, ymax=hi66), colour=INK2, linewidth=0.8, alpha=0.9) +
  geom_point(aes(x=xi, y=med, fill=Dose), shape=21, size=1.6, colour=INK, stroke=0.3) +
  scale_fill_manual(values=DOSE_RAMP, limits=dose_levels, guide="none") +
  scale_x_continuous(breaks=1:8, labels=dose_levels, limits=XLIM) +
  labs(x=expression("Prothioconazole (mg L"^-1*")"),
       y=expression(Delta*italic(E)[a]*" (eV)")) +
  theme_pub()

fig <- (pA + pB) / (pC + pD) / (pE + pF) / pG / pGs +
  plot_layout(guides="collect", heights=c(1,1,1,0.78,0.42)) &
  theme(legend.position="bottom",
        legend.title=element_text(size=BASE),
        plot.tag=element_text(size=LABELPT, face="bold", family=FONT))
ggsave(file.path(OUT,"Figure1_withG.png"), fig, width=180, height=238,
       units="mm", dpi=450, bg="white")
ok <- tryCatch({ ggsave(file.path(OUT,"Figure1_withG.pdf"), fig, width=180, height=238,
         units="mm", bg="white", device=cairo_pdf); TRUE }, error = function(e) FALSE)
if (!ok) message("PDF export failed in this R build - PNG written; PDF available in repo.")
cat("Figure1_withG done\n")
