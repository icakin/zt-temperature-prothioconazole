# fig_style.R — shared publication style (R/ggplot side).
# Mirror of fig_style.py — SAME hex values, SAME font, SAME sizes.
suppressPackageStartupMessages({library(ggplot2)})

FONT <- "DejaVu Sans"

# temperature: ordered cool->warm ramp (7 levels)
TEMPS     <- c(15, 18, 21, 24, 26, 27, 28)
TEMP_RAMP <- c("#27519E", "#457BBE", "#64A5DE", "#C9A227",
               "#E08214", "#D64B21", "#A31621")
names(TEMP_RAMP) <- as.character(TEMPS)
TEMP3 <- c("15" = "#27519E", "21" = "#64A5DE", "27" = "#D64B21")  # exact ramp members, CVD-validated

# dose: control gray + purple ramp light->dark (ordered)
DOSES     <- c("Control", "0.06", "0.12", "0.25", "0.5", "1", "2", "4")
DOSE_RAMP <- c("#8C8C8C", "#A1DA38", "#4FC36B", "#21A585", "#25848E",
               "#33638D", "#433D84", "#471063")
names(DOSE_RAMP) <- DOSES

CTRL <- "#9A9A9A"; DRUG <- "#433D84"
DIV_LOW <- "#2166AC"; DIV_MID <- "#F7F7F7"; DIV_HIGH <- "#B2182B"
GRAY_PT <- "#C8C8C8"; INK <- "#1A1A1A"; INK2 <- "#555555"

BASE <- 7; LABELPT <- 10
W2COL_MM <- 183

theme_pub <- function(base_size = BASE) {
  theme_classic(base_size = base_size, base_family = FONT) +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = INK2),
      axis.ticks = element_line(linewidth = 0.3, colour = INK2),
      axis.text = element_text(size = base_size - 0.5, colour = INK2),
      axis.title = element_text(size = base_size, colour = INK),
      legend.text = element_text(size = base_size - 0.5, colour = INK),
      legend.title = element_text(size = base_size, colour = INK),
      legend.key.size = unit(3.2, "mm"),
      legend.background = element_blank(),
      plot.title = element_text(size = base_size + 1, face = "bold", colour = INK),
      plot.tag = element_text(size = LABELPT, face = "bold", colour = INK),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold", colour = INK)
    )
}
