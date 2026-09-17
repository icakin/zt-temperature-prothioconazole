#!/usr/bin/env Rscript
# =============================================================================
# 49_ptc_sham_review.R -- inspect every fitted curve of the prothioconazole x
# SHAM factorial, one at a time, and record any curve you reject or re-window.
#
#   Rscript scripts/49_ptc_sham_review.R          # run from the repository root
#
# The app is a reviewer, not a trim selector. The interval it shows is the one
# the rule in 46_ptc_sham_prepared_traces.R chose; you can compare rules from the
# dropdown to see how the interval moves, and you can override an interval by
# hand, but an override is written to a file of its own, carries a reason, and
# is counted in the output of 47_ptc_sham_rates.R, which reports the contrast
# both with and without overrides. Nothing you do here is invisible.
#
# Panels, per curve: grey points, raw readings; black points, the analysis-ready
# series (step offset removed, denoised); shaded band, the fitting interval;
# blue line, the fitted model over that interval, dashed where projected; the
# fitted r, K, R2 and the oxygen consumed.
#
# Reads : tables/aox/ptc_sham_prepared_traces.csv.gz, ptc_sham_windows.csv (48)
# Writes: tables/aox/ptc_sham/plot_exclude_points.csv   curves marked failed
#         tables/aox/ptc_sham/manual_overrides.csv      hand-set intervals
# =============================================================================
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
source(file.path(ROOT, "scripts/aox_common.R"))
need <- c("shiny", "ggplot2"); miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("install first: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({library(shiny); library(ggplot2)})

PREP <- file.path(ROOT, "tables/aox/ptc_sham_prepared_traces.csv.gz")
WINS <- file.path(ROOT, "tables/aox/ptc_sham_windows.csv")
if (!file.exists(PREP)) stop("run scripts/46_ptc_sham_prepared_traces.R first")
L  <- read.csv(PREP, stringsAsFactors = FALSE)
Wn <- read.csv(WINS, stringsAsFactors = FALSE)
OUT <- file.path(ROOT, "tables/aox/ptc_sham"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
EXC_CSV <- file.path(OUT, "plot_exclude_points.csv"); OVR_CSV <- file.path(OUT, "manual_overrides.csv")

RULES <- list("peak + 36 h (primary)" = list(type = "peak", dur = 36*60),
              "peak + 30 h"           = list(type = "peak", dur = 30*60),
              "peak + 24 h"           = list(type = "peak", dur = 24*60),
              "fixed 3-40 h"          = list(type = "fixed", a = 180, b = 40*60),
              "drawdown 5-80%"        = list(type = "draw",  a = 0.05, b = 0.80),
              "stable r"              = list(type = "stable_r"))

Wn$key   <- sprintf("%d_%s_R%d%s", Wn$temp, Wn$condition, Wn$replicate, Wn$well)
Wn$label <- sprintf("%d C  %-4s  culture %d  %s", Wn$temp, Wn$condition, Wn$replicate, Wn$well)
Wn <- Wn[order(Wn$temp, Wn$replicate, match(Wn$condition, PS_ORDER), Wn$well), ]
L$key <- sprintf("%d_%s_R%d%s", L$temp, L$condition, L$replicate, L$well)

rd <- function(f, cols) if (file.exists(f)) { x <- read.csv(f, stringsAsFactors = FALSE)
  if (all(cols %in% names(x))) x else NULL } else NULL
exc0 <- rd(EXC_CSV, c("T","Dose","Replicate"))
excl_init <- if (is.null(exc0)) character(0) else sprintf("%s_%s_%s", exc0$T, exc0$Dose, toupper(exc0$Replicate))
ovr_init <- rd(OVR_CSV, c("T","Dose","Replicate","fit_start","fit_end"))
if (is.null(ovr_init)) ovr_init <- data.frame(T = numeric(0), Dose = character(0), Replicate = character(0),
  fit_start = numeric(0), fit_end = numeric(0), reason = character(0), stringsAsFactors = FALSE)
if (!"reason" %in% names(ovr_init)) ovr_init$reason <- ""

ui <- fluidPage(
  titlePanel("Prothioconazole x SHAM - review every fitted curve"),
  sidebarLayout(
    sidebarPanel(width = 4,
      selectInput("curve", "Curve:", setNames(Wn$key, Wn$label)),
      fluidRow(column(6, actionButton("prev", "< Prev", width = "100%")),
               column(6, actionButton("nxt", "Next >", width = "100%"))),
      tags$hr(),
      selectInput("rule", "Interval rule (view):", names(RULES)),
      helpText("The primary rule is what 46 fits. Switching rules here only changes what you see."),
      tags$hr(),
      checkboxInput("exclude", "Reject this curve (failed well)", FALSE),
      tags$hr(),
      strong("Override the interval by hand"),
      helpText("Click the plot to set the start, then the end. Every override is saved with a reason and reported separately from the rule-based result."),
      radioButtons("mode", "Click sets:", c("Start", "End"), "Start", inline = TRUE),
      fluidRow(column(6, numericInput("ostart", "Start (min)", NA)),
               column(6, numericInput("oend", "End (min)", NA))),
      textInput("reason", "Reason for the override", ""),
      fluidRow(column(6, actionButton("save_ovr", "Save override", width = "100%")),
               column(6, actionButton("clear_ovr", "Back to rule", width = "100%"))),
      tags$hr(),
      strong(textOutput("counts")),
      div(style = "max-height:220px; overflow-y:auto;", tableOutput("ovr_tbl"))),
    mainPanel(width = 8, plotOutput("plot", height = "600px", click = "click"), verbatimTextOutput("info"))))

server <- function(input, output, session) {
  excl <- reactiveVal(excl_init); ovr <- reactiveVal(ovr_init)
  row_of <- function(k) Wn[Wn$key == k, ][1, ]
  ovr_of <- function(k) { r <- row_of(k); o <- ovr()
    o[o$T == r$temp & o$Dose == r$condition & toupper(o$Replicate) == sprintf("R%d%s", r$replicate, r$well), ][1, ] }

  observeEvent(input$curve, {
    updateCheckboxInput(session, "exclude", value = input$curve %in% excl())
    o <- ovr_of(input$curve)
    updateNumericInput(session, "ostart", value = if (nrow(o) && is.finite(o$fit_start)) o$fit_start else NA)
    updateNumericInput(session, "oend",   value = if (nrow(o) && is.finite(o$fit_end))   o$fit_end   else NA)
    updateTextInput(session, "reason", value = if (nrow(o) && !is.na(o$reason)) o$reason else "")
    updateRadioButtons(session, "mode", selected = "Start")
  }, ignoreInit = FALSE)
  observeEvent(input$prev, { i <- match(input$curve, Wn$key); if (i > 1) updateSelectInput(session, "curve", selected = Wn$key[i-1]) })
  observeEvent(input$nxt,  { i <- match(input$curve, Wn$key); if (i < nrow(Wn)) updateSelectInput(session, "curve", selected = Wn$key[i+1]) })
  observeEvent(input$exclude, { k <- input$curve
    excl(if (isTRUE(input$exclude)) union(excl(), k) else setdiff(excl(), k)) }, ignoreInit = TRUE)
  observeEvent(input$click, { v <- round(input$click$x, 1)
    if (identical(input$mode, "Start")) { updateNumericInput(session, "ostart", value = v); updateRadioButtons(session, "mode", selected = "End") }
    else updateNumericInput(session, "oend", value = v) })

  observeEvent(input$save_ovr, {
    r <- row_of(input$curve); a <- input$ostart; b <- input$oend
    if (!is.finite(a) || !is.finite(b) || b <= a) { showNotification("Set a start and a larger end first.", type = "warning"); return() }
    if (!nzchar(trimws(input$reason))) { showNotification("Give a reason: every override is reported.", type = "warning"); return() }
    o <- ovr(); rp <- sprintf("R%d%s", r$replicate, r$well)
    o <- o[!(o$T == r$temp & o$Dose == r$condition & toupper(o$Replicate) == rp), , drop = FALSE]
    ovr(rbind(o, data.frame(T = r$temp, Dose = r$condition, Replicate = rp, fit_start = a, fit_end = b,
                            reason = trimws(input$reason), stringsAsFactors = FALSE)))
    showNotification("Override saved and logged.", type = "message")
  })
  observeEvent(input$clear_ovr, { r <- row_of(input$curve); rp <- sprintf("R%d%s", r$replicate, r$well); o <- ovr()
    ovr(o[!(o$T == r$temp & o$Dose == r$condition & toupper(o$Replicate) == rp), , drop = FALSE])
    updateNumericInput(session, "ostart", value = NA); updateNumericInput(session, "oend", value = NA)
    updateTextInput(session, "reason", value = "") })

  observe({ e <- excl()
    d <- if (!length(e)) Wn[0, c("temp","condition","replicate","well")] else Wn[Wn$key %in% e, c("temp","condition","replicate","well")]
    out <- data.frame(T = d$temp, Dose = d$condition, Replicate = sprintf("R%d%s", d$replicate, d$well))
    tryCatch(write.csv(out, EXC_CSV, row.names = FALSE), error = function(x) NULL) })
  observe({ tryCatch(write.csv(ovr(), OVR_CSV, row.names = FALSE), error = function(x) NULL) })

  cur <- reactive({
    k <- input$curve; r <- row_of(k); g <- L[L$key == k, ]
    w <- ps_window(g$time_min, g$o2_prepared, RULES[[input$rule]])
    o <- ovr_of(k); manual <- nrow(o) > 0 && is.finite(o$fit_start)
    a <- if (manual) o$fit_start else w[1]; b <- if (manual) o$fit_end else w[2]
    m <- g$time_min >= a & g$time_min <= b & is.finite(g$o2_prepared)
    f <- if (sum(m) > 30) fit_o2_model(g$time_min[m], g$o2_prepared[m]) else NULL
    list(r = r, g = g, a = a, b = b, m = m, f = f, manual = manual, rule_w = w)
  })

  output$plot <- renderPlot({
    z <- cur(); g <- z$g
    p <- ggplot() +
      annotate("rect", xmin = z$a, xmax = z$b, ymin = -Inf, ymax = Inf, fill = "#EDF2F7") +
      geom_point(data = g, aes(time_min, o2_raw), colour = "grey78", size = 0.7) +
      geom_point(data = g, aes(time_min, o2_prepared), colour = "grey20", size = 0.8) +
      geom_vline(xintercept = c(z$a, z$b), colour = c("forestgreen", "red"), linewidth = 1) +
      labs(title = z$r$label, x = "Time (min)", y = expression("O"[2]*" (mg/L)"),
           subtitle = paste(c(if (z$manual) "MANUAL OVERRIDE" else sprintf("rule: %s", input$rule),
                              if (input$curve %in% excl()) "REJECTED"), collapse = "   |   ")) +
      theme_classic(13) + theme(plot.subtitle = element_text(face = "bold",
        colour = if (z$manual || input$curve %in% excl()) "red" else "grey30"))
    if (!is.null(z$f)) {
      tt <- seq(min(g$time_min), max(g$time_min), length.out = 400)
      gy <- resp_model(tt - z$a, z$f$r_per_h/60, z$f$K, z$f$O2_0)
      gd <- data.frame(tt, gy); gd <- gd[is.finite(gd$gy) & gd$gy > min(g$o2_prepared, na.rm = TRUE) - 1 &
                                          gd$gy < max(g$o2_prepared, na.rm = TRUE) + 1, ]
      p <- p + geom_line(data = gd[gd$tt < z$a | gd$tt > z$b, ], aes(tt, gy), colour = "blue", linetype = "dashed", linewidth = 0.7) +
               geom_line(data = gd[gd$tt >= z$a & gd$tt <= z$b, ], aes(tt, gy), colour = "blue", linewidth = 1.1)
    }
    p
  })
  output$info <- renderText({ z <- cur()
    if (is.null(z$f)) return("no fit for this interval")
    sprintf("interval %.0f - %.0f min (%.1f h, %d points)   r = %.4f h-1   K = %.5f   R2 = %.4f   O2 consumed = %.2f mg/L\nrule would give %.0f - %.0f min",
            z$a, z$b, (z$b - z$a)/60, sum(z$m), z$f$r_per_h, z$f$K, z$f$R2,
            diff(range(z$g$o2_prepared[z$m], na.rm = TRUE)), z$rule_w[1], z$rule_w[2]) })
  output$counts <- renderText(sprintf("rejected: %d    overridden: %d of %d", length(excl()), nrow(ovr()), nrow(Wn)))
  output$ovr_tbl <- renderTable(ovr())
}
shiny::runApp(shinyApp(ui, server), launch.browser = TRUE)
