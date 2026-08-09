# =============================================================================
# 03_trim_selector.R
# Click on each O2 curve to set its fit-window START and END (manual trimming).
# -----------------------------------------------------------------------------
# - Shows one curve at a time (full O2 vs Time trace).
# - Grey dashed lines  = the AUTOMATIC start/end (from 02 metadata), for reference.
# - Choose "Click sets: Start" or "End", then click on the plot at the time you
#   want. Green line = your manual start, red line = your manual end.
# - Navigate curves with the dropdown or Prev/Next.
# - The set of manual windows is written to tables/manual_fit_windows.csv and
#   printed as a ready-to-paste MANUAL_FIT_WINDOWS block.
# - SELECT / DESELECT SAMPLES: tick "Don't include this sample" to discard a
#   curve (or untick to keep it). The discard list is written to
#   tables/plot_exclude_points.csv and printed as a ready-to-paste EXCLUDE_POINTS
#   block for 04_oxygen_fits.R. The currently-discarded samples (the EXCLUDE_POINTS
#   already in 04_oxygen_fits.R) are RESTORED on startup from that CSV.
# - Curves are keyed by T + Dose + Replicate (this project has no "Clade").
# - Persists: reloads tables/manual_fit_windows.csv + plot_exclude_points.csv.
#
# RUN:
#   RStudio: open this file -> "Run App".
#   Terminal: Rscript scripts/03_trim_selector.R
#
# Needs 01_longdata.R + 02_trimming.R to have run (for the data + metadata).
#
# PATHS: this app anchors strictly to THIS project. base_dir is the PARENT of
# the scripts/ folder (identical rule to 00_config.R), so every file it reads
# or writes stays inside this project's tables/ folder. There are no hard-coded
# external paths, so nothing is ever generated outside this project.
# =============================================================================

# ---- Locate data ------------------------------------------------------------
# Robustly find the folder THIS script lives in — whether launched via RStudio
# "Run App", source(), or Rscript — then derive the project root from it.
.script_dir <- local({
  d <- tryCatch({
    if (requireNamespace("rstudioapi", quietly = TRUE) &&
        rstudioapi::isAvailable() &&
        nzchar(rstudioapi::getActiveDocumentContext()$path)) {
      dirname(rstudioapi::getActiveDocumentContext()$path)
    } else NA_character_
  }, error = function(e) NA_character_)
  if (length(d) == 0 || is.na(d) || !nzchar(d)) {
    d <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
  }
  if (length(d) == 0 || is.na(d) || !nzchar(d)) d <- getwd()
  normalizePath(d, mustWork = FALSE)
})

# base_dir = parent of scripts/ (same rule as 00_config.R); tables/ lives there.
base_dir   <- dirname(.script_dir)
tables_dir <- file.path(base_dir, "tables")

# Fallback ONLY to tables/ dirs relative to the working directory — never to any
# hard-coded external location — so files can only ever land inside this project.
if (!dir.exists(tables_dir)) {
  cand <- c(file.path(getwd(), "tables"), file.path(getwd(), "..", "tables"))
  hit  <- cand[dir.exists(cand)]
  if (length(hit)) tables_dir <- normalizePath(hit[1], mustWork = FALSE)
}

find_file <- function(fname) {
  p <- file.path(tables_dir, fname)
  if (file.exists(p)) normalizePath(p, mustWork = FALSE) else NA_character_
}

long_path <- find_file("Oxygen_All_Long.csv")
meta_path <- find_file("Oxygen_Trimmed_Series_Metadata.csv")
if (is.na(long_path)) stop("Could not find Oxygen_All_Long.csv in ", tables_dir,
                           ". Run 01_longdata.R first.")
if (is.na(meta_path)) stop("Could not find Oxygen_Trimmed_Series_Metadata.csv in ",
                           tables_dir, ". Run 02_trimming.R first.")
message("Trim selector reading/writing tables in: ", tables_dir)
out_csv <- file.path(tables_dir, "manual_fit_windows.csv")

# ---- Packages ---------------------------------------------------------------
need <- c("shiny", "ggplot2", "readr", "dplyr")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0) {
  stop("Install missing package(s) first:\n  install.packages(c(",
       paste(sprintf('"%s"', miss), collapse = ", "), "))")
}
library(shiny); library(ggplot2); suppressPackageStartupMessages(library(dplyr))

# Model used for the guide curve (identical to the pipeline's resp_model).
# IMPORTANT: the guide fit here (start values + bounds in fit_resp / winfit below)
# is kept IDENTICAL to 04_oxygen_fits.R's fit_one, so the blue guide line is the
# exact model 03 will fit. If you change bounds here, change them in 03 too.
resp_model  <- function(r, K, t, O2_0) O2_0 + (K / r) * (1 - exp(r * t))
has_minpack <- requireNamespace("minpack.lm", quietly = TRUE)

# Suggested exponential-phase bounds from the RAW trace (reference markers only —
# they do NOT constrain your selection). start = time of the O2 peak (where the
# draw-down begins); end = where O2 has fallen to within 5% of its post-peak
# minimum (draw-down essentially complete). No smoothing is used.
exp_phase_bounds <- function(tt, yy) {
  ok <- is.finite(tt) & is.finite(yy)
  tt <- tt[ok]; yy <- yy[ok]
  if (length(tt) < 5) return(c(NA_real_, NA_real_))
  ord <- order(tt); tt <- tt[ord]; yy <- yy[ord]
  win <- which(tt <= min(tt) + 400)          # peak is always early; avoid late noise
  if (!length(win)) win <- seq_along(tt)
  pk      <- win[which.max(yy[win])]
  peak_t  <- tt[pk]; o2pk <- yy[pk]
  seg     <- which(tt >= peak_t)
  o2min   <- min(yy[seg])
  thr     <- o2min + 0.05 * (o2pk - o2min)
  hit     <- seg[yy[seg] <= thr]
  end_t   <- if (length(hit)) tt[hit[1]] else tt[length(tt)]
  c(peak_t, end_t)
}

# 95%-down endpoint measured from a GIVEN start time: take O2 at that start,
# find the floor after it, and return the first time O2 falls to within 5% of
# that floor. Lets the end suggestion follow wherever you click the start.
suggest_end_95 <- function(tt, yy, start_t, frac = 0.95) {
  ok <- is.finite(tt) & is.finite(yy)
  tt <- tt[ok]; yy <- yy[ok]
  ord <- order(tt); tt <- tt[ord]; yy <- yy[ord]
  seg <- which(tt >= start_t)
  if (length(seg) < 3) return(NA_real_)
  o2_start <- yy[seg[1]]
  o2min    <- min(yy[seg])
  thr      <- o2min + (1 - frac) * (o2_start - o2min)   # frac of the way down
  hit      <- seg[yy[seg] <= thr]
  if (length(hit)) tt[hit[1]] else tt[length(tt)]
}

# Fit resp_model to a (time, oxygen) vector; return coefs + R2 + RMSE.
fit_resp <- function(tt, yy) {
  out <- list(ok = FALSE, co = NULL, r2 = NA_real_, rmse = NA_real_, n = length(tt))
  if (!has_minpack || length(tt) < 6) return(out)
  t0     <- tt - min(tt)
  slope0 <- suppressWarnings(stats::median(diff(yy) / diff(t0), na.rm = TRUE))
  K0     <- if (is.finite(slope0)) max(abs(slope0), 1e-6) else 1e-3
  ft <- try(minpack.lm::nlsLM(
    Oxygen ~ resp_model(r, K, Time0, O2_0),
    data  = data.frame(Time0 = t0, Oxygen = yy),
    start = list(r = 1e-3, K = K0, O2_0 = yy[1]),
    lower = c(r = 1e-6, K = 1e-10, O2_0 = min(yy) - 1),
    upper = c(r = 0.15, K = 1,     O2_0 = max(yy) + 1),
    control = minpack.lm::nls.lm.control(maxiter = 200)
  ), silent = TRUE)
  if (inherits(ft, "try-error")) return(out)
  co    <- coef(ft)
  preds <- resp_model(co[["r"]], co[["K"]], t0, co[["O2_0"]])
  resid <- yy - preds
  sstot <- sum((yy - mean(yy))^2)
  out$ok <- TRUE; out$co <- co
  out$rmse <- sqrt(mean(resid^2))
  out$r2   <- if (sstot > 1e-12) 1 - sum(resid^2) / sstot else NA_real_
  out
}

# Auto-detect a window. Search START times from the O2 peak forward (up to
# peak + start_search minutes, never before the peak), and for each start take
# the LONGEST window that keeps R2 >= r2_target and RMSE <= rmse_max (>= min_pts
# points). Pick the overall longest such window (tie-break: higher R2). Falls
# back to a peak-anchored best-R2 window if none qualify. Returns c(start, end).
auto_detect_window <- function(tt, yy, r2_target, rmse_max, min_pts,
                               start_search = 60, start_mode = "peak",
                               start_drawdown_frac = 0.05,
                               t_peak_min = 10, t_peak_max = 200,
                               start_step = 15, end_step = 6) {
  ok <- is.finite(tt) & is.finite(yy)
  tt <- tt[ok]; yy <- yy[ok]
  ord <- order(tt); tt <- tt[ord]; yy <- yy[ord]
  if (length(tt) < min_pts + 1) return(c(NA_real_, NA_real_))
  win <- which(tt >= t_peak_min & tt <= t_peak_max)
  if (!length(win)) win <- seq_along(tt)
  pk_idx <- win[which.max(yy[win])]
  peak   <- tt[pk_idx]
  o2_pk  <- yy[pk_idx]

  if (identical(start_mode, "drawdown")) {
    # Fixed rule: start where O2 has fallen start_drawdown_frac of its peak->min
    # drop. Same point for every curve (no per-curve floating).
    seg    <- which(tt >= peak)
    o2_min <- min(yy[seg])
    target <- o2_pk - start_drawdown_frac * (o2_pk - o2_min)
    hit    <- seg[yy[seg] <= target]
    starts <- if (length(hit)) tt[hit[1]] else peak
  } else {
    # Peak mode: candidate starts from peak .. peak + start_search.
    s_pool <- tt[tt >= peak & tt <= peak + max(0, start_search)]
    if (!length(s_pool)) s_pool <- peak
    starts <- peak; last <- peak
    for (s in s_pool) if (s - last >= start_step) { starts <- c(starts, s); last <- s }
  }

  best <- list(start = NA_real_, end = NA_real_, len = -Inf, r2 = -Inf)
  fb   <- list(start = NA_real_, end = NA_real_, r2 = -Inf)   # best-R2 fallback
  for (s in starts) {
    after <- which(tt > s)
    if (length(after) < min_pts) next
    cand <- after[seq(min_pts, length(after), by = end_step)]
    for (ei in cand) {
      idx <- which(tt >= s & tt <= tt[ei])
      if (length(idx) < min_pts) next
      fr <- fit_resp(tt[idx], yy[idx])
      if (!isTRUE(fr$ok) || !is.finite(fr$r2)) next
      pass <- fr$r2 >= r2_target &&
        (!is.finite(rmse_max) || (is.finite(fr$rmse) && fr$rmse <= rmse_max))
      len <- tt[ei] - s
      if (pass && (len > best$len || (len == best$len && fr$r2 > best$r2))) {
        best <- list(start = s, end = tt[ei], len = len, r2 = fr$r2)
      }
      if (fr$r2 > fb$r2) { fb$r2 <- fr$r2; fb$start <- s; fb$end <- tt[ei] }
    }
  }
  if (is.finite(best$len) && best$len > 0) return(c(best$start, best$end))
  c(fb$start, fb$end)
}

# ---- Load -------------------------------------------------------------------
# NOTE: in THIS project the second experimental dimension is Dose (a factor with
# levels like "Control", "0.06", ... "4"), NOT Clade. Curves are therefore keyed
# by T + Dose + Replicate, matching Oxygen_All_Long.csv and the EXCLUDE_POINTS
# list in 04_oxygen_fits.R.
long <- readr::read_csv(long_path, show_col_types = FALSE) %>%
  dplyr::mutate(T = as.numeric(T), Dose = as.character(Dose),
                Replicate = toupper(as.character(Replicate)),
                key = paste(T, Dose, Replicate, sep = "_"))

meta <- readr::read_csv(meta_path, show_col_types = FALSE) %>%
  dplyr::mutate(T = as.numeric(T), Dose = as.character(Dose),
                Replicate = toupper(as.character(Replicate)),
                key = paste(T, Dose, Replicate, sep = "_"),
                main_run_start_time = as.numeric(main_run_start_time),
                steepest_drop_time  = as.numeric(steepest_drop_time)) %>%
  dplyr::select(key, auto_start = main_run_start_time, auto_end = steepest_drop_time)

# Order doses sensibly: Control first, then ascending numeric.
.dose_key <- function(x) ifelse(tolower(x) == "control", -Inf,
                                suppressWarnings(as.numeric(x)))
curves <- long %>%
  dplyr::distinct(T, Dose, Replicate, key) %>%
  dplyr::arrange(T, .dose_key(Dose), Replicate) %>%
  dplyr::mutate(label = sprintf("T=%g  Dose=%s  %s", T, Dose, Replicate))

# preload existing manual windows
init <- data.frame(key = character(0), fit_start = numeric(0), fit_end = numeric(0),
                   stringsAsFactors = FALSE)
if (file.exists(out_csv)) {
  prev <- tryCatch(readr::read_csv(out_csv, show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(prev) && all(c("T", "Dose", "Replicate") %in% names(prev))) {
    init <- prev %>%
      dplyr::mutate(key = paste(as.numeric(T), as.character(Dose), toupper(Replicate), sep = "_"),
                    fit_start = suppressWarnings(as.numeric(fit_start)),
                    fit_end   = suppressWarnings(as.numeric(fit_end))) %>%
      dplyr::select(key, fit_start, fit_end)
  }
}

# preload existing DISCARDED samples so they are restored on startup. Read from
# tables/plot_exclude_points.csv (seeded from EXCLUDE_POINTS in 04_oxygen_fits.R).
excl_csv <- file.path(tables_dir, "plot_exclude_points.csv")
excl_init <- character(0)
if (file.exists(excl_csv)) {
  pe <- tryCatch(readr::read_csv(excl_csv, show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(pe) && all(c("T", "Dose", "Replicate") %in% names(pe))) {
    excl_init <- paste(as.numeric(pe$T), as.character(pe$Dose), toupper(pe$Replicate), sep = "_")
  }
}

# ---- UI ---------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Trim selector — click to set each curve's fit window"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      selectInput("curve", "Curve:", choices = setNames(curves$key, curves$label)),
      fluidRow(column(6, actionButton("prev", "◀ Prev", width = "100%")),
               column(6, actionButton("nxt",  "Next ▶", width = "100%"))),
      br(),
      radioButtons("mode", "Click sets:", c("Start", "End"), selected = "Start", inline = TRUE),
      checkboxInput("show_guide", "Show model guide curve (blue)", value = TRUE),
      checkboxInput("show_exp", "Show suggested exponential phase (orange dotted)", value = TRUE),
      numericInput("exp_frac", "Suggested end at % of draw-down (from your start)",
                   value = 95, min = 50, max = 99.9, step = 1),
      actionButton("clear_one", "Clear this curve (back to auto)"),
      helpText("Click sets Start first (resets to 'Start' on each new curve), then switch to End. Green = your start, Red = your end, orange dashed = suggested end (computed from your start, at the % below). It's a guide only — you can still select any window you like. The long low-O2 plateau is cropped from view."),
      hr(),
      strong("Auto-detect window"),
      fluidRow(
        column(3, numericInput("r2_target", "R2 >=", value = 0.99,
                               min = 0.5, max = 1, step = 0.005)),
        column(3, numericInput("rmse_max", "RMSE <=", value = 0.08,
                               min = 0, step = 0.01)),
        column(3, numericInput("min_pts", "Min pts", value = 50, min = 6, step = 5)),
        column(3, numericInput("start_search", "Start +min", value = 60,
                               min = 0, step = 15))
      ),
      radioButtons("start_mode", "Start rule:",
                   c("Peak (+ optional forward search)" = "peak",
                     "Fixed % drawdown from peak (consistent)" = "drawdown"),
                   selected = "peak"),
      conditionalPanel(
        condition = "input.start_mode == 'drawdown'",
        numericInput("start_drawdown", "Start drawdown %", value = 5,
                     min = 0, max = 50, step = 1)
      ),
      fluidRow(
        column(6, actionButton("auto_one", "Auto THIS curve", width = "100%")),
        column(6, actionButton("auto_all", "Auto ALL curves", width = "100%"))
      ),
      helpText("Start = O2 peak; end = longest window keeping R2 >= target. You can still nudge afterwards."),
      hr(),
      strong("Trim ALL by time interval"),
      fluidRow(
        column(6, numericInput("all_start", "All start (min)", value = NA)),
        column(6, numericInput("all_end",   "All end (min)",   value = NA))
      ),
      actionButton("apply_all_time", "Apply interval to ALL curves"),
      helpText("Sets the same start/end (min) on every curve. Leave a box blank to keep that side unchanged."),
      checkboxInput("exclude_this",
                    "Don't include this sample (exclude from plots/SS fits)",
                    value = FALSE),
      hr(),
      strong(textOutput("count")),
      div(style = "max-height:200px; overflow-y:auto;", tableOutput("tbl")),
      downloadButton("dl", "Download manual_fit_windows.csv"),
      hr(),
      strong(textOutput("excl_count")),
      strong("COPY EXCLUDE_POINTS into 04_oxygen_fits.R; MANUAL_FIT_WINDOWS is saved to tables/:"),
      verbatimTextOutput("combo"),
      downloadButton("excl_dl", "Download plot_exclude_points.csv")
    ),
    mainPanel(
      width = 8,
      plotOutput("plot", height = "640px", click = "plot_click"),
      textOutput("curinfo")
    )
  )
)

# ---- Server -----------------------------------------------------------------
server <- function(input, output, session) {
  wins <- reactiveVal(init)
  excl <- reactiveVal(excl_init)

  # On each new curve: reset the click mode to "Start" and sync the exclude box.
  observeEvent(input$curve, {
    updateRadioButtons(session, "mode", selected = "Start")
    updateCheckboxInput(session, "exclude_this", value = input$curve %in% excl())
  }, ignoreInit = FALSE)

  observeEvent(input$exclude_this, {
    k <- input$curve; cur <- excl()
    excl(if (isTRUE(input$exclude_this)) union(cur, k) else setdiff(cur, k))
  }, ignoreInit = TRUE)

  get_row <- function(k) {
    w <- wins(); r <- w[w$key == k, , drop = FALSE]
    if (nrow(r) == 0) list(fit_start = NA_real_, fit_end = NA_real_)
    else list(fit_start = r$fit_start[1], fit_end = r$fit_end[1])
  }
  set_val <- function(k, side, val) {
    w <- wins()
    if (!(k %in% w$key)) w <- rbind(w, data.frame(key = k, fit_start = NA_real_, fit_end = NA_real_))
    w[w$key == k, side] <- val
    # drop fully-empty rows
    w <- w[!(is.na(w$fit_start) & is.na(w$fit_end)), , drop = FALSE]
    wins(w)
  }

  observeEvent(input$prev, {
    i <- match(input$curve, curves$key); if (i > 1) updateSelectInput(session, "curve", selected = curves$key[i - 1])
  })
  observeEvent(input$nxt, {
    i <- match(input$curve, curves$key); if (i < nrow(curves)) updateSelectInput(session, "curve", selected = curves$key[i + 1])
  })

  observeEvent(input$plot_click, {
    k <- input$curve; x <- input$plot_click$x
    if (is.null(x)) return()
    set_val(k, if (input$mode == "Start") "fit_start" else "fit_end", round(x, 1))
  })
  observeEvent(input$clear_one, {
    w <- wins(); wins(w[w$key != input$curve, , drop = FALSE])
  })

  observeEvent(input$apply_all_time, {
    s <- suppressWarnings(as.numeric(input$all_start))
    e <- suppressWarnings(as.numeric(input$all_end))
    if (!is.finite(s) && !is.finite(e)) {
      showNotification("Enter a start and/or end time first.", type = "warning")
      return()
    }
    w <- wins()
    for (k in curves$key) {
      if (!(k %in% w$key))
        w <- rbind(w, data.frame(key = k, fit_start = NA_real_, fit_end = NA_real_))
      if (is.finite(s)) w[w$key == k, "fit_start"] <- s
      if (is.finite(e)) w[w$key == k, "fit_end"]   <- e
    }
    w <- w[!(is.na(w$fit_start) & is.na(w$fit_end)), , drop = FALSE]
    wins(w)
    showNotification("Applied time interval to all curves.", type = "message")
  })

  observeEvent(input$auto_one, {
    k <- input$curve
    d <- long %>% dplyr::filter(key == k, is.finite(Time), is.finite(Oxygen)) %>%
      dplyr::arrange(Time)
    sd_frac <- (if (is.null(input$start_drawdown)) 5 else input$start_drawdown) / 100
    se <- auto_detect_window(d$Time, d$Oxygen, input$r2_target, input$rmse_max,
                             input$min_pts, input$start_search,
                             start_mode = input$start_mode, start_drawdown_frac = sd_frac)
    if (all(is.finite(se))) {
      set_val(k, "fit_start", round(se[1], 1))
      set_val(k, "fit_end",   round(se[2], 1))
    } else {
      showNotification("Auto-detect could not find a valid window for this curve.",
                       type = "warning")
    }
  })

  observeEvent(input$auto_all, {
    n <- nrow(curves)
    w <- wins()
    sd_frac <- (if (is.null(input$start_drawdown)) 5 else input$start_drawdown) / 100
    withProgress(message = "Auto-detecting all curves...", value = 0, {
      for (i in seq_len(n)) {
        k <- curves$key[i]
        d <- long %>% dplyr::filter(key == k, is.finite(Time), is.finite(Oxygen)) %>%
          dplyr::arrange(Time)
        se <- auto_detect_window(d$Time, d$Oxygen, input$r2_target, input$rmse_max,
                                 input$min_pts, input$start_search,
                                 start_mode = input$start_mode, start_drawdown_frac = sd_frac)
        if (all(is.finite(se))) {
          if (!(k %in% w$key))
            w <- rbind(w, data.frame(key = k, fit_start = NA_real_, fit_end = NA_real_))
          w[w$key == k, "fit_start"] <- round(se[1], 1)
          w[w$key == k, "fit_end"]   <- round(se[2], 1)
        }
        incProgress(1 / n)
      }
    })
    w <- w[!(is.na(w$fit_start) & is.na(w$fit_end)), , drop = FALSE]
    wins(w)
    showNotification("Auto-detect ALL done. Review and tweak as needed.", type = "message")
  })

  # Fit resp_model to the CURRENT window and report how well the raw data matches
  # it (R2 and RMSE). Shared by the plot guide and the readout below the plot.
  winfit <- reactive({
    k <- input$curve
    d <- long %>% dplyr::filter(key == k, is.finite(Time), is.finite(Oxygen)) %>%
      dplyr::arrange(Time)
    a <- meta[meta$key == k, , drop = FALSE]
    r <- get_row(k)
    auto_start <- if (nrow(a)) a$auto_start[1] else NA_real_
    auto_end   <- if (nrow(a)) a$auto_end[1]   else NA_real_
    eff_start  <- if (is.finite(r$fit_start)) r$fit_start else auto_start
    eff_end    <- if (is.finite(r$fit_end))   r$fit_end   else auto_end
    res <- list(ok = FALSE, eff_start = eff_start, eff_end = eff_end,
                co = NULL, rmse = NA_real_, r2 = NA_real_, n = 0L)
    if (has_minpack && is.finite(eff_start) && is.finite(eff_end) && eff_end > eff_start) {
      dw <- d[d$Time >= eff_start & d$Time <= eff_end, , drop = FALSE]
      if (nrow(dw) >= 6) {
        t0     <- dw$Time - min(dw$Time)
        slope0 <- suppressWarnings(stats::median(diff(dw$Oxygen) / diff(t0), na.rm = TRUE))
        K0     <- if (is.finite(slope0)) max(abs(slope0), 1e-6) else 1e-3
        ft <- try(minpack.lm::nlsLM(
          Oxygen ~ resp_model(r, K, Time0, O2_0),
          data  = data.frame(Time0 = t0, Oxygen = dw$Oxygen),
          start = list(r = 1e-3, K = K0, O2_0 = dw$Oxygen[1]),
          lower = c(r = 1e-6, K = 1e-10, O2_0 = min(dw$Oxygen) - 1),
          upper = c(r = 0.15, K = 1,     O2_0 = max(dw$Oxygen) + 1),
          control = minpack.lm::nls.lm.control(maxiter = 200)
        ), silent = TRUE)
        if (!inherits(ft, "try-error")) {
          co    <- coef(ft)
          preds <- resp_model(co[["r"]], co[["K"]], t0, co[["O2_0"]])
          resid <- dw$Oxygen - preds
          sstot <- sum((dw$Oxygen - mean(dw$Oxygen))^2)
          res$ok   <- TRUE
          res$co   <- co
          res$rmse <- sqrt(mean(resid^2))
          res$r2   <- if (sstot > 1e-12) 1 - sum(resid^2) / sstot else NA_real_
          res$n    <- nrow(dw)
        }
      }
    }
    res
  })

  output$plot <- renderPlot({
    k <- input$curve
    d <- long %>% dplyr::filter(key == k, is.finite(Time), is.finite(Oxygen)) %>%
      dplyr::arrange(Time)
    a <- meta[meta$key == k, , drop = FALSE]
    r <- get_row(k)
    is_excl <- k %in% excl()

    auto_start <- if (nrow(a)) a$auto_start[1] else NA_real_
    auto_end   <- if (nrow(a)) a$auto_end[1]   else NA_real_
    start_ref  <- if (is.finite(r$fit_start)) r$fit_start else auto_start
    eff_start  <- start_ref
    eff_end    <- if (is.finite(r$fit_end)) r$fit_end else auto_end

    # Suggested END — recomputed LIVE from wherever Start currently is (your
    # clicked start if set, else the peak), at the % below. Moves on every click.
    frac     <- (if (is.null(input$exp_frac) || !is.finite(input$exp_frac)) 95 else input$exp_frac) / 100
    peak_t   <- exp_phase_bounds(d$Time, d$Oxygen)[1]
    anchor_s <- if (is.finite(r$fit_start)) r$fit_start else peak_t
    sugg_end <- if (isTRUE(input$show_exp)) suggest_end_95(d$Time, d$Oxygen, anchor_s, frac) else NA_real_

    # Crop the long low-O2 plateau, but keep room past the end so the projected
    # guide curve AND the suggested-end guide stay on screen.
    o2min <- min(d$Oxygen); o2max <- max(d$Oxygen); yrng <- o2max - o2min
    thr   <- o2min + 0.03 * yrng
    below <- d$Time[d$Oxygen <= thr]
    t_bot <- if (length(below)) min(below) else max(d$Time)
    x_hi  <- min(max(d$Time),
                 max(t_bot + 60, eff_start + 60, auto_end + 60, eff_end + 90,
                     if (is.finite(sugg_end)) sugg_end + 90 else -Inf, na.rm = TRUE),
                 na.rm = TRUE)

    # Guide curve: use the shared window fit (winfit). Draw it SOLID within the
    # window and DASHED beyond it (back- and forward-projection).
    wf <- winfit()
    guide <- NULL
    if (isTRUE(input$show_guide) && isTRUE(wf$ok)) {
      co <- wf$co
      gx <- seq(min(d$Time), x_hi, length.out = 400)
      gy <- resp_model(co[["r"]], co[["K"]], gx - eff_start, co[["O2_0"]])
      guide <- data.frame(Time = gx, Oxygen = gy)
      guide <- guide[is.finite(guide$Oxygen) &
                       guide$Oxygen >= o2min - 0.12 * yrng &
                       guide$Oxygen <= o2max + 0.12 * yrng, , drop = FALSE]
    }

    p <- ggplot(d, aes(Time, Oxygen)) +
      geom_point(size = 1.1, alpha = 0.7, colour = if (is_excl) "grey65" else "black") +
      labs(title = curves$label[match(k, curves$key)],
           subtitle = if (is_excl) "EXCLUDED - not in plots / SS fits" else NULL,
           x = "Time (min)", y = expression("O"[2]*" (mg/L)")) +
      coord_cartesian(xlim = c(min(d$Time), x_hi),
                      ylim = c(o2min - 0.12 * yrng, o2max + 0.05 * yrng)) +
      theme_classic(13) +
      theme(plot.subtitle = ggplot2::element_text(colour = "red", face = "bold"))
    if (!is.null(guide) && nrow(guide) > 1) {
      g_in   <- guide[guide$Time >= eff_start & guide$Time <= eff_end, , drop = FALSE]  # solid
      g_pre  <- guide[guide$Time <= eff_start, , drop = FALSE]   # back-projection (dashed)
      g_post <- guide[guide$Time >= eff_end,   , drop = FALSE]   # forward projection (dashed)
      if (nrow(g_pre) > 1)
        p <- p + geom_line(data = g_pre,  aes(Time, Oxygen), colour = "blue",
                           linewidth = 0.8, linetype = "dashed")
      if (nrow(g_post) > 1)
        p <- p + geom_line(data = g_post, aes(Time, Oxygen), colour = "blue",
                           linewidth = 0.8, linetype = "dashed")
      if (nrow(g_in) > 1)
        p <- p + geom_line(data = g_in,   aes(Time, Oxygen), colour = "blue", linewidth = 1)
    }
    # Only guide drawn: the suggested END (orange), live from the current Start,
    # with a label marking it.
    if (is.finite(sugg_end)) {
      p <- p +
        geom_vline(xintercept = sugg_end, colour = "darkorange",
                   linetype = "dashed", linewidth = 1.1) +
        annotate("label", x = sugg_end, y = o2max, hjust = 0.5, vjust = 1,
                 label = sprintf("%g%% end\n%.0f min", frac * 100, sugg_end),
                 colour = "darkorange", fill = "white", label.size = 0, size = 3.1)
    }
    if (is.finite(r$fit_start))
      p <- p + geom_vline(xintercept = r$fit_start, colour = "forestgreen", linewidth = 1.1)
    if (is.finite(r$fit_end))
      p <- p + geom_vline(xintercept = r$fit_end, colour = "red", linewidth = 1.1)
    p
  })

  output$curinfo <- renderText({
    r  <- get_row(input$curve)
    wf <- winfit()
    qual <- if (isTRUE(wf$ok)) {
      sprintf("FIT MATCH:  R2 = %.4f   |   RMSE = %.4f   (n = %d points)",
              wf$r2, wf$rmse, wf$n)
    } else {
      "FIT MATCH:  (need >= 6 points between start and end)"
    }
    # Suggested end at frac% of draw-down measured from your start (else the peak).
    d <- long %>% dplyr::filter(key == input$curve, is.finite(Time), is.finite(Oxygen))
    frac   <- (if (is.null(input$exp_frac) || !is.finite(input$exp_frac)) 95 else input$exp_frac) / 100
    peak_t <- exp_phase_bounds(d$Time, d$Oxygen)[1]
    anchor <- if (is.finite(r$fit_start)) r$fit_start else peak_t
    sugg   <- suggest_end_95(d$Time, d$Oxygen, anchor, frac)
    sugg_txt <- if (is.finite(sugg))
      sprintf("SUGGESTED END (%g%% down from %s): %.1f min",
              frac * 100, if (is.finite(r$fit_start)) "your start" else "peak", sugg)
    else "SUGGESTED END: (n/a)"
    sprintf("Start: %s    End: %s    (blank = automatic)\n%s\n%s",
            ifelse(is.finite(r$fit_start), r$fit_start, "auto"),
            ifelse(is.finite(r$fit_end),   r$fit_end,   "auto"),
            sugg_txt, qual)
  })

  win_df <- reactive({
    w <- wins()
    if (nrow(w) == 0) return(curves[0, c("T", "Dose", "Replicate")] %>%
                               dplyr::mutate(fit_start = numeric(0), fit_end = numeric(0)))
    curves %>% dplyr::inner_join(w, by = "key") %>%
      dplyr::arrange(T, .dose_key(Dose), Replicate) %>%
      dplyr::select(T, Dose, Replicate, fit_start, fit_end)
  })

  output$count <- renderText(sprintf("Manual windows set: %d curve(s)", nrow(win_df())))
  output$tbl <- renderTable(win_df(), digits = 1, na = "auto")

  observe({ tryCatch(readr::write_csv(win_df(), out_csv), error = function(e) NULL) })

  # Discarded samples (excluded from plots / SS fits). Shared with EXCLUDE_POINTS.
  excl_df <- reactive({
    keys <- excl()
    if (length(keys) == 0) return(curves[0, c("T", "Dose", "Replicate")])
    curves %>% dplyr::filter(key %in% keys) %>%
      dplyr::distinct(T, Dose, Replicate) %>%
      dplyr::arrange(T, .dose_key(Dose), Replicate)
  })
  output$excl_count <- renderText(sprintf("Excluded samples: %d", nrow(excl_df())))
  observe({ tryCatch(readr::write_csv(excl_df(), excl_csv), error = function(e) NULL) })
  output$excl_dl <- downloadHandler(
    filename = function() "plot_exclude_points.csv",
    content = function(file) readr::write_csv(excl_df(), file)
  )

  output$combo <- renderText({
    w <- win_df(); e <- excl_df()
    fmt <- function(v) ifelse(is.na(v), "NA", format(v, trim = TRUE))

    # EXCLUDE_POINTS block — same shape as the tibble in 04_oxygen_fits.R.
    if (nrow(e) == 0) {
      excl_txt <- paste0("EXCLUDE_POINTS <- tibble(\n",
        "  T = numeric(0), Dose = character(0), Replicate = character(0)\n)")
    } else {
      excl_txt <- paste0(
        "EXCLUDE_POINTS <- tibble(\n",
        "  T         = c(", paste(e$T, collapse = ", "), "),\n",
        "  Dose      = c(", paste(sprintf('"%s"', e$Dose), collapse = ", "), "),\n",
        "  Replicate = c(", paste(sprintf('"%s"', e$Replicate), collapse = ", "), ")\n)")
    }

    if (nrow(w) == 0) {
      win_txt <- paste0("MANUAL_FIT_WINDOWS <- data.frame(\n",
        "  T = numeric(0), Dose = character(0), Replicate = character(0),\n",
        "  fit_start = numeric(0), fit_end = numeric(0),\n",
        "  stringsAsFactors = FALSE\n)")
    } else {
      win_txt <- paste0(
        "MANUAL_FIT_WINDOWS <- data.frame(\n",
        "  T         = c(", paste(w$T, collapse = ", "), "),\n",
        "  Dose      = c(", paste(sprintf('"%s"', w$Dose), collapse = ", "), "),\n",
        "  Replicate = c(", paste(sprintf('"%s"', w$Replicate), collapse = ", "), "),\n",
        "  fit_start = c(", paste(vapply(w$fit_start, fmt, ""), collapse = ", "), "),\n",
        "  fit_end   = c(", paste(vapply(w$fit_end,   fmt, ""), collapse = ", "), "),\n",
        "  stringsAsFactors = FALSE\n)")
    }

    paste0("# --- Paste into 04_oxygen_fits.R (replaces EXCLUDE_POINTS) ---\n",
           excl_txt, "\n\n", win_txt)
  })

  output$dl <- downloadHandler(
    filename = function() "manual_fit_windows.csv",
    content = function(file) readr::write_csv(win_df(), file)
  )
}

shinyApp(ui, server)
