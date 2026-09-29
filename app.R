###############################################################################
# MALO Risk Prediction — Shiny App
# Fine-Gray subdistribution hazard model (competing risks) for Major Adverse
# Liver Outcomes in at-risk SLD with low fibrosis burden.
#
# HOW A SHINY APP WORKS — three-paragraph orientation for new users
# ─────────────────────────────────────────────────────────────────
# A Shiny app is split into two pieces that talk to each other:
#
# 1.  ui  (User Interface) — defines the *layout and widgets* the browser
#     shows.  Think of it as the HTML template: input controls on the left,
#     output placeholders on the right.  It runs once, when the page loads.
#
# 2.  server — contains the *R logic* that runs in response to user actions.
#     Every time the user clicks a button or changes a slider, Shiny re-runs
#     the relevant reactive blocks inside server and pushes new values into
#     the output placeholders defined in ui.
#
# The call shinyApp(ui, server) at the bottom ties them together.
#
# WHAT CHANGED FROM THE PREVIOUS VERSION
# ──────────────────────────────────────
# This app previously called a ~1 GB random survival forest (randomForestSRC)
# downloaded from Dropbox, and explained each prediction with Kernel SHAP
# against a background sample of the training data.  It now uses the linear
# Fine-Gray model exported by src/02_FineGray_main_analysis.R, which ships as a
# 5 KB file inside this repository.  Consequences:
#
#   * no download, no background dataset, no individual-level UK Biobank data
#     anywhere in the deployment — the model is coefficients plus a baseline
#     cumulative subdistribution hazard;
#   * prediction is a closed-form formula, so it is instantaneous (no progress
#     bar or button-locking needed);
#   * dependencies drop to shiny + bslib;
#   * SHAP is replaced by an exact counterfactual attribution — see the
#     "RISK ATTRIBUTION" note in section 4 below.
#
# REQUIRED FILES
# ──────────────
#   model/FG_clinical_model_deploy.rds  — written by src/02_FineGray_main_analysis.R
###############################################################################


# ── 1.  Packages ──────────────────────────────────────────────────────────────
# shiny: the web-app framework
# bslib: Bootstrap 5 theme helpers (cards, layouts, colour themes)
library(shiny)
library(bslib)


# ── 2.  Load the model at start-up (OUTSIDE the server function) ──────────────
# Placing readRDS() here means R loads the model into memory *once* when the
# app process starts and shares that single copy across all user sessions.
#
# The saved object is a plain list assembled in Section 6 of the analysis
# script.  The pieces this app uses:
#
#   $predictors    chr[6]      predictor names, in model-matrix order
#   $coefficients  num[6]      named log-subdistribution-hazard ratios
#   $baseline      data.frame  H0(t) on a fixed 0.05-year grid (t = 0 … 19)
#   $preprocessing list        BMI winsorising bounds, alcohol cap
#   $xlevels       list        factor codings for sex / has_t2dm / smoking
#   $input_range   list        training 1st–99th percentiles, for range warnings
#   $vcov          num[6,6]    covariance of the coefficients (for risk CIs)
#   $risk_cut_10y  num         median predicted 10-year risk in the cohort
#   $cohort_summary, $cv_performance   aggregate metadata shown in the footer
#   $predict_fn    function    CIF(t|x) = 1 - exp(-H0(t) * exp(x %*% beta))
#
# $predict_fn does its own factor coercion, BMI winsorising and alcohol capping,
# but it does NOT clamp age and cannot floor BMI, so the app applies the input
# policy in section 3b (prepare_inputs) *before* calling it.  The fitted model
# object itself is unchanged.
MODEL_PATH <- "model/FG_clinical_model_deploy.rds"
if (!file.exists(MODEL_PATH))
  stop("Model file not found: ", normalizePath(MODEL_PATH, mustWork = FALSE))
model <- readRDS(MODEL_PATH)
message("Loaded ", model$model_name, " (created ", model$date_created, ")")


# ── 3.  Clinical constants and display metadata ───────────────────────────────

# Surveillance decision threshold, as an absolute 10-year MALO risk in percent.
# NOTE: the model object also carries model$risk_cut_10y — the median predicted
# 10-year risk in the training cohort (0.47%), i.e. the data-driven cut-point
# between standard and aggressive surveillance used in the analysis.  The app
# deliberately keeps the round guideline-facing 0.5% figure; change the line
# below to 100 * model$risk_cut_10y to switch to the cohort median instead.
RISK_THRESHOLD_10Y <- 0.5

# Alcohol-assessment banner thresholds (g/week), independent of the model.
# NOT a surveillance override: a sex*alcohol interaction was tested as a
# sensitivity analysis (src/07_sex_alcohol_interaction.R) and did not support
# raising risk more steeply in women, so the model and its recommendation are
# unchanged. These thresholds match has_excess_alcohol in
# src/01_preprocess_clean_ukbb_cohort.R and simply flag a patient for a
# clinical alcohol assessment; the calculated risk and the surveillance
# recommendation below are always the model's, never overridden by this flag.
ALCOHOL_ASSESSMENT_CUT <- c(Male = 210, Female = 140)

# BMI is a linear term in the model, so a BMI below 25 would lower the predicted
# risk (BMI 20 vs 25: x0.71).  Lean SLD is not low-risk, so BMI is floored at
# this value in the linear predictor; patients at or below it are "at target".
# This must equal OPTIMAL$bmi below so that the BMI counterfactual is exactly x1.
BMI_FLOOR <- 25

# Optimal profile for the modifiable factors.  Each modifiable factor's
# contribution is measured by moving it to its value here while the patient's
# other factors — including age and sex, which are not modifiable and are held
# at the patient's own values — stay where they are.
OPTIMAL <- list(
  bmi                = BMI_FLOOR,
  has_t2dm           = "0",
  alcohol_grams_week = 0,
  smoking_binary     = "0"
)
MODIFIABLE <- names(OPTIMAL)

VAR_LABELS <- c(
  age                = "Age",
  sex                = "Sex",
  bmi                = "BMI",
  has_t2dm           = "Type 2 diabetes",
  alcohol_grams_week = "Alcohol intake",
  smoking_binary     = "Smoking status"
)

# Human-readable rendering of a predictor value, for the contributions table.
fmt_value <- function(var, x) {
  switch(var,
    age                = paste0(round(as.numeric(x)), " years"),
    sex                = if (as.character(x) == "1") "Male" else "Female",
    bmi                = paste0(format(round(as.numeric(x), 1), nsmall = 1), " kg/m²"),
    has_t2dm           = if (as.character(x) == "1") "Yes" else "No",
    alcohol_grams_week = paste0(round(as.numeric(x)), " g/week"),
    smoking_binary     = if (as.character(x) == "1") "Current" else "Non-current",
    as.character(x)
  )
}

# Risks here are small (often well under 1%), so a flat round(x, 2) would print
# "0.08%" as "0.08%" but "0.004%" as "0%".  Widen the precision for small values.
fmt_pct <- function(x, suffix = "%") {
  if (is.na(x)) return("—")
  paste0(sprintf(if (x == 0) "%.2f" else if (abs(x) < 0.1) "%.3f" else "%.2f", x), suffix)
}


# ── 3b.  Input policy: clamp to observed bounds, floor BMI ────────────────────
#
# One policy for every continuous input (age, BMI, alcohol): clamp to the bounds
# observed in training and predict from the clamped value, never extrapolate.
#   age      model$input_range$age            (training 1st–99th percentiles)
#   bmi      model$preprocessing$bmi_winsor   (training 1st–99th percentiles)
#   alcohol  0 … model$preprocessing$alcohol_cap
# The single exception is BMI, which is then floored at BMI_FLOOR (see above).
#
# This is applied once per prediction, before anything else touches the patient,
# so the risk, its CI and the counterfactual table all use the same values.
# Returns the prepared patient plus what was done, for the notes and table.
prepare_inputs <- function(model, patient) {
  bounds <- list(
    age                = model$input_range$age,
    bmi                = unname(model$preprocessing$bmi_winsor),
    alcohol_grams_week = c(0, model$preprocessing$alcohol_cap)
  )
  out     <- patient[1L, , drop = FALSE]
  clamped <- setNames(rep(FALSE, length(bounds)), names(bounds))
  for (v in names(bounds)) {
    raw       <- as.numeric(patient[[v]][1])
    out[[v]]  <- min(max(raw, bounds[[v]][1]), bounds[[v]][2])
    clamped[v] <- out[[v]] != raw
  }

  # BMI: the floor supersedes the lower clamp (both bounds sit below it), so
  # only report an upper-bound clamp; a low BMI is reported as floored instead.
  floored <- out$bmi < BMI_FLOOR
  out$bmi <- max(out$bmi, BMI_FLOOR)
  clamped[["bmi"]] <- as.numeric(patient$bmi[1]) > bounds$bmi[2]

  list(patient = out, raw = patient[1L, , drop = FALSE],
       clamped = clamped, bounds = bounds, bmi_floored = floored)
}


# ── 4.  Risk attribution (the replacement for SHAP) ───────────────────────────
#
# RISK ATTRIBUTION — why this is not SHAP, and why it does not need to be
# ──────────────────────────────────────────────────────────────────────
# Kernel SHAP was needed for the random forest because a forest has no closed
# form: the only way to learn what a feature did was to resample a background
# dataset and watch the prediction move.  The Fine-Gray model is linear in the
# log-subdistribution-hazard, so the same question has an exact answer:
#
#     CIF(t | x) = 1 - exp(-H0(t) * exp(sum(beta_j * x_j)))
#
# Setting modifiable predictor j to its optimal value opt_j and leaving the rest
# alone (age and sex stay at the patient's own values) multiplies the
# subdistribution hazard by exactly exp(beta_j * (x_j - opt_j)), independently
# of the other predictors and of t.  So for each modifiable factor we compute
# two numbers, both exact and both from one call to $predict_fn:
#
#   * risk multiple  — exp(beta_j * (x_j - opt_j)): the factor by which this
#     patient's value of factor j multiplies their risk relative to the same
#     patient at the optimal value.  These multiply exactly: the product across
#     the four modifiable factors is the patient's total risk multiple versus
#     the same patient at the fully optimal profile.  Recovered from the
#     predictions as log(1-r_patient) / log(1-r_j), because
#     -log(1 - CIF) = H0(t) * exp(lp) and H0(t) cancels.
#
#   * delta  — r_patient - r_j: how many percentage points of this patient's
#     absolute 10-year risk are attributable to factor j sitting where it does
#     rather than at its optimal value, i.e. the risk reduction that reaching
#     the optimal value would achieve, holding everything else fixed.
#
# The deltas are not additive (absolute risk is a nonlinear function of the
# linear predictor), which is why the table below reports them alongside the
# multiples rather than as a decomposition that sums to the total.
#
# Everything routes through model$predict_fn, so there is exactly one
# implementation of the risk equation in the deployment.  `prep` is the output
# of prepare_inputs(): the patient row and the counterfactual rows all start
# from the already-clamped, BMI-floored values.  Because the floored BMI equals
# OPTIMAL$bmi for any patient at or below 25, that row's multiple is exactly x1
# and its delta exactly 0 ("at target").
attribute_risk <- function(model, prep, time = 10) {
  patient <- prep$patient
  preds   <- MODIFIABLE

  # Row 1 is the patient; rows 2..5 are the patient with a single modifiable
  # predictor moved to its optimal value.  One predict_fn call covers all of them.
  nd <- patient[rep(1L, 1L + length(preds)), , drop = FALSE]
  for (i in seq_along(preds)) nd[[preds[i]]][i + 1L] <- OPTIMAL[[preds[i]]]

  r     <- model$predict_fn(model, nd, times = time)[, 1]
  r_pat <- r[1]
  r_ref <- setNames(r[-1], preds)

  # The table must show the value the model actually used, not the raw input:
  # a raw 900 g/week would sit next to a multiple computed at 500 g/week.  The
  # exception is a floored BMI, where the entered value is shown (the patient is
  # at target) and marked with a double dagger.
  show_value <- function(v) {
    if (v == "bmi" && prep$bmi_floored)
      return(paste0(fmt_value(v, prep$raw$bmi), "\u2021"))
    txt <- fmt_value(v, patient[[v]][1])
    if (isTRUE(unname(prep$clamped[v]))) paste0(txt, "\u2020") else txt
  }

  data.frame(
    var       = preds,
    label     = unname(VAR_LABELS[preds]),
    value     = vapply(preds, show_value, character(1)),
    opt_value = vapply(preds, function(v) fmt_value(v, OPTIMAL[[v]]), character(1)),
    # exp(beta_j * (x_j - opt_j)); > 1 raises risk, < 1 lowers it
    multiple  = log1p(-r_pat) / log1p(-r_ref),
    delta_pp  = 100 * (r_pat - r_ref),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


# ── 4b.  Confidence intervals for the predicted risk ──────────────────────────
#
# The 95% CI propagates the sampling uncertainty of the coefficient vector
# (model$vcov, the covariance matrix from the Fine-Gray fit) through the linear
# predictor:  se(lp) = sqrt(x' V x),  lp ± z * se(lp),  then the same monotone
# transform to the probability scale.  Since 1 - CIF = exp(-H0 * exp(lp)),
# shifting lp by d gives  1 - CIF_new = (1 - CIF)^exp(d), so the bounds come
# straight from the point estimate and H0 never has to be touched again.
#
# What it does NOT include: uncertainty in the baseline cumulative hazard H0(t).
# The Fine-Gray fit (cmprsk::crr) returns H0 as a point estimate with no
# variance, so it cannot be recovered from the deployed model.  With ~1,800
# events H0 is estimated far more tightly than the coefficients, so this is a
# small omission, but the interval is best read as a slightly optimistic one.
# Nor does it capture model misspecification or performance in other cohorts.
#
# `patient` must already have passed through prepare_inputs().
risk_ci <- function(model, patient, times = c(5, 10), level = 0.95) {
  nd <- patient[1L, , drop = FALSE]
  for (v in names(model$xlevels))
    nd[[v]] <- factor(as.character(nd[[v]]), levels = model$xlevels[[v]])

  X <- model.matrix(reformulate(model$predictors), nd)[, names(model$coefficients),
                                                       drop = FALSE]
  stopifnot(identical(dim(model$vcov), rep(length(model$coefficients), 2L)))
  se_lp <- sqrt(drop(X %*% model$vcov %*% t(X)))
  z     <- qnorm(1 - (1 - level) / 2)

  r <- model$predict_fn(model, patient, times = times)[1, ]
  list(lower = 1 - (1 - r)^exp(-z * se_lp),
       upper = 1 - (1 - r)^exp( z * se_lp),
       se_lp = se_lp)
}

# ── 5.  UI ────────────────────────────────────────────────────────────────────
# page_sidebar() gives a fixed left sidebar and a scrollable main area.
# All widget functions follow the pattern:
#   widgetType(inputId, label, ...)
# The inputId is the name you use in server to read the widget's value
# as input$<inputId>.

ui <- page_sidebar(
  title = "Major Adverse Liver Outcomes (MALO) Risk Calculator: Early At-Risk SLD with Low Fibrosis Burden",
  theme = bs_theme(bootswatch = "flatly"),

  # page_sidebar() always wraps the page in bslib's page_fillable(), which
  # pins html/body to exactly one screen's height regardless of this
  # `fillable` argument (that argument only affects the main content area's
  # internal flex behavior, per bslib's own source) -- so the sidebar-layout
  # grid inherits a fixed height and its main panel (div.main) scrolls inside
  # its own cramped little box the moment results don't fit. fillable = FALSE
  # turns off the fill-item flex treatment on the main panel (a smaller
  # improvement); the CSS overrides below (in tags$head) do the rest, letting
  # the main panel grow to its natural content height so the page itself gets
  # longer and scrolls normally. The sidebar keeps its own width/behavior.
  fillable = FALSE,

  # Extra vertical breathing room in the MAIN panel only (the sidebar is left
  # alone). bslib's "flatly" theme is fairly compact by default, which crowds
  # the risk numbers, the contributions table and the "About this model"
  # accordion together on taller/wider screens. This loosens spacing between
  # and inside those blocks without changing the sidebar's layout.
  tags$head(tags$style(HTML("
    /* Let the page grow to its content's natural height and scroll normally,
       instead of bslib's page_fillable() pinning it to one screen and making
       div.main scroll inside its own short internal box. The sidebar's own
       column/width is untouched by any of this. */
    html, body.bslib-page-sidebar          { height: auto !important; min-height: 100vh !important; overflow-y: visible !important; }
    main.bslib-page-main                   { height: auto !important; }
    .bslib-sidebar-layout                  { height: auto !important; }
    .bslib-sidebar-layout > .main          { height: auto !important; max-height: none !important; overflow-y: visible !important; }

    #results_ui .card       { margin-bottom: 1.75rem !important; }
    #results_ui .card-body  { padding: 1.5rem 1.75rem !important; }
    #results_ui .card-body table.table td,
    #results_ui .card-body table.table th { padding-top: .6rem !important; padding-bottom: .6rem !important; }
    #results_ui hr          { margin-top: 1.5rem !important; margin-bottom: 1.5rem !important; }
    #results_ui .alert      { padding: 1rem 1.25rem !important; margin-bottom: 1.5rem !important; }
    #results_ui .accordion  { margin-top: 1.75rem !important; }
    #results_ui p           { margin-bottom: .75rem !important; }
  "))),

  # ── Left sidebar: input controls ──────────────────────────────────────────
  sidebar = sidebar(
    width = 310,

    h5("Patient Characteristics", class = "mt-1 mb-3"),

    numericInput(
      inputId = "age",
      label   = "Age (years)",
      value   = 55, min = 18, max = 100, step = 1
    ),

    selectInput(
      inputId = "sex",
      label   = "Sex",
      choices = c("Male", "Female")
    ),

    numericInput(
      inputId = "bmi",
      label   = "BMI (kg/m²)",
      value   = 27.0, min = 10, max = 70, step = 0.1
    ),

    selectInput(
      inputId = "t2dm",
      label   = "Type 2 Diabetes",
      # The *value* sent to server is "0" or "1"; the *label* is what the
      # user sees in the drop-down.
      choices = c("No" = "0", "Yes" = "1")
    ),

    radioButtons(
      inputId  = "alcohol_mode",
      label    = "Alcohol intake — input unit",
      choices  = c("Standard drinks / week" = "drinks",
                   "Grams / week"           = "grams"),
      selected = "drinks",
      inline   = TRUE
    ),
    conditionalPanel(
      condition = "input.alcohol_mode == 'drinks'",
      numericInput(
        inputId = "alcohol_drinks",
        label   = "Alcohol (standard drinks / week)",
        value   = 5, min = 0, max = 250, step = 1
      )
    ),
    conditionalPanel(
      condition = "input.alcohol_mode == 'grams'",
      numericInput(
        inputId = "alcohol_grams",
        label   = "Alcohol (grams / week)",
        value   = 70, min = 0, max = 1000, step = 1
      )
    ),
    helpText(
      "1 drink (one 12oz/355mL beer, 5oz/150mL wine, or 1.5oz/45mL spirits) ≈ 14g of alcohol.",
      tags$br(),
      "Values above 500g/week are capped per training data."
    ),

    selectInput(
      inputId = "smoking",
      label   = "Smoking status",
      choices = c(
        "Non-current smoker (never / ex)" = "0",
        "Current smoker"                  = "1"
      )
    ),

    hr(),

    # actionButton() does nothing by itself — the server listens for clicks
    # via observeEvent(input$predict_btn, { … })
    actionButton(
      inputId = "predict_btn",
      label   = "Calculate Risk",
      class   = "btn-primary w-100"
    )
  ),

  # ── Right main panel: output placeholder ──────────────────────────────────
  # uiOutput() is an empty container that server fills in dynamically.
  # Before the button is clicked this area renders nothing.
  uiOutput("results_ui"),

  div(
    class = "mt-4 pt-3 border-top text-muted",
    style = "font-size: 0.85em;",
    tags$em(
      "For use in at-risk SLD patients in whom advanced fibrosis has been excluded (FIB-4 <2.67)",
      "and after exclusion of competing liver disease. Intended to individualize fibrosis surveillance",
      "intervals within the existing guideline framework. Not validated for use outside this population."
    )
  )
)


# ── 6.  Server ────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # observeEvent() runs its code block whenever the named event fires.
  # Here: run the whole pipeline whenever the button is clicked.
  observeEvent(input$predict_btn, {

    # ── 6a. Input validation ──────────────────────────────────────────────────
    # Collect any problems into a character vector; if any exist, show them
    # inside the results panel and abort the rest of the pipeline with return().
    errors <- character(0)

    if (is.na(input$age) || input$age < 18 || input$age > 100)
      errors <- c(errors, "Age must be between 18 and 100 years.")

    if (is.na(input$bmi) || input$bmi < 10 || input$bmi > 70)
      errors <- c(errors, "BMI must be between 10 and 70 kg/m².")

    if (input$alcohol_mode == "drinks") {
      if (is.na(input$alcohol_drinks) || input$alcohol_drinks < 0)
        errors <- c(errors, "Alcohol intake cannot be negative.")
    } else {
      if (is.na(input$alcohol_grams) || input$alcohol_grams < 0)
        errors <- c(errors, "Alcohol intake cannot be negative.")
    }

    if (length(errors) > 0) {
      output$results_ui <- renderUI({
        div(
          class = "alert alert-danger mt-3",
          strong("Please fix the following:"),
          tags$ul(lapply(errors, tags$li))
        )
      })
      return()   # stop here; do not run the model
    }

    # ── 6b. Build new-patient data frame ─────────────────────────────────────
    # The Fine-Gray model was fitted with these exact codings, recorded in
    # model$xlevels and model$input_coding:
    #   sex:            "1" = male, "2" = female
    #   has_t2dm:       "0" = no,   "1" = yes
    #   smoking_binary: "0" = never/previous, "1" = current
    #   alcohol_grams_week: numeric g/week
    # model$predict_fn coerces these character codes to the trained factor
    # levels itself.  Continuous inputs are clamped and BMI floored by
    # prepare_inputs() just below (section 3b).
    alcohol_g <- if (input$alcohol_mode == "drinks") {
      as.numeric(input$alcohol_drinks) * 14
    } else {
      as.numeric(input$alcohol_grams)
    }

    # Alcohol-assessment banner: evaluated on the raw reported intake (not the
    # clamped value predict_fn uses), so it reflects what the patient actually
    # reported. This only ever adds a banner; it never changes `patient`,
    # `risk`, or the surveillance recommendation computed below.
    alcohol_assessment_flag <- alcohol_g > ALCOHOL_ASSESSMENT_CUT[[input$sex]]

    patient <- data.frame(
      age                = as.numeric(input$age),
      sex                = if (input$sex == "Male") "1" else "2",
      bmi                = as.numeric(input$bmi),
      has_t2dm           = as.character(input$t2dm),
      alcohol_grams_week = alcohol_g,
      smoking_binary     = as.character(input$smoking),
      stringsAsFactors   = FALSE
    )

    # Apply the input policy once; everything downstream uses `patient`
    # (clamped, BMI-floored).  `prep` keeps the raw values and what was changed.
    prep    <- prepare_inputs(model, patient)
    patient <- prep$patient

    # ── 6c. Predict cumulative incidence at 5 and 10 years ────────────────────
    # predict_fn returns a 1-row matrix with columns risk_5y and risk_10y, on
    # the probability scale.  H0(t) is tabulated on a 0.05-year grid on which
    # 5 and 10 are exact knots, so no interpolation error enters here.
    risk <- model$predict_fn(model, patient, times = c(5, 10))
    risk_5  <- 100 * risk[1, "risk_5y"]
    risk_10 <- 100 * risk[1, "risk_10y"]
    ratio_10 <- round(risk_10 / RISK_THRESHOLD_10Y, 1)

    # 95% CI (coefficient uncertainty only — see risk_ci() in section 4b)
    ci    <- risk_ci(model, patient, times = c(5, 10))
    ci_5  <- 100 * ci$lower[["risk_5y"]];  ci_5u  <- 100 * ci$upper[["risk_5y"]]
    ci_10 <- 100 * ci$lower[["risk_10y"]]; ci_10u <- 100 * ci$upper[["risk_10y"]]

    # ── 6d. Attribute the risk across the four modifiable predictors ─────────
    contrib <- attribute_risk(model, prep, time = 10)
    contrib <- contrib[order(-abs(contrib$delta_pp)), ]

    # ── 6e. Clamping notice ──────────────────────────────────────────────────
    # Inputs outside the training bounds are clamped (never extrapolated), so
    # say which values the estimate is actually based on.
    notes <- character(0)
    if (prep$clamped[["age"]])
      notes <- c(notes, sprintf("Age %g is outside the training range (%g–%g years); the estimate uses %g and is not extrapolated.",
                                prep$raw$age, prep$bounds$age[1], prep$bounds$age[2],
                                patient$age))
    if (prep$clamped[["bmi"]])
      notes <- c(notes, sprintf("BMI %g is above the training range; the estimate uses %.1f kg/m² and is not extrapolated.",
                                prep$raw$bmi, patient$bmi))
    if (prep$clamped[["alcohol_grams_week"]])
      notes <- c(notes, sprintf("Alcohol intake %g g/week was capped at %g g/week before prediction.",
                                prep$raw$alcohol_grams_week, patient$alcohol_grams_week))

    # ── 6f. Surveillance recommendation ──────────────────────────────────────
    high_risk     <- risk_10 >= RISK_THRESHOLD_10Y
    surveillance  <- if (high_risk) {
      "Repeat FIB-4 assessment every 1–2 years (intensified surveillance)"
    } else {
      "Repeat FIB-4 assessment every 3 years (standard surveillance)"
    }

    # ── 6g. Build the output UI ──────────────────────────────────────────────
    # renderUI() lets server construct arbitrary HTML/Shiny widgets at
    # run-time and push them into the uiOutput("results_ui") placeholder.
    # layout_columns() / card() / card_header() / card_body() come from bslib.
    output$results_ui <- renderUI({
      tagList(

        hr(),

        # ── Alcohol-assessment banner ─────────────────────────────────────────
        # A flag, not a surveillance override: the recommendation below is
        # always the model's own, computed from the continuous risk. See
        # ALCOHOL_ASSESSMENT_CUT above and src/07_sex_alcohol_interaction.R.
        if (isTRUE(alcohol_assessment_flag)) div(
          class = "alert alert-danger",
          strong("⚠ Alcohol assessment recommended. "),
          sprintf(
            "Reported intake (%s g/week) exceeds the level associated with alcohol-related harm for %s (>%s g/week).",
            round(alcohol_g), tolower(input$sex), ALCOHOL_ASSESSMENT_CUT[[input$sex]]
          ),
          " This flag does not change the calculated risk or the surveillance recommendation below; ",
          "it is a separate prompt to assess alcohol use clinically.",
          if (!high_risk) tagList(
            tags$br(),
            tags$em(
              "Despite this flag, this patient's calculated 10-year risk is below the surveillance ",
              "threshold, so the recommendation below remains standard surveillance, based on the model."
            )
          )
        ),

        if (length(notes) > 0) div(
          class = "alert alert-warning py-2",
          style = "font-size: 0.9em;",
          tags$ul(class = "mb-0", lapply(notes, tags$li))
        ),

        # ── Card 1: predicted risk (with CI) and surveillance recommendation ─
        card(
          card_header(
            class = "fw-semibold",
            "Predicted Risk of Major Adverse Liver Outcomes"
          ),
          card_body(
            tags$table(
              class = "table table-sm table-borderless mb-0",
              tags$tbody(
                tags$tr(
                  tags$th("5-year risk:"),
                  tags$td(strong(fmt_pct(risk_5)),
                          tags$span(class = "text-muted",
                                    sprintf(" (95%% CI %s\u2013%s)",
                                            fmt_pct(ci_5), fmt_pct(ci_5u))))
                ),
                tags$tr(
                  tags$th("10-year risk:"),
                  tags$td(strong(fmt_pct(risk_10)),
                          tags$span(class = "text-muted",
                                    sprintf(" (95%% CI %s\u2013%s)",
                                            fmt_pct(ci_10), fmt_pct(ci_10u))))
                )
              )
            ),
            p(
              class = "mt-2 mb-0 text-muted",
              style = "font-size: 0.9em;",
              sprintf("This patient's risk is %sx the surveillance threshold (%s at 10 years).",
                      ratio_10, fmt_pct(RISK_THRESHOLD_10Y))
            ),
            p(
              class = "mt-1 mb-0 text-muted",
              style = "font-size: 0.8em;",
              "95% CI reflects sampling uncertainty in the model coefficients; ",
              "uncertainty in the baseline hazard is not included."
            ),
            hr(class = "my-2"),
            p(class = "mb-0", strong("Recommendation: "), surveillance)
          )
        ),

        # ── Card 2: modifiable-factor contribution breakdown ─────────────────
        card(
          class = "mt-3",
          card_header(class = "fw-semibold", "Modifiable Risk Factor Contributions"),
          card_body(
            p(
              class = "text-muted mb-2",
              style = "font-size: 0.9em;",
              "Each modifiable factor is compared with the same patient at an optimal profile (BMI 25, no diabetes, no alcohol, non-smoking); age and sex are held at the patient's own values."
            ),
            tags$table(
              class = "table table-sm align-middle mb-2",
              tags$thead(tags$tr(
                tags$th("Factor"),
                tags$th("Patient"),
                tags$th("Optimal"),
                tags$th(class = "text-end", "Risk multiple"),
                tags$th(class = "text-end", "Δ 10-year risk")
              )),
              tags$tbody(lapply(seq_len(nrow(contrib)), function(i) {
                row <- contrib[i, ]
                tags$tr(
                  tags$td(row$label),
                  tags$td(row$value),
                  tags$td(class = "text-muted", row$opt_value),
                  tags$td(class = "text-end", sprintf("×%.2f", row$multiple)),
                  tags$td(
                    class = paste("text-end",
                                  if (row$delta_pp > 0) "text-danger"
                                  else if (row$delta_pp < 0) "text-success"
                                  else "text-muted"),
                    paste0(if (row$delta_pp > 0) "+" else "", fmt_pct(row$delta_pp))
                  )
                )
              }))
            ),
            p(
              class = "text-muted mb-0",
              style = "font-size: 0.85em;",
              tags$strong("Reading this table. "),
              "Risk multiples are exact and multiply together: their product is ",
              sprintf("this patient's total risk multiple versus the same patient at the optimal profile (×%.2f). ",
                      prod(contrib$multiple)),
              "Δ 10-year risk is the absolute risk attributable to that factor, ",
              "holding the other factors fixed; because absolute risk is a nonlinear ",
              "function of the linear predictor, the Δ column does not sum to the total. ",
              "These are model counterfactuals, not estimates of treatment effect.",
              if (any(grepl("†", contrib$value, fixed = TRUE)))
                tagList(tags$br(),
                        "† Value clamped to the training bounds before prediction; ",
                        "the risk estimate uses the clamped value shown here."),
              if (any(grepl("\u2021", contrib$value, fixed = TRUE)))
                tagList(tags$br(),
                        "\u2021 BMI at or below 25 is floored at 25 in the model, so this patient ",
                        "is at target (\u00d71.00, \u0394 0%).")
            )
          )
        ),

        # ── Model provenance ─────────────────────────────────────────────────
        accordion(
          class = "mt-3",
          open = FALSE,
          accordion_panel(
            "About this model",
            tags$ul(
              class = "mb-0",
              style = "font-size: 0.9em;",
              # Drop the internal "BFA" project prefix and re-capitalise the first word.
              tags$li(sub("^(.)", "\\U\\1", sub("^BFA\\s+", "", model$model_name), perl = TRUE)),
              tags$li(sprintf("Fine-Gray subdistribution hazard model; %s", model$equation)),
              tags$li(sprintf(
                "Fitted on %s participants (%s MALO events, %s non-liver deaths; median follow-up %.1f years).",
                format(model$cohort_summary$n, big.mark = ","),
                format(model$cohort_summary$n_malo, big.mark = ","),
                format(model$cohort_summary$n_nonliver_death, big.mark = ","),
                model$cohort_summary$median_followup_years)),
              tags$li(sprintf(
                "Cross-validated performance: C-index %.3f; time-dependent AUC %.3f at 5 years, %.3f at 10 years.",
                model$cv_performance$mean_cindex,
                model$cv_performance$mean_auc_t5,
                model$cv_performance$mean_auc_t10)),
              tags$li("BMI is floored at 25 to avoid treating low BMI as protective."),
              tags$li(sprintf(
                "An alcohol-assessment banner is shown above %s g/week (men) or %s g/week (women); it flags a patient for clinical alcohol assessment but does not change the calculated risk or recommendation.",
                ALCOHOL_ASSESSMENT_CUT[["Male"]], ALCOHOL_ASSESSMENT_CUT[["Female"]])),
              tags$li(sprintf("Model exported %s under R %s.",
                              model$date_created, model$r_version))
            )
          )
        )
      )   # end tagList
    })    # end renderUI
  })      # end observeEvent
}


# ── 7.  Launch ────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
