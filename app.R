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
#   $risk_cut_10y  num         median predicted 10-year risk in the cohort
#   $cohort_summary, $cv_performance   aggregate metadata shown in the footer
#   $predict_fn    function    CIF(t|x) = 1 - exp(-H0(t) * exp(x %*% beta))
#
# $predict_fn does its own preprocessing (BMI winsorising, alcohol capping,
# factor coercion), so the app passes raw user input straight through and never
# reimplements the training transforms.
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

# Reference ("counterfactual target") profile.  Every factor's contribution is
# measured against this patient.  Values: age = midpoint of the training
# 1st–99th percentile range; sex = the model's baseline level (male); and for
# the four modifiable factors, the clinical target a clinician would aim for.
REFERENCE <- list(
  age                = 55,
  sex                = "1",   # male — the reference level of the sex term
  bmi                = 25,
  has_t2dm           = "0",
  alcohol_grams_week = 0,
  smoking_binary     = "0"
)

REFERENCE_TEXT <- "a 55-year-old man with BMI 25 kg/m², no type 2 diabetes, no alcohol intake and no current smoking"

# Factors a clinician can act on.  Age and sex are excluded: they contribute to
# the risk estimate but cannot be targets of advice.
MODIFIABLE <- c("bmi", "has_t2dm", "alcohol_grams_week", "smoking_binary")

VAR_LABELS <- c(
  age                = "Age",
  sex                = "Sex",
  bmi                = "BMI",
  has_t2dm           = "Type 2 diabetes",
  alcohol_grams_week = "Alcohol intake",
  smoking_binary     = "Smoking status"
)

# How to describe moving each factor to its reference value, in a sentence of
# the form "<phrase> is projected to lower the 10-year risk from A% to B%".
ACTION_PHRASE <- c(
  bmi                = "Reducing BMI to 25 kg/m²",
  has_t2dm           = "Absence of type 2 diabetes",
  alcohol_grams_week = "Abstaining from alcohol",
  smoking_binary     = "Stopping smoking"
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
# Setting predictor j to a reference value ref_j and leaving the rest alone
# multiplies the subdistribution hazard by exactly exp(beta_j * (x_j - ref_j)),
# independently of the other predictors and of t.  So for each factor we
# compute two numbers, both exact and both from one call to $predict_fn:
#
#   * risk multiple  — exp(beta_j * (x_j - ref_j)): the factor by which this
#     patient's value of factor j multiplies their risk relative to the
#     reference patient.  These multiply exactly: the product across all six
#     factors is the patient's total risk multiple versus the reference.
#     Recovered from the predictions as log(1-r_patient) / log(1-r_j), because
#     -log(1 - CIF) = H0(t) * exp(lp) and H0(t) cancels.
#
#   * delta  — r_patient - r_j: how many percentage points of this patient's
#     absolute 10-year risk are attributable to factor j sitting where it does
#     rather than at its reference value.  For the four modifiable factors the
#     reference IS the clinical target, so this is the risk reduction that
#     reaching that target would achieve, holding everything else fixed.
#
# The deltas are not additive (absolute risk is a nonlinear function of the
# linear predictor), which is why the table below reports them alongside the
# multiples rather than as a decomposition that sums to the total.  Ranking by
# delta is what identifies the principal driver: it answers "which single
# factor is contributing the most absolute risk for this patient", which is the
# question the SHAP panel was there to answer.
#
# Everything routes through model$predict_fn, so the app inherits the training
# preprocessing (BMI winsorising, alcohol capping) for the counterfactual rows
# as well as the patient row, and there is exactly one implementation of the
# risk equation in the deployment.
# Reproduce the clamping that model$predict_fn applies internally, so the UI can
# display (and flag) the values the risk estimate is actually based on.
effective_values <- function(model, patient) {
  vals <- lapply(model$predictors, function(v) patient[[v]][1])
  names(vals) <- model$predictors
  adj <- setNames(as.list(rep(FALSE, length(vals))), names(vals))

  bw <- model$preprocessing$bmi_winsor
  clamped <- max(min(as.numeric(vals$bmi), bw[2]), bw[1])
  adj$bmi  <- !isTRUE(all.equal(clamped, as.numeric(vals$bmi)))
  vals$bmi <- clamped

  capped <- min(as.numeric(vals$alcohol_grams_week), model$preprocessing$alcohol_cap)
  adj$alcohol_grams_week  <- !isTRUE(all.equal(capped, as.numeric(vals$alcohol_grams_week)))
  vals$alcohol_grams_week <- capped

  list(values = vals, adjusted = adj)
}

attribute_risk <- function(model, patient, time = 10) {
  preds <- model$predictors

  # Row 1 is the patient; rows 2..7 are the patient with a single predictor
  # moved to its reference value.  One predict_fn call covers all of them.
  nd <- patient[rep(1L, 1L + length(preds)), , drop = FALSE]
  for (i in seq_along(preds)) nd[[preds[i]]][i + 1L] <- REFERENCE[[preds[i]]]

  r     <- model$predict_fn(model, nd, times = time)[, 1]
  r_pat <- r[1]
  r_ref <- setNames(r[-1], preds)

  # The table must show the value the model actually used, not the raw input:
  # predict_fn winsorises BMI and caps alcohol, so a raw 900 g/week sits next to
  # a multiple computed at 500 g/week unless we substitute the effective value.
  eff <- effective_values(model, patient)

  data.frame(
    var       = preds,
    label     = unname(VAR_LABELS[preds]),
    value     = vapply(preds, function(v) {
                  txt <- fmt_value(v, eff$values[[v]])
                  if (isTRUE(eff$adjusted[[v]])) paste0(txt, "\u2020") else txt
                }, character(1)),
    ref_value = vapply(preds, function(v) fmt_value(v, REFERENCE[[v]]), character(1)),
    # exp(beta_j * (x_j - ref_j)); > 1 raises risk, < 1 lowers it
    multiple  = log1p(-r_pat) / log1p(-r_ref),
    risk_if_ref = 100 * r_ref,
    delta_pp    = 100 * (r_pat - r_ref),
    modifiable  = preds %in% MODIFIABLE,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


# ── 5.  UI ────────────────────────────────────────────────────────────────────
# page_sidebar() gives a fixed left sidebar and a scrollable main area.
# All widget functions follow the pattern:
#   widgetType(inputId, label, ...)
# The inputId is the name you use in server to read the widget's value
# as input$<inputId>.

ui <- page_sidebar(
  title = "MALO Risk Calculator - SLD At-Risk, Low Fibrosis Burden",
  theme = bs_theme(bootswatch = "flatly"),

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
      "For use in at-risk SLD patients with low fibrosis burden (FIB-4 < 2.67) after exclusion",
      "of competing liver disease. This tool is intended to individualize fibrosis surveillance",
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
    # levels itself, winsorises BMI to the training bounds and caps alcohol, so
    # the raw values go through unmodified here.
    alcohol_g <- if (input$alcohol_mode == "drinks") {
      as.numeric(input$alcohol_drinks) * 14
    } else {
      as.numeric(input$alcohol_grams)
    }

    patient <- data.frame(
      age                = as.numeric(input$age),
      sex                = if (input$sex == "Male") "1" else "2",
      bmi                = as.numeric(input$bmi),
      has_t2dm           = as.character(input$t2dm),
      alcohol_grams_week = alcohol_g,
      smoking_binary     = as.character(input$smoking),
      stringsAsFactors   = FALSE
    )

    # ── 6c. Predict cumulative incidence at 5 and 10 years ────────────────────
    # predict_fn returns a 1-row matrix with columns risk_5y and risk_10y, on
    # the probability scale.  H0(t) is tabulated on a 0.05-year grid on which
    # 5 and 10 are exact knots, so no interpolation error enters here.
    risk <- model$predict_fn(model, patient, times = c(5, 10))
    risk_5  <- 100 * risk[1, "risk_5y"]
    risk_10 <- 100 * risk[1, "risk_10y"]
    ratio_10 <- round(risk_10 / RISK_THRESHOLD_10Y, 1)

    # ── 6d. Attribute the risk across the six predictors ─────────────────────
    contrib <- attribute_risk(model, patient, time = 10)
    contrib <- contrib[order(-abs(contrib$delta_pp)), ]

    # Principal modifiable driver: the modifiable factor carrying the most
    # absolute risk relative to its clinical target.  A patient already at or
    # better than target on all four has no positive delta — say so rather than
    # naming whichever factor happens to be least negative.
    mod <- contrib[contrib$modifiable, ]
    top <- if (any(mod$delta_pp > 0)) mod[which.max(mod$delta_pp), ] else NULL

    # ── 6e. Out-of-range notice ──────────────────────────────────────────────
    # model$input_range holds the training 1st–99th percentiles.  Inputs outside
    # them are extrapolation (or, for BMI and alcohol, are actively winsorised
    # by predict_fn), so the estimate deserves a caveat.
    rng   <- model$input_range
    notes <- character(0)
    if (patient$age < rng$age[1] || patient$age > rng$age[2])
      notes <- c(notes, sprintf("Age %g is outside the training range (%g–%g years); the estimate is an extrapolation.",
                                patient$age, rng$age[1], rng$age[2]))
    if (patient$bmi < model$preprocessing$bmi_winsor[1] ||
        patient$bmi > model$preprocessing$bmi_winsor[2])
      notes <- c(notes, sprintf("BMI %g was winsorised to the training bounds (%.1f–%.1f kg/m²) before prediction.",
                                patient$bmi, model$preprocessing$bmi_winsor[1],
                                model$preprocessing$bmi_winsor[2]))
    if (patient$alcohol_grams_week > model$preprocessing$alcohol_cap)
      notes <- c(notes, sprintf("Alcohol intake %g g/week was capped at %g g/week before prediction.",
                                patient$alcohol_grams_week, model$preprocessing$alcohol_cap))

    # ── 6f. Risk category and recommendation ─────────────────────────────────
    high_risk     <- risk_10 >= RISK_THRESHOLD_10Y
    risk_category <- if (high_risk) "High Risk" else "Low Risk"
    threshold_txt <- paste0(if (high_risk) "≥ " else "< ", fmt_pct(RISK_THRESHOLD_10Y))
    card_colour   <- if (high_risk) "danger" else "success"
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

        if (length(notes) > 0) div(
          class = "alert alert-warning py-2",
          style = "font-size: 0.9em;",
          tags$ul(class = "mb-0", lapply(notes, tags$li))
        ),

        layout_columns(
          col_widths = c(5, 7),

          # ── Card 1: predicted risk numbers ──────────────────────────────
          card(
            height = "100%",
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
                    tags$td(strong(fmt_pct(risk_5)))
                  ),
                  tags$tr(
                    tags$th("10-year risk:"),
                    tags$td(strong(fmt_pct(risk_10)))
                  )
                )
              ),
              p(
                class = "mt-2 mb-0 text-muted",
                style = "font-size: 0.9em;",
                sprintf("This patient's risk is %sx the surveillance threshold (%s at 10 years).",
                        ratio_10, fmt_pct(RISK_THRESHOLD_10Y))
              )
            )
          ),

          # ── Card 2: risk category + principal modifiable driver ──────────
          card(
            height = "100%",
            card_header(
              class = paste0("bg-", card_colour, " text-white fw-semibold"),
              paste0("Risk Category: ", risk_category,
                     "  (10-year risk ", threshold_txt, ")")
            ),
            card_body(
              if (is.null(top)) {
                p(
                  strong("Principal Modifiable Risk Factor: "), "none — ",
                  "this patient is already at or better than target on BMI, ",
                  "diabetes status, alcohol intake and smoking. The remaining ",
                  "risk is driven by age and sex."
                )
              } else {
                tagList(
                  p(
                    strong("Principal Modifiable Risk Factor: "), top$label,
                    sprintf(" (%s)", top$value)
                  ),
                  p(
                    class = "mb-0",
                    sprintf(
                      "It multiplies this patient's 10-year risk by %.2f. %s is projected to lower that risk from %s to %s (%s percentage points).",
                      top$multiple, ACTION_PHRASE[[top$var]],
                      fmt_pct(risk_10), fmt_pct(top$risk_if_ref),
                      fmt_pct(-top$delta_pp, suffix = "")
                    )
                  )
                )
              },
              hr(class = "my-2"),
              p(class = "mb-0", strong("Recommendation: "), surveillance)
            )
          )
        ), # end layout_columns

        # ── Card 3: full contribution breakdown ──────────────────────────────
        card(
          class = "mt-3",
          card_header(class = "fw-semibold", "Risk Factor Contributions"),
          card_body(
            p(
              class = "text-muted mb-2",
              style = "font-size: 0.9em;",
              "Each factor is compared with a reference patient: ", REFERENCE_TEXT,
              sprintf(" (10-year risk %s).",
                      fmt_pct(100 * model$predict_fn(
                        model, as.data.frame(REFERENCE, stringsAsFactors = FALSE),
                        times = 10)[, 1]))
            ),
            tags$table(
              class = "table table-sm align-middle mb-2",
              tags$thead(tags$tr(
                tags$th("Factor"),
                tags$th("Patient"),
                tags$th("Reference"),
                tags$th(class = "text-end", "Risk multiple"),
                tags$th(class = "text-end", "Δ 10-year risk")
              )),
              tags$tbody(lapply(seq_len(nrow(contrib)), function(i) {
                row <- contrib[i, ]
                tags$tr(
                  tags$td(row$label,
                          if (!row$modifiable)
                            tags$span(class = "text-muted",
                                      style = "font-size: 0.85em;",
                                      " (not modifiable)")),
                  tags$td(row$value),
                  tags$td(class = "text-muted", row$ref_value),
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
              sprintf("this patient's total risk multiple versus the reference patient (×%.2f). ",
                      prod(contrib$multiple)),
              "Δ 10-year risk is the absolute risk attributable to that factor, ",
              "holding the other five fixed; because absolute risk is a nonlinear ",
              "function of the linear predictor, the Δ column does not sum to the total. ",
              "These are model counterfactuals, not estimates of treatment effect.",
              if (any(grepl("†", contrib$value, fixed = TRUE)))
                tagList(tags$br(),
                        "† Value clamped to the training bounds before prediction; ",
                        "the risk estimate uses the clamped value shown here.")
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
              tags$li(model$model_name),
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
