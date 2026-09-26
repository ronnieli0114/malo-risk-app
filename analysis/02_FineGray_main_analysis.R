###############################################################################
# Main analysis using Fine-Gray competing-risks regression
# Ronnie Li
#
# Outcome: MALO (cause 1); competing event: non-liver death (cause 2).
# All predictors enter the Fine-Gray models as linear terms. Whether allowing
# nonlinearity improves model fit and prediction is assessed with the random
# survival forest in 03_RSF_sensitivity_analysis.R.
#
# All cross-validated results use the stratified 5-fold split (seed 1234).
# The out-of-fold (OOF) predictions and fold assignment are saved at the end
# and re-used by 03_RSF_sensitivity_analysis.R, so that the RSF is evaluated
# on exactly the same folds.
#
# 1. 5-fold CV: Fine-Gray (clinical, precision) vs CLivD vs FIB-4
#    time-dependent AUC, C-index, Brier score
# 2. Risk stratification by median predicted risk: above median (aggressive
#    surveillance) vs below median (standard surveillance) -- MALO and competing
#    non-liver death (numbers, Aalen-Johansen cumulative incidence, Gray's
#    test), KM curves
# 3. Decision curve analysis vs the standard-of-care surveillance rules
#    (type 2 diabetes, >= 2 CMRFs, and the guideline "T2DM or >= 2 CMRFs")
# 4. Calibration plots (decile-based, AJ-CIF)
# 5. Full-data Fine-Gray fits: SHR forest plots/tables and model equation
# 6. Export the deployable clinical Fine-Gray model (coefficients + baseline
#    cumulative subdistribution hazard only; no individual-level data)
# 7. Save all tables and the OOF predictions
#
# Outputs: tables  -> results/02_FineGray_results.xlsx (sheets 01_, 02_, ...)
#          figures -> results/plots/02.01_*.png, 02.02_*.png, ...
#          model   -> results/FG_clinical_model_deploy.rds (for the Shiny app)
###############################################################################

library(riskRegression)
library(timeROC)
library(survival)
library(pec)       # for ipcw()
library(cmprsk)
library(prodlim)
library(tidyverse)
library(arrow)
library(writexl)
library(patchwork)
library(dcurves)
library(ggsci)

# ── Directories ───────────────────────────────────────────────────────────────
project_dir <- "/mnt/d/Projects/BFA"
source(file.path(project_dir, "src/clivd_scores.R"))
data_dir   <- file.path(project_dir, "data")
result_dir <- file.path(project_dir, "results")
plot_dir   <- file.path(result_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
oof_path   <- file.path(data_dir, "BFA_FineGray_OOF_predictions.feather")
results_xlsx <- file.path(result_dir, "02_FineGray_results.xlsx")
deploy_rds   <- file.path(result_dir, "FG_clinical_model_deploy.rds")

tables <- list()   # every output table; written to one workbook at the end

# ── Helpers ───────────────────────────────────────────────────────────────────
# Winsorizing function caps the variable between the 1st and 99th percentiles
winsorise <- function(x, probs = c(0.01, 0.99)) {
  q <- quantile(x, probs, na.rm = TRUE)
  pmax(pmin(x, q[2]), q[1])
}

make_fg_formula <- function(preds) {
  as.formula(paste("Hist(time_to_event, status) ~", paste(preds, collapse = " + ")))
}

fit_fg <- function(preds, data) FGR(make_fg_formula(preds), data = data, cause = 1)

# Absolute risk of MALO (CIF) at time t
fg_risk <- function(fg_fit, newdata, t) {
  as.numeric(predictRisk(fg_fit, newdata = newdata, times = t, cause = 1))
}

# Recalibrate a single score (CLivD LP, FIB-4 or one predictor) to absolute
# MALO risk with a one-covariate Fine-Gray model fitted in the training data
fit_fg_score <- function(train, score) {
  FGR(Hist(time_to_event, status) ~ score, cause = 1,
      data = data.frame(time_to_event = train$time_to_event,
                        status = train$status, score = score))
}

# Time-dependent AUC for cause 1 (controls: event-free at t)
td_auc <- function(df, marker, t) {
  roc <- timeROC(T = df$time_to_event, delta = df$status, marker = marker,
                 cause = 1, times = t, iid = FALSE)
  unname(roc$AUC_1[match(t, roc$times)])
}

# Harrell's C-index for MALO (competing events treated as censored)
c_index <- function(df, marker) {
  concordance(Surv(df$time_to_event, df$status == 1) ~ marker, reverse = TRUE)$concordance
}

# IPCW Brier score for competing risks (cause 1)
compute_brier <- function(pred_risk, t, df) {
  iw <- pec::ipcw(Surv(time_to_event, status != 0) ~ 1, data = df,
                  method = "marginal", times = t, subjectTimes = df$time_to_event)
  # pec::ipcw() is an internal helper of pec() and assumes the data are sorted
  # by time (its own help-page example sorts first). It returns
  # IPCW.subjectTimes in THAT sorted order, not in the row order of `df`, so the
  # values must be mapped back onto the original rows. Using them as returned
  # pairs G(T_i-) with the wrong subject and can hand a participant who had an
  # event well before t a weight of 100+ instead of ~1.
  G_subject <- numeric(nrow(df))
  G_subject[order(df$time_to_event)] <- as.numeric(iw$IPCW.subjectTimes)
  y <- as.numeric(df$time_to_event <= t & df$status == 1)
  w <- ifelse(df$time_to_event <= t & df$status != 0, 1 / G_subject,            # any event by t
       ifelse(df$time_to_event > t, 1 / iw$IPCW.times, 0))                     # censored before t: 0
  mean(w * (y - pred_risk)^2)
}

# ── Data preparation ──────────────────────────────────────────────────────────
adata <- read_feather(file.path(data_dir, "BFA_principal_data.feather"))

model_data <- adata %>%
  filter(included_in_cohort == TRUE) %>%
  mutate(
    bmi                = winsorise(bmi),
    hip_circ_raw       = hip_circ,
    waist_hip_ratio    = waist_circ / hip_circ,          # unwinsorised WHR for CLivD
    waist_circ         = winsorise(waist_circ),          # winsorised WC for precision model
    triglycerides      = winsorise(triglycerides),
    hdl                = winsorise(hdl),
    trig_hdl_ratio     = triglycerides / hdl,
    alcohol_grams_week = pmin(alcohol_grams_week, 500),  # CLivD-aligned cap
    hba1c              = winsorise(pmin(hba1c, 70)),
    smoking_binary   = as.factor(case_when(
      smoking %in% c("Never", "Previous") ~ 0L,
      smoking == "Current"                ~ 1L
    )),
    has_t2dm         = as.factor(as.integer(has_t2dm)),
    has_hypertension = as.factor(as.integer(has_hypertension)),
    sex              = as.factor(case_when(sex == "Male" ~ 1L, sex == "Female" ~ 2L)),
    smoking          = factor(smoking, levels = c("Never", "Previous", "Current")),
    status           = case_when(
      event_malo == TRUE            ~ 1L,
      event_non_liver_death == TRUE ~ 2L,
      TRUE                         ~ 0L
    )
  ) %>%
  drop_na(age, sex, bmi, waist_circ, alcohol_grams_week, smoking_binary, 
          PNPLA3_rs738409_G, TM6SF2_rs58542926_T, HSD17B13_rs9992651_A, 
          time_to_event, status)

# Predictor sets
preds_clin <- c("age", "sex", "bmi", "has_t2dm", "alcohol_grams_week", "smoking_binary")
preds_prec <- c(preds_clin, "waist_circ", "PNPLA3_rs738409_G", "TM6SF2_rs58542926_T", "HSD17B13_rs9992651_A")
stopifnot(all(complete.cases(model_data[, preds_clin])))
eval_times <- c(5, 10)

# Stratified 5-fold CV (same seed as the RSF sensitivity analysis)
set.seed(1234)
n_folds   <- 5
malo_idx  <- which(model_data$status == 1)
rest_idx  <- which(model_data$status != 1)
model_data$fold <- NA_integer_
model_data$fold[malo_idx] <- sample(rep(1:n_folds, length.out = length(malo_idx)))
model_data$fold[rest_idx] <- sample(rep(1:n_folds, length.out = length(rest_idx)))

# Multivariable Fine-Gray models fitted in each fold
fg_pred_sets <- list(fg_clin = preds_clin, fg_prec = preds_prec)

# Single scores recalibrated to absolute risk in each fold (comparators + DCA)
score_fns <- list(
  clivd_nl   = compute_clivd_nonlab,
  clivd_l    = compute_clivd,
  fib4_fg    = function(d) d$fib4,
  # Single predictors recalibrated to absolute risk. Section 3 now benchmarks
  # against the standard-of-care decision rules instead, but these are still
  # computed and saved in the OOF file for reference.
  sv_alcohol = function(d) d$alcohol_grams_week,
  sv_bmi     = function(d) d$bmi,
  sv_t2dm    = function(d) as.numeric(as.character(d$has_t2dm))
)

# Models evaluated in Section 1 (label -> column prefix)
perf_models <- c("FG (clinical)"   = "fg_clin",
                 "FG (precision)"  = "fg_prec",
                 "CLivD (non-lab)" = "clivd_nl",
                 "CLivD (lab)"     = "clivd_l",
                 "FIB-4"           = "fib4_fg")
risk_cols <- as.vector(outer(c(names(fg_pred_sets), names(score_fns)), eval_times, paste, sep = "_"))

###############################################################################
# SECTION 1: 5-fold CV — Fine-Gray vs CLivD vs FIB-4
###############################################################################
cat("\n=== Section 1: CV loop ===\n")

fold_metric_list <- list()
oof_list         <- list()

for (fold in 1:n_folds) {
  cat("  Fold", fold, "\n")
  tr <- model_data[model_data$fold != fold, ]
  te <- model_data[model_data$fold == fold, ]
  cat("    Train events:", sum(tr$status == 1), "| Test events:", sum(te$status == 1), "\n")

  # ── Multivariable Fine-Gray models ──
  for (m in names(fg_pred_sets)) {
    cat("    Fitting", m, "...\n")
    fit <- fit_fg(fg_pred_sets[[m]], tr)
    for (t in eval_times) te[[paste0(m, "_", t)]] <- fg_risk(fit, te, t)
  }

  # ── CLivD, FIB-4 and single predictors: one-covariate Fine-Gray recalibration ──
  for (m in names(score_fns)) {
    fit  <- fit_fg_score(tr, score_fns[[m]](tr))
    newd <- data.frame(score = score_fns[[m]](te))
    for (t in eval_times) te[[paste0(m, "_", t)]] <- fg_risk(fit, newd, t)
  }

  # ── Per-fold AUC (marker = risk at t), C-index (risk at 10 y), Brier ──
  fold_metric_list[[fold]] <- bind_rows(lapply(names(perf_models), function(m) {
    k <- perf_models[[m]]
    tibble(
      fold      = fold,
      model     = m,
      auc_t5    = td_auc(te, te[[paste0(k, "_5")]],  5),
      auc_t10   = td_auc(te, te[[paste0(k, "_10")]], 10),
      cindex    = c_index(te, te[[paste0(k, "_10")]]),
      brier_t5  = compute_brier(te[[paste0(k, "_5")]],  5,  te),
      brier_t10 = compute_brier(te[[paste0(k, "_10")]], 10, te)
    )
  }))

  oof_list[[fold]] <- te %>% dplyr::select(eid, fold, time_to_event, status, all_of(risk_cols))
  gc()
}

# Standard-of-care surveillance rules (binary decisions, not risk scores).
# Guideline practice is to screen for advanced fibrosis when a patient has
# type 2 diabetes or >= 2 cardiometabolic risk factors; n_cmrfs is the count of
# the five CMRFs defined in 01_preprocess_clean_ukbb_cohort.R (large waist,
# prediabetes/T2DM, hypertension, hypertriglyceridaemia, low HDL).
soc_flags <- model_data %>%
  transmute(eid,
            soc_t2dm      = as.integer(as.character(has_t2dm)),
            soc_cmrf2     = as.integer(n_cmrfs >= 2),
            soc_guideline = as.integer(soc_t2dm == 1L | n_cmrfs >= 2))

all_oof      <- bind_rows(oof_list) %>% left_join(soc_flags, by = "eid")
stopifnot(!anyNA(all_oof[, c("soc_t2dm", "soc_cmrf2", "soc_guideline")]))
fold_metrics <- bind_rows(fold_metric_list)

# ── Summarize Section 1 ──
perf_summary <- fold_metrics %>%
  group_by(model) %>%
  summarise(across(c(auc_t5, auc_t10, cindex, brier_t5, brier_t10),
                   list(mean = mean, sd = sd), .names = "{.fn}_{.col}"),
            .groups = "drop") %>%
  mutate(model = factor(model, levels = names(perf_models))) %>%
  arrange(model)

tables[["01_Performance"]]          <- perf_summary
tables[["02_Performance_per_fold"]] <- fold_metrics
print(perf_summary)

# AUC plot
auc_models <- names(perf_models)
p_auc <- fold_metrics %>%
  pivot_longer(c(auc_t5, auc_t10), names_to = "time", values_to = "auc") %>%
  group_by(model, time) %>%
  summarise(mean_auc = mean(auc), sd_auc = sd(auc), .groups = "drop") %>%
  mutate(model_label = factor(gsub(" ", "\n", model), levels = gsub(" ", "\n", auc_models)),
         time_label  = factor(if_else(time == "auc_t5", "t = 5 years", "t = 10 years"),
                              levels = c("t = 5 years", "t = 10 years"))) %>%
  ggplot(aes(x = model_label, y = mean_auc)) +
  geom_point(size = 2.5, color = "steelblue") +
  geom_errorbar(aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc),
                color = "steelblue", width = 0.3) +
  geom_hline(yintercept = 0.50, color = "darkorange",
             linewidth = 0.6, linetype = "dashed") +
  facet_wrap(~ time_label) +
  scale_y_continuous(limits = c(0.5, 0.85),
                     labels = scales::number_format(accuracy = 0.01)) +
  labs(x = "Model", y = "Time-dependent AUC",
       title = "Model comparison: Fine-Gray vs. CLivD vs. FIB-4",
       subtitle = "Mean ± SD across 5 CV folds") +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey90"))

ggsave(file.path(plot_dir, "02.01_FG_AUC_comparison.png"),
       p_auc, width = 8, height = 4.5, dpi = 400, bg = "white")

###############################################################################
# SECTION 2: Risk stratification by median predicted risk
# Two surveillance strata rather than tertiles: participants above the median
# predicted 10-year MALO risk would be triaged to aggressive surveillance,
# those below the median to standard surveillance. The cut-point is the median
# of the out-of-fold Fine-Gray (clinical) CIF at t = 10 years, and is exported
# with the deployable model in Section 6 so the calculator can apply it.
# MALO and competing non-liver death: numbers, Aalen-Johansen CIF, Gray's test
###############################################################################
cat("\n=== Section 2: Risk stratification ===\n")

risk_group_levels <- c("Below median", "Above median")
risk_group_labels <- c("Below median (standard surveillance)",
                       "Above median (aggressive surveillance)")

risk_cut_10y <- median(all_oof$fg_clin_10)
cat(sprintf("Median OOF predicted 10-year MALO risk (cut-point): %.4f%%\n",
            100 * risk_cut_10y))

all_oof <- all_oof %>%
  mutate(risk_group = factor(if_else(fg_clin_10 > risk_cut_10y,
                                     "Above median", "Below median"),
                             levels = risk_group_levels))

ci_fit     <- cmprsk::cuminc(ftime = all_oof$time_to_event, fstatus = all_oof$status,
                             group = all_oof$risk_group)
ci_overall <- cmprsk::cuminc(ftime = all_oof$time_to_event, fstatus = all_oof$status)

# ── 3a. Numbers and cumulative incidence (95% CI, log-log) at 5 and 10 years ──
cif_at_times <- function(ci_obj, times = eval_times) {
  tp  <- cmprsk::timepoints(ci_obj, times)
  est <- tp$est
  a   <- qnorm(0.975) * sqrt(tp$var) / abs(est * log(est))
  tibble(curve = rep(rownames(est), times = length(times)),
         time  = rep(times, each = nrow(est)),
         est   = as.vector(est),
         lower = as.vector(est^exp(a)),
         upper = as.vector(est^exp(-a))) %>%
    mutate(group = sub(" [0-9]+$", "", curve),
           cause = sub("^.* ", "", curve))
}

cif_long <- bind_rows(
  cif_at_times(ci_overall) %>% mutate(group = "Overall"),
  cif_at_times(ci_fit)
) %>%
  mutate(event = if_else(cause == "1", "MALO", "Non-liver death"),
         across(c(est, lower, upper), ~ round(100 * .x, 2))) %>%
  dplyr::select(group, event, time, cif_pct = est, lower_pct = lower, upper_pct = upper)

group_counts <- bind_rows(
  all_oof %>% mutate(group = "Overall"),
  all_oof %>% mutate(group = as.character(risk_group))
) %>%
  group_by(group) %>%
  summarise(n                   = n(),
            pred_risk_10y_min   = round(100 * min(fg_clin_10), 3),
            pred_risk_10y_max   = round(100 * max(fg_clin_10), 3),
            n_malo              = sum(status == 1),
            pct_malo            = round(100 * mean(status == 1), 2),
            n_nonliver_death    = sum(status == 2),
            pct_nonliver_death  = round(100 * mean(status == 2), 2),
            n_censored          = sum(status == 0),
            .groups = "drop")

cif_wide <- cif_long %>%
  mutate(value = sprintf("%.2f (%.2f–%.2f)", cif_pct, lower_pct, upper_pct),
         key   = paste0("CIF_", if_else(event == "MALO", "malo", "nonliver_death"),
                        "_", time, "y_pct_95CI")) %>%
  dplyr::select(group, key, value) %>%
  pivot_wider(names_from = key, values_from = value)

competing_summary <- group_counts %>%
  left_join(cif_wide, by = "group") %>%
  mutate(group = factor(group, levels = c("Overall", risk_group_levels))) %>%
  arrange(group)

grays_test <- tibble(
  event     = c("MALO", "Non-liver death"),
  statistic = ci_fit$Tests[c("1", "2"), "stat"],
  df        = ci_fit$Tests[c("1", "2"), "df"],
  p_value   = ci_fit$Tests[c("1", "2"), "pv"]
)

cat("Events and cumulative incidence by surveillance risk group:\n")
print(competing_summary, width = Inf)
print(grays_test)

tables[["03_Competing_events_by_risk_group"]] <- competing_summary
tables[["04_CIF_by_risk_group_long"]]         <- cif_long
tables[["05_Grays_test"]]                     <- grays_test

# ── 3b. Aalen-Johansen CIF curves by risk group (MALO and non-liver death) ──
cif_curves <- function(ci_obj, cause) {
  bind_rows(lapply(names(ci_obj), function(nm) {
    if (!is.list(ci_obj[[nm]]) || !endsWith(nm, paste0(" ", cause))) return(NULL)
    tibble(time  = ci_obj[[nm]]$time,
           cif   = ci_obj[[nm]]$est,
           group = sub(" [0-9]+$", "", nm))
  })) %>%
    mutate(group = factor(group, levels = risk_group_levels))
}

p_label <- function(p, test) {
  paste0(test, " p ", if (p < 0.001) "< 0.001" else paste0("= ", formatC(p, digits = 3, format = "f")))
}

plot_cif <- function(cause, ylab, title) {
  df <- cif_curves(ci_fit, cause)
  ggplot(df, aes(x = time, y = cif, colour = group)) +
    geom_step(linewidth = 0.8) +
    annotate("label", x = 0, y = max(df$cif), hjust = 0, vjust = 1, size = 3,
             label = p_label(ci_fit$Tests[as.character(cause), "pv"], "Gray's")) +
    scale_color_d3(labels = risk_group_labels) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    labs(x = "Time (years)", y = ylab, colour = "Surveillance strategy", title = title,
         subtitle = sprintf("Split at the median Fine-Gray clinical OOF CIF at t = 10 years (%.2f%%)",
                            100 * risk_cut_10y)) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom")
}

p_cif_malo  <- plot_cif(1, "Cumulative incidence (MALO)", "Aalen-Johansen CIF: MALO")
p_cif_death <- plot_cif(2, "Cumulative incidence (non-liver death)",
                        "Aalen-Johansen CIF: competing non-liver death")

# ── 3c. Kaplan-Meier curves ──
km_fit  <- survfit(Surv(time_to_event, status == 1) ~ risk_group, data = all_oof)
lr_test <- survdiff(Surv(time_to_event, status == 1) ~ risk_group, data = all_oof)
lr_p    <- 1 - pchisq(lr_test$chisq, df = length(lr_test$obs) - 1)

km_df <- broom::tidy(km_fit) %>%
  mutate(strata = gsub("risk_group=", "", strata),
         strata = factor(strata, levels = risk_group_levels))

p_km <- ggplot(km_df, aes(x = time, y = estimate, colour = strata)) +
  geom_step(linewidth = 0.8) +
  annotate("label", x = 0, y = min(km_df$estimate), hjust = 0, vjust = 0, size = 3,
           label = p_label(lr_p, "Log-rank")) +
  scale_color_d3(labels = risk_group_labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(x = "Time (years)", y = "Kaplan-Meier survival (MALO)",
       colour = "Surveillance strategy",
       title = "Kaplan-Meier curves by median predicted risk",
       subtitle = "Note: competing events treated as censored") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")

p_risk_strat <- p_cif_malo + p_cif_death + p_km + plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10))

ggsave(file.path(plot_dir, "02.02_FG_risk_stratification.png"),
       p_risk_strat, width = 15, height = 4.5, dpi = 400, bg = "white")

###############################################################################
# SECTION 3: Decision curve analysis
# Note: DCA net benefit is small in absolute terms because MALO is rare in this
# general-population cohort (< ~1% 5-year event rate). "Treat All" dominates
# at very low thresholds because the FP penalty is negligible there. The
# clinically relevant window is ~0.5–3% for t=5 and ~1–5% for t=10 — zoom
# in on that range to see meaningful separation from "Treat None".
#
# Panel B benchmarks the model against the standard-of-care surveillance
# decisions: refer a patient for fibrosis assessment if they have type 2
# diabetes, if they have >= 2 cardiometabolic risk factors, or (the guideline
# rule) if either holds. These are all-or-nothing decision rules rather than
# risk scores, so they enter the DCA as 0/1 indicators and refer exactly the
# same patients at every threshold. Their net benefit still declines with the
# threshold (the false-positive weight pt/(1-pt) grows), but they cannot trade
# sensitivity for specificity as the threshold moves — which is the point of
# the comparison: the model can tighten or loosen its referral rate to match
# the threshold, the standard-of-care rules cannot.
###############################################################################
cat("\n=== Section 3: Decision curve analysis ===\n")

thresholds <- list("5" = seq(0, 0.03, by = 0.001), "10" = seq(0, 0.03, by = 0.001))
dca_tables <- list()

for (t in eval_times) {
  dca_dat <- all_oof %>%
    transmute(
      time_to_event,
      status,
      `Fine-Gray (clinical)`        = .data[[paste0("fg_clin_", t)]],
      `Type 2 diabetes`             = soc_t2dm,
      `Two or more CMRFs`           = soc_cmrf2,
      `T2DM or two or more CMRFs`   = soc_guideline
    )

  thres <- thresholds[[as.character(t)]]

  # dca_A <- dca(Surv(time_to_event, status == 1) ~ `Fine-Gray (clinical)`,
  #              data = dca_dat, time = t, thresholds = thres)
  # p_dca_A <- dca_A %>%
  #   plot(smooth = TRUE) +
  #   ggsci::scale_color_d3() +
  #   labs(title = paste0("DCA at t = ", t, " years: Fine-Gray net benefit"),
  #        x = "Threshold probability", y = "Net benefit") +
  #   theme_bw(base_size = 8) +
  #   theme(legend.position = "bottom")

  dca_B <- dca(Surv(time_to_event, status == 1) ~ `Fine-Gray (clinical)` +
                 `Type 2 diabetes` + `Two or more CMRFs`,
               data = dca_dat, time = t, thresholds = thres)
  p_dca_B <- dca_B %>%
    plot(smooth = TRUE) +
    ggsci::scale_color_d3() +
    labs(title = paste0("DCA at t = ", t, " years: Fine-Gray vs standard of care"),
         x = "Threshold probability", y = "Net benefit") +
    theme_bw(base_size = 8) +
    theme(legend.position = "bottom")

  # p_dca_combined <- wrap_plots(p_dca_A, p_dca_B, ncol = 2) +
  #   plot_annotation(title = paste0("Decision Curve Analysis at t = ", t, " years"))
  p_dca_combined <- p_dca_B
  
  fig_no <- 2 + match(t, eval_times)   # 02.03 (t = 5), 02.04 (t = 10)
  ggsave(file.path(plot_dir, sprintf("02.%02d_FG_DCA_t%d.png", fig_no, t)),
         p_dca_combined, width = 6, height = 4.5, dpi = 400, bg = "white")

  # Panel B contains every curve shown in Panel A
  dca_tables[[as.character(t)]] <- as_tibble(dca_B) %>% mutate(eval_time = t)
}
tables[["06_DCA_net_benefit"]] <- bind_rows(dca_tables)

# How many patients each standard-of-care rule sends to surveillance, and how
# many MALO cases they capture — the referral burden behind the DCA curves
soc_rule_summary <- bind_rows(
  lapply(list(`Type 2 diabetes`            = all_oof$soc_t2dm == 1,
              `Two or more CMRFs`          = all_oof$soc_cmrf2 == 1,
              `T2DM or two or more CMRFs`  = all_oof$soc_guideline == 1,
              `Above median predicted risk` = all_oof$risk_group == "Above median"),
         function(flag) {
           tibble(n_referred        = sum(flag),
                  pct_referred      = round(100 * mean(flag), 2),
                  n_malo_captured   = sum(flag & all_oof$status == 1),
                  pct_malo_captured = round(100 * sum(flag & all_oof$status == 1) /
                                              sum(all_oof$status == 1), 2))
         }),
  .id = "rule")
print(soc_rule_summary, width = Inf)
tables[["07_SOC_rule_referral_burden"]] <- soc_rule_summary

###############################################################################
# SECTION 4: Calibration plots — Fine-Gray clinical and precision
###############################################################################
cat("\n=== Section 4: Calibration plots ===\n")

# Observed (Aalen-Johansen) vs mean predicted CIF across predicted-risk deciles
# within each fold, pooled across folds
compute_aj_cif <- function(oof, pred_prefix, model_label) {
  bind_rows(lapply(eval_times, function(t) {
    oof %>%
      mutate(pred = .data[[paste0(pred_prefix, "_", t)]]) %>%
      group_by(fold) %>%
      mutate(bin = ntile(pred, 10)) %>%
      group_by(bin) %>%
      summarise(
        mean_pred = mean(pred),
        obs_cif   = {
          sub <- pick(everything())
          if (sum(sub$status == 1) == 0) 0 else
            as.numeric(predict(prodlim(Hist(time_to_event, status) ~ 1, data = sub),
                               times = t, cause = 1))
        },
        .groups = "drop"
      ) %>%
      mutate(eval_time = t)
  })) %>%
    mutate(model      = model_label,
           time_label = factor(paste0("t = ", eval_time, " years"),
                               levels = paste0("t = ", eval_times, " years")))
}

cal_all_df <- bind_rows(
  compute_aj_cif(all_oof, "fg_clin", "FG (clinical)"),
  compute_aj_cif(all_oof, "fg_prec", "FG (precision)")
)
tables[["08_Calibration_deciles"]] <- cal_all_df %>% dplyr::select(-time_label)

p_cal <- ggplot(cal_all_df, aes(x = mean_pred, y = obs_cif)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2, colour = "#2166AC") +
  geom_smooth(method = "loess", se = TRUE,
              colour = "#2166AC", fill = "#2166AC", alpha = 0.15, linewidth = 0.7) +
  facet_grid(model ~ time_label, scales = "free") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(x = "Mean predicted CIF",
       y = "Observed CIF (Aalen-Johansen)",
       title = "Calibration plots: Fine-Gray clinical and precision models",
       subtitle = "Observed vs. predicted CIF across predicted risk deciles (pooled across folds)") +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey90"))

ggsave(file.path(plot_dir, "02.05_FG_calibration.png"),
       p_cal, width = 8, height = 6, dpi = 400, bg = "white")

###############################################################################
# SECTION 5: Full-data Fine-Gray fits — SHR forest plots and model equation
###############################################################################
cat("\n=== Section 5: Full-data Fine-Gray fits ===\n")

# Scale data to have more reasonable SHRs
model_data_scaled <- model_data %>%
  mutate(
    age        = age / 3,                         # SHR per 3 years
    bmi        = bmi / 3,                         # SHR per 3 kg/m²
    waist_circ = waist_circ / 3,                  # SHR per 3 cm
    alcohol_grams_week = alcohol_grams_week / 20  # SHR per 20 g/week
  )

fg_clinical_full  <- fit_fg(preds_clin, model_data_scaled)
fg_precision_full <- fit_fg(preds_prec, model_data_scaled)

# Unscaled clinical fit: the one the model equation and the deployable model
# are built from, so that the calculator takes predictors in natural units
fg_clin_full_raw <- fit_fg(preds_clin, model_data)

# Tidy Fine-Gray results
tidy_fg <- function(fg_model, model_label) {
  s  <- summary(fg_model$crrFit)
  cf <- as.data.frame(s$coef)
  colnames(cf) <- c("coef", "SHR", "se_coef", "z", "p")
  cf$p_adj        <- p.adjust(cf$p, method = "BH")
  cf$significance <- factor(ifelse(cf$p_adj < 0.05, "Yes", "No"), levels = c("No", "Yes"))
  cf$term         <- rownames(cf)
  cf$label        <- model_label
  cf$SHR_lo       <- exp(cf$coef - 1.96 * cf$se_coef)
  cf$SHR_hi       <- exp(cf$coef + 1.96 * cf$se_coef)
  cf
}

# Pretty labels
nice_names <- c(
  age                  = "Age (per 3 years)",
  sex2                 = "Sex (Female vs. Male)",
  bmi                  = "BMI (per 3 kg/m²)",
  has_t2dm1            = "Type 2 diabetes",
  alcohol_grams_week   = "Alcohol intake (per 20 g/week)",
  smoking_binary1      = "Current smoker",
  waist_circ           = "Waist circumference (per 3 cm)",
  waist_hip_ratio      = "Waist-hip ratio",
  trig_hdl_ratio       = "Triglyceride:HDL ratio",
  PNPLA3_rs738409_G    = "PNPLA3 (rs738409-G)",
  TM6SF2_rs58542926_T  = "TM6SF2 (rs58542926-T)",
  HSD17B13_rs9992651_A = "HSD17B13 (rs9992651-A)"
)

relabel <- function(df) {
  df$display <- ifelse(df$term %in% names(nice_names), nice_names[df$term], df$term)
  df
}

df_clinical  <- relabel(tidy_fg(fg_clinical_full,  "Clinical"))
df_precision <- relabel(tidy_fg(fg_precision_full, "Precision"))

make_forest_plot <- function(df, title) {
  df <- df %>% arrange(SHR)
  df$display <- factor(df$display, levels = df$display)

  ggplot(df, aes(x = SHR, y = display, color = significance)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_errorbar(aes(xmin = SHR_lo, xmax = SHR_hi), orientation = "y",
                  width = 0.25, colour = "grey30", linewidth = 0.5) +
    geom_point(size = 2.5) +
    scale_colour_manual(values = c("No" = "firebrick", "Yes" = "orange")) +
    scale_x_log10() +
    labs(
      title  = title,
      colour = "Significant?\n(BH-adj p < 0.05)",
      x      = "Subdistribution Hazard Ratio (95% CI)",
      y      = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 9))
}

p_forest_clinical  <- make_forest_plot(df_clinical,  "Clinical model: Fine-Gray SHRs")
p_forest_precision <- make_forest_plot(df_precision, "Precision model: Fine-Gray SHRs")

ggsave(file.path(plot_dir, "02.06_FG_forest_clinical.png"),
       p_forest_clinical,  width = 7, height = 3.5, dpi = 400, bg = "white")
ggsave(file.path(plot_dir, "02.07_FG_forest_precision.png"),
       p_forest_precision, width = 7, height = 4.5, dpi = 400, bg = "white")

prettify_hr <- function(df) {
  df %>%
    mutate(Variable = display, Model = label, SHR = round(SHR, 3),
           SHR_lo = round(SHR_lo, 2), SHR_hi = round(SHR_hi, 2),
           SHR = sprintf("%s [%s-%s]", SHR, SHR_lo, SHR_hi),
           z_value = round(z, 3),
           p_value = ifelse(p < 0.001, "<0.001", p),
           p_adj = ifelse(p_adj < 0.001, "<0.001", p_adj)) %>%
    dplyr::select(Variable, Model, SHR, z_value, p_value, p_adj)
}

# Clinical model equation (unscaled predictors):
#   CIF(t | x) = 1 - exp(-H0(t) * exp(sum(beta * x)))
# where H0(t) is the baseline cumulative subdistribution hazard at x = 0
crr_clin <- fg_clin_full_raw$crrFit
H0       <- sapply(eval_times, function(t) sum(crr_clin$bfitj[crr_clin$uftime <= t]))
equation_tbl <- tibble(
  term  = c(names(crr_clin$coef), paste0("H0(", eval_times, " years)")),
  value = c(unname(crr_clin$coef), H0)
)

# Check the equation reproduces predictRisk()
chk_rows <- model_data[1:200, ]
X_chk    <- model.matrix(reformulate(preds_clin), chk_rows)[, names(crr_clin$coef)]
risk_eq  <- 1 - exp(-H0[eval_times == 10] * exp(drop(X_chk %*% crr_clin$coef)))
stopifnot(max(abs(risk_eq - fg_risk(fg_clin_full_raw, chk_rows, 10))) < 1e-8)

equation_notes <- tibble(note = c(
  "CIF(t | x) = 1 - exp(-H0(t) * exp(sum(beta * x)))",
  "age in years; bmi in kg/m2; alcohol_grams_week in g/week",
  "sex2 = female (vs. male); has_t2dm1 = type 2 diabetes; smoking_binary1 = current smoker",
  sprintf("BMI winsorised to %.2f-%.2f kg/m2; alcohol capped at 500 g/week",
          min(model_data$bmi), max(model_data$bmi))
))

tables[["09_SHR_clinical"]]     <- prettify_hr(df_clinical)
tables[["10_SHR_precision"]]    <- prettify_hr(df_precision)
tables[["11_Clinical_equation"]] <- equation_tbl   # read by 03_RSF_sensitivity_analysis.R
tables[["12_Equation_notes"]]   <- equation_notes

###############################################################################
# SECTION 6: Export the deployable clinical Fine-Gray model
#
# The Shiny calculator needs to turn six predictor values into an absolute MALO
# risk. That requires only the coefficient vector and the baseline cumulative
# subdistribution hazard H0(t) -- nothing else. We therefore do NOT save the
# fitted FGR/crr object: it carries the model response (every participant's
# follow-up time and event status) and the per-subject score residuals, i.e.
# individual-level UK Biobank data that must not leave the approved
# environment. Instead we assemble a small, self-contained list:
#
#   * coefficients and their covariance matrix (aggregate quantities)
#   * H0(t) tabulated on a FIXED 0.05-year grid rather than at the observed
#     event times -- a set of unique event times would itself be individual-
#     level information, and the grid removes it at negligible accuracy cost
#   * the preprocessing constants needed to reproduce the training transforms
#     (BMI winsorising bounds, alcohol cap) and the factor coding
#   * aggregate cohort/performance metadata for display in the app
#   * a predict function whose environment is set to the stats namespace, so
#     that serialising it cannot drag the analysis workspace along with it
#
# Everything is then audited: the object must contain no element whose length
# matches the cohort size, and must reproduce predictRisk() exactly.
###############################################################################
cat("\n=== Section 6: Export deployable clinical model ===\n")

# ── 7a. Baseline cumulative subdistribution hazard on a fixed time grid ──
h0_grid_by  <- 0.05                                     # ~18 days
h0_grid_max <- ceiling(max(model_data$time_to_event))   # whole years, not an
                                                        # individual's max follow-up
t_grid <- round(seq(0, h0_grid_max, by = h0_grid_by), 2)
stopifnot(all(eval_times %in% t_grid))                  # 5 and 10 are exact knots

H0_grid <- vapply(t_grid,
                  function(tt) sum(crr_clin$bfitj[crr_clin$uftime <= tt]),
                  numeric(1))

# ── 7b. Prediction function (no captured data; see note above) ──
predict_malo_risk <- function(model, newdata, times = c(5, 10)) {
  newdata <- as.data.frame(newdata)
  need    <- model$predictors
  missing_vars <- setdiff(need, names(newdata))
  if (length(missing_vars))
    stop("Missing predictor(s): ", paste(missing_vars, collapse = ", "))

  # Reproduce the training-cohort preprocessing exactly
  bw <- model$preprocessing$bmi_winsor
  newdata$age <- as.numeric(newdata$age)
  newdata$bmi <- pmax(pmin(as.numeric(newdata$bmi), bw[2]), bw[1])
  newdata$alcohol_grams_week <- pmin(as.numeric(newdata$alcohol_grams_week),
                                     model$preprocessing$alcohol_cap)
  for (v in names(model$xlevels))
    newdata[[v]] <- factor(as.character(newdata[[v]]), levels = model$xlevels[[v]])
  if (anyNA(newdata[need]))
    stop("Predictors contain NA or a value outside the coded factor levels.")

  X  <- model.matrix(reformulate(need), newdata)[, names(model$coefficients), drop = FALSE]
  lp <- drop(X %*% model$coefficients)

  # H0 is a step function: hold the last grid value (constant interpolation)
  H0 <- approx(model$baseline$time, model$baseline$H0, xout = times,
               method = "constant", rule = 2, f = 0)$y

  out <- 1 - exp(-outer(exp(lp), H0))
  dimnames(out) <- list(NULL, paste0("risk_", times, "y"))
  out
}
# Detach the closure from the analysis workspace: it can now see only the
# stats namespace and base, never model_data or anything else defined here.
environment(predict_malo_risk) <- asNamespace("stats")

# ── 7c. Assemble the export object (aggregate quantities only) ──
q01_99 <- function(x) unname(round(quantile(x, c(0.01, 0.99), na.rm = TRUE), 2))

fg_clinical_deploy <- list(
  model_name  = "BFA clinical Fine-Gray model for MALO in low-fibrosis SLD",
  outcome     = "Major adverse liver outcome (cause 1); competing event: non-liver death (cause 2)",
  fitted_on   = "Full UK Biobank at-risk cohort (no train/test split)",
  date_created = as.character(Sys.Date()),
  r_version   = paste(R.version$major, R.version$minor, sep = "."),

  predictors   = preds_clin,
  coefficients = crr_clin$coef,                       # named, model.matrix scale
  vcov         = crr_clin$var,                        # for linear-predictor CIs only
  # plain data.frame, not a tibble: the app needs only base R + stats to use this
  baseline     = data.frame(time = t_grid, H0 = H0_grid), # step function, fixed grid

  preprocessing = list(
    bmi_winsor  = c(lower = min(model_data$bmi), upper = max(model_data$bmi)),
    alcohol_cap = 500,
    note = "BMI winsorised to the training 1st-99th percentiles; alcohol capped at 500 g/week"
  ),

  xlevels = list(
    sex            = levels(model_data$sex),            # "1" = male, "2" = female
    has_t2dm       = levels(model_data$has_t2dm),       # "0" = no, "1" = yes
    smoking_binary = levels(model_data$smoking_binary)  # "0" = never/previous, "1" = current
  ),

  input_coding = c(
    age                = "years",
    sex                = "1 = male, 2 = female",
    bmi                = "kg/m2",
    has_t2dm           = "0 = no, 1 = yes",
    alcohol_grams_week = "g/week",
    smoking_binary     = "0 = never or previous, 1 = current"
  ),

  # Plausible input ranges for the app's sliders: training 1st-99th percentiles,
  # NOT observed minima/maxima, so no individual's extreme value is exposed
  input_range = list(
    age                = q01_99(model_data$age),
    bmi                = q01_99(model_data$bmi),
    alcohol_grams_week = q01_99(model_data$alcohol_grams_week)
  ),

  # Median out-of-fold predicted 10-year risk: the Section 2 cut-point between
  # standard and aggressive surveillance
  risk_cut_10y = risk_cut_10y,

  cohort_summary = list(
    n                     = nrow(model_data),
    n_malo                = sum(model_data$status == 1),
    n_nonliver_death      = sum(model_data$status == 2),
    median_followup_years = round(median(model_data$time_to_event), 2)
  ),

  cv_performance = as.data.frame(perf_summary %>% filter(model == "FG (clinical)")),

  equation = "CIF(t | x) = 1 - exp(-H0(t) * exp(sum(beta * x)))",
  predict_fn = predict_malo_risk,
  usage = paste(
    'm <- readRDS("FG_clinical_model_deploy.rds")',
    'm$predict_fn(m, data.frame(age = 60, sex = "1", bmi = 31,',
    '             has_t2dm = "1", alcohol_grams_week = 120,',
    '             smoking_binary = "0"), times = c(5, 10))',
    sep = "\n"
  )
)

# ── 7d. Disclosure audit: nothing of cohort length may be in the object ──
cohort_n <- nrow(model_data)
audit_lengths <- function(x, n) {
  if (is.function(x)) return(0L)
  if (is.data.frame(x))
    return(as.integer(nrow(x) == n) + sum(vapply(x, audit_lengths, integer(1), n = n)))
  if (is.list(x)) return(sum(vapply(x, audit_lengths, integer(1), n = n)))
  as.integer(length(x) == n)
}
stopifnot(audit_lengths(fg_clinical_deploy, cohort_n) == 0L)
stopifnot(!any(grepl("eid", names(unlist(fg_clinical_deploy[
  setdiff(names(fg_clinical_deploy), "predict_fn")])), ignore.case = TRUE)))
# The predict function must not have captured the analysis workspace
stopifnot(identical(environment(fg_clinical_deploy$predict_fn), asNamespace("stats")))

# ── 7e. Correctness check: the export must reproduce predictRisk() ──
deploy_risk <- fg_clinical_deploy$predict_fn(fg_clinical_deploy, chk_rows, times = eval_times)
for (t in eval_times) {
  stopifnot(max(abs(deploy_risk[, paste0("risk_", t, "y")] -
                      fg_risk(fg_clin_full_raw, chk_rows, t))) < 1e-8)
}

saveRDS(fg_clinical_deploy, deploy_rds, compress = "xz")
cat("Deployable clinical model saved to:", deploy_rds,
    sprintf("(%.1f KB)\n", file.size(deploy_rds) / 1024))

###############################################################################
# SECTION 7: Save all tables and the OOF predictions
###############################################################################
write_xlsx(tables, results_xlsx)
cat("Tables saved to:", results_xlsx, "\n")

write_feather(all_oof, oof_path)
cat("OOF predictions saved to:", oof_path, "\n")

cat("\n=== All analyses complete ===\n")
cat("Results:", result_dir, "\n")
cat("Plots:  ", plot_dir,   "\n")
cat("Model:  ", deploy_rds, "\n")
