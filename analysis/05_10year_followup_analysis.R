###############################################################################
# 10-year follow-up recompute of surveillance-rule capture and flag rates
# Ronnie Li
#
# The Results currently report capture/flag rates over FULL follow-up (any
# event, at any time -- median follow-up is ~15-18 years; see
# src/02_FineGray_main_analysis.R, "SOC rule referral burden" table):
#   model (>median predicted 10y risk)         74.9% captured / 50.0% flagged
#   composite (T2DM or >= 2 CMRFs)             88.2% captured / 79.0% flagged
#   T2DM only                                  21.9% captured / 11.5% flagged
#
# This script recomputes the same quantities restricted to events observed by
# 10 years, using the existing out-of-fold (OOF) Fine-Gray predictions -- no
# refitting. It reuses fg_clin_10 (already a 10-year predicted risk) and the
# soc_t2dm / soc_cmrf2 / soc_guideline flags exported alongside it.
#
# CASE DEFINITION (10-year capture):
#   A "10-year case" is time_to_event <= 10 & status == 1 (MALO within 10
#   years). This mirrors how the full-follow-up figures above are computed
#   (status == 1, unrestricted by time) but adds the 10-year cutoff. Capture
#   (sensitivity) only conditions on realized cases, so it needs no IPCW
#   weighting for censoring: an event that has already happened is observed
#   exactly, regardless of what happens to other participants after year 10.
#
#   Participants censored before 10 years without an event (administrative
#   censoring or dropout) have unknown 10-year MALO status. They cannot be
#   cases (no event was observed), but a handful may have gone on to have
#   MALO had they been followed to year 10 -- an unavoidable limitation of a
#   fixed time-window analysis, shared with the full-follow-up figures above
#   (which do not have this issue because follow-up there is not truncated).
#   The count is reported in the cohort-summary sheet for transparency.
#
#   Competing events (non-liver death) by 10 years are, correctly, not cases.
#
# FLAG DEFINITIONS (unchanged from the full-follow-up figures; flags are
# static, baseline-covariate rules and do not depend on the follow-up window):
#   model (0.5%)      fg_clin_10 >= 0.005 (the guideline-facing surveillance
#                      threshold used in the Shiny app)
#   model (median)     fg_clin_10 > median(fg_clin_10) (the cohort-median cut
#                      used in Section 2 of 02_FineGray_main_analysis.R, ~0.47%,
#                      which is what the cited full-follow-up 74.9%/50.0%
#                      figures actually used, reported loosely as "0.5%" in
#                      the Results text -- both are reported below so the two
#                      can be compared directly against each other and against
#                      the original full-follow-up figures)
#   composite  soc_guideline == 1  (T2DM or >= 2 CMRFs; requested explicitly
#              because this column is exported to all_oof but is NOT one of
#              the curves fed into decision curve analysis in Section 3 of
#              02_FineGray_main_analysis.R, which uses only the T2DM-alone and
#              CMRFs-alone columns)
#   T2D only   soc_t2dm == 1
#
# Outputs:
#   results/05_10year_followup.xlsx
#     01_cohort_summary          n, 10y cases, competing events, censored-early
#     02_10yr_capture_flag       capture/flag at the four rules above, 10y
#     03_matched_flag_rates      model thresholds matched to the composite's
#                                and T2D-only's flag rates, with 10y capture
#     04_crosstab_model_x_composite   2x2 N and 10y events per cell
#     05_full_vs_10yr_comparison      each rule's capture at full follow-up
#                                (any time) side by side with 10-year capture
#     06_DCA_net_benefit         net benefit at 0.25/0.5/1/2% thresholds, 10y,
#                                for the model, treat-all, treat-none and the
#                                composite guideline rule
#     07_DCA_net_interv_avoided_05pct   interventions avoided per 100 vs
#                                treat-all, at the 0.5% threshold
#     08_DCA_bootstrap_diffs     bootstrap 95% CI for net-benefit differences
#                                (model - composite, model - treat-all) at
#                                each threshold (B = 1000 resamples)
###############################################################################

library(tidyverse)
library(arrow)
library(writexl)
library(dcurves)
library(survival)

# ── Directories ───────────────────────────────────────────────────────────────
project_dir <- "/mnt/d/Projects/BFA"
data_dir    <- file.path(project_dir, "data")
result_dir  <- file.path(project_dir, "results")

# ── Load the existing out-of-fold predictions (no refitting) ──────────────────
all_oof <- read_feather(file.path(data_dir, "BFA_FineGray_OOF_predictions.feather"))

MODEL_THRESHOLD_10Y <- 0.005                  # 0.5%, matches app.R's RISK_THRESHOLD_10Y
MODEL_MEDIAN_10Y    <- median(all_oof$fg_clin_10)   # cohort-median cut, as in Section 2

all_oof <- all_oof %>%
  mutate(
    case_10y        = time_to_event <= 10 & status == 1,
    case_any        = status == 1,                      # for the full-follow-up comparison
    censored_10y    = time_to_event <  10 & status == 0, # unknown 10y status
    flag_model      = fg_clin_10 >= MODEL_THRESHOLD_10Y,
    flag_model_med  = fg_clin_10 >  MODEL_MEDIAN_10Y,
    flag_composite  = soc_guideline == 1,
    flag_t2dm       = soc_t2dm == 1
  )
cat(sprintf("Cohort-median predicted 10y risk cut-point: %.4f%%\n", 100 * MODEL_MEDIAN_10Y))

# ── 01. Cohort summary ─────────────────────────────────────────────────────────
cohort_summary <- tibble(
  metric = c(
    "N (total, OOF cohort)",
    "N 10-year MALO cases (time_to_event <= 10 & status == 1)",
    "N 10-year competing events (time_to_event <= 10 & status == 2)",
    "N censored before 10y with unknown 10y status (time_to_event < 10 & status == 0)",
    "N any-time MALO cases (status == 1, unrestricted -- for comparison to current Results)"
  ),
  value = c(
    nrow(all_oof),
    sum(all_oof$case_10y),
    sum(all_oof$time_to_event <= 10 & all_oof$status == 2),
    sum(all_oof$censored_10y),
    sum(all_oof$status == 1)
  )
)
print(cohort_summary, n = Inf)

# ── 02. 10-year capture and flag rate for each rule ────────────────────────────
capture_flag <- function(flag, case, rule) {
  tibble(
    rule           = rule,
    n_flagged      = sum(flag),
    pct_flagged    = round(100 * mean(flag), 1),
    n_10y_cases    = sum(case),
    n_captured     = sum(flag & case),
    pct_capture_10y = round(100 * sum(flag & case) / sum(case), 1)
  )
}

tbl_capture_flag <- bind_rows(
  capture_flag(all_oof$flag_model,      all_oof$case_10y, "Model (predicted 10y risk >= 0.5%)"),
  capture_flag(all_oof$flag_model_med,  all_oof$case_10y, "Model (predicted 10y risk > cohort median, ~0.47%)"),
  capture_flag(all_oof$flag_composite,  all_oof$case_10y, "Composite (T2DM or >= 2 CMRFs)"),
  capture_flag(all_oof$flag_t2dm,       all_oof$case_10y, "T2DM only")
)
print(tbl_capture_flag, width = Inf)

# ── 03. Matched flag rates: model thresholds matched to composite / T2DM-only ─
# Choose the model threshold that flags exactly as many participants (by count,
# for precision) as each comparator rule, then report 10-year capture there.
match_threshold_capture <- function(score, n_target, case, label) {
  thr    <- sort(score, decreasing = TRUE)[n_target]
  flag   <- score >= thr
  tibble(
    matched_to        = label,
    model_threshold    = thr,
    n_flagged          = sum(flag),
    pct_flagged        = round(100 * mean(flag), 1),
    n_10y_cases        = sum(case),
    n_captured         = sum(flag & case),
    pct_capture_10y    = round(100 * sum(flag & case) / sum(case), 1)
  )
}

tbl_matched <- bind_rows(
  match_threshold_capture(all_oof$fg_clin_10, sum(all_oof$flag_composite), all_oof$case_10y,
                           "Composite flag rate (T2DM or >= 2 CMRFs)"),
  match_threshold_capture(all_oof$fg_clin_10, sum(all_oof$flag_t2dm), all_oof$case_10y,
                           "T2DM-only flag rate")
)
print(tbl_matched, width = Inf)

# ── 04. 2x2 cross-tab: model flag x composite flag ─────────────────────────────
tbl_crosstab <- all_oof %>%
  mutate(
    model_flag     = if_else(flag_model, "Model: flagged (>=0.5%)", "Model: not flagged (<0.5%)"),
    composite_flag = if_else(flag_composite, "Composite: flagged", "Composite: not flagged")
  ) %>%
  count(model_flag, composite_flag, name = "n") %>%
  left_join(
    all_oof %>%
      mutate(
        model_flag     = if_else(flag_model, "Model: flagged (>=0.5%)", "Model: not flagged (<0.5%)"),
        composite_flag = if_else(flag_composite, "Composite: flagged", "Composite: not flagged")
      ) %>%
      group_by(model_flag, composite_flag) %>%
      summarise(n_10y_events = sum(case_10y), .groups = "drop"),
    by = c("model_flag", "composite_flag")
  ) %>%
  arrange(desc(model_flag), desc(composite_flag))
print(tbl_crosstab, width = Inf)

stopifnot(sum(tbl_crosstab$n) == nrow(all_oof))
stopifnot(sum(tbl_crosstab$n_10y_events) == sum(all_oof$case_10y))

# ── 05. Full follow-up vs 10-year, side by side ────────────────────────────────
# Same rules and flags as sheet 02, but each rule's capture is reported both
# over full follow-up (case_any, matching the cited Results figures exactly)
# and restricted to 10 years (case_10y), so the two are directly comparable.
# Flag rate is identical in both columns because flags are static, baseline
# rules -- only the case definition (and hence the capture denominator/
# numerator) changes with the follow-up window.
capture_pct <- function(flag, case) round(100 * sum(flag & case) / sum(case), 1)

tbl_comparison <- tibble(
  rule = c(
    "Model (predicted 10y risk >= 0.5%)",
    "Model (predicted 10y risk > cohort median, ~0.47%)",
    "Composite (T2DM or >= 2 CMRFs)",
    "T2DM only"
  ),
  pct_flagged = c(
    round(100 * mean(all_oof$flag_model), 1),
    round(100 * mean(all_oof$flag_model_med), 1),
    round(100 * mean(all_oof$flag_composite), 1),
    round(100 * mean(all_oof$flag_t2dm), 1)
  ),
  pct_capture_full_followup = c(
    capture_pct(all_oof$flag_model,     all_oof$case_any),
    capture_pct(all_oof$flag_model_med, all_oof$case_any),
    capture_pct(all_oof$flag_composite, all_oof$case_any),
    capture_pct(all_oof$flag_t2dm,      all_oof$case_any)
  ),
  pct_capture_10yr = c(
    capture_pct(all_oof$flag_model,     all_oof$case_10y),
    capture_pct(all_oof$flag_model_med, all_oof$case_10y),
    capture_pct(all_oof$flag_composite, all_oof$case_10y),
    capture_pct(all_oof$flag_t2dm,      all_oof$case_10y)
  )
)
print(tbl_comparison, width = Inf)

# The median-cut row must reproduce the cited full-follow-up Results exactly
# (74.9% captured / 50.0% flagged) -- a check that this script's case/flag
# definitions match what Section 2/3 of 02_FineGray_main_analysis.R used.
stopifnot(tbl_comparison$pct_flagged[2] == 50.0,
          tbl_comparison$pct_capture_full_followup[2] == 74.9,
          tbl_comparison$pct_flagged[3] == 79.0,
          tbl_comparison$pct_capture_full_followup[3] == 88.2,
          tbl_comparison$pct_flagged[4] == 11.5,
          tbl_comparison$pct_capture_full_followup[4] == 21.9)

# ── Write output ────────────────────────────────────────────────────────────
tables <- list(
  "01_cohort_summary"                = cohort_summary,
  "02_10yr_capture_flag"             = tbl_capture_flag,
  "03_matched_flag_rates"            = tbl_matched,
  "04_crosstab_model_x_composite"    = tbl_crosstab,
  "05_full_vs_10yr_comparison"       = tbl_comparison
)
out_path <- file.path(result_dir, "05_10year_followup.xlsx")
write_xlsx(tables, out_path)
cat("\nSaved:", out_path, "\n")

###############################################################################
# DECISION CURVE ANALYSIS (10 years): net benefit at fixed thresholds, with
# bootstrap CIs for the model-vs-standard-of-care differences
#
# Reuses the Fine-Gray (clinical) 10-year OOF predictions above; no
# refitting. Uses the composite guideline rule (T2DM or >= 2 CMRFs), for the
# same reason as the capture table above and the fix now applied to
# 02_FineGray_main_analysis.R Section 3: the composite is the rule guidelines
# actually specify, and using it here keeps this table, the capture table,
# and the corrected DCA figure all consistent with each other.
#
# Net benefit follows the Vickers decision-curve formula, extended to
# time-to-event data via dcurves::dca()'s default Kaplan-Meier estimator at
# t = 10 years -- the same method 02_FineGray_main_analysis.R's DCA uses.
# As there, competing events (non-liver death) are treated as censored for
# this calculation, a known limitation of the standard survival-DCA method
# used throughout this project (Kaplan-Meier overestimates cause-1 incidence
# under competing risks relative to the Aalen-Johansen CIF used in Section 2
# of that script), not something introduced here.
#
# Bootstrap: 1000 resamples of the full OOF cohort (with replacement),
# matching the n_boot = 1000 convention in 03_RSF_sensitivity_analysis.R.
# Percentile 95% CIs for the two requested differences: model minus
# composite, and model minus treat-all (dcurves' built-in reference
# strategy, not a separate column).
###############################################################################
cat("\n=== Decision curve analysis (10 years) ===\n")

DCA_THRESHOLDS <- c(0.0025, 0.005, 0.01, 0.02)   # 0.25%, 0.5%, 1%, 2%
RULE_MODEL     <- "Fine-Gray (clinical)"
RULE_COMPOSITE <- "Composite (T2DM or >=2 CMRFs)"

dca_dat <- all_oof %>%
  transmute(
    time_to_event,
    status,
    `Fine-Gray (clinical)`          = fg_clin_10,
    `Composite (T2DM or >=2 CMRFs)` = soc_guideline
  )

dca_fit <- function(data) {
  dca(Surv(time_to_event, status == 1) ~ `Fine-Gray (clinical)` +
        `Composite (T2DM or >=2 CMRFs)`,
      data = data, time = 10, thresholds = DCA_THRESHOLDS)
}

dca_point <- dca_fit(dca_dat)

# ── 06. Net benefit at each threshold, for model / treat-all / treat-none / composite ─
tbl_dca_nb <- as_tibble(dca_point) %>%
  filter(threshold %in% DCA_THRESHOLDS) %>%
  transmute(
    rule          = as.character(label),
    threshold_pct = 100 * threshold,
    n             = n,
    tp_rate       = tp_rate,
    fp_rate       = fp_rate,
    net_benefit   = net_benefit
  ) %>%
  arrange(threshold_pct, rule)
print(tbl_dca_nb, n = Inf)

# ── 07. Net interventions avoided vs treat-all, per 100, at 0.5% ─────────────
tbl_nia_05 <- as_tibble(net_intervention_avoided(dca_point, nper = 100)) %>%
  filter(threshold == 0.005) %>%
  transmute(
    rule                              = as.character(label),
    threshold_pct                     = 100 * threshold,
    net_interventions_avoided_per_100 = net_intervention_avoided
  ) %>%
  arrange(desc(net_interventions_avoided_per_100))
print(tbl_nia_05, n = Inf)

# ── 08. Bootstrap CIs for the net-benefit differences ────────────────────────
n_boot  <- 1000
n_cores <- if (.Platform$OS.type == "windows") 1L else
  max(1L, min(8L, parallel::detectCores() - 1L))
n_obs   <- nrow(dca_dat)
set.seed(1234)

cat("  Bootstrapping DCA net-benefit differences (B =", n_boot, ", cores =", n_cores, ")...\n")

boot_mat <- do.call(rbind, parallel::mclapply(seq_len(n_boot), function(b) {
  idx <- sample.int(n_obs, n_obs, replace = TRUE)
  nb  <- as_tibble(dca_fit(dca_dat[idx, ])) %>%
    filter(threshold %in% DCA_THRESHOLDS) %>%
    select(label, threshold, net_benefit) %>%
    arrange(threshold) %>%
    pivot_wider(names_from = label, values_from = net_benefit) %>%
    arrange(threshold)
  setNames(
    c(nb[[RULE_MODEL]] - nb[[RULE_COMPOSITE]], nb[[RULE_MODEL]] - nb[["Treat All"]]),
    c(paste0("model_minus_composite_", DCA_THRESHOLDS),
      paste0("model_minus_treatall_",  DCA_THRESHOLDS))
  )
}, mc.cores = n_cores))
stopifnot(nrow(boot_mat) == n_boot, !anyNA(boot_mat))

point_nb <- tbl_dca_nb %>% select(rule, threshold_pct, net_benefit) %>%
  pivot_wider(names_from = rule, values_from = net_benefit) %>%
  arrange(threshold_pct)

tbl_boot_diffs <- bind_rows(
  tibble(
    comparison    = "Model - Composite",
    threshold_pct = 100 * DCA_THRESHOLDS,
    diff          = point_nb[[RULE_MODEL]] - point_nb[[RULE_COMPOSITE]],
    boot_se       = apply(boot_mat[, paste0("model_minus_composite_", DCA_THRESHOLDS)], 2, sd),
    ci_lower      = apply(boot_mat[, paste0("model_minus_composite_", DCA_THRESHOLDS)], 2, quantile, 0.025),
    ci_upper      = apply(boot_mat[, paste0("model_minus_composite_", DCA_THRESHOLDS)], 2, quantile, 0.975)
  ),
  tibble(
    comparison    = "Model - Treat All",
    threshold_pct = 100 * DCA_THRESHOLDS,
    diff          = point_nb[[RULE_MODEL]] - point_nb[["Treat All"]],
    boot_se       = apply(boot_mat[, paste0("model_minus_treatall_", DCA_THRESHOLDS)], 2, sd),
    ci_lower      = apply(boot_mat[, paste0("model_minus_treatall_", DCA_THRESHOLDS)], 2, quantile, 0.025),
    ci_upper      = apply(boot_mat[, paste0("model_minus_treatall_", DCA_THRESHOLDS)], 2, quantile, 0.975)
  )
) %>%
  mutate(boot_p_value = 2 * pnorm(-abs(diff) / boot_se)) %>%
  arrange(comparison, threshold_pct)
print(tbl_boot_diffs, n = Inf)

# ── Append the DCA sheets to the same workbook ───────────────────────────────
tables_dca <- c(tables, list(
  "06_DCA_net_benefit"                     = tbl_dca_nb,
  "07_DCA_net_interv_avoided_05pct" = tbl_nia_05,
  "08_DCA_bootstrap_diffs"                  = tbl_boot_diffs
))
write_xlsx(tables_dca, out_path)
cat("\nSaved (with DCA sheets):", out_path, "\n")
