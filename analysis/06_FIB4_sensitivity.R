###############################################################################
# FIB-4 <1.3 sensitivity analysis
# Ronnie Li
#
# The primary cohort (and the deployed calculator) targets low fibrosis
# burden, FIB-4 < 2.67. This checks whether the Fine-Gray clinical model's
# discrimination and calibration hold up in the stricter FIB-4 < 1.3 subset
# (indeterminate/low-risk fibrosis by the standard two-cutoff rule), versus
# the 1.3-2.67 stratum immediately above it.
#
# Reuses the existing out-of-fold (OOF) Fine-Gray (clinical) predictions and
# the fold assignment from 02_FineGray_main_analysis.R -- no refitting. FIB-4
# itself is not in the OOF file (which carries predictions and outcomes, not
# raw covariates), so it is joined in from BFA_principal_data.feather by eid.
#
# 1. FIB-4 < 1.3 subset only (excludes 1.3-2.67): N, events, time-dependent
#    AUC at 5 and 10 years, and decile-based AJ-CIF calibration -- the same
#    metrics and methods as 02_FineGray_main_analysis.R Sections 1 and 4,
#    just restricted to this subset's OOF rows.
# 2. Events and 10-year cumulative incidence (Aalen-Johansen) by FIB-4
#    stratum (<1.3 vs 1.3-2.67), with Gray's test between them.
#
# Outputs: tables  -> results/06_FIB4_sensitivity.xlsx (sheets 01_, 02_, ...)
#          figures -> results/plots/06.01_*.png, 06.02_*.png
###############################################################################

library(timeROC)
library(survival)
library(cmprsk)
library(prodlim)
library(tidyverse)
library(arrow)
library(writexl)
library(ggsci)

# ── Directories ───────────────────────────────────────────────────────────────
project_dir <- "/mnt/d/Projects/BFA"
data_dir    <- file.path(project_dir, "data")
result_dir  <- file.path(project_dir, "results")
plot_dir    <- file.path(result_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

eval_times <- c(5, 10)
tables     <- list()

# ── Load the existing OOF predictions and join FIB-4 (no refitting) ───────────
all_oof <- read_feather(file.path(data_dir, "BFA_FineGray_OOF_predictions.feather"))

fib4_lookup <- read_feather(file.path(data_dir, "BFA_principal_data.feather")) %>%
  filter(included_in_cohort == TRUE) %>%
  transmute(eid, fib4)

all_oof <- all_oof %>% left_join(fib4_lookup, by = "eid")
stopifnot(!anyNA(all_oof$fib4))
stopifnot(all(all_oof$fib4 < 2.67))   # cohort inclusion criterion; sanity check

FIB4_LOW_CUT <- 1.3
strat_levels <- c("<1.3", "1.3-2.67")

all_oof <- all_oof %>%
  mutate(fib4_stratum = factor(if_else(fib4 < FIB4_LOW_CUT, "<1.3", "1.3-2.67"),
                               levels = strat_levels))

oof_low <- all_oof %>% filter(fib4_stratum == "<1.3")
cat(sprintf("FIB-4 < 1.3 subset: N = %d (of %d, %.1f%% of cohort)\n",
            nrow(oof_low), nrow(all_oof), 100 * nrow(oof_low) / nrow(all_oof)))

###############################################################################
# SECTION 1: FIB-4 < 1.3 subset -- N, events, time-dependent AUC, calibration
###############################################################################
cat("\n=== Section 1: FIB-4 < 1.3 subset performance ===\n")

# Time-dependent AUC for cause 1 (same helper/method as 02_FineGray_main_analysis.R)
td_auc <- function(df, marker, t) {
  roc <- timeROC(T = df$time_to_event, delta = df$status, marker = marker,
                 cause = 1, times = t, iid = FALSE)
  unname(roc$AUC_1[match(t, roc$times)])
}

subset_summary <- function(df, label) {
  tibble(
    subset             = label,
    n                  = nrow(df),
    n_malo             = sum(df$status == 1),
    n_nonliver_death   = sum(df$status == 2),
    n_malo_10y         = sum(df$time_to_event <= 10 & df$status == 1),
    median_followup_yr = round(median(df$time_to_event), 2),
    td_auc_5y          = round(td_auc(df, df$fg_clin_5, 5), 3),
    td_auc_10y         = round(td_auc(df, df$fg_clin_10, 10), 3)
  )
}

tbl_subset_summary <- bind_rows(
  subset_summary(oof_low, "FIB-4 < 1.3"),
  subset_summary(all_oof, "Full cohort (for comparison)")
)
print(tbl_subset_summary, width = Inf)
tables[["01_FIB4lt13_subset_summary"]] <- tbl_subset_summary

# ── Calibration: decile-based AJ-CIF, Fine-Gray (clinical), FIB-4 < 1.3 only ──
# Same method as 02_FineGray_main_analysis.R Section 4 (compute_aj_cif()):
# bin predicted risk into deciles within each fold, pool across folds, compare
# mean predicted CIF against observed Aalen-Johansen CIF per bin.
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

cal_low_df <- compute_aj_cif(oof_low, "fg_clin", "FG (clinical), FIB-4 < 1.3")
tables[["02_FIB4lt13_calibration_deciles"]] <- cal_low_df %>% dplyr::select(-time_label)

p_cal_low <- ggplot(cal_low_df, aes(x = mean_pred, y = obs_cif)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2, colour = "#2166AC") +
  geom_smooth(method = "loess", se = TRUE,
              colour = "#2166AC", fill = "#2166AC", alpha = 0.15, linewidth = 0.7) +
  facet_wrap(~ time_label, scales = "free") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(x = "Mean predicted CIF",
       y = "Observed CIF (Aalen-Johansen)",
       title = "Calibration: Fine-Gray clinical model, FIB-4 < 1.3 subset",
       subtitle = "Observed vs. predicted CIF across predicted risk deciles (pooled across folds)") +
  theme_bw(base_size = 12)

ggsave(file.path(plot_dir, "06.01_FIB4lt13_calibration.png"),
       p_cal_low, width = 8, height = 4.5, dpi = 400, bg = "white")

###############################################################################
# SECTION 2: Events and 10-year cumulative incidence by FIB-4 stratum
###############################################################################
cat("\n=== Section 2: Events and 10-year incidence by FIB-4 stratum ===\n")

ci_fib4 <- cmprsk::cuminc(ftime = all_oof$time_to_event, fstatus = all_oof$status,
                          group = all_oof$fib4_stratum)

# Same log-log-CI helper as 02_FineGray_main_analysis.R's cif_at_times()
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

cif_long_fib4 <- cif_at_times(ci_fib4) %>%
  mutate(event = if_else(cause == "1", "MALO", "Non-liver death"),
         across(c(est, lower, upper), ~ round(100 * .x, 2))) %>%
  dplyr::select(group, event, time, cif_pct = est, lower_pct = lower, upper_pct = upper)

counts_by_stratum <- all_oof %>%
  group_by(group = as.character(fib4_stratum)) %>%
  summarise(n                  = n(),
            n_malo             = sum(status == 1),
            pct_malo           = round(100 * mean(status == 1), 2),
            n_nonliver_death   = sum(status == 2),
            pct_nonliver_death = round(100 * mean(status == 2), 2),
            n_censored         = sum(status == 0),
            n_malo_10y         = sum(time_to_event <= 10 & status == 1),
            .groups = "drop")

cif_wide_fib4 <- cif_long_fib4 %>%
  mutate(value = sprintf("%.2f (%.2f-%.2f)", cif_pct, lower_pct, upper_pct),
         key   = paste0("CIF_", if_else(event == "MALO", "malo", "nonliver_death"),
                        "_", time, "y_pct_95CI")) %>%
  dplyr::select(group, key, value) %>%
  pivot_wider(names_from = key, values_from = value)

events_by_stratum <- counts_by_stratum %>%
  left_join(cif_wide_fib4, by = "group") %>%
  mutate(group = factor(group, levels = strat_levels)) %>%
  arrange(group)

grays_test_fib4 <- tibble(
  event     = c("MALO", "Non-liver death"),
  statistic = ci_fib4$Tests[c("1", "2"), "stat"],
  df        = ci_fib4$Tests[c("1", "2"), "df"],
  p_value   = ci_fib4$Tests[c("1", "2"), "pv"]
)

cat("Events and cumulative incidence by FIB-4 stratum:\n")
print(events_by_stratum, width = Inf)
print(grays_test_fib4)

tables[["03_events_by_FIB4_stratum"]]      <- events_by_stratum
tables[["04_CIF_by_FIB4_stratum_long"]]    <- cif_long_fib4
tables[["05_Grays_test_FIB4_stratum"]]     <- grays_test_fib4

# ── AJ-CIF curves by FIB-4 stratum (MALO), same style as 02's plot_cif() ──
cif_curves <- function(ci_obj, cause, group_levels) {
  bind_rows(lapply(names(ci_obj), function(nm) {
    if (!is.list(ci_obj[[nm]]) || !endsWith(nm, paste0(" ", cause))) return(NULL)
    tibble(time  = ci_obj[[nm]]$time,
           cif   = ci_obj[[nm]]$est,
           group = sub(" [0-9]+$", "", nm))
  })) %>%
    mutate(group = factor(group, levels = group_levels))
}

p_label <- function(p, test) {
  paste0(test, " p ", if (p < 0.001) "< 0.001" else paste0("= ", formatC(p, digits = 3, format = "f")))
}

df_cif_malo <- cif_curves(ci_fib4, 1, strat_levels)

p_cif_fib4 <- ggplot(df_cif_malo, aes(x = time, y = cif, colour = group)) +
  geom_step(linewidth = 0.8) +
  annotate("label", x = 0, y = max(df_cif_malo$cif), hjust = 0, vjust = 1, size = 3,
           label = p_label(ci_fib4$Tests["1", "pv"], "Gray's")) +
  scale_color_d3(labels = strat_levels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(x = "Time (years)", y = "Cumulative incidence (MALO)", colour = "FIB-4 stratum",
       title = "Aalen-Johansen CIF: MALO by FIB-4 stratum") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(plot_dir, "06.02_FIB4_stratum_AJ_CIF.png"),
       p_cif_fib4, width = 7, height = 5, dpi = 400, bg = "white")

# ── Write output ────────────────────────────────────────────────────────────
out_path <- file.path(result_dir, "06_FIB4_sensitivity.xlsx")
write_xlsx(tables, out_path)
cat("\nSaved:", out_path, "\n")
