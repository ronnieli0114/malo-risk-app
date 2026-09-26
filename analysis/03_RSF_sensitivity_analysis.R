###############################################################################
# Sensitivity analysis using random survival forest (RSF)
# Ronnie Li
#
# The Fine-Gray models in 02_FineGray_main_analysis.R use linear predictors
# only. The RSF allows nonlinear effects and interactions, and is used here to
# check whether this meaningfully improves discrimination and overall accuracy.
#
# Run after 02_FineGray_main_analysis.R: re-uses its fold assignment, its
# out-of-fold (OOF) Fine-Gray predictions and its clinical model equation, so
# both models are evaluated on exactly the same folds.
#
# Forests are grown with rfsrc.anonymous() and then stripped of all remaining
# participant-level training data (see strip_training_data()), so the saved
# model objects contain no training data.
#
# 1. 5-fold CV: RSF (clinical, precision) — AUC, C-index, Brier
# 2. Fine-Gray vs RSF: per-fold differences and bootstrap test of the
#    difference in pooled OOF time-dependent AUC, C-index and Brier score
# 3. Decision curve analysis: Fine-Gray vs RSF
# 4. Calibration plots (decile-based, AJ-CIF) for RSF clinical and precision
# 5. Full-data anonymous RSF training and saving (clinical + precision)
# 6. Partial dependence: linear Fine-Gray vs RSF
# 7. Variable importance: VIMP + SHAP beeswarm/bar plots
#
# Outputs: tables  -> results/03_RSF_results.xlsx (sheets 01_, 02_, ...)
#          figures -> results/plots/03.01_*.png, 03.02_*.png, ...
###############################################################################

library(randomForestSRC)
library(timeROC)
library(survival)
library(pec)       # for ipcw()
library(prodlim)
library(tidyverse)
library(arrow)
library(writexl)
library(readxl)
library(patchwork)
library(dcurves)
library(ggsci)
library(fastshap)
library(shapviz)

# ── Directories ───────────────────────────────────────────────────────────────
project_dir <- "/mnt/d/Projects/BFA"
data_dir   <- file.path(project_dir, "data")
result_dir <- file.path(project_dir, "results")
plot_dir   <- file.path(result_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
oof_path   <- file.path(data_dir, "BFA_FineGray_OOF_predictions.feather")
fg_results_xlsx <- file.path(result_dir, "02_FineGray_results.xlsx")
results_xlsx    <- file.path(result_dir, "03_RSF_results.xlsx")

tables <- list()   # every output table; written to one workbook at the end

# Bootstrap settings for the Fine-Gray vs RSF comparison
n_boot    <- 1000
boot_seed <- 20260914
n_cores   <- if (.Platform$OS.type == "windows") 1L else
  max(1L, min(8L, parallel::detectCores() - 1L))

# ── Helpers ───────────────────────────────────────────────────────────────────
winsorise <- function(x, probs = c(0.01, 0.99)) {
  q <- quantile(x, probs, na.rm = TRUE)
  pmax(pmin(x, q[2]), q[1])
}

extract_cif <- function(pred_obj, t, cause = 1) {
  ti <- which.min(abs(pred_obj$time.interest - t))
  pred_obj$cif[, ti, cause]
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

# IPCW weights and cause-1 indicator for the Brier score at time t
# (computed once per data set and shared across models)
brier_weights <- function(df, t) {
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
  list(y = as.numeric(df$time_to_event <= t & df$status == 1),
       w = ifelse(df$time_to_event <= t & df$status != 0, 1 / G_subject,            # any event by t
           ifelse(df$time_to_event > t, 1 / iw$IPCW.times, 0)))                     # censored before t: 0
}

# IPCW Brier score for competing risks (cause 1)
compute_brier <- function(pred_risk, t, df, bw = brier_weights(df, t)) {
  mean(bw$w * (bw$y - pred_risk)^2)
}

# rfsrc.anonymous() removes the predictors (xvar) but still keeps the training
# outcomes (yvar, event.info$time/cens), per-subject in-bag/OOB ensembles,
# censoring weights and — through the model call — the entire training data
# frame. None of these are needed by predict(), so remove them all.
strip_training_data <- function(fit) {
  fit[intersect(c("call", "yvar", "predicted", "predicted.oob", "chf", "chf.oob",
                  "cif", "cif.oob", "survival", "survival.oob", "hazard", "hazard.oob",
                  "membership", "inbag", "proximity", "distance", "imputed.data",
                  "imputed.indv", "case.depth", "forest.wt", "subj", "uno.weights"),
                names(fit))] <- NULL
  fit$forest[intersect(c("yvar", "case.wt", "subj", "uno.weights"), names(fit$forest))] <- NULL
  slim_event_info <- function(ei) {
    ei$time  <- NULL
    ei$cens  <- sort(unique(ei$cens))   # predict() only checks the set of event codes
    ei$event <- sort(unique(ei$event))
    ei
  }
  fit$event.info         <- slim_event_info(fit$event.info)
  fit$forest$event.info  <- slim_event_info(fit$forest$event.info)
  fit$forest$impute.mean <- fit$forest$impute.mean[fit$xvar.names]
  # exact name checks (`$` would partially match e.g. yvar -> yvar.names)
  stopifnot(!any(c("xvar", "yvar", "call") %in% names(fit)),
            !any(c("xvar", "yvar") %in% names(fit$forest)))
  fit
}

# Anonymous competing-risks RSF with the subsampling settings of rfsrc.fast()
# (used in the original analysis). Only the model columns are passed so that
# no other variables reach the stored imputation summaries.
fit_rsf <- function(preds, data, importance = FALSE, ...) {
  fit <- rfsrc.anonymous(
    as.formula(paste("Surv(time_to_event, status) ~", paste(preds, collapse = "+"))),
    data = as.data.frame(data[, c("time_to_event", "status", preds)]),
    ntree = 300, nodesize = 50, nsplit = 5, ntime = 50,
    bootstrap = "by.root", samptype = "swor",
    sampsize = function(x) min(x * .632, max(150, x^(3/4))),
    splitrule = "logrankCR", cause = 1, importance = importance, seed = -1, ...
  )
  strip_training_data(fit)
}

predict_rsf <- function(fit, newdata) {
  predict(fit, newdata = as.data.frame(newdata[, fit$xvar.names]))
}

# ── Data preparation (identical to 02_FineGray_main_analysis.R) ───────────────
adata <- read_feather(file.path(data_dir, "BFA_principal_data.feather"))

model_data <- adata %>%
  filter(included_in_cohort == TRUE) %>%
  mutate(
    bmi                = winsorise(bmi),
    hip_circ_raw       = hip_circ,
    waist_hip_ratio    = waist_circ / hip_circ,
    waist_circ         = winsorise(waist_circ),
    triglycerides      = winsorise(triglycerides),
    hdl                = winsorise(hdl),
    trig_hdl_ratio     = triglycerides / hdl,
    alcohol_grams_week = pmin(alcohol_grams_week, 500),
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
  drop_na(age, sex, bmi, waist_circ, trig_hdl_ratio, alcohol_grams_week,
          smoking, time_to_event, status)

# Predictor sets
preds_clin <- c("age", "sex", "bmi", "has_t2dm", "alcohol_grams_week", "smoking_binary")
preds_prec <- c(preds_clin, "waist_circ", "PNPLA3_rs738409_G", "TM6SF2_rs58542926_T", "HSD17B13_rs9992651_A")
eval_times <- c(5, 10)

# Fold assignment and Fine-Gray OOF predictions from the main analysis
fg_oof <- read_feather(oof_path)
stopifnot(setequal(model_data$eid, fg_oof$eid))
model_data <- model_data %>%
  inner_join(fg_oof %>% dplyr::select(eid, fold, starts_with("fg_"),
                                      oof_time = time_to_event, oof_status = status),
             by = "eid")
stopifnot(all(model_data$time_to_event == model_data$oof_time),
          all(model_data$status == model_data$oof_status))
n_folds <- max(model_data$fold)

rsf_pred_sets <- list(rsf_clin = preds_clin, rsf_prec = preds_prec)

# Models evaluated (label -> column prefix)
perf_models <- c("FG (clinical)"   = "fg_clin",  "RSF (clinical)"  = "rsf_clin",
                 "FG (precision)"  = "fg_prec",  "RSF (precision)" = "rsf_prec")

###############################################################################
# SECTION 1: 5-fold CV — RSF (same folds as Fine-Gray)
###############################################################################
cat("\n=== Section 1: RSF CV loop ===\n")

fold_metric_list <- list()
oof_list         <- list()

for (fold in 1:n_folds) {
  cat("  Fold", fold, "\n")
  tr <- model_data[model_data$fold != fold, ]
  te <- model_data[model_data$fold == fold, ]

  for (m in names(rsf_pred_sets)) {
    fit <- fit_rsf(rsf_pred_sets[[m]], tr)
    pr  <- predict_rsf(fit, te)
    for (t in eval_times) te[[paste0(m, "_", t)]] <- extract_cif(pr, t)
    rm(fit, pr)
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

  oof_list[[fold]] <- te %>%
    dplyr::select(eid, fold, time_to_event, status,
                  all_of(as.vector(outer(unname(perf_models), eval_times, paste, sep = "_"))))
  gc()
}

all_oof      <- bind_rows(oof_list)
fold_metrics <- bind_rows(fold_metric_list)

perf_summary <- fold_metrics %>%
  group_by(model) %>%
  summarise(across(c(auc_t5, auc_t10, cindex, brier_t5, brier_t10),
                   list(mean = mean, sd = sd), .names = "{.fn}_{.col}"),
            .groups = "drop") %>%
  mutate(model = factor(model, levels = names(perf_models))) %>%
  arrange(model)
print(perf_summary)

###############################################################################
# SECTION 2: Fine-Gray (linear) vs RSF (nonlinear) — statistical comparison
# (a) per-fold paired differences (RSF − FG), (b) nonparametric bootstrap
#     (resampling participants) of the pooled OOF predictions: difference in
#     time-dependent AUC and Brier score at 5 and 10 years and in C-index, with
#     95% percentile CI and a two-sided Wald p-value using the bootstrap SE.
#     Higher AUC/C-index and lower Brier score are better.
###############################################################################
cat("\n=== Section 2: Fine-Gray vs RSF ===\n")

comparisons <- list(Clinical  = c(FG = "fg_clin", RSF = "rsf_clin"),
                    Precision = c(FG = "fg_prec", RSF = "rsf_prec"))
metric_labels <- c(auc_t5   = "AUC (t = 5 years)",   auc_t10   = "AUC (t = 10 years)",
                   cindex   = "C-index",
                   brier_t5 = "Brier (t = 5 years)", brier_t10 = "Brier (t = 10 years)")

# (a) Per-fold paired differences
fold_delta <- fold_metrics %>%
  mutate(comparison = str_extract(model, "(?<=\\()[a-z]+(?=\\))") %>% str_to_title(),
         method     = str_extract(model, "^[A-Z]+")) %>%
  dplyr::select(fold, comparison, method, all_of(names(metric_labels))) %>%
  pivot_longer(all_of(names(metric_labels)), names_to = "metric") %>%
  pivot_wider(names_from = method, values_from = value) %>%
  mutate(delta = RSF - FG)

fold_delta_summary <- fold_delta %>%
  group_by(comparison, metric) %>%
  summarise(mean_FG = mean(FG), mean_RSF = mean(RSF),
            mean_delta = mean(delta), sd_delta = sd(delta), .groups = "drop")

# (b) Bootstrap of pooled OOF predictions
pooled_metrics <- function(df) {
  bw  <- lapply(setNames(eval_times, eval_times), function(t) brier_weights(df, t))
  out <- c()
  for (cmp in names(comparisons)) for (mod in names(comparisons[[cmp]])) {
    k <- comparisons[[cmp]][[mod]]
    key <- function(met) paste(cmp, mod, met, sep = "|")
    out[key("auc_t5")]    <- td_auc(df, df[[paste0(k, "_5")]],  5)
    out[key("auc_t10")]   <- td_auc(df, df[[paste0(k, "_10")]], 10)
    out[key("cindex")]    <- c_index(df, df[[paste0(k, "_10")]])
    out[key("brier_t5")]  <- compute_brier(df[[paste0(k, "_5")]],  5,  df, bw[["5"]])
    out[key("brier_t10")] <- compute_brier(df[[paste0(k, "_10")]], 10, df, bw[["10"]])
  }
  out
}

boot_df   <- as.data.frame(all_oof)
point_est <- pooled_metrics(boot_df)

cat("  Bootstrapping pooled OOF metrics (B =", n_boot, ", cores =", n_cores, ")...\n")
boot_mat <- do.call(rbind, parallel::mclapply(seq_len(n_boot), function(b) {
  set.seed(boot_seed + b)
  pooled_metrics(boot_df[sample.int(nrow(boot_df), replace = TRUE), ])
}, mc.cores = n_cores))

boot_tests <- bind_rows(lapply(names(comparisons), function(cmp) {
  bind_rows(lapply(names(metric_labels), function(met) {
    fg_key  <- paste(cmp, "FG",  met, sep = "|")
    rsf_key <- paste(cmp, "RSF", met, sep = "|")
    d_boot  <- boot_mat[, rsf_key] - boot_mat[, fg_key]
    delta   <- point_est[[rsf_key]] - point_est[[fg_key]]
    p_value <- 2 * pnorm(-abs(delta / sd(d_boot)))
    tibble(
      comparison   = cmp,
      metric       = metric_labels[[met]],
      FG           = point_est[[fg_key]],
      FG_lower     = quantile(boot_mat[, fg_key],  0.025, names = FALSE),
      FG_upper     = quantile(boot_mat[, fg_key],  0.975, names = FALSE),
      RSF          = point_est[[rsf_key]],
      RSF_lower    = quantile(boot_mat[, rsf_key], 0.025, names = FALSE),
      RSF_upper    = quantile(boot_mat[, rsf_key], 0.975, names = FALSE),
      delta_RSF_minus_FG = delta,
      delta_lower  = quantile(d_boot, 0.025, names = FALSE),
      delta_upper  = quantile(d_boot, 0.975, names = FALSE),
      boot_se      = sd(d_boot),
      p_value      = p_value,
      significant  = if_else(p_value < 0.05, "Yes", "No")
    )
  }))
})) %>%
  mutate(n_boot = n_boot)

cat("Fine-Gray vs RSF (pooled OOF, bootstrap):\n")
print(boot_tests, width = Inf)

tables[["01_FG_vs_RSF_bootstrap"]]      <- boot_tests
tables[["02_FG_vs_RSF_per_fold_delta"]] <- fold_delta_summary
tables[["03_Performance"]]              <- perf_summary
tables[["04_Performance_per_fold"]]     <- fold_metrics

# Pooled OOF estimates with bootstrap 95% CI
p_cmp <- boot_tests %>%
  dplyr::select(comparison, metric, FG_est = FG, FG_lower, FG_upper,
                RSF_est = RSF, RSF_lower, RSF_upper) %>%
  pivot_longer(-c(comparison, metric), names_to = c("model", ".value"), names_sep = "_") %>%
  mutate(metric = factor(metric, levels = metric_labels),
         comparison = factor(comparison, levels = names(comparisons))) %>%
  ggplot(aes(x = comparison, y = est, colour = model)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5)) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_color_d3(labels = c(FG = "Fine-Gray (linear)", RSF = "RSF")) +
  labs(x = "Predictor set", y = "Estimate (pooled OOF, bootstrap 95% CI)", colour = NULL,
       title = "Linear Fine-Gray vs. random survival forest") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "grey90"),
        axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(file.path(plot_dir, "03.01_FG_vs_RSF_performance.png"),
       p_cmp, width = 12, height = 4.2, dpi = 400, bg = "white")

###############################################################################
# SECTION 3: Decision curve analysis — Fine-Gray vs RSF (clinical)
###############################################################################
cat("\n=== Section 3: Decision curve analysis ===\n")

thresholds <- seq(0, 0.03, by = 0.001)
dca_tables <- list()

for (t in eval_times) {
  dca_dat <- all_oof %>%
    transmute(time_to_event, status,
              `Fine-Gray (clinical)` = .data[[paste0("fg_clin_",  t)]],
              `RSF (clinical)`       = .data[[paste0("rsf_clin_", t)]])

  dca_fr <- dca(Surv(time_to_event, status == 1) ~ `Fine-Gray (clinical)` + `RSF (clinical)`,
                data = dca_dat, time = t, thresholds = thresholds)
  p_dca <- dca_fr %>%
    plot(smooth = TRUE) +
    scale_color_d3() +
    labs(title = paste0("DCA at t = ", t, " years: Fine-Gray vs RSF"),
         x = "Threshold probability", y = "Net benefit") +
    theme_bw(base_size = 8) +
    theme(legend.position = "bottom")

  fig_no <- 1 + match(t, eval_times)   # 03.02 (t = 5), 03.03 (t = 10)
  ggsave(file.path(plot_dir, sprintf("03.%02d_FG_vs_RSF_DCA_t%d.png", fig_no, t)),
         p_dca, width = 5.5, height = 4.5, dpi = 400, bg = "white")
  dca_tables[[as.character(t)]] <- as_tibble(dca_fr) %>% mutate(eval_time = t)
}
tables[["05_DCA_net_benefit"]] <- bind_rows(dca_tables)

###############################################################################
# SECTION 4: Calibration plots — RSF clinical and precision
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
  compute_aj_cif(all_oof, "rsf_clin", "RSF (clinical)"),
  compute_aj_cif(all_oof, "rsf_prec", "RSF (precision)")
)
tables[["06_Calibration_deciles"]] <- cal_all_df %>% dplyr::select(-time_label)

p_cal <- ggplot(cal_all_df, aes(x = mean_pred, y = obs_cif)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2, colour = "steelblue") +
  geom_smooth(method = "loess", se = TRUE,
              colour = "steelblue", fill = "steelblue", alpha = 0.15, linewidth = 0.7) +
  facet_grid(model ~ time_label, scales = "free") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(x = "Mean predicted CIF",
       y = "Observed CIF (Aalen-Johansen)",
       title = "Calibration plots: RSF clinical and precision models",
       subtitle = "Observed vs. predicted CIF across predicted risk deciles (pooled across folds)") +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey90"))

ggsave(file.path(plot_dir, "03.04_RSF_calibration.png"),
       p_cal, width = 8, height = 6, dpi = 400, bg = "white")

###############################################################################
# SECTION 5: Full-data anonymous RSF training and saving
###############################################################################
cat("\n=== Section 5: Full-dataset RSF training ===\n")

rsf_clin_path <- file.path(result_dir, "RSF_final_clinical_model.rds")
cat("  Fitting full RSF clinical model...\n")
final_rsf_clin <- fit_rsf(preds_clin, model_data, importance = TRUE, do.trace = 60)
saveRDS(final_rsf_clin, rsf_clin_path)
cat("  Saved to:", rsf_clin_path, "\n")

rsf_prec_path <- file.path(result_dir, "RSF_final_precision_model.rds")
cat("  Fitting full RSF precision model (waist_circ)...\n")
final_rsf_prec <- fit_rsf(preds_prec, model_data, importance = TRUE, do.trace = 60)
saveRDS(final_rsf_prec, rsf_prec_path)
cat("  Saved to:", rsf_prec_path, "\n")

###############################################################################
# SECTION 6: Partial dependence — linear Fine-Gray vs RSF (clinical model)
# Both curves use the same subsample and grid, on the same scale (10-year
# MALO CIF), so any curvature in the RSF shows what the linear model omits.
###############################################################################
cat("\n=== Section 6: Partial dependence ===\n")

# Average predicted 10-year CIF over a random subsample with `var` fixed at each
# grid value
partial_dependence <- function(pred_fun, data, var, grid_n = 30, n_sub = 1000) {
  set.seed(1234)
  sub  <- data[sample(nrow(data), min(n_sub, nrow(data))), ]
  grid <- seq(quantile(data[[var]], 0.02, na.rm = TRUE),
              quantile(data[[var]], 0.98, na.rm = TRUE), length.out = grid_n)
  cif_vals <- sapply(grid, function(val) {
    sub[[var]] <- val
    mean(pred_fun(sub))
  })
  tibble(variable = var, x = grid, cif = cif_vals, cif_centered = cif_vals - mean(cif_vals))
}

rsf_cif10 <- function(newdata) extract_cif(predict_rsf(final_rsf_clin, newdata), 10)

# Full-data Fine-Gray clinical model, from the equation saved by script 02
fg_eq    <- read_xlsx(fg_results_xlsx, sheet = "11_Clinical_equation")
fg_beta  <- with(filter(fg_eq, !startsWith(term, "H0")), setNames(value, term))
fg_H0_10 <- fg_eq$value[fg_eq$term == "H0(10 years)"]
fg_cif10 <- function(newdata) {
  X <- model.matrix(reformulate(preds_clin), newdata)[, names(fg_beta), drop = FALSE]
  1 - exp(-fg_H0_10 * exp(drop(X %*% fg_beta)))
}

var_labels <- c(age = "Age (years)",
                bmi = "BMI (kg/m²)",
                alcohol_grams_week = "Alcohol (g/week)")

cat("  Computing partial dependence (might take a while)...\n")
pdp_data <- bind_rows(lapply(names(var_labels), function(v) bind_rows(
  partial_dependence(fg_cif10,  model_data, v) %>% mutate(model = "Fine-Gray (linear)"),
  partial_dependence(rsf_cif10, model_data, v) %>% mutate(model = "RSF")
)))
tables[["07_Partial_dependence"]] <- pdp_data

p_pdp <- pdp_data %>%
  mutate(variable = var_labels[variable]) %>%
  ggplot(aes(x = x, y = cif_centered, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ variable, scales = "free_x", nrow = 1) +
  scale_color_d3() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(x = NULL, y = "Partial 10-year CIF\n(centered)", colour = NULL,
       title = "Partial dependence of 10-year MALO risk: linear Fine-Gray vs RSF") +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey90"), legend.position = "bottom")

ggsave(file.path(plot_dir, "03.05_FG_vs_RSF_partial_dependence.png"),
       p_pdp, width = 9, height = 3.8, dpi = 400, bg = "white")

###############################################################################
# SECTION 7: Variable importance — VIMP + SHAP
###############################################################################
cat("\n=== Section 7: Variable importance ===\n")

variable_importance_analysis <- function(rsf_fit, predictor_list, model_label) {

  # ── VIMP plot ──────────────────────────────────────────────────────────────
  vimp_df <- data.frame(
    variable   = rownames(rsf_fit$importance),
    importance = as.numeric(rsf_fit$importance[, 1])
  ) %>%
    arrange(importance) %>%
    mutate(
      variable = case_when(
        variable == "bmi"                  ~ "BMI",
        variable == "alcohol_grams_week"   ~ "Alcohol (g/week)",
        variable == "smoking_binary"       ~ "Smoking",
        variable == "sex"                  ~ "Sex",
        variable == "has_t2dm"             ~ "Type 2 diabetes",
        variable == "age"                  ~ "Age",
        variable == "waist_circ"           ~ "Waist circumference (cm)",
        variable == "waist_hip_ratio"      ~ "Waist:Hip ratio",
        variable == "trig_hdl_ratio"       ~ "Triglyceride:HDL ratio",
        variable == "PNPLA3_rs738409_G"    ~ "PNPLA3 rs738409:G",
        variable == "TM6SF2_rs58542926_T"  ~ "TM6SF2 rs58542926:T",
        variable == "HSD17B13_rs9992651_A" ~ "HSD17B13 rs9992651:A",
        TRUE ~ variable
      ),
      variable = factor(variable, levels = variable)
    )

  p_vimp <- ggplot(vimp_df, aes(x = importance, y = variable)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    labs(x = "Variable importance (MALO)", y = NULL,
         title = paste("Variable importance:", model_label)) +
    theme_bw(base_size = 9)

  # ── SHAP values ────────────────────────────────────────────────────────────
  pred_fun <- function(object, newdata) {
    p        <- predict(object, newdata = newdata)
    time_idx <- which.min(abs(p$time.interest - 10))
    p$cif[, time_idx, 1]
  }

  X <- model_data %>% dplyr::select(all_of(predictor_list))
  X[] <- lapply(X, function(x) as.numeric(as.character(x)))

  set.seed(1234)
  idx <- sample(nrow(X), min(1000, nrow(X)))

  cat("  Computing SHAP for", model_label, "(", length(idx), "rows, nsim = 20)...\n")
  shap_values <- fastshap::explain(
    rsf_fit,
    X            = X,
    pred_wrapper = pred_fun,
    nsim         = 20,
    newdata      = as.data.frame(X[idx, ])
  )

  column_mapping <- c(
    "age"                  = "Age",
    "bmi"                  = "BMI",
    "alcohol_grams_week"   = "Alcohol (g/week)",
    "smoking_binary"       = "Smoking",
    "sex"                  = "Sex",
    "has_t2dm"             = "Type 2 diabetes",
    "waist_circ"           = "Waist circumference (cm)",
    "waist_hip_ratio"      = "Waist:Hip ratio",
    "trig_hdl_ratio"       = "Triglyceride:HDL ratio",
    "PNPLA3_rs738409_G"    = "PNPLA3 rs738409:G",
    "TM6SF2_rs58542926_T"  = "TM6SF2 rs58542926:T",
    "HSD17B13_rs9992651_A" = "HSD17B13 rs9992651:A"
  )

  sv <- shapviz(shap_values, X = as.data.frame(X[idx, ]))
  # Rename columns to human-readable labels
  shared_cols <- intersect(colnames(sv$X), names(column_mapping))
  colnames(sv$X)[match(shared_cols, colnames(sv$X))] <- column_mapping[shared_cols]
  colnames(sv$S)[match(shared_cols, colnames(sv$S))] <- column_mapping[shared_cols]

  p_beeswarm <- sv_importance(sv, kind = "bee") +
    theme_bw(base_size = 9) +
    labs(title = paste("SHAP beeswarm:", model_label))

  p_barplot <- sv_importance(sv, kind = "bar") +
    theme_bw(base_size = 9) +
    labs(title = paste("Mean |SHAP|:", model_label))

  color_var <- if ("PNPLA3 rs738409:G" %in% colnames(sv$S)) "PNPLA3 rs738409:G" else "Smoking"
  p_dependence <- sv_dependence(sv, v = "Alcohol (g/week)", color_var = color_var) +
    theme_bw(base_size = 9) +
    labs(title = "SHAP dependence: Alcohol (g/week)",
         x = "Alcohol (g/week)", color = color_var)

  list(plots = list(VIMP_plot       = p_vimp,
                    beeswarm_plot   = p_beeswarm,
                    bar_plot        = p_barplot,
                    dependence_plot = p_dependence),
       vimp  = vimp_df %>% mutate(model = model_label, variable = as.character(variable)))
}

vimp_clin <- variable_importance_analysis(final_rsf_clin, preds_clin, "RSF clinical")
vimp_prec <- variable_importance_analysis(final_rsf_prec, preds_prec, "RSF precision")
tables[["08_VIMP"]] <- bind_rows(vimp_clin$vimp, vimp_prec$vimp)

# Save plots: 03.06_ onwards
fig_no <- 6
for (spec in list(list(res = vimp_clin, tag = "clinical"),
                  list(res = vimp_prec, tag = "precision"))) {
  for (nm in names(spec$res$plots)) {
    fname <- sprintf("03.%02d_RSF_%s_%s.png", fig_no, spec$tag, nm)
    ggsave(file.path(plot_dir, fname),
           plot = spec$res$plots[[nm]],
           width = 6, height = 5, dpi = 400, bg = "white")
    fig_no <- fig_no + 1
  }
}

write_xlsx(tables, results_xlsx)
cat("Tables saved to:", results_xlsx, "\n")

cat("\n=== All analyses complete ===\n")
cat("Results:", result_dir, "\n")
cat("Plots:  ", plot_dir,   "\n")
