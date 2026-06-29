# =============================================================================
# Fine-Gray Competing Risks Baseline Models using 5-fold Cross-Validation
# Outputs: SHR forest plots (full-data fit), CV AUC, C-index, Brier score,
#          calibration plots
# Ronnie Li
# =============================================================================

library(timeROC)
library(tidyverse)
library(arrow)
library(writexl)
library(prodlim)
library(survival)
library(riskRegression)
library(patchwork)
library(ggsci)
library(pec) # for ipcw()

# -----------------------------------------------------------------
# Directories
# -----------------------------------------------------------------
project_dir <- "/mnt/d/Projects/BFA"
source(file.path(project_dir, "src/BFA_CLivD_scores.R"))
data_dir    <- file.path(project_dir, "data")
result_dir <- file.path(project_dir, "results")
plot_dir <- file.path(result_dir, "plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# -----------------------------------------------------------------
# Winsorizing function caps the variable between the 1st and 99th percentiles
# -----------------------------------------------------------------
winsorise <- function(x, probs = c(0.01, 0.99)) {
  q <- quantile(x, probs, na.rm = TRUE)
  pmax(pmin(x, q[2]), q[1])
}

# -----------------------------------------------------------------
# Read and prepare data
# -----------------------------------------------------------------
adata <- read_feather(file.path(data_dir, "BFA_principal_data.feather"))

model_data <- adata %>%
  filter(included_in_cohort == TRUE) %>%
  mutate(
    bmi                = winsorise(bmi),
    hip_circ_raw       = hip_circ,
    waist_hip_ratio    = waist_circ / hip_circ,          # unwinsorised WHR for CLivD
    waist_circ         = winsorise(waist_circ),          # winsorised WC for precision Fine-Gray
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
  drop_na(age, sex, bmi, waist_circ, trig_hdl_ratio, alcohol_grams_week,
          smoking, time_to_event, status)


preds_clin <- c("age", "sex", "bmi", "has_t2dm", "alcohol_grams_week", "smoking_binary")
preds_prec <- c(preds_clin, "waist_circ", "PNPLA3_rs738409_G", "TM6SF2_rs58542926_T", "HSD17B13_rs9992651_A")
eval_times <- c(5, 10)

# -----------------------------------------------------------------
# Stratified folds
# -----------------------------------------------------------------
set.seed(1234)
n_folds    <- 5
malo_cases <- which(model_data$status == 1)
non_cases  <- which(model_data$status != 1)
model_data$fold <- NA_integer_
model_data$fold[malo_cases] <- sample(rep(1:n_folds, length.out = length(malo_cases)))
model_data$fold[non_cases]  <- sample(rep(1:n_folds, length.out = length(non_cases)))

# -----------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------

make_fg_formula <- function(preds) {
  as.formula(paste("Hist(time_to_event, status) ~", paste(preds, collapse = " + ")))
}

# Predict absolute CIF at a given time from an FGR model fit on training data,
# applied to test data.  Returns a numeric vector of length nrow(test_df).
predict_fg_cif <- function(fg_fit, test_df, t) {
  p <- predictRisk(fg_fit, newdata = test_df, times = t, cause = 1)
  as.numeric(p)
}

# IPCW Brier score for competing risks (cause 1), mirrors RSF script
compute_brier <- function(pred_risk, t, df) {
  ipcw_obj <- pec::ipcw(
    Surv(time_to_event, status != 0) ~ 1,
    data         = df,
    method       = "marginal",
    times        = t,
    subjectTimes = df$time_to_event
  )
  event_ind <- as.numeric(df$time_to_event <= t & df$status == 1)
  w <- numeric(nrow(df))
  for (i in seq_len(nrow(df))) {
    if (df$time_to_event[i] <= t && df$status[i] == 0) next   # censored before t
    if (df$time_to_event[i] <= t && df$status[i] != 0) {
      w[i] <- (event_ind[i] - pred_risk[i])^2 / ipcw_obj$IPCW.subjectTimes[i]
    } else {
      w[i] <- (0 - pred_risk[i])^2 / ipcw_obj$IPCW.times
    }
  }
  mean(w)
}

# -----------------------------------------------------------------
# 5-fold CV loop
# -----------------------------------------------------------------
fold_aucs    <- list()
fold_metrics <- list()
cal_data_clinical_list  <- list()
cal_data_precision_list <- list()

for (fold in 1:n_folds) {
  cat("\n=== Fold", fold, "===\n")
  
  train_df <- model_data %>% filter(fold != !!fold)
  test_df  <- model_data %>% filter(fold == !!fold)
  cat("  Train events:", sum(train_df$status == 1),
      "| Test events:", sum(test_df$status == 1), "\n")
  
  # --- Fine-Gray: clinical ---
  cat("  Fitting Fine-Gray (clinical)...\n")
  fg_clinical_fit <- FGR(make_fg_formula(preds_clin), data = train_df, cause = 1)
  
  # --- Fine-Gray: precision ---
  cat("  Fitting Fine-Gray (precision)...\n")
  fg_precision_fit <- FGR(make_fg_formula(preds_prec), data = train_df, cause = 1)
  
  # Predict CIF on test set at each eval time
  for (t in eval_times) {
    test_df[[paste0("fg_clinical_cif_",  t)]] <- predict_fg_cif(fg_clinical_fit,  test_df, t)
    test_df[[paste0("fg_precision_cif_", t)]] <- predict_fg_cif(fg_precision_fit, test_df, t)
  }
  
  # --- Time-dependent AUC ---
  cat("  Computing time-dependent AUCs...\n")
  compute_timeroc_fold <- function(marker) {
    timeROC(T = test_df$time_to_event, delta = test_df$status,
            marker = marker, cause = 1, times = eval_times, iid = FALSE)
  }
  
  roc_fg_clinical  <- compute_timeroc_fold(test_df$fg_clinical_cif_10)
  roc_fg_precision <- compute_timeroc_fold(test_df$fg_precision_cif_10)
  
  fold_aucs[[fold]] <- tibble(
    fold  = fold,
    model = c("FG (clinical)", "FG (clinical)",
              "FG (precision)", "FG (precision)"),
    time  = rep(eval_times, 2),
    auc   = c(roc_fg_clinical$AUC_1[1],  roc_fg_clinical$AUC_1[2],
              roc_fg_precision$AUC_1[1], roc_fg_precision$AUC_1[2])
  )
  
  # --- C-index (cause-specific, Harrell) ---
  test_surv <- Surv(test_df$time_to_event, test_df$status == 1)
  
  cindex_fg_clinical  <- concordance(test_surv ~ test_df$fg_clinical_cif_10,  reverse = TRUE)$concordance
  cindex_fg_precision <- concordance(test_surv ~ test_df$fg_precision_cif_10, reverse = TRUE)$concordance
  
  # --- Brier score ---
  cat("  Computing Brier scores...\n")
  brier_results <- tibble(
    model = c("FG (clinical)", "FG (precision)")
  )
  for (t in eval_times) {
    brier_results[[paste0("brier_t", t)]] <- c(
      compute_brier(test_df[[paste0("fg_clinical_cif_",    t)]], t, test_df),
      compute_brier(test_df[[paste0("fg_precision_cif_",   t)]], t, test_df)
    )
  }
  
  fold_metrics[[fold]] <- brier_results %>%
    mutate(
      fold   = fold,
      cindex = c(cindex_fg_clinical, cindex_fg_precision)
    )
  
  # --- Calibration data (collect across folds) ---
  for (t in eval_times) {
    for (mod in c("clinical", "precision")) {
      pred_col <- paste0("fg_", mod, "_cif_", t)
      cal_df <- tibble(
        pred      = test_df[[pred_col]],
        status    = test_df$status,
        time      = test_df$time_to_event,
        bin       = ntile(test_df[[pred_col]], 10),
        eval_time = t,
        fold      = fold
      )
      if (mod == "clinical") {
        cal_data_clinical_list[[paste0(fold, "_", t)]]  <- cal_df
      } else {
        cal_data_precision_list[[paste0(fold, "_", t)]] <- cal_df
      }
    }
  }
  
  gc()
}

# -----------------------------------------------------------------
# Summarize AUC
# -----------------------------------------------------------------
auc_summary <- bind_rows(fold_aucs) %>%
  group_by(model, time) %>%
  summarise(mean_auc = mean(auc), sd_auc = sd(auc), .groups = "drop")

# Summarize C-index and Brier
metrics_summary <- bind_rows(fold_metrics) %>%
  group_by(model) %>%
  summarise(
    mean_cindex    = mean(cindex),    sd_cindex    = sd(cindex),
    mean_brier_t5  = mean(brier_t5),  sd_brier_t5  = sd(brier_t5),
    mean_brier_t10 = mean(brier_t10), sd_brier_t10 = sd(brier_t10),
    .groups = "drop"
  )

print(auc_summary)
print(metrics_summary)

write_xlsx(
  list(AUC = auc_summary, Metrics = metrics_summary),
  file.path(result_dir, "BFA01.1.3_FineGray_CV_performance.xlsx")
)

# -----------------------------------------------------------------
# AUC comparison plot
# -----------------------------------------------------------------
model_order <- c("FG (clinical)", "FG (precision)")

auc_plot_data <- auc_summary %>%
  mutate(
    time_label = factor(paste0("t = ", time, " years"),
                        levels = paste0("t = ", eval_times, " years")),
    model = factor(model, levels = model_order),
    lower = mean_auc - sd_auc,
    upper = mean_auc + sd_auc
  )

p_auc <- ggplot(auc_plot_data, aes(x = model, y = mean_auc, colour = model)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  facet_wrap(~ time_label) +
  scale_color_d3() +
  scale_y_continuous(limits = c(0.5, 0.80),
                     labels = scales::number_format(accuracy = 0.01)) +
  labs(
    title    = "Time-dependent AUC: Fine-Gray regression",
    subtitle = "Mean \u00b1 SD across 5 cross-validation folds",
    x = NULL, y = "AUC (cause-specific)", colour = "Model"
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position  = "bottom",
    axis.text.x      = element_text(angle = 0, hjust = 0.5),
    strip.background = element_rect(fill = "grey90"))

print(p_auc)

ggsave(file.path(plot_dir, "BFA01.1.3_FineGray_AUC_comparison.png"),
       p_auc, width = 6, height = 4, dpi = 400, bg = "white")

# -----------------------------------------------------------------
# Calibration plot function (Aalen-Johansen within decile bins)
# -----------------------------------------------------------------
fg_calibration_plot <- function(cal_data_list, title = "Fine-Gray calibration") {
  
  cal_all <- bind_rows(cal_data_list)
  
  cal_summary <- cal_all %>%
    group_by(eval_time, bin) %>%
    summarise(
      mean_pred = mean(pred),
      obs_cif   = {
        sub  <- pick(everything())
        t_pt <- eval_time[1]
        if (sum(sub$status == 1) == 0) {
          0
        } else {
          fit          <- prodlim(Hist(time, status) ~ 1, data = sub)
          pred_prodlim <- predict(fit, times = t_pt, cause = 1)
          as.numeric(pred_prodlim)
        }
      },
      .groups = "drop"
    ) %>%
    mutate(
      time_label = factor(paste0("t = ", eval_time, " years"),
                          levels = paste0("t = ", eval_times, " years"))
    )
  
  ggplot(cal_summary, aes(x = mean_pred, y = obs_cif)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(size = 2.5, colour = "#2166AC") +
    geom_smooth(method = "loess", se = TRUE, colour = "#2166AC",
                fill = "#2166AC", alpha = 0.15, linewidth = 0.8) +
    facet_wrap(~ time_label, scales = "free") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    labs(
      title    = title,
      subtitle = "Observed vs predicted CIF (cause 1: MALO) across deciles",
      x = "Mean predicted CIF", y = "Observed CIF (Aalen-Johansen)"
    ) +
    theme_bw(base_size = 10) +
    theme(strip.background = element_rect(fill = "grey90"))
}

p_cal_clinical  <- fg_calibration_plot(cal_data_clinical_list,
                                       "Fine-Gray clinical model calibration")
p_cal_precision <- fg_calibration_plot(cal_data_precision_list,
                                       "Fine-Gray precision model calibration")

ggsave(file.path(plot_dir, "BFA01.1.3_FineGray_calibration_clinical.png"),
       p_cal_clinical,  width = 7, height = 4, dpi = 400, bg = "white")
ggsave(file.path(plot_dir, "BFA01.1.3_FineGray_calibration_precision.png"),
       p_cal_precision, width = 7, height = 4, dpi = 400, bg = "white")

# -----------------------------------------------------------------
# Full-data Fine-Gray fit --> SHR forest plots
# (fitted once on entire dataset for coefficient reporting)
# -----------------------------------------------------------------
cat("\nFitting full-data Fine-Gray models for SHR reporting...\n")

# Scale data to have more reasonable SHRs
model_data_scaled <- model_data %>%
  mutate(
    age        = age / 3,          # coef: SHR per 3 years
    bmi        = bmi / 3,          # coef: SHR per 3 kg/m²
    waist_circ = waist_circ / 3,    # coef: SHR per 3 cm
    alcohol_grams_week = alcohol_grams_week / 20   # coef: alcohol per 20 g
  )

# Fit full models
fg_clinical_full  <- FGR(make_fg_formula(preds_clin), data = model_data_scaled, cause = 1)
fg_precision_full <- FGR(make_fg_formula(preds_prec), data = model_data_scaled, cause = 1)

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
  bmi                  = "BMI (per 3 kg/m\u00b2)",
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

# Create data frames of results for plotting
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
print(p_forest_clinical)
print(p_forest_precision)

ggsave(file.path(plot_dir, "BFA01.1.3_FineGray_forest_clinical.png"),
       p_forest_clinical,  width = 7, height = 3.5, dpi = 400, bg = "white")
ggsave(file.path(plot_dir, "BFA01.1.3_FineGray_forest_precision.png"),
       p_forest_precision, width = 7, height = 4.5, dpi = 400, bg = "white")

# Save scaling parameters
# saveRDS(scale_params, file.path(result_dir, "BFA01.1.3_Fine-Gray_scaling_parameters.rds"))

# Save hazard ratios 
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

write_xlsx(list("FG_clinical" = prettify_hr(df_clinical), 
                "FG_precision" = prettify_hr(df_precision)),
           file.path(result_dir, "BFA01.1.3_Fine-Gray_hazard_ratios.xlsx"))

cat("\nDone. All outputs written to:", outdir, "\n")