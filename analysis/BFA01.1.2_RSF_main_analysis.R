###############################################################################
# Main analysis using random survival forest
# Ronnie Li
#
# 1. RSF vs CLivD vs FIB-4 time-dependent AUC, Brier score, C-index (precision uses waist_circ)
# 2. Extended clinical RSF (+ HTN, HDL, TG, HbA1c): C-index comparison
# 3. Nonlinearity justification: LRT + RSF PDPs vs cause-specific Cox splines
# 4. Risk-stratified CIF (Aalen-Johansen) and KM curves by predicted tertile
# 5. Decision curve analysis — Panel A (RSF vs FG) and Panel B (RSF vs single vars)
# 6. Calibration plots (decile-based, AJ-CIF) for RSF clinical and precision
# 7. Full-dataset model training and saving (clinical + precision)
# 8. Variable importance: VIMP + SHAP beeswarm/bar plots
###############################################################################

library(randomForestSRC)
library(riskRegression)
library(timeROC)
library(survival)
library(pec)
library(cmprsk)
library(prodlim)
library(tidyverse)
library(arrow)
library(writexl)
library(patchwork)
library(splines)
library(broom)
library(dcurves)
library(ggsci)
library(fastshap)
library(shapviz)

# ── Directories ───────────────────────────────────────────────────────────────
project_dir <- "/mnt/d/Projects/BFA"
source(file.path(project_dir, "src/BFA_CLivD_scores.R"))
data_dir   <- file.path(project_dir, "data")
result_dir <- file.path(project_dir, "results")
plot_dir   <- file.path(result_dir, "plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ── Helpers ───────────────────────────────────────────────────────────────────
winsorise <- function(x, probs = c(0.01, 0.99)) {
  q <- quantile(x, probs, na.rm = TRUE)
  pmax(pmin(x, q[2]), q[1])
}

extract_cif <- function(pred_obj, t, cause = 1) {
  ti <- which.min(abs(pred_obj$time.interest - t))
  pred_obj$cif[, ti, cause]
}

compute_brier <- function(pred_risk, t, df) {
  iw <- pec::ipcw(Surv(time_to_event, status != 0) ~ 1, data = df,
                  method = "marginal", times = t, subjectTimes = df$time_to_event)
  ev <- as.numeric(df$time_to_event <= t & df$status == 1)
  w  <- numeric(nrow(df))
  for (i in seq_len(nrow(df))) {
    if (df$time_to_event[i] <= t && df$status[i] == 0) next
    if (df$time_to_event[i] <= t)
      w[i] <- (ev[i] - pred_risk[i])^2 / iw$IPCW.subjectTimes[i]
    else
      w[i] <- (0 - pred_risk[i])^2 / iw$IPCW.times
  }
  mean(w)
}

# Cause-specific Cox linear predictor (LP) -> probability at time t
lp_to_prob <- function(cox_fit, lp_vec, t) {
  bh <- basehaz(cox_fit, centered = FALSE)
  h0 <- bh$hazard[which.min(abs(bh$time - t))]
  1 - exp(-h0 * exp(lp_vec))
}

# ── Data preparation ──────────────────────────────────────────────────────────
adata <- read_feather(file.path(data_dir, "BFA_principal_data.feather"))

model_data <- adata %>%
  filter(included_in_cohort == TRUE) %>%
  mutate(
    bmi                = winsorise(bmi),
    hip_circ_raw       = hip_circ,
    waist_hip_ratio    = waist_circ / hip_circ,          # unwinsorised WHR for CLivD
    waist_circ         = winsorise(waist_circ),          # winsorised WC for precision RSF
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

# Predictor sets
preds_clin <- c("age", "sex", "bmi", "has_t2dm", "alcohol_grams_week", "smoking_binary")
preds_prec <- c(preds_clin, "waist_circ", "PNPLA3_rs738409_G", "TM6SF2_rs58542926_T", "HSD17B13_rs9992651_A")
preds_ext  <- c(preds_clin, "has_hypertension", "trig_hdl_ratio", "hba1c")
eval_times <- c(5, 10)

# Stratified 5-fold CV (same seed as existing scripts)
set.seed(1234)
n_folds   <- 5
malo_idx  <- which(model_data$status == 1)
rest_idx  <- which(model_data$status != 1)
model_data$fold <- NA_integer_
model_data$fold[malo_idx] <- sample(rep(1:n_folds, length.out = length(malo_idx)))
model_data$fold[rest_idx] <- sample(rep(1:n_folds, length.out = length(rest_idx)))

###############################################################################
# SECTION 1 & DCA DATA: 5-fold CV
###############################################################################
cat("\n=== Section 1 CV loop ===\n")

fold_aucs        <- list()
fold_cidx        <- list()
fold_brier       <- list()
cal_clin_list    <- list()
cal_prec_list    <- list()
all_test         <- data.frame()

for (fold in 1:n_folds) {
  cat("  Fold", fold, "\n")
  tr    <- model_data[model_data$fold != fold, ]
  te    <- model_data[model_data$fold == fold, ]
  tr_cs <- mutate(tr, ev1 = as.integer(status == 1))
  
  # ── RSF clinical ──
  fit_rc <- rfsrc.fast(
    as.formula(paste("Surv(time_to_event, status) ~", paste(preds_clin, collapse = "+"))),
    data = tr, ntree = 300, nodesize = 50, nsplit = 5, forest = TRUE, seed = -1,
    save.memory = TRUE, cause = 1, splitrule = "logrankCR", importance = FALSE
  )
  pr_rc <- predict(fit_rc, newdata = te)
  for (t in eval_times) te[[paste0("rsf_clin_", t)]] <- extract_cif(pr_rc, t)
  
  # ── RSF precision (waist_circ) ──
  fit_rp <- rfsrc.fast(
    as.formula(paste("Surv(time_to_event, status) ~", paste(preds_prec, collapse = "+"))),
    data = tr, ntree = 300, nodesize = 50, nsplit = 5, forest = TRUE, seed = -1,
    save.memory = TRUE, cause = 1, splitrule = "logrankCR", importance = FALSE
  )
  pr_rp <- predict(fit_rp, newdata = te)
  for (t in eval_times) te[[paste0("rsf_prec_", t)]] <- extract_cif(pr_rp, t)
  
  # ── Fine-Gray clinical (for decision curve analysis, DCA) ──
  fit_fg <- FGR(
    as.formula(paste("Hist(time_to_event, status) ~", paste(preds_clin, collapse = "+"))),
    data = tr, cause = 1
  )
  for (t in eval_times) {
    te[[paste0("fg_clin_", t)]] <- as.numeric(
      predictRisk(fit_fg, newdata = te, times = t, cause = 1)
    )
  }
  
  # ── CLivD linear predictors ──
  te$clivd_nl_lp <- compute_clivd_nonlab(te)
  te$clivd_l_lp  <- compute_clivd(te)
  
  # ── FIB-4 cause-specific Cox ──
  fit_fib4   <- coxph(Surv(time_to_event, ev1) ~ fib4, data = tr_cs)
  te$fib4_lp <- predict(fit_fib4, newdata = te, type = "lp")
  
  # ── CLivD / FIB-4 linear predictors to probabilities at eval_times ──
  tr_nl  <- mutate(tr_cs, lp = compute_clivd_nonlab(tr_cs))
  tr_l   <- mutate(tr_cs, lp = compute_clivd(tr_cs))
  cs_nl  <- coxph(Surv(time_to_event, ev1) ~ lp, data = tr_nl)
  cs_l   <- coxph(Surv(time_to_event, ev1) ~ lp, data = tr_l)
  nl_lp  <- compute_clivd_nonlab(te)
  l_lp   <- compute_clivd(te)
  for (t in eval_times) {
    te[[paste0("clivd_nl_", t)]]  <- lp_to_prob(cs_nl,   nl_lp,        t)
    te[[paste0("clivd_l_", t)]]   <- lp_to_prob(cs_l,    l_lp,         t)
    te[[paste0("fib4_", t)]]      <- lp_to_prob(fit_fib4, te$fib4_lp,  t)
  }
  
  # ── Single-variable models for DCA Panel B ──
  for (sv in c("alcohol_grams_week", "bmi")) {
    fit_sv <- coxph(as.formula(paste("Surv(time_to_event, ev1) ~", sv)), data = tr_cs)
    lp_sv  <- predict(fit_sv, newdata = te, type = "lp")
    for (t in eval_times)
      te[[paste0(sv, "_sv_", t)]] <- lp_to_prob(fit_sv, lp_sv, t)
  }
  fit_tdm   <- coxph(Surv(time_to_event, ev1) ~ has_t2dm, data = tr_cs)
  lp_tdm    <- predict(fit_tdm, newdata = te, type = "lp")
  for (t in eval_times)
    te[[paste0("t2dm_sv_", t)]] <- lp_to_prob(fit_tdm, lp_tdm, t)
  
  # ── Per-fold AUC ──
  auc_markers <- list(
    "RSF (clinical)"      = te$rsf_clin_10,
    "RSF (precision)"     = te$rsf_prec_10,
    "CLivD (non-lab)"     = te$clivd_nl_lp,
    "CLivD (lab)"         = te$clivd_l_lp,
    "FIB-4"              = te$fib4_lp
  )
  fold_auc_rows <- bind_rows(lapply(names(auc_markers), function(m) {
    roc <- timeROC(T = te$time_to_event, delta = te$status,
                   marker = auc_markers[[m]], cause = 1,
                   times = eval_times, iid = FALSE)
    tibble(fold = fold, model = m, time = eval_times, auc = roc$AUC_1)
  }))
  fold_aucs[[fold]] <- fold_auc_rows
  
  # ── Per-fold C-index ──
  ts <- Surv(te$time_to_event, te$status == 1)
  fold_cidx[[fold]] <- tibble(
    fold   = fold,
    model  = names(auc_markers),
    cindex = sapply(auc_markers, function(m)
      concordance(ts ~ m, reverse = TRUE)$concordance)
  )
  
  # ── Per-fold Brier score ──
  brier_probs <- list(
    "RSF (clinical)"  = list("5" = te$rsf_clin_5,  "10" = te$rsf_clin_10),
    "RSF (precision)" = list("5" = te$rsf_prec_5,  "10" = te$rsf_prec_10),
    "CLivD (non-lab)" = list("5" = te$clivd_nl_5,  "10" = te$clivd_nl_10),
    "CLivD (lab)"     = list("5" = te$clivd_l_5,   "10" = te$clivd_l_10),
    "FIB-4"          = list("5" = te$fib4_5,       "10" = te$fib4_10)
  )
  fold_brier[[fold]] <- tibble(
    fold  = fold,
    model = names(brier_probs),
    brier_t5  = sapply(brier_probs, function(p) compute_brier(p[["5"]],  5,  te)),
    brier_t10 = sapply(brier_probs, function(p) compute_brier(p[["10"]], 10, te))
  )
  
  # ── Calibration data (RSF clinical and precision) ──
  for (t in eval_times) {
    for (mod in c("clin", "prec")) {
      pc <- paste0("rsf_", mod, "_", t)
      cal_df <- tibble(
        pred      = te[[pc]],
        status    = te$status,
        time      = te$time_to_event,
        bin       = ntile(te[[pc]], 10),
        eval_time = t,
        fold      = fold
      )
      key <- paste0(fold, "_", t)
      if (mod == "clin") cal_clin_list[[key]] <- cal_df
      else               cal_prec_list[[key]] <- cal_df
    }
  }
  
  all_test <- bind_rows(all_test, te)
  rm(fit_rc, fit_rp, fit_fg, pr_rc, pr_rp)
  gc()
}

# ── Summarize Section 1 ──
auc_summary <- bind_rows(fold_aucs) %>%
  group_by(model, time) %>%
  summarise(mean_auc = mean(auc), sd_auc = sd(auc), .groups = "drop") %>%
  mutate(model = factor(model, levels = c("RSF (clinical)", "RSF (precision)",
                                          "CLivD (non-lab)", "CLivD (lab)", "FIB-4")))

cindex_summary <- bind_rows(fold_cidx) %>%
  group_by(model) %>%
  summarise(mean_cindex = mean(cindex), sd_cindex = sd(cindex), .groups = "drop")

brier_summary <- bind_rows(fold_brier) %>%
  group_by(model) %>%
  summarise(
    mean_brier_t5  = mean(brier_t5),  sd_brier_t5  = sd(brier_t5),
    mean_brier_t10 = mean(brier_t10), sd_brier_t10 = sd(brier_t10),
    .groups = "drop"
  )

# Combined performance table
perf_combined <- auc_summary %>%
  pivot_wider(names_from = time, values_from = c(mean_auc, sd_auc),
              names_glue = "{.value}_t{time}") %>%
  left_join(cindex_summary, by = "model") %>%
  left_join(brier_summary,  by = "model")

write_xlsx(
  list(
    Performance = perf_combined,
    AUC_by_time = auc_summary,
    C_index     = cindex_summary,
    Brier_score = brier_summary
  ),
  file.path(result_dir, "BFA01.1.2_sec1_model_comparison.xlsx")
)
cat("Section 1 output written.\n")
print(perf_combined)

# AUC plot
p_auc <- auc_summary %>%
  mutate(model_label = gsub(" ", "\n", model),
         model_label = factor(model_label, levels = c("RSF\n(clinical)", "RSF\n(precision)",
                                                      "CLivD\n(non-lab)", "CLivD\n(lab)", "FIB-4")),
         time_label = factor(paste0("t = ", time, " years"), levels = c("t = 5 years", "t = 10 years")),
         lower = mean_auc - sd_auc,
         upper = mean_auc + sd_auc) %>%
  ggplot(aes(x = model_label, y = mean_auc)) +
  geom_point(size = 2.5, color = "steelblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                color = "steelblue", width = 0.3) +
  geom_hline(yintercept = 0.50, color = "darkorange",
             linewidth = 0.6, linetype = "dashed") +
  facet_wrap(~ time_label) +
  scale_y_continuous(limits = c(0.5, 0.85),
                     labels = scales::number_format(accuracy = 0.01)) +
  labs(x = "Model", y = "Time-dependent AUC",
       title = "Model comparison: RSF vs. CLivD vs. FIB-4",
       subtitle = "Mean \u00b1 SD across 5 CV folds") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey90"))

ggsave(file.path(plot_dir, "BFA01.1.2A_AUC_comparison.png"),
       p_auc, width = 8, height = 4.5, dpi = 400, bg = "white")

###############################################################################
# SECTION 2: Extended clinical RSF — C-index / AUC comparison
# Note: AIC is not applicable to RSF (non-likelihood-based method)
###############################################################################
cat("\n=== Section 2: Extended RSF model ===\n")

model_data_ext <- adata %>%
  filter(included_in_cohort == TRUE) %>%
  mutate(
    bmi                = winsorise(bmi),
    waist_circ         = winsorise(waist_circ),
    triglycerides      = winsorise(triglycerides),
    hdl                = winsorise(hdl),
    trig_hdl_ratio     = triglycerides / hdl,
    alcohol_grams_week = pmin(alcohol_grams_week, 500),
    hba1c              = winsorise(pmin(hba1c, 70)),
    smoking_binary     = as.factor(case_when(
      smoking %in% c("Never", "Previous") ~ 0L,
      smoking == "Current"                ~ 1L
    )),
    has_t2dm         = as.factor(as.integer(has_t2dm)),
    has_hypertension = as.factor(as.integer(has_hypertension)),
    sex              = case_when(sex == "Male" ~ 1L, sex == "Female" ~ 2L),
    status           = case_when(
      event_malo == TRUE            ~ 1L,
      event_non_liver_death == TRUE ~ 2L,
      TRUE                         ~ 0L
    )
  ) %>%
  drop_na(age, sex, bmi, waist_circ, trig_hdl_ratio, alcohol_grams_week,
          smoking, time_to_event, status,
          has_hypertension, hdl, triglycerides, hba1c)

set.seed(1234)
malo_ext <- which(model_data_ext$status == 1)
rest_ext <- which(model_data_ext$status != 1)
model_data_ext$fold <- NA_integer_
model_data_ext$fold[malo_ext] <- sample(rep(1:n_folds, length.out = length(malo_ext)))
model_data_ext$fold[rest_ext] <- sample(rep(1:n_folds, length.out = length(rest_ext)))

fold_sec2 <- list()
for (fold in 1:n_folds) {
  cat("  Fold", fold, "\n")
  tr_e <- model_data_ext[model_data_ext$fold != fold, ]
  te_e <- model_data_ext[model_data_ext$fold == fold, ]
  
  fit_base <- rfsrc.fast(
    as.formula(paste("Surv(time_to_event, status) ~", paste(preds_clin, collapse = "+"))),
    data = tr_e, ntree = 300, nodesize = 50, nsplit = 5, forest = TRUE, seed = -1,
    save.memory = TRUE, cause = 1, splitrule = "logrankCR", importance = FALSE
  )
  fit_extd <- rfsrc.fast(
    as.formula(paste("Surv(time_to_event, status) ~", paste(preds_ext, collapse = "+"))),
    data = tr_e, ntree = 300, nodesize = 50, nsplit = 5, forest = TRUE, seed = -1,
    save.memory = TRUE, cause = 1, splitrule = "logrankCR", importance = FALSE
  )
  
  pr_b <- predict(fit_base, newdata = te_e)
  pr_e <- predict(fit_extd, newdata = te_e)
  cif_b5  <- extract_cif(pr_b, 5);  cif_b10 <- extract_cif(pr_b, 10)
  cif_e5  <- extract_cif(pr_e, 5);  cif_e10 <- extract_cif(pr_e, 10)
  
  ts_e <- Surv(te_e$time_to_event, te_e$status == 1)
  roc_b5  <- timeROC(T = te_e$time_to_event, delta = te_e$status,
                     marker = cif_b5,  cause = 1, times = 5,  iid = FALSE)
  roc_b10 <- timeROC(T = te_e$time_to_event, delta = te_e$status,
                     marker = cif_b10, cause = 1, times = 10, iid = FALSE)
  roc_e5  <- timeROC(T = te_e$time_to_event, delta = te_e$status,
                     marker = cif_e5,  cause = 1, times = 5,  iid = FALSE)
  roc_e10 <- timeROC(T = te_e$time_to_event, delta = te_e$status,
                     marker = cif_e10, cause = 1, times = 10, iid = FALSE)
  
  fold_sec2[[fold]] <- tibble(
    fold   = fold,
    model  = c("RSF clinical (base)", "RSF clinical (extended)"),
    cindex = c(concordance(ts_e ~ cif_b10, reverse = TRUE)$concordance,
               concordance(ts_e ~ cif_e10, reverse = TRUE)$concordance),
    auc_t5  = c(roc_b5$AUC_1[2],  roc_e5$AUC_1[2]),
    auc_t10 = c(roc_b10$AUC_1[2], roc_e10$AUC_1[2])
  )
  rm(fit_base, fit_extd, pr_b, pr_e)
  gc()
}

sec2_summary <- bind_rows(fold_sec2) %>%
  group_by(model) %>%
  summarise(
    mean_cindex = mean(cindex), sd_cindex = sd(cindex),
    mean_auc_t5 = mean(auc_t5), sd_auc_t5 = sd(auc_t5),
    mean_auc_t10 = mean(auc_t10), sd_auc_t10 = sd(auc_t10),
    .groups = "drop"
  )

sec2_delta <- sec2_summary %>%
  dplyr::select(model, mean_cindex, mean_auc_t5, mean_auc_t10) %>%
  pivot_wider(names_from = model,
              values_from = c(mean_cindex, mean_auc_t5, mean_auc_t10)) %>%
  transmute(
    delta_cindex = `mean_cindex_RSF clinical (extended)` - `mean_cindex_RSF clinical (base)`,
    delta_auc_t5 = `mean_auc_t5_RSF clinical (extended)` - `mean_auc_t5_RSF clinical (base)`,
    delta_auc_t10= `mean_auc_t10_RSF clinical (extended)`- `mean_auc_t10_RSF clinical (base)`
  )

cat("Section 2: AIC is not applicable to RSF (tree-based, non-likelihood method).\n")
print(sec2_summary)
print(sec2_delta)

write_xlsx(
  list(Extended_vs_Base = sec2_summary, Delta = sec2_delta),
  file.path(result_dir, "BFA01.1.2_sec2_extended_model.xlsx")
)

###############################################################################
# SECTION 3: Nonlinearity — LRT + RSF PDPs vs CS-Cox splines
###############################################################################
cat("\n=== Section 3: Nonlinearity ===\n")

model_cs <- model_data %>%
  mutate(ev1 = as.integer(status == 1))

# ── 3a. LRT: linear vs natural spline (df=3) for each continuous predictor ──
cont_preds_test <- list(
  age                = "age",
  bmi                = "bmi",
  alcohol_grams_week = "alcohol_grams_week"
)
other_preds <- c("sex", "has_t2dm", "smoking_binary")

lrt_rows <- lapply(names(cont_preds_test), function(vname) {
  v       <- cont_preds_test[[vname]]
  others  <- setdiff(preds_clin, c(v, other_preds))
  all_lin <- c(others, other_preds, v)
  all_ns  <- c(others, other_preds, paste0("ns(", v, ", df = 3)"))
  
  f_lin <- as.formula(paste("Surv(time_to_event, ev1) ~", paste(all_lin, collapse = "+")))
  f_ns  <- as.formula(paste("Surv(time_to_event, ev1) ~", paste(all_ns,  collapse = "+")))
  
  m_lin <- coxph(f_lin, data = model_cs)
  m_ns  <- coxph(f_ns,  data = model_cs)
  lrt   <- anova(m_lin, m_ns)
  
  tibble(
    predictor = vname,
    chi2      = round(lrt[["Chisq"]][2], 2),
    df        = lrt[["Df"]][2],
    p_value   = lrt[["Pr(>|Chi|)"]][2],
    p_adj_BH  = NA_real_
  )
})
lrt_table <- bind_rows(lrt_rows) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

# ── 3b. LRT: key interaction terms ──
cox_main    <- coxph(Surv(time_to_event, ev1) ~ age + sex + bmi + has_t2dm +
                       alcohol_grams_week + smoking_binary, data = model_cs)
cox_int_AxS <- coxph(Surv(time_to_event, ev1) ~ age + sex + bmi + has_t2dm +
                       alcohol_grams_week + smoking_binary +
                       alcohol_grams_week:sex, data = model_cs)
cox_int_BxD <- coxph(Surv(time_to_event, ev1) ~ age + sex + bmi + has_t2dm +
                       alcohol_grams_week + smoking_binary +
                       bmi:has_t2dm, data = model_cs)
cox_int_AxD <- coxph(Surv(time_to_event, ev1) ~ age + sex + bmi + has_t2dm +
                       alcohol_grams_week + smoking_binary +
                       age:sex, data = model_cs)

lrt_int <- bind_rows(
  { lrt <- anova(cox_main, cox_int_AxS)
  tibble(interaction = "Alcohol \u00d7 Sex", chi2 = round(lrt$Chisq[2], 2),
         df = lrt$Df[2], p_value = lrt$`Pr(>|Chi|)`[2]) },
  { lrt <- anova(cox_main, cox_int_BxD)
  tibble(interaction = "BMI \u00d7 T2DM",   chi2 = round(lrt$Chisq[2], 2),
         df = lrt$Df[2], p_value = lrt$`Pr(>|Chi|)`[2]) },
  { lrt <- anova(cox_main, cox_int_AxD)
  tibble(interaction = "Age \u00d7 Sex",     chi2 = round(lrt$Chisq[2], 2),
         df = lrt$Df[2], p_value = lrt$`Pr(>|Chi|)`[2]) }
) %>% mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

write_xlsx(
  list(Nonlinearity_LRT = lrt_table, Interaction_LRT = lrt_int),
  file.path(result_dir, "BFA01.1.2_sec3_nonlinearity_LRT.xlsx")
)
cat("Nonlinearity LRT table:\n"); print(lrt_table)
cat("Interaction LRT table:\n");  print(lrt_int)

# ── 3c. CS-Cox with natural splines — partial effects ──
cox_ns_full <- coxph(
  Surv(time_to_event, ev1) ~ ns(age, df = 3) + sex + ns(bmi, df = 3) + has_t2dm +
    ns(alcohol_grams_week, df = 3) + smoking_binary,
  data = model_cs
)

spline_partial <- function(cox_fit, data, var, grid_n = 60) {
  grid <- seq(quantile(data[[var]], 0.02, na.rm = TRUE),
              quantile(data[[var]], 0.98, na.rm = TRUE), length.out = grid_n)
  ref <- as.data.frame(lapply(names(data), function(cn) {
    x <- data[[cn]]
    if (is.factor(x))  return(factor(levels(x)[1], levels = levels(x)))
    if (is.numeric(x)) return(mean(x, na.rm = TRUE))
    return(x[1])
  }))
  names(ref) <- names(data)
  newdat <- ref[rep(1, grid_n), ]
  newdat[[var]] <- grid
  lp <- predict(cox_fit, newdata = newdat, type = "lp")
  tibble(x = grid, lp_centered = lp - mean(lp), variable = var)
}

spline_data <- bind_rows(lapply(
  c("age", "bmi", "alcohol_grams_week"),
  spline_partial, cox_fit = cox_ns_full, data = model_cs
))

# ── 3d. RSF partial dependence plots ──
rsf_full_path <- file.path(result_dir, "RSF_final_clinical_model.rds")
if (file.exists(rsf_full_path)) {
  cat("  Loading saved RSF clinical model for PDPs...\n")
  final_rsf_clin <- readRDS(rsf_full_path)
} else {
  cat("  Full model not found; will be fitted in Section 7.\n")
  final_rsf_clin <- NULL
}

rsf_pdp <- function(rsf_fit, data, var, t = 10, grid_n = 30, n_sub = 1000) {
  set.seed(1234)
  sub  <- data[sample(nrow(data), min(n_sub, nrow(data))), ]
  grid <- seq(quantile(data[[var]], 0.02, na.rm = TRUE),
              quantile(data[[var]], 0.98, na.rm = TRUE), length.out = grid_n)
  ti   <- which.min(abs(rsf_fit$time.interest - t))
  cif_vals <- sapply(grid, function(val) {
    sub[[var]] <- val
    mean(predict(rsf_fit, newdata = sub)$cif[, ti, 1])
  })
  tibble(x = grid, cif_centered = cif_vals - mean(cif_vals), variable = var)
}

var_labels <- c(age = "Age (years)",
                bmi = "BMI (kg/m\u00b2)",
                alcohol_grams_week = "Alcohol (g/week)")

p_spline <- spline_data %>%
  mutate(variable = var_labels[variable]) %>%
  ggplot(aes(x = x, y = lp_centered)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(colour = "#D62728", linewidth = 0.8) +
  facet_wrap(~ variable, scales = "free_x", nrow = 1) +
  labs(x = NULL, y = "Partial log-hazard\n(centered)",
       title = "Cause-specific Cox with natural splines") +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey90"))

if (!is.null(final_rsf_clin)) {
  cat("  Computing RSF PDPs (might take a while)...\n")
  pdp_data <- bind_rows(lapply(
    c("age", "bmi", "alcohol_grams_week"),
    rsf_pdp, rsf_fit = final_rsf_clin, data = model_data
  ))
  
  p_pdp <- pdp_data %>%
    mutate(variable = var_labels[variable]) %>%
    ggplot(aes(x = x, y = cif_centered)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_line(colour = "#1F77B4", linewidth = 0.8) +
    facet_wrap(~ variable, scales = "free_x", nrow = 1) +
    labs(x = NULL, y = "Partial CIF at t=10\n(centered)",
         title = "RSF partial dependence plots") +
    theme_bw(base_size = 10) +
    theme(strip.background = element_rect(fill = "grey90"))
  
  p_nonlin <- p_spline / p_pdp +
    plot_annotation(title = "Nonlinear predictor effects: CS-Cox splines vs RSF",
                    tag_levels = "A")
} else {
  p_nonlin <- p_spline +
    plot_annotation(title = "CS-Cox spline effects (RSF PDPs computed after Section 7)")
}

# print(p_nonlin)
# ggsave(file.path(plot_dir, "BFA01.1.2B_nonlinearity.png"),
#        p_nonlin, width = 9, height = 5, dpi = 400, bg = "white")

###############################################################################
# SECTION 4: Risk-stratified CIF and KM curves (predicted tertiles)
###############################################################################
cat("\n=== Section 4: Risk stratification ===\n")

all_test <- all_test %>%
  mutate(risk_tertile = factor(ntile(rsf_clin_10, 3),
                               labels = c("Low", "Medium", "High")))

# ── 4a. Aalen-Johansen CIF by tertile (cause 1 = MALO) ──
ci_fit <- cmprsk::cuminc(
  ftime   = all_test$time_to_event,
  fstatus = all_test$status,
  group   = all_test$risk_tertile
)

ci_df <- bind_rows(lapply(names(ci_fit), function(nm) {
  if (!is.list(ci_fit[[nm]])) return(NULL)
  if (endsWith(nm, " 1")) {
    tibble(time  = ci_fit[[nm]]$time,
           cif   = ci_fit[[nm]]$est,
           group = sub(" 1$", "", nm))
  }
})) %>%
  mutate(group = factor(group, levels = c("Low", "Medium", "High")))

gray_p     <- ci_fit$Tests["1", "pv"]
gray_label <- paste0("Gray's p ", if (gray_p < 0.001) "< 0.001" else
  paste0("= ", formatC(gray_p, digits = 3, format = "f")))

p_cif <- ggplot(ci_df, aes(x = time, y = cif, colour = group)) +
  geom_step(linewidth = 0.8) +
  geom_label(data = tibble(x = 15, y = 0.0001),
             aes(x = x, y = y, label = gray_label),
             colour = "black", size = 3, inherit.aes = FALSE) +
  scale_color_d3() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(x = "Time (years)", y = "Cumulative incidence (MALO)",
       colour = "Risk tertile",
       title = "Aalen-Johansen CIF by predicted risk tertile",
       subtitle = "Tertiles based on RSF clinical model CIF at t = 10 years") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")

# ── 4b. Kaplan-Meier curves ──
km_fit  <- survfit(Surv(time_to_event, status == 1) ~ risk_tertile, data = all_test)
lr_test <- survdiff(Surv(time_to_event, status == 1) ~ risk_tertile, data = all_test)
lr_p    <- 1 - pchisq(lr_test$chisq, df = length(lr_test$obs) - 1)
lr_label <- paste0("Log-rank p ", if (lr_p < 0.001) "< 0.001" else
  paste0("= ", formatC(lr_p, digits = 3, format = "f")))

km_df <- tidy(km_fit) %>%
  mutate(strata = gsub("risk_tertile=", "", strata),
         strata = factor(strata, levels = c("Low", "Medium", "High")))

p_km <- ggplot(km_df, aes(x = time, y = estimate, colour = strata)) +
  geom_step(linewidth = 0.8) +
  geom_label(data = tibble(x = 3, y = 0.985),
             aes(x = x, y = y, label = lr_label),
             colour = "black", size = 3, inherit.aes = FALSE) +
  scale_color_d3() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(x = "Time (years)", y = "Kaplan-Meier survival (MALO)",
       colour = "Risk tertile",
       title = "Kaplan-Meier curves by predicted risk tertile",
       subtitle = "Note: competing events treated as censored") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")

p_risk_strat <- p_cif + p_km + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(plot_dir, "BFA01.1.2C_risk_stratification.png"),
       p_risk_strat, width = 9, height = 4.5, dpi = 400, bg = "white")

counts_tertile <- all_test %>%
  group_by(risk_tertile) %>%
  summarise(n = n(), n_malo = sum(status == 1),
            n_comp_death = sum(status == 2), n_censored = sum(status == 0))
write_xlsx(list(Risk_tertile_counts = counts_tertile),
           file.path(result_dir, "BFA01.1.2_sec4_risk_tertile_counts.xlsx"))

###############################################################################
# SECTION 5: Decision curve analysis
# Note: DCA net benefit is small in absolute terms because MALO is rare in this
# general-population cohort (< ~1% 5-year event rate). "Treat All" dominates
# at very low thresholds because the FP penalty is negligible there. The
# clinically relevant window is ~0.5–3% for t=5 and ~1–5% for t=10 — zoom
# in on that range to see meaningful separation from "Treat None".
###############################################################################
cat("\n=== Section 5: Decision curve analysis ===\n")

thresholds <- list("5" = seq(0, 0.03, by = 0.001), "10" = seq(0, 0.03, by = 0.001))

for (t in c(5, 10)) {
  dca_dat <- all_test %>%
    transmute(
      time_to_event,
      status,
      rsf_clinical = .data[[paste0("rsf_clin_", t)]],
      fg_clinical  = .data[[paste0("fg_clin_",  t)]],
      alcohol_only = .data[[paste0("alcohol_grams_week_sv_", t)]],
      bmi_only     = .data[[paste0("bmi_sv_",   t)]],
      t2dm_only    = .data[[paste0("t2dm_sv_",  t)]]
    ) %>%
    dplyr::rename("RSF (clinical)" = "rsf_clinical",
                  "Fine-Gray (clinical)" = "fg_clinical",
                  "Alcohol (g/week)" = "alcohol_only",
                  "BMI" = "bmi_only", "Type 2 diabetes mellitus" = "t2dm_only")
  
  thres <- thresholds[[as.character(t)]]
  
  dca_A <- dca(Surv(time_to_event, status == 1) ~ `RSF (clinical)`,
               data = dca_dat, time = t, thresholds = thres)
  p_dca_A <- dca_A %>%
    plot(smooth = TRUE) +
    # coord_cartesian(xlim = xlim) +
    scale_color_d3() +
    labs(title    = paste0("DCA at t = ", t, " years: RSF net benefit"),
         x = "Threshold probability", y = "Net benefit") +
    theme_bw(base_size = 8) +
    theme(legend.position = "bottom")
  
  dca_B <- dca(Surv(time_to_event, status == 1) ~ `RSF (clinical)` + `Alcohol (g/week)` + BMI + `Type 2 diabetes mellitus`,
               data = dca_dat, time = t, thresholds = thres)
  p_dca_B <- dca_B %>%
    plot(smooth = TRUE) +
    viridis::scale_color_viridis(discrete = TRUE) +
    labs(title    = paste0("DCA at t = ", t, " years: RSF vs single predictors"),
         x = "Threshold probability", y = "Net benefit") +
    theme_bw(base_size = 8) +
    theme(legend.position = "bottom")
  
  p_dca_combined <- wrap_plots(p_dca_A, p_dca_B, ncol = 2) +
    plot_annotation(title = paste0("Decision Curve Analysis at t = ", t, " years"))
  
  ggsave(file.path(plot_dir, paste0("BFA01.1.2D_DCA_t", t, ".png")),
         p_dca_combined, width = 10, height = 4.5, dpi = 400, bg = "white")
  
  write_xlsx(
    list(Panel_A = as_tibble(dca_A) %>% mutate(eval_time = t),
         Panel_B = as_tibble(dca_B) %>% mutate(eval_time = t)),
    file.path(result_dir, paste0("BFA01.1.2_sec5_DCA_netbenefit_t", t, ".xlsx"))
  )
}

###############################################################################
# SECTION 6: Calibration plots — RSF clinical and precision
###############################################################################
cat("\n=== Section 6: Calibration plots ===\n")

compute_aj_cif <- function(cal_list, model_label) {
  bind_rows(cal_list) %>%
    group_by(eval_time, bin) %>%
    summarise(
      mean_pred = mean(pred),
      obs_cif   = {
        sub  <- pick(everything())
        t_pt <- eval_time[1]
        if (sum(sub$status == 1) == 0) {
          0
        } else {
          fit  <- prodlim(Hist(time, status) ~ 1, data = sub)
          as.numeric(predict(fit, times = t_pt, cause = 1))
        }
      },
      .groups = "drop"
    ) %>%
    mutate(
      model      = model_label,
      time_label = factor(paste0("t = ", eval_time, " years"),
                          levels = paste0("t = ", eval_times, " years"))
    )
}

cal_clin_df <- compute_aj_cif(cal_clin_list, "RSF (clinical)")
cal_prec_df <- compute_aj_cif(cal_prec_list, "RSF (precision)")
cal_all_df  <- bind_rows(cal_clin_df, cal_prec_df)

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

ggsave(file.path(plot_dir, "BFA01.1.2E_RSF_calibration.png"),
       p_cal, width = 8, height = 6, dpi = 400, bg = "white")

###############################################################################
# SECTION 7: Full-dataset model training and saving
###############################################################################
cat("\n=== Section 7: Full-dataset model training ===\n")

# ── Clinical model ──
rsf_clin_path <- file.path(result_dir, "RSF_final_clinical_model.rds")
cat("  Fitting full RSF clinical model...\n")
final_rsf_clin <- rfsrc.fast(
  as.formula(paste("Surv(time_to_event, status) ~", paste(preds_clin, collapse = "+"))),
  data = model_data, ntree = 300, nodesize = 50, nsplit = 5,
  save.memory = FALSE, forest = TRUE, cause = 1, splitrule = "logrankCR",
  importance = TRUE, seed = -1, do.trace = 60
)
saveRDS(final_rsf_clin, rsf_clin_path)
cat("  Saved to:", rsf_clin_path, "\n")

# ── Precision model (waist_circ) ──
rsf_prec_path <- file.path(result_dir, "RSF_final_precision_model.rds")
cat("  Fitting full RSF precision model (waist_circ)...\n")
final_rsf_prec <- rfsrc.fast(
  as.formula(paste("Surv(time_to_event, status) ~", paste(preds_prec, collapse = "+"))),
  data = model_data, ntree = 300, nodesize = 50, nsplit = 5,
  save.memory = FALSE, forest = TRUE, cause = 1, splitrule = "logrankCR",
  importance = TRUE, seed = -1, do.trace = 60
)
saveRDS(final_rsf_prec, rsf_prec_path)
cat("  Saved to:", rsf_prec_path, "\n")

# Now recompute PDPs with the freshly trained model and overwrite Section 3 figure
cat("  Recomputing RSF PDPs with full model...\n")
pdp_data <- bind_rows(lapply(
  c("age", "bmi", "alcohol_grams_week"),
  rsf_pdp, rsf_fit = final_rsf_clin, data = model_data
))
p_pdp <- pdp_data %>%
  mutate(variable = var_labels[variable]) %>%
  ggplot(aes(x = x, y = cif_centered)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(colour = "#1F77B4", linewidth = 0.8) +
  facet_wrap(~ variable, scales = "free_x", nrow = 1) +
  labs(x = NULL, y = "Partial CIF at t=10\n(centered)",
       title = "RSF partial dependence plots") +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey90"))

p_nonlin_final <- p_spline / p_pdp +
  plot_annotation(title = "Nonlinear predictor effects: CS-Cox splines vs RSF",
                  tag_levels = "A")
ggsave(file.path(plot_dir, "BFA01.1.2B_nonlinearity.png"),
       p_nonlin_final, width = 9, height = 5, dpi = 400, bg = "white")

###############################################################################
# SECTION 8: Variable importance — VIMP + SHAP
###############################################################################
cat("\n=== Section 8: Variable importance ===\n")

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
  
  cat("  Computing SHAP for", model_label, "(", length(idx), "rows, nsim = 50)...\n")
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
  
  list(VIMP_plot       = p_vimp,
       beeswarm_plot   = p_beeswarm,
       bar_plot        = p_barplot,
       dependence_plot = p_dependence)
}

vimp_clin <- variable_importance_analysis(final_rsf_clin, preds_clin, "RSF clinical")
vimp_prec <- variable_importance_analysis(final_rsf_prec, preds_prec, "RSF precision")

# Save plots
plot_specs <- list(
  list(plots = vimp_clin, tag = "clinical"),
  list(plots = vimp_prec, tag = "precision")
)
for (spec in plot_specs) {
  for (nm in names(spec$plots)) {
    fname <- sprintf("BFA01.1.5F_RSF_%s_%s.png", spec$tag, nm)
    ggsave(file.path(plot_dir, fname),
           plot = spec$plots[[nm]],
           width = 6, height = 5, dpi = 400, bg = "white")
  }
}

cat("\n=== All analyses complete ===\n")
cat("Results:", result_dir, "\n")
cat("Plots:  ", plot_dir,   "\n")