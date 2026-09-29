###############################################################################
# Sex x alcohol interaction: sensitivity analysis
# Ronnie Li
#
# Prompted by a stress test on the deployed calculator: at the same age, BMI,
# diabetes and smoking status, a woman drinking 21 drinks/week (294 g/week)
# gets a lower calculated risk than a man at the same intake, because the
# fitted sex term (SHR 0.63, female vs. male) is large relative to the
# alcohol term. The question raised was whether this reflects a missing
# sex*alcohol interaction (alcohol raising risk more steeply per gram in
# women) rather than a genuine, uniform sex difference.
#
# This refits the full cohort (same model_data construction, same clinical
# predictors, as 02_FineGray_main_analysis.R) with a sex*alcohol interaction
# added, in three specifications:
#   1. Cause-specific Cox model, continuous alcohol (per 100 g/week) x female
#   2. Fine-Gray (crr), continuous alcohol (per 100 g/week) x female --
#      the same model family as the deployed calculator
#   3. Cause-specific Cox model, a sex-specific "excess alcohol" indicator
#      (>210 g/week men, >140 g/week women -- the app's banner thresholds)
#      x female, as a categorical check that doesn't assume linearity
#
# CONCLUSION (see results/07_sex_alcohol_interaction.xlsx for the numbers):
# none of the three specifications supports alcohol raising risk more
# steeply in women; if anything the point estimates go the other way, and
# all are compatible with no interaction. The deployed model and its
# recommendation are UNCHANGED by this analysis -- see app.R, which instead
# adds a non-overriding alcohol-assessment banner.
#
# Outputs: tables  -> results/07_sex_alcohol_interaction.xlsx
#          figures -> results/plots/07.01_sex_alcohol_interaction_forest.png
###############################################################################

library(survival)
library(cmprsk)
library(tidyverse)
library(arrow)
library(writexl)

# ── Directories ───────────────────────────────────────────────────────────────
project_dir <- "/mnt/d/Projects/BFA"
data_dir    <- file.path(project_dir, "data")
result_dir  <- file.path(project_dir, "results")
plot_dir    <- file.path(result_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

tables <- list()

# ── Data preparation -- identical to 02_FineGray_main_analysis.R, so this
# sensitivity analysis uses exactly the same modeling population (N, exclusions,
# preprocessing) as the deployed model ──────────────────────────────────────────
winsorise <- function(x, probs = c(0.01, 0.99)) {
  q <- quantile(x, probs, na.rm = TRUE)
  pmax(pmin(x, q[2]), q[1])
}

adata <- read_feather(file.path(data_dir, "BFA_principal_data.feather"))

model_data <- adata %>%
  filter(included_in_cohort == TRUE) %>%
  mutate(
    bmi                = winsorise(bmi),
    waist_circ         = winsorise(waist_circ),
    triglycerides      = winsorise(triglycerides),
    hdl                = winsorise(hdl),
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

preds_clin <- c("age", "sex", "bmi", "has_t2dm", "alcohol_grams_week", "smoking_binary")
stopifnot(all(complete.cases(model_data[, preds_clin])))

cat(sprintf("N = %s, MALO events = %s (matches 02_FineGray_main_analysis.R's model_data)\n",
            format(nrow(model_data), big.mark = ","), format(sum(model_data$status == 1), big.mark = ",")))

# Working copy with convenience numeric columns for the interaction models
d <- model_data %>%
  mutate(
    sex_female = as.integer(as.character(sex) == "2"),
    smk        = as.integer(as.character(smoking_binary) == "1"),
    t2         = as.integer(as.character(has_t2dm) == "1"),
    alc100     = alcohol_grams_week / 100,   # per 100 g/week, ~7 drinks
    excess_alcohol = as.integer(
      (as.character(sex) == "1" & alcohol_grams_week > 210) |   # male
      (as.character(sex) == "2" & alcohol_grams_week > 140)     # female
    )   # same thresholds as app.R's ALCOHOL_ASSESSMENT_CUT
  )

tables[["01_cohort_summary"]] <- tibble(
  metric = c("N", "N MALO events", "N female", "N male",
             "N female excess alcohol (>140 g/wk)", "N male excess alcohol (>210 g/wk)"),
  value  = c(nrow(d), sum(d$status == 1), sum(d$sex_female == 1), sum(d$sex_female == 0),
             sum(d$sex_female == 1 & d$excess_alcohol == 1),
             sum(d$sex_female == 0 & d$excess_alcohol == 1))
)
print(tables[["01_cohort_summary"]], n = Inf)

###############################################################################
# 1. Cause-specific Cox: continuous alcohol x female
###############################################################################
cat("\n=== Cause-specific Cox: alcohol (per 100 g/wk) x female ===\n")

cox0 <- coxph(Surv(time_to_event, status == 1) ~ age + sex_female + bmi + t2 + alc100 + smk, data = d)
cox1 <- update(cox0, . ~ . + sex_female:alc100)

cox1_sum <- summary(cox1)$coefficients
tbl_cox_interaction <- as_tibble(cox1_sum, rownames = "term") %>%
  rename(HR_log = coef, HR = `exp(coef)`, se = `se(coef)`, z = z, p_value = `Pr(>|z|)`) %>%
  mutate(HR_lower = exp(HR_log - 1.96 * se), HR_upper = exp(HR_log + 1.96 * se))
lrt_cox <- anova(cox0, cox1)
cox_lrt_p <- lrt_cox$`Pr(>|Chi|)`[2]

print(tbl_cox_interaction, width = Inf)
cat("LRT p (interaction term):", cox_lrt_p, "\n")

tables[["02_coxph_continuous_interact"]] <- tbl_cox_interaction
tables[["03_coxph_continuous_LRT"]] <- tibble(
  model = c("Main effects", "+ sex*alcohol interaction"),
  loglik = c(cox0$loglik[2], cox1$loglik[2]),
  df     = c(length(coef(cox0)), length(coef(cox1))),
  lrt_p  = c(NA, cox_lrt_p)
)

###############################################################################
# 2. Fine-Gray (crr): continuous alcohol x female -- same model family as the
#    deployed calculator
###############################################################################
cat("\n=== Fine-Gray (crr): alcohol (per 100 g/wk) x female ===\n")

X0 <- model.matrix(~ age + sex_female + bmi + t2 + alc100 + smk, d)[, -1]
X1 <- cbind(X0, `sex_female:alc100` = d$sex_female * d$alc100)

t0 <- Sys.time()
fg0 <- crr(d$time_to_event, d$status, cov1 = X0, failcode = 1, cencode = 0)
fg1 <- crr(d$time_to_event, d$status, cov1 = X1, failcode = 1, cencode = 0)
cat("  crr fit time:", format(Sys.time() - t0), "\n")

fg1_sum <- summary(fg1)$coef
tbl_fg_interaction <- as_tibble(fg1_sum, rownames = "term") %>%
  rename(coef = coef, SHR = `exp(coef)`, se = `se(coef)`, z = z, p_value = `p-value`) %>%
  mutate(SHR_lower = exp(coef - 1.96 * se), SHR_upper = exp(coef + 1.96 * se))
fg_pseudo_lrt_p <- pchisq(2 * (fg1$loglik - fg0$loglik), df = 1, lower.tail = FALSE)

print(tbl_fg_interaction, width = Inf)
cat("Pseudo-LRT p (interaction term):", fg_pseudo_lrt_p, "\n")

tables[["04_FineGray_continuous_interact"]] <- tbl_fg_interaction
tables[["05_FineGray_continuous_LRT"]] <- tibble(
  model = c("Main effects", "+ sex*alcohol interaction"),
  loglik = c(fg0$loglik, fg1$loglik),
  df     = c(length(fg0$coef), length(fg1$coef)),
  pseudo_lrt_p = c(NA, fg_pseudo_lrt_p)
)

###############################################################################
# 3. Cause-specific Cox: categorical excess-alcohol indicator x female
#    (app.R's own sex-specific thresholds, so this doesn't assume linearity)
###############################################################################
cat("\n=== Cause-specific Cox: excess-alcohol indicator x female ===\n")

cox2 <- update(cox0, . ~ . + excess_alcohol)
cox3 <- update(cox2, . ~ . + excess_alcohol:sex_female)

cox3_sum <- summary(cox3)$coefficients
tbl_cox_categorical <- as_tibble(cox3_sum, rownames = "term") %>%
  rename(HR_log = coef, HR = `exp(coef)`, se = `se(coef)`, z = z, p_value = `Pr(>|z|)`) %>%
  mutate(HR_lower = exp(HR_log - 1.96 * se), HR_upper = exp(HR_log + 1.96 * se))
lrt_cox_cat <- anova(cox2, cox3)
cox_cat_lrt_p <- lrt_cox_cat$`Pr(>|Chi|)`[2]

print(tbl_cox_categorical, width = Inf)
cat("LRT p (interaction term):", cox_cat_lrt_p, "\n")

tables[["06_coxph_categorical_interaction"]] <- tbl_cox_categorical
tables[["07_coxph_categorical_LRT"]] <- tibble(
  model = c("+ excess-alcohol indicator", "+ indicator*sex interaction"),
  loglik = c(cox2$loglik[2], cox3$loglik[2]),
  df     = c(length(coef(cox2)), length(coef(cox3))),
  lrt_p  = c(NA, cox_cat_lrt_p)
)

###############################################################################
# 4. The deployed model is unchanged: predicted risk for the stress-test case
#    (age 55, BMI 27, no T2DM, no smoking, 294 g/week = 21 drinks/week) under
#    the deployed model (main effects only) vs. what the interaction model
#    would have predicted -- to show how little the interaction moves it
###############################################################################
cat("\n=== Deployed model vs. interaction model: stress-test case ===\n")

deploy <- readRDS(file.path(result_dir, "FG_clinical_model_deploy.rds"))
case <- data.frame(age = 55, bmi = 27, has_t2dm = "0", alcohol_grams_week = 294, smoking_binary = "0")

deployed_risk <- function(sex_code) {
  nd <- cbind(sex = sex_code, case)
  100 * deploy$predict_fn(deploy, nd, times = 10)[, "risk_10y"]
}

# Fine-Gray + interaction model's 10y risk needs its own baseline hazard;
# recompute H0(10) from fg1's residuals the same way 02_FineGray_main_analysis.R
# does for the deployed model (sum of increments up to t = 10).
H0_10_fg1 <- sum(fg1$bfitj[fg1$uftime <= 10])
lp_fg1 <- function(sex_female, alc_g) {
  x <- c(age = 55, sex_female = sex_female, bmi = 27, t2 = 0, alc100 = alc_g / 100, smk = 0,
         `sex_female:alc100` = sex_female * alc_g / 100)
  sum(x * fg1$coef)
}
risk_fg1 <- function(sex_female, alc_g) 100 * (1 - exp(-H0_10_fg1 * exp(lp_fg1(sex_female, alc_g))))

tbl_stress_test <- tibble(
  sex = c("Male", "Female"),
  deployed_model_risk_10y_pct     = c(deployed_risk("1"), deployed_risk("2")),
  with_interaction_risk_10y_pct   = c(risk_fg1(0, 294), risk_fg1(1, 294))
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

print(tbl_stress_test, width = Inf)
tables[["08_stress_test_case_294g_week"]] <- tbl_stress_test

# ── Forest plot of the three interaction estimates ─────────────────────────────
forest_df <- bind_rows(
  tibble(spec = "Cox: alcohol (per 100g) x female",
         estimate = tbl_cox_interaction$HR[tbl_cox_interaction$term == "sex_female:alc100"],
         lower    = tbl_cox_interaction$HR_lower[tbl_cox_interaction$term == "sex_female:alc100"],
         upper    = tbl_cox_interaction$HR_upper[tbl_cox_interaction$term == "sex_female:alc100"],
         p_value  = tbl_cox_interaction$p_value[tbl_cox_interaction$term == "sex_female:alc100"]),
  tibble(spec = "Fine-Gray: alcohol (per 100g) x female",
         estimate = tbl_fg_interaction$SHR[tbl_fg_interaction$term == "sex_female:alc100"],
         lower    = tbl_fg_interaction$SHR_lower[tbl_fg_interaction$term == "sex_female:alc100"],
         upper    = tbl_fg_interaction$SHR_upper[tbl_fg_interaction$term == "sex_female:alc100"],
         p_value  = tbl_fg_interaction$p_value[tbl_fg_interaction$term == "sex_female:alc100"]),
  tibble(spec = "Cox: excess-alcohol indicator x female",
         # NB: base R names interaction terms by the order the main effects
         # first appear in the full formula, not the order written in the
         # interaction spec -- sex_female precedes excess_alcohol here.
         estimate = tbl_cox_categorical$HR[tbl_cox_categorical$term == "sex_female:excess_alcohol"],
         lower    = tbl_cox_categorical$HR_lower[tbl_cox_categorical$term == "sex_female:excess_alcohol"],
         upper    = tbl_cox_categorical$HR_upper[tbl_cox_categorical$term == "sex_female:excess_alcohol"],
         p_value  = tbl_cox_categorical$p_value[tbl_cox_categorical$term == "sex_female:excess_alcohol"])
) %>%
  mutate(spec = factor(spec, levels = rev(spec)))

tables[["09_interaction_forest_data"]] <- forest_df

p_forest <- ggplot(forest_df, aes(x = estimate, y = spec)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.2, colour = "grey30") +
  geom_point(size = 3, colour = "#2166AC") +
  scale_x_log10() +
  labs(
    x = "Interaction hazard/subdistribution-hazard ratio (95% CI, log scale)",
    y = NULL,
    title = "Sex x alcohol interaction: three specifications",
    subtitle = "Ratio > 1 = alcohol raises risk more steeply in women; all three CIs cross 1"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(plot_dir, "07.01_sex_alcohol_interaction_forest.png"),
       p_forest, width = 10, height = 4, dpi = 400, bg = "white")

# ── Write output ────────────────────────────────────────────────────────────
out_path <- file.path(result_dir, "07_sex_alcohol_interaction.xlsx")
write_xlsx(tables, out_path)
cat("\nSaved:", out_path, "\n")
