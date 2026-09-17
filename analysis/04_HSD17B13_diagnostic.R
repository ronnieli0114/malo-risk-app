# =============================================================================
# HSD17B13 diagnostic: check why the HSD17B13 variant ALT/ALT genotype
# seems deleterious rather than protective.
# =============================================================================

library(tidyverse)
library(survival)
library(prodlim)
library(patchwork)
library(ggstance)
library(arrow)
library(ggsci)
library(writexl)
library(car)

# -----------------------------------------------------------------
# Directories
# -----------------------------------------------------------------
project_dir <- "/mnt/d/Projects/BFA"
data_dir    <- file.path(project_dir, "data")
outdir      <- file.path(project_dir, "plots")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# -----------------------------------------------------------------
# Read and prepare data
# -----------------------------------------------------------------
adata <- read_feather(file.path(data_dir, "BFA_principal_data.feather"))

model_data <- adata %>%
  filter(included_in_cohort == TRUE) %>%
  mutate(
    alcohol_grams_week = pmin(alcohol_grams_week, 500),
    hba1c              = pmin(hba1c, 70),
    waist_hip_ratio    = waist_circ / hip_circ,
    smoking_binary     = case_when(
      smoking == "Never" | smoking == "Previous" ~ 0,
      smoking == "Current"                       ~ 1
    ),
    smoking_binary   = as.factor(smoking_binary),
    trig_hdl_ratio   = triglycerides / hdl,
    alcohol_category = case_when(
      sex == "Male"   & alcohol_grams_week < 210                              ~ 0,
      sex == "Female" & alcohol_grams_week < 140                              ~ 0,
      sex == "Male"   & alcohol_grams_week >= 210 & alcohol_grams_week <= 420 ~ 1,
      sex == "Female" & alcohol_grams_week >= 140 & alcohol_grams_week <= 350 ~ 1,
      sex == "Male"   & alcohol_grams_week > 420                              ~ 2,
      sex == "Female" & alcohol_grams_week > 350                              ~ 2,
      TRUE ~ NA_integer_
    ),
    has_t2dm         = as.factor(as.numeric(has_t2dm)),
    alcohol_cmrf     = as.factor(alcohol_category * n_cmrfs),
    sex              = case_when(sex == "Male" ~ 1, sex == "Female" ~ 2),
    smoking          = factor(smoking, levels = c("Never", "Previous", "Current")),
    alcohol_category = as.factor(alcohol_category),
    status           = case_when(
      event_malo            == TRUE ~ 1,
      event_non_liver_death == TRUE ~ 2,
      TRUE                          ~ 0
    )
  ) %>%
  drop_na(age, sex, bmi, waist_circ, trig_hdl_ratio, alcohol_grams_week,
          smoking, time_to_event, status)

predictors           <- c("age", "sex", "bmi", "has_t2dm", "alcohol_grams_week", "smoking_binary")
predictors_precision <- c(predictors, "waist_hip_ratio", "trig_hdl_ratio",
                          "PNPLA3_rs738409_G", "TM6SF2_rs58542926_T", "HSD17B13_rs9992651_A")
eval_times <- c(5,10)

# -----------------------------------------------------------------
# 1. Basic frequency and event rate by genotype
# -----------------------------------------------------------------
hsd_summary <- model_data %>%
  mutate(genotype = factor(HSD17B13_rs9992651_A,
                           levels = c(0, 1, 2),
                           labels = c("0 (ref/ref)", "1 (ref/alt)", "2 (alt/alt)"))) %>%
  group_by(genotype) %>%
  summarise(
    n          = n(),
    n_malo     = sum(status == 1),
    n_comp     = sum(status == 2),
    n_censored = sum(status == 0),
    event_rate = mean(status == 1),
    .groups    = "drop"
  )

cat("=== HSD17B13 genotype distribution and raw event rates ===\n")
print(hsd_summary)
write_xlsx(hsd_summary, file.path(outdir, "BFA01.1.4_HSD17B13_event_summary.xlsx"))

# -----------------------------------------------------------------
# 2. Cumulative incidence (Aalen-Johansen) by genotype
#    This is the proper competing-risks estimate
# -----------------------------------------------------------------
aj_data <- model_data %>%
  mutate(genotype = factor(HSD17B13_rs9992651_A,
                           levels = c(0, 1, 2),
                           labels = c("0 (ref/ref)", "1 (ref/alt)", "2 (alt/alt)")))

aj_fit <- prodlim(Hist(time_to_event, status) ~ genotype, data = aj_data)

# Extract CIF at eval_times for cause 1
cif_by_genotype <- lapply(eval_times, function(t) {
  pred <- predict(aj_fit, times = t, cause = 1,
                  newdata = data.frame(genotype = factor(
                    c("0 (ref/ref)", "1 (ref/alt)", "2 (alt/alt)"),
                    levels = c("0 (ref/ref)", "1 (ref/alt)", "2 (alt/alt)")
                  )))
  tibble(
    genotype = c("0 (ref/ref)", "1 (ref/alt)", "2 (alt/alt)"),
    cif      = as.numeric(pred),
    time     = t
  )
}) %>% bind_rows()

cat("\n=== Aalen-Johansen CIF by genotype ===\n")
print(cif_by_genotype)

# -----------------------------------------------------------------
# 3. Kaplan-Meier-style cumulative incidence curves (ggplot)
# -----------------------------------------------------------------

# Extract full AJ curve data manually
aj_plot_data <- lapply(c("0 (ref/ref)", "1 (ref/alt)", "2 (alt/alt)"), function(g) {
  sub <- aj_data %>% filter(genotype == g)
  fit <- prodlim(Hist(time_to_event, status) ~ 1, data = sub)
  # time points and CIF for cause 1
  tibble(
    time     = fit$time,
    cif      = fit$cuminc[[1]],   # cause 1
    genotype = g
  )
}) %>% bind_rows()

p_cif <- ggplot(aj_plot_data, aes(x = time, y = cif, colour = genotype)) +
  geom_step(linewidth = 0.8) +
  geom_vline(xintercept = eval_times, linetype = "dotted",
             colour = "grey50", linewidth = 0.4) +
  scale_color_d3() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(
    title    = "Cumulative incidence of MALO by HSD17B13 rs9992651 genotype",
    subtitle = "Aalen-Johansen estimator; non-liver death as competing event",
    x        = "Time (years)", y = "Cumulative incidence (MALO)",
    colour   = "Genotype\n(copies of A allele)"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

print(p_cif)

# -----------------------------------------------------------------
# 4. Cause-specific Cox regression: unadjusted and age/sex adjusted
#    Treats genotype as both numeric (additive) and factor (per-genotype)
# -----------------------------------------------------------------
cs_data <- model_data %>%
  mutate(status_cs = as.integer(status == 1))

# Additive (per-allele) — tests linear dose-response
cox_additive <- coxph(Surv(time_to_event, status_cs) ~ HSD17B13_rs9992651_A,
                      data = cs_data)

cox_additive_adj <- coxph(Surv(time_to_event, status_cs) ~
                            HSD17B13_rs9992651_A + age + sex,
                          data = cs_data)

# Factor — estimates HR per genotype vs ref (0)
cox_factor <- coxph(Surv(time_to_event, status_cs) ~
                      factor(HSD17B13_rs9992651_A),
                    data = cs_data)

cox_factor_adj <- coxph(Surv(time_to_event, status_cs) ~
                          factor(HSD17B13_rs9992651_A) + age + sex,
                        data = cs_data)

cat("\n=== Cox: additive (per A allele), unadjusted ===\n")
print(summary(cox_additive)$conf.int)

cat("\n=== Cox: additive (per A allele), age+sex adjusted ===\n")
print(summary(cox_additive_adj)$conf.int)

cat("\n=== Cox: per-genotype vs ref=0, unadjusted ===\n")
print(summary(cox_factor)$conf.int)

cat("\n=== Cox: per-genotype vs ref=0, age+sex adjusted ===\n")
print(summary(cox_factor_adj)$conf.int)

# -----------------------------------------------------------------
# 5. Forest plot of per-genotype HRs
# -----------------------------------------------------------------
tidy_cox <- function(fit, label) {
  s  <- summary(fit)$conf.int
  tibble(
    term   = gsub("_"," ", rownames(s)),
    HR     = s[, 1],
    HR_lo  = s[, 3],
    HR_hi  = s[, 4],
    p      = summary(fit)$coefficients[, 5],
    model  = label
  ) 
}

hr_data <- bind_rows(
  tidy_cox(cox_factor,     "Unadjusted"),
  tidy_cox(cox_factor_adj, "Age+sex adjusted")
) %>%
  filter(grepl("HSD", term)) %>%
  mutate(
    genotype = case_when(
      grepl(")1", term, fixed=TRUE) ~ "1 (ref/alt)",
      grepl(")2", term, fixed=TRUE) ~ "2 (alt/alt)"
    ),
    model = factor(model, levels = c("Unadjusted", "Age+sex adjusted")),
    genotype = factor(genotype, levels = c("2 (alt/alt)", "1 (ref/alt)"))
  )

p_hr <- ggplot(hr_data, aes(x = HR, y = genotype, color = model)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_errorbar(aes(xmin = HR_lo, xmax = HR_hi),
                position = position_dodgev(height = 0.5),
                orientation = "y", width = 0.2, linewidth = 0.5) +
  geom_point(size = 2.5,
             position = position_dodgev(height = 0.5)) +
  geom_label(aes(x = HR_hi * 1.1,
                 label = paste0("p=", formatC(p, digits=4, format="f"))),
             size = 3, color = "black", show.legend=FALSE,
             position = position_dodge2v(height = 0.5)) +
  scale_color_d3() + scale_x_log10(limits = c(0.7, 1.5)) +
  labs(
    title    = "Cause-specific HR for MALO",
    subtitle = "(additive HSD17B13 effect)",
    x = "Hazard Ratio (95% CI)", y = NULL, colour = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

print(p_hr)

# -----------------------------------------------------------------
# 6. Combined output
# -----------------------------------------------------------------
p_combined <- p_cif / p_hr +
  plot_annotation(tag_levels = "A")
print(p_combined)

ggsave(file.path(outdir, "BFA01.1.4_HSD17B13_diagnostic.png"),
       p_combined, width = 7, height = 9, dpi = 400, bg = "white")

cat("\nPlot saved. Check event rates and HRs to determine allele direction.\n")





# =============================================================================
# Multicollinearity diagnostics for precision model predictors
# VIF + correlation matrix + genotype frequency check
# =============================================================================

# -----------------------------------------------------------------
# 1. Genotype frequency and event rates for HSD17B13
#    Check if genotype 2 is simply too rare for stable SHAP estimates
# -----------------------------------------------------------------
geno_summary <- model_data %>%
  mutate(genotype = factor(HSD17B13_rs9992651_A,
                           levels = c(0, 1, 2),
                           labels = c("0 (ref/ref)", "1 (ref/alt)", "2 (alt/alt)"))) %>%
  group_by(genotype) %>%
  summarise(
    n           = n(),
    pct         = n() / nrow(model_data) * 100,
    n_malo      = sum(status == 1),
    malo_rate   = mean(status == 1) * 100,
    .groups     = "drop"
  )

cat("=== HSD17B13 rs9992651 genotype frequencies ===\n")
print(geno_summary)

# -----------------------------------------------------------------
# 2. VIF via linear regression
#    For each predictor, regress it on all others and compute R²
#    VIF = 1 / (1 - R²)
#    Use numeric versions of all predictors
# -----------------------------------------------------------------

# Convert all predictors to numeric for VIF calculation
vif_data <- model_data %>%
  dplyr::select(all_of(predictors_precision)) %>%
  mutate(across(everything(), ~ as.numeric(as.character(.x)))) %>%
  drop_na()

cat("\nN observations for VIF calculation:", nrow(vif_data), "\n")

# Fit a multivariate linear model (outcome irrelevant for VIF)
# Use a dummy continuous outcome — any works, we only care about predictor matrix
vif_lm <- lm(
  age ~ sex + bmi + has_t2dm + alcohol_grams_week + smoking_binary + waist_hip_ratio + trig_hdl_ratio +
    PNPLA3_rs738409_G + TM6SF2_rs58542926_T + HSD17B13_rs9992651_A,
  data = vif_data
)

# car::vif() on this model gives VIF for all predictors except the "outcome" (age)
# To get VIF for ALL predictors including age, compute manually via R²
compute_all_vif <- function(data, predictors) {
  vif_vals <- sapply(predictors, function(y) {
    others <- setdiff(predictors, y)
    f      <- as.formula(paste(y, "~", paste(others, collapse = " + ")))
    fit    <- lm(f, data = data)
    r2     <- summary(fit)$r.squared
    1 / (1 - r2)
  })
  tibble(
    predictor = predictors,
    VIF       = round(vif_vals, 3)
  ) %>% arrange(desc(VIF))
}

vif_results <- compute_all_vif(vif_data, predictors_precision)

cat("\n=== Variance Inflation Factors (precision model) ===\n")
cat("Rule of thumb: VIF > 5 = moderate concern, VIF > 10 = serious concern\n\n")
print(vif_results)

# -----------------------------------------------------------------
# 3. Pairwise Spearman correlation matrix
#    Spearman is appropriate because genotype is ordinal
# -----------------------------------------------------------------
cor_matrix <- cor(vif_data, method = "spearman", use = "pairwise.complete.obs")

nice_names <- c(
  age                  = "Age",
  sex                  = "Sex",
  bmi                  = "BMI",
  has_t2dm             = "T2DM",
  alcohol_grams_week   = "Alcohol",
  smoking_binary       = "Smoking",
  waist_hip_ratio      = "WHR",
  trig_hdl_ratio       = "Trig:HDL",
  PNPLA3_rs738409_G    = "PNPLA3",
  TM6SF2_rs58542926_T  = "TM6SF2",
  HSD17B13_rs9992651_A = "HSD17B13"
)

rownames(cor_matrix) <- nice_names[rownames(cor_matrix)]
colnames(cor_matrix) <- nice_names[colnames(cor_matrix)]

# Melt for ggplot
cor_long <- as.data.frame(cor_matrix) %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "rho") %>%
  mutate(
    var1 = factor(var1, levels = rev(colnames(cor_matrix))),
    var2 = factor(var2, levels = colnames(cor_matrix))
  )

p_cor <- ggplot(cor_long, aes(x = var2, y = var1, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D62728",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Spearman\n\u03c1") +
  labs(
    title = "Spearman correlation: precision model predictors",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid       = element_blank()
  )

print(p_cor)

# -----------------------------------------------------------------
# 4. VIF bar plot
# -----------------------------------------------------------------
p_vif <- ggplot(vif_results,
                aes(x = VIF, y = reorder(predictor, VIF))) +
  geom_col(aes(fill = VIF > 5), alpha = 0.8) +
  geom_vline(xintercept = c(5, 10), linetype = "dashed",
             colour = c("orange", "red"), linewidth = 0.5) +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "darkorange"),
                    labels = c("FALSE" = "VIF \u2264 5", "TRUE" = "VIF > 5"),
                    name   = NULL) +
  scale_y_discrete(labels = function(x) nice_names[x]) +
  labs(
    title = "Variance Inflation Factors: precision model",
    x = "VIF", y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

print(p_vif)

# -----------------------------------------------------------------
# 5. HSD17B13 correlation with each other predictor (focus plot)
# -----------------------------------------------------------------
hsd_cor <- cor_long %>%
  filter(var1 == "HSD17B13" | var2 == "HSD17B13") %>%
  filter(var1 != var2) %>%
  mutate(other = ifelse(var1 == "HSD17B13", as.character(var2), as.character(var1))) %>%
  distinct(other, .keep_all = TRUE) %>%
  arrange(desc(abs(rho)))

cat("\n=== Spearman correlations with HSD17B13 ===\n")
print(hsd_cor %>% dplyr::select(other, rho))

# -----------------------------------------------------------------
# 6. Combined output
# -----------------------------------------------------------------
p_combined <- (p_vif | p_cor) +
  plot_annotation(tag_levels = "A")
print(p_combined)

ggsave(file.path(outdir, "BFA01.1.4_HSD17B13_multicollinearity_diagnostics.png"),
       p_combined, width = 10, height = 5, dpi = 400, bg = "white")

# Export VIF table
write_xlsx(
  list(VIF = vif_results,
       Correlations = as.data.frame(cor_matrix) %>% rownames_to_column("predictor")),
  file.path(outdir, "BFA01.1.4_HSD17B13_multicollinearity_diagnostics.xlsx")
)

cat("\nDone. Outputs written to:", outdir, "\n")