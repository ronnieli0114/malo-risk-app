# MALO Risk Calculator

Individualized prediction of **major adverse liver outcomes (MALO)** in at-risk steatotic liver disease (SLD) patients with low fibrosis burden, using a random survival forest (RSF) trained on UK Biobank data.

---

## Background

Current guidelines recommend fibrosis surveillance for all at-risk SLD patients with FIB-4 ≤ 2.67, but apply a one-size-fits-all interval. This tool uses a competing-risks survival model to estimate each patient's absolute 10-year MALO risk and compare it to the 0.5% surveillance threshold, enabling individualized surveillance planning.

**Target population:** At-risk SLD patients with low fibrosis burden (FIB-4 ≤ 2.67), after exclusion of competing liver disease.

---

## Analytical Pipeline

### 1. Data preparation (`BFA01.1.1_preprocess_and_clean.R`)

- Extracted and cleaned UK Biobank cohort data, including baseline clinical variables and individual genotypes for three MASLD-associated variants (including *HSD17B13*)
- Defined the at-risk SLD subcohort and applied inclusion/exclusion criteria
- Derived the composite MALO endpoint and competing event (non-liver death)
- Output: `BFA_principal_data.feather`

### 2. Main RSF analysis (`BFA01.1.2_RSF_main_analysis.R`)

- Trained a **Random Survival Forest** (competing risks, cause 1 = MALO, cause 2 = non-liver death) using `randomForestSRC`
- Predictors: age, sex, BMI, type 2 diabetes, alcohol intake (g/week), smoking status
- Model evaluation: time-dependent AUC, Brier score, and C-index compared against FIB-4 alone and the CLivD score
- Benchmarked an extended "precision" model adding HTN, HDL, triglycerides, and HbA1c
- Assessed nonlinearity via LRT and partial dependence plots vs. cause-specific Cox splines
- Generated risk-stratified cumulative incidence curves (Aalen-Johansen), decision curve analysis, and calibration plots (decile-based)
- Computed variable importance (VIMP) and SHAP values (beeswarm and bar plots)
- Saved the final full-dataset clinical model as `RSF_final_clinical_model.rds`

### 3. Fine-Gray comparison (`BFA01.1.3_Fine_Gray_main_model.R`)

- Fit a **Fine-Gray subdistribution hazard model** as a parametric benchmark
- Compared discrimination and calibration against the RSF

### 4. *HSD17B13* diagnostic analysis (`BFA01.1.4_HSD17B13_diagnostic.R`)

- Evaluated whether *HSD17B13* rs72613567 genotype adds incremental prognostic value beyond the clinical model

---

### **Note**: The section below has been developed but not uploaded to GitHub due to UK Biobank rules on data governance. We are actively working with UK Biobank to see whether this tool can be hosted publicly.

## Shiny App (`app.R`)

An interactive risk calculator hosted as a Shiny web app.

### Inputs

| Variable | Details |
|---|---|
| Age | years |
| Sex | Male / Female |
| BMI | kg/m² |
| Type 2 diabetes | Yes / No |
| Alcohol intake | Standard drinks/week **or** grams/week (toggle between modes) |
| Smoking status | Current / Non-current |

### Outputs

- **Predicted MALO risk** at 5 and 10 years (cumulative incidence from competing-risks RSF)
- **Threshold comparison:** fold-difference relative to the 0.5% at 10 years surveillance threshold
- **Risk category** (High ≥ 0.5% / Low < 0.5%) with suggested surveillance interval
- **Top modifiable risk factor** identified via Kernel SHAP, with direction of effect

### Technical notes

- The RSF model (~1 GB) is downloaded from Dropbox on first launch and cached for the process lifetime
- SHAP values are computed at runtime using `kernelshap` against a pre-computed background dataset (`background_data.rds`)
- Alcohol values above 500 g/week are capped to match training data bounds
- The Calculate Risk button is disabled during computation to prevent duplicate submissions

---

## Disclaimer

For use in at-risk SLD patients with low fibrosis burden (FIB-4 ≤ 1.67) after exclusion of competing liver disease. This tool is intended to individualize fibrosis surveillance intervals within the existing guideline framework. Not validated for use outside this population.
