# MALO Risk Calculator

Individualized prediction of **major adverse liver outcomes (MALO)** in at-risk steatotic liver disease (SLD) patients with low fibrosis burden, using a competing-risks Fine-Gray model trained on UK Biobank data.

---

## Background

Current guidelines recommend fibrosis surveillance for all at-risk SLD patients with FIB-4 below 2.67, but apply a one-size-fits-all interval. This tool uses a competing-risks survival model to estimate each patient's absolute 10-year MALO risk and compare it to the 0.5% surveillance threshold, enabling individualized surveillance planning.

**Target population:** At-risk SLD patients (type 2 diabetes, medically complicated obesity, metabolic syndrome, or excess alcohol intake) with low fibrosis burden (FIB-4 < 2.67), after exclusion of competing liver disease.

---

## Analytical Pipeline

Scripts run in numbered order; each depends on the outputs of the previous one. Paths are set by the `project_dir` variable at the top of each script and will need adjusting outside the original environment.

### 1. Data preparation (`analysis/01_preprocess_clean_ukbb_cohort.R`)

- Cleans and renames the raw UK Biobank cohort extract, including baseline clinical variables and individual genotypes for three MASLD-associated variants (*PNPLA3* rs738409, *TM6SF2* rs58542926, *HSD17B13* rs9992651)
- Defines the at-risk SLD subcohort (type 2 diabetes, medically complicated obesity, metabolic syndrome, or excess alcohol intake) and reports completeness of the six clinical predictors before any exclusion
- Applies exclusions: missing FIB-4 components, FIB-4 ≥ 2.67, and competing liver disease at baseline; reports a full exclusion cascade
- Derives the composite MALO endpoint — cirrhosis, hepatic decompensation, HCC, liver transplant, or liver-related death — with a 183-day washout, and the competing event (non-liver death)
- Outputs: `data/BFA_principal_data.feather`, `results/01_cohort_summaries.xlsx`

### 2. Main analysis — Fine-Gray (`analysis/02_FineGray_main_analysis.R`)

The primary model. All predictors enter as **linear** terms; whether nonlinearity helps is deferred to the RSF sensitivity analysis in script 03.

- **Clinical model predictors:** age, sex, BMI, type 2 diabetes, alcohol intake (g/week), smoking status
- **Precision model:** the six clinical predictors plus waist circumference and the three genotypes
- 5-fold cross-validation (stratified on MALO, seed 1234): time-dependent AUC, C-index and Brier score at 5 and 10 years, versus the CLivD score and FIB-4 alone
- Risk stratification at the median predicted 10-year risk (aggressive vs. standard surveillance): event counts, Aalen-Johansen cumulative incidence with Gray's test, and KM curves
- Decision curve analysis against the standard-of-care surveillance rules (type 2 diabetes; ≥ 2 CMRFs; the guideline "T2DM or ≥ 2 CMRFs")
- Decile-based calibration plots against the Aalen-Johansen CIF
- Full-data fits: subdistribution hazard ratio forest plots/tables and the model equation
- Exports the deployable clinical model (see below)
- Outputs: `results/02_FineGray_results.xlsx`, `results/plots/02.*.png`, `results/FG_clinical_model_deploy.rds`, and the out-of-fold predictions and fold assignment at `data/BFA_FineGray_OOF_predictions.feather`

### 3. Sensitivity analysis — random survival forest (`analysis/03_RSF_sensitivity_analysis.R`)

Checks whether allowing nonlinear effects and interactions meaningfully improves on the linear Fine-Gray model. Re-uses script 02's fold assignment, out-of-fold predictions and model equation, so both models are evaluated on exactly the same folds.

- 5-fold CV of the RSF (clinical and precision) using `randomForestSRC`
- Fine-Gray vs. RSF: per-fold differences and a 1000-replicate bootstrap test of the difference in pooled out-of-fold time-dependent AUC, C-index and Brier score
- Decision curve analysis and decile-based calibration plots for both models
- Partial dependence: linear Fine-Gray vs. RSF
- Variable importance (VIMP) plus SHAP beeswarm and bar plots (`fastshap`, `shapviz`)
- Forests are grown with `rfsrc.anonymous()` and then stripped of all remaining participant-level training data, so the saved model objects contain none
- Outputs: `results/03_RSF_results.xlsx`, `results/plots/03.*.png`

### 4. *HSD17B13* diagnostic (`analysis/04_HSD17B13_diagnostic.R`)

Investigates why the *HSD17B13* rs9992651 ALT/ALT genotype appears deleterious rather than protective in this cohort.

- Genotype frequencies and event rates; Aalen-Johansen cumulative incidence and curves by genotype
- Cause-specific Cox regression, unadjusted and age/sex adjusted, with a per-genotype HR forest plot
- Multicollinearity diagnostics: VIF, pairwise Spearman correlation matrix, and correlation of *HSD17B13* with each other predictor

### Supporting code (`analysis/clivd_scores.R`)

`compute_clivd()` and `compute_clivd_nonlab()` — the published CLivD score (with and without the laboratory GGT term), including the restricted cubic spline in alcohol intake. Sourced by script 02 as a comparator.

---

### **Note on data governance**

`app.R` and `model/FG_clinical_model_deploy.rds` are included in this repository. The deployable model is an export of aggregate quantities only — coefficients, a baseline cumulative subdistribution hazard on a fixed time grid, preprocessing constants and cohort-level summaries — assembled in Section 6 of `02_FineGray_main_analysis.R` and audited there to contain no element of cohort length and no individual-level UK Biobank data. The upstream analysis scripts read UK Biobank extracts that are **not** distributed here.

Public hosting of the calculator remains subject to UK Biobank approval.

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

- **Predicted MALO risk** at 5 and 10 years (cumulative incidence from the competing-risks Fine-Gray model)
- **Threshold comparison:** fold-difference relative to the 0.5% at 10 years surveillance threshold
- **Risk category** (High ≥ 0.5% / Low < 0.5%) with suggested surveillance interval
- **Principal modifiable risk factor**, with the absolute risk reduction that reaching its clinical target would achieve
- **Risk factor contributions** table covering all six predictors, against a fixed reference patient
- **Range warnings** when an input falls outside the training 1st–99th percentiles or is winsorised/capped before prediction

### Technical notes

- The app serves the linear Fine-Gray model exported by `02_FineGray_main_analysis.R` as
  `model/FG_clinical_model_deploy.rds` (5 KB, versioned in this repository). The file holds only
  coefficients, the baseline cumulative subdistribution hazard `H0(t)` tabulated on a fixed
  0.05-year grid, the preprocessing constants and aggregate metadata — no individual-level
  UK Biobank data. There is no model download at launch.
- Risk is the closed-form `CIF(t | x) = 1 - exp(-H0(t) * exp(sum(beta * x)))`, evaluated by the
  `predict_fn` shipped inside the `.rds`, so the app never reimplements the training transforms.
- **Attribution replaces SHAP.** Because the model is linear in the log-subdistribution-hazard,
  moving predictor *j* to a reference value multiplies risk by exactly `exp(beta_j * (x_j - ref_j))`,
  independently of the other predictors. For each factor the app reports that exact multiple
  (these multiply to the patient's total risk multiple versus the reference patient) and the
  absolute percentage points of 10-year risk attributable to it. The principal modifiable factor
  is the one carrying the most absolute risk relative to its clinical target. Reference patient:
  55-year-old man, BMI 25 kg/m², no type 2 diabetes, no alcohol, non-smoker.
- Alcohol values above 500 g/week are capped and BMI is winsorised to the training bounds, by the
  model's own `predict_fn`; the UI flags when this has happened.
- Dependencies are `shiny` and `bslib` only. Prediction is closed-form and instantaneous, so the
  progress bar and button-locking of the RSF version are gone.

---

## Disclaimer

For use in at-risk SLD patients with low fibrosis burden (FIB-4 < 2.67) after exclusion of competing liver disease. This tool is intended to individualize fibrosis surveillance intervals within the existing guideline framework. Not validated for use outside this population.
