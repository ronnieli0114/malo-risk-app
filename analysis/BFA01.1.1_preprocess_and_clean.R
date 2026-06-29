### This script cleans and preprocesses the data extracted from the UK Biobank cohort
### along with the individual genotypes for three MASLD-associated variants,
### and outputs a clean data frame (in .feather) format readable by arrow (R) or pyarrow (Python).

library(tidyverse)
library(arrow)
library(writexl)

project_dir <- "/mnt/d/Projects/BFA"
data_dir <- file.path(project_dir, "data")

# Path to raw data
path_to_raw_csv <- file.path(data_dir, "2026-02-19_cohort_extract.csv")

# Read in raw data
raw <- read.csv(path_to_raw_csv, header = TRUE)
cat("Read in raw data. Cleaning...\n")

# Rename columns to interpretable names
raw_renamed <- raw %>%
  dplyr::rename(
    
    # Participant ID
    participant_id = eid,
    
    # Assessment dates
    assessment_date_i0 = p53_i0,
    assessment_date_i1 = p53_i1,
    assessment_date_i2 = p53_i2,
    assessment_date_i3 = p53_i3,
    
    # Age at assessment
    age_i0 = p21003_i0,
    age_i1 = p21003_i1,
    age_i2 = p21003_i2,
    age_i3 = p21003_i3,
    
    # Sex
    sex = p31,
    
    # Ethnic background
    ethnic_background_i0 = p21000_i0,
    ethnic_background_i1 = p21000_i1,
    ethnic_background_i2 = p21000_i2,
    ethnic_background_i3 = p21000_i3,
    
    # BMI (kg/m²)
    bmi_i0 = p21001_i0,
    bmi_i1 = p21001_i1,
    bmi_i2 = p21001_i2,
    bmi_i3 = p21001_i3,
    
    # Smoking status
    smoking_status_i0 = p20116_i0,
    smoking_status_i1 = p20116_i1,
    smoking_status_i2 = p20116_i2,
    smoking_status_i3 = p20116_i3,
    
    # Waist circumference (cm)
    waist_circ_i0 = p48_i0,
    
    # Hip circumference (cm)
    hip_circ_i0 = p49_i0,
    
    # Fasting duration
    fasting_i0 = p74_i0,
    fasting_i1 = p74_i1,
    fasting_i2 = p74_i2,
    fasting_i3 = p74_i3,
    
    # Serum glucose
    serum_glucose_i0 = p30740_i0,
    serum_glucose_i1 = p30740_i1,
    
    # Diabetes diagnosed by doctor
    diabetes_doctor_i0 = p2443_i0,
    diabetes_doctor_i1 = p2443_i1,
    diabetes_doctor_i2 = p2443_i2,
    diabetes_doctor_i3 = p2443_i3,
    
    # Age at diabetes diagnosis
    age_diabetes_dx_i0 = p2976_i0,
    age_diabetes_dx_i1 = p2976_i1,
    age_diabetes_dx_i2 = p2976_i2,
    age_diabetes_dx_i3 = p2976_i3,
    
    # HbA1c (mmol/mol)
    hba1c_i0 = p30750_i0,
    hba1c_i1 = p30750_i1,
    
    # GGT (U/L)
    ggt_i0 = p30730_i0,
    
    # Systolic blood pressure (mmHg)
    sbp_i0_a0 = p4080_i0_a0,
    sbp_i0_a1 = p4080_i0_a1,
    sbp_i1_a0 = p4080_i1_a0,
    sbp_i1_a1 = p4080_i1_a1,
    sbp_i2_a0 = p4080_i2_a0,
    sbp_i2_a1 = p4080_i2_a1,
    sbp_i3_a0 = p4080_i3_a0,
    sbp_i3_a1 = p4080_i3_a1,
    
    # Diastolic blood pressure (mmHg)
    dbp_i0_a0 = p4079_i0_a0,
    dbp_i0_a1 = p4079_i0_a1,
    dbp_i1_a0 = p4079_i1_a0,
    dbp_i1_a1 = p4079_i1_a1,
    dbp_i2_a0 = p4079_i2_a0,
    dbp_i2_a1 = p4079_i2_a1,
    dbp_i3_a0 = p4079_i3_a0,
    dbp_i3_a1 = p4079_i3_a1,
    
    # Triglycerides (mmol/L)
    triglycerides_i0 = p30870_i0,
    triglycerides_i1 = p30870_i1,
    
    # HDL cholesterol (mmol/L)
    hdl_i0 = p30760_i0,
    hdl_i1 = p30760_i1,
    
    # Alcohol intake - Red wine (glasses/week)
    red_wine_i0 = p1568_i0,
    red_wine_i1 = p1568_i1,
    red_wine_i2 = p1568_i2,
    red_wine_i3 = p1568_i3,
    
    # Alcohol intake - White wine/champagne (glasses/week)
    white_wine_i0 = p1578_i0,
    white_wine_i1 = p1578_i1,
    white_wine_i2 = p1578_i2,
    white_wine_i3 = p1578_i3,
    
    # Alcohol intake - Beer/cider (pints/week)
    beer_cider_i0 = p1588_i0,
    beer_cider_i1 = p1588_i1,
    beer_cider_i2 = p1588_i2,
    beer_cider_i3 = p1588_i3,
    
    # Alcohol intake - Spirits (measures/week)
    spirits_i0 = p1598_i0,
    spirits_i1 = p1598_i1,
    spirits_i2 = p1598_i2,
    spirits_i3 = p1598_i3,
    
    # Alcohol intake - Fortified wine (glasses/week)
    fortified_wine_i0 = p1608_i0,
    fortified_wine_i1 = p1608_i1,
    fortified_wine_i2 = p1608_i2,
    fortified_wine_i3 = p1608_i3,
    
    # AST (U/L)
    ast_i0 = p30650_i0,
    ast_i1 = p30650_i1,
    
    # ALT (U/L)
    alt_i0 = p30620_i0,
    alt_i1 = p30620_i1,
    
    # Platelet count (10⁹/L)
    platelets_i0 = p30080_i0,
    platelets_i1 = p30080_i1,
    platelets_i2 = p30080_i2,
    
    # Date of death
    date_of_death_i0 = p40000_i0,
    date_of_death_i1 = p40000_i1,
    
    # Primary cause of death (ICD-10)
    death_cause_primary_i0 = p40001_i0,
    death_cause_primary_i1 = p40001_i1,
    
    # OPCS procedure code
    opcs_procedure = p41272,
    
    # ICD-10 diagnosis codes (from hesin)
    icd10_codes = p41270,
    
    # ICD-9 diagnosis codes (from hesin)
    icd9_codes = p41271
  )

# Rename array columns with a loop for efficiency
# Self-reported illness codes (field 20002)
for(i in 0:3) {
  for(j in 0:36) {
    old_name <- paste0("p20002_i", i, "_a", j)
    new_name <- paste0("illness_code_i", i, "_a", j)
    if(old_name %in% names(raw_renamed)) {
      raw_renamed <- raw_renamed %>% rename(!!new_name := !!old_name)
    }
  }
}

# Self-reported medications - men (field 6177)
for(i in 0:3) {
  old_name <- paste0("p6177_i", i)
  new_name <- paste0("medication_men_i", i)
  if(old_name %in% names(raw_renamed)) {
    raw_renamed <- raw_renamed %>% rename(!!new_name := !!old_name)
  }
}

# Self-reported medications - women (field 6153)
for(i in 0:3) {
  old_name <- paste0("p6153_i", i)
  new_name <- paste0("medication_women_i", i)
  if(old_name %in% names(raw_renamed)) {
    raw_renamed <- raw_renamed %>% rename(!!new_name := !!old_name)
  }
}

# Prescription medication codes (field 20003)
for(i in 0:3) {
  for(j in 0:47) {
    old_name <- paste0("p20003_i", i, "_a", j)
    new_name <- paste0("prescription_code_i", i, "_a", j)
    if(old_name %in% names(raw_renamed)) {
      raw_renamed <- raw_renamed %>% rename(!!new_name := !!old_name)
    }
  }
}

# ICD-10 diagnosis dates (field 41280)
for(j in 0:258) {
  old_name <- paste0("p41280_a", j)
  new_name <- paste0("icd10_date_a", j)
  if(old_name %in% names(raw_renamed)) {
    raw_renamed <- raw_renamed %>% rename(!!new_name := !!old_name)
  }
}

# ICD-9 diagnosis dates (field 41281)
for(j in 0:46) {
  old_name <- paste0("p41281_a", j)
  new_name <- paste0("icd9_date_a", j)
  if(old_name %in% names(raw_renamed)) {
    raw_renamed <- raw_renamed %>% rename(!!new_name := !!old_name)
  }
}

# Cancer diagnosis dates (field 40005)
for(j in 0:21) {
  old_name <- paste0("p40005_i", j)
  new_name <- paste0("cancer_date_i", j)
  if(old_name %in% names(raw_renamed)) {
    raw_renamed <- raw_renamed %>% rename(!!new_name := !!old_name)
  }
}

# Cancer ICD-10 codes (field 40006)
for(j in 0:21) {
  old_name <- paste0("p40006_i", j)
  new_name <- paste0("cancer_icd10_i", j)
  if(old_name %in% names(raw_renamed)) {
    raw_renamed <- raw_renamed %>% rename(!!new_name := !!old_name)
  }
}

# Secondary causes of death (field 40002)
for(i in 0:1) {
  for(j in 0:14) {
    old_name <- paste0("p40002_i", i, "_a", j)
    new_name <- paste0("death_cause_secondary_i", i, "_a", j)
    if(old_name %in% names(raw_renamed)) {
      raw_renamed <- raw_renamed %>% rename(!!new_name := !!old_name)
    }
  }
}

raw_renamed[raw_renamed == ""] <- NA
all_na_cols <- sapply(colnames(raw_renamed), function(col) sum(is.na(raw_renamed[[col]])) == nrow(raw_renamed))
raw_renamed <- raw_renamed[,!all_na_cols]

# Clean up the data
data_clean <- raw_renamed %>%
  
  # ============================================================================
  # 1. Collapse array fields into comma-separated strings
  # ============================================================================
  
  # Combine illness codes (all instances and arrays)
  unite("illnesses_all", starts_with("illness_code_"), 
        sep = ",", na.rm = TRUE, remove = TRUE) %>%
  
  # Combine prescription codes
  unite("prescriptions_all", starts_with("prescription_code_"), 
        sep = ",", na.rm = TRUE, remove = TRUE) %>%
  
  # Combine ICD-10 diagnosis dates
  unite("icd10_dates_all", starts_with("icd10_date_"), 
        sep = ",", na.rm = TRUE, remove = TRUE) %>%
  
  # Combine ICD-9 diagnosis dates
  unite("icd9_dates_all", starts_with("icd9_date_"), 
        sep = ",", na.rm = TRUE, remove = TRUE) %>%
  
  # Combine cancer dates
  unite("cancer_dates_all", starts_with("cancer_date_"), 
        sep = ",", na.rm = TRUE, remove = TRUE) %>%
  
  # Combine cancer ICD-10 codes
  unite("cancer_icd10_all", starts_with("cancer_icd10_"), 
        sep = ",", na.rm = TRUE, remove = TRUE) %>%
  
  # Combine secondary death causes
  unite("death_cause_secondary_all", starts_with("death_cause_secondary_"), 
        sep = ",", na.rm = TRUE, remove = TRUE) %>%
  
  # ============================================================================
  # 2. Calculate derived variables for baseline (instance 0)
  # ============================================================================
  
  mutate(
    
    # Alcohol intake (replace NA with 0)
    red_wine_i0 = ifelse(is.na(red_wine_i0), 0, red_wine_i0),
    white_wine_i0 = ifelse(is.na(white_wine_i0), 0, white_wine_i0),
    beer_cider_i0 = ifelse(is.na(beer_cider_i0), 0, beer_cider_i0),
    spirits_i0 = ifelse(is.na(spirits_i0), 0, spirits_i0),
    fortified_wine_i0 = ifelse(is.na(fortified_wine_i0), 0, fortified_wine_i0),
    
    # Convert -1 (don't know) and -3 (prefer not to answer) to NA for alcohol
    across(starts_with("red_wine_") | starts_with("white_wine_") | 
             starts_with("beer_cider_") | starts_with("spirits_") | 
             starts_with("fortified_wine_"),
           ~ifelse(. %in% c("Do not know", "Prefer not to answer"), NA, .)),
    
    # Average blood pressure readings (instance 0)
    sbp_avg_i0 = rowMeans(select(., sbp_i0_a0, sbp_i0_a1), na.rm = TRUE),
    dbp_avg_i0 = rowMeans(select(., dbp_i0_a0, dbp_i0_a1), na.rm = TRUE),
    
    # FIB-4 score (baseline only)
    fib4_i0 = (age_i0 * ast_i0) / (platelets_i0 * sqrt(alt_i0)),
    
    # Alcohol intake in grams per week (baseline)
    # Formula: (drinks/week × units/drink) × 8g
    alcohol_grams_week_i0 = tryCatch({(
        (as.numeric(red_wine_i0) * 2) +        # 2 units per glass
        (as.numeric(white_wine_i0) * 2) +      # 2 units per glass
        (as.numeric(beer_cider_i0) * 2) +      # 2 units per pint
        (as.numeric(spirits_i0) * 1) +         # 1 unit per measure
        (as.numeric(fortified_wine_i0) * 1)    # 1 unit per glass
    ) * 8}, error = function(e) NA)
  
  ) %>%
  
  # ============================================================================
  # 3. Keep only baseline (instance 0) for key variables
  # ============================================================================

  dplyr::select(
    # Identifiers
    participant_id,
    
    # Demographics (baseline)
    assessment_date_i0,
    age_i0,
    sex,
    ethnic_background_i0,
    
    # Anthropometrics (baseline)
    bmi_i0,
    waist_circ_i0,
    hip_circ_i0,
    
    # Clinical labs (baseline)
    ast_i0,
    alt_i0,
    ggt_i0,
    platelets_i0,
    hba1c_i0,
    fasting_i0,
    serum_glucose_i0,
    hdl_i0,
    triglycerides_i0,
    
    # Blood pressure (baseline - averaged)
    sbp_avg_i0,
    dbp_avg_i0,
    
    # Risk factors (baseline)
    smoking_status_i0,
    diabetes_doctor_i0,
    age_diabetes_dx_i0,
    
    # Alcohol (baseline)
    alcohol_grams_week_i0,
    
    # Calculated scores
    fib4_i0,
    
    # Medications (baseline)
    medication_men_i0,
    medication_women_i0,
    
    # Collapsed arrays
    illnesses_all,
    prescriptions_all,
    
    # Death data
    date_of_death_i0,
    death_cause_primary_i0,
    death_cause_secondary_all,
    
    # ICD codes and dates
    icd10_codes,
    icd9_codes,
    icd10_dates_all,
    icd9_dates_all,
    
    # Cancer registry
    cancer_dates_all,
    cancer_icd10_all,
    
    # OPCS procedures
    opcs_procedure

  ) %>%
  
  # ============================================================================
  # 4. Clean up empty strings in collapsed fields
  # ============================================================================

  mutate(
    across(ends_with("_all"), ~ifelse(. == "", NA, .))
  ) %>%
  
  # ============================================================================
  # 5. Rename to shorter, cleaner names
  # ============================================================================

  dplyr::rename(
    eid = participant_id,
    assessment_date = assessment_date_i0,
    age = age_i0,
    ethnicity = ethnic_background_i0,
    bmi = bmi_i0,
    waist_circ = waist_circ_i0,
    hip_circ = hip_circ_i0,
    ast = ast_i0,
    alt = alt_i0,
    ggt = ggt_i0,
    platelets = platelets_i0,
    hba1c = hba1c_i0,
    fasting_hours = fasting_i0,
    serum_glucose = serum_glucose_i0,
    hdl = hdl_i0,
    triglycerides = triglycerides_i0,
    systolic_avg = sbp_avg_i0,
    diastolic_avg = dbp_avg_i0,
    smoking = smoking_status_i0,
    diabetes_dx = diabetes_doctor_i0,
    age_diabetes = age_diabetes_dx_i0,
    alcohol_grams_week = alcohol_grams_week_i0,
    medication_men = medication_men_i0,
    medication_women = medication_women_i0,
    fib4 = fib4_i0
  )

# ============================================================================
# 6. Create helper columns for filtering
# ============================================================================

data_clean <- data_clean %>%

  mutate(
    
    # Waist-hip ratio
    waist_hip_ratio = waist_circ / hip_circ,
    
    # Has prediabetes (from supplement)
    has_prediabetes = case_when(
      hba1c >= 39 & hba1c < 48 ~ TRUE, # in mmol/mol (5.7-6.5%)
      fasting_hours >= 8 & serum_glucose >= 5.6 & serum_glucose <= 6.9 ~ TRUE,   # in mmol/L (100-125 mg/dL)
      TRUE ~ FALSE
    ),

    # Has Type 2 Diabetes (from supplement criteria)
    has_t2dm = case_when(
      diabetes_dx == "Yes" ~ TRUE,  # Self-reported diagnosis
      grepl("diabetes", illnesses_all, ignore.case = TRUE) ~ TRUE,  # Illness codes
      hba1c >= 48 ~ TRUE,  # Lab evidence (HbA1c >= 6.5%)
      grepl("insulin|metformin|rosiglitazone|        # Medications
            gliclazide|glimepiride|glipizide|
            tolbutamide|chlorpropamide|pioglitazone|
            rosiglitazone|repaglinide|nateglinide|acarbose", prescriptions_all,
            ignore.case = TRUE) ~ TRUE,  
      TRUE ~ FALSE
    ),
    
    # Prediabetes OR Type 2 diabetes
    prediabetes_or_t2dm = has_prediabetes | has_t2dm,
    
    # Obese (BMI >= 30)
    is_obese = bmi >= 30,
    
    # Has hypertension
    has_hypertension = case_when(
      grepl("hypertension", illnesses_all) ~ TRUE,  # Self-reported
      systolic_avg >= 130 | diastolic_avg >= 80 ~ TRUE,  # Lab evidence
      grepl("Blood pressure medication", medication_men, ignore.case = TRUE) ~ TRUE,
      grepl("Blood pressure medication", medication_women, ignore.case = TRUE) ~ TRUE,
      grepl("lisinopril|ramipril|enalapril|perindopril|losartan|candesartan|
            valsartan|irbesartan|amlodipine|felodipine|nifedipine|bendroflumethiazide|
            indapamide|hydrochlorothiazide|chlortalidone", prescriptions_all,
            ignore.case = TRUE) ~ TRUE,
      TRUE ~ FALSE
    ),
    
    # Has hypertriglyceridemia
    has_hypertrig = case_when(
      triglycerides > 1.7 ~ TRUE,  # > 150 mg/dL
      grepl("Cholesterol lowering medication", medication_men, ignore.case = TRUE) ~ TRUE,
      grepl("Cholesterol lowering medication", medication_women, ignore.case = TRUE) ~ TRUE,
      grepl("simvastatin|atorvastatin|pravastatin|rosuvastatin|fluvastatin|
            bezafibrate|fenofibrate|gemfibrozil|ciprofibrate|ezetimibe", 
            prescriptions_all, ignore.case = TRUE) ~ TRUE,
      TRUE ~ FALSE
    ),
      
    # Has low HDL
    has_low_hdl = case_when(
      sex == "Male" & hdl < 1.03 ~ TRUE,  # Males < 40 mg/dL
      sex == "Female" & hdl < 1.29 ~ TRUE,  # Females < 50 mg/dL
      TRUE ~ FALSE
    ),
    
    # Large waist circumference
    has_large_waist = case_when(
      sex == "Male" & waist_circ > 102 ~ TRUE,  # Males
      sex == "Female" & waist_circ > 88 ~ TRUE,   # Females
      TRUE ~ FALSE
    ),

    # Excessive alcohol consumption
    has_excess_alcohol = case_when(
      sex == "Male" & alcohol_grams_week > 210 ~ TRUE,  # Males
      sex == "Female" & alcohol_grams_week > 140 ~ TRUE,  # Females
      TRUE ~ FALSE
    ),
    
    # FIB-4 > 2.67
    fib4_high = fib4 >= 2.67
  ) %>%
  
  # Calculate n_cmrfs after the flags are created
  mutate(
    # Count cardiometabolic risk factors (CMRFs)
    n_cmrfs = rowSums(cbind(has_large_waist, prediabetes_or_t2dm,
                            has_hypertension, has_hypertrig, 
                            has_low_hdl), na.rm = TRUE),
    
    # Has metabolic syndrome (>= 3 CMRFs)
    has_metabolic_syndrome = n_cmrfs >= 3,
    
    # Has medically complicated obesity (obesity + at least 1 CMRF)
    has_complicated_obesity = is_obese & n_cmrfs >= 1,
    
    # Has dyslipidemia (HDL < 1.0 or triglycerides > 1.7)
    has_dyslipidemia = has_low_hdl | has_hypertrig
    
  ) %>%
  mutate(
    eid = paste0("UKBB_", eid),
    assessment_date = as.Date(assessment_date),
    date_of_death = as.Date(date_of_death_i0),
    age_diabetes = ifelse(age_diabetes=='NaN', NA_real_, as.numeric(age_diabetes))
  ) %>%
  dplyr::select(-date_of_death_i0)

#===============================================================================
# Incorporate genotype information (.raw files)
#===============================================================================

df_list <- list()
for (f in list.files(data_dir, pattern = ".raw")) {
  
  df <- read.table(file.path(data_dir, f), sep = "\t", header = TRUE)
  df <- df[,c(2,7)] # select IID and genotype columns
  df <- df %>% 
    mutate(eid = paste0("UKBB_", IID)) %>%
    column_to_rownames("eid") %>%
    dplyr::select(-IID)
  df[,1] <- as.integer(2 - df[,1]) # convert to minor allele dosage
  df_list[[f]] <- df
}
gt <- do.call(cbind, df_list) %>%
  rownames_to_column("eid")
colnames(gt) <- c("eid","TM6SF2_rs58542926_T","PNPLA3_rs738409_G","HSD17B13_rs9992651_A")

data_clean_with_gt <- inner_join(data_clean, gt, by = "eid")


#===============================================================================
# Helper function to check if diagnosis occurred before assessment date
#===============================================================================

has_diagnosis_before_baseline <- function(icd_dates_all, icd_codes, pattern, assessment_date) {
  if(is.na(icd_dates_all) || is.na(icd_codes) || is.na(assessment_date)) return(FALSE)
  
  dates <- strsplit(as.character(icd_dates_all), ",")[[1]]
  codes <- strsplit(as.character(icd_codes), ",")[[1]]
  
  if(length(dates) == 0 || length(codes) == 0) return(FALSE)
  
  # Find matching codes
  matches <- grepl(pattern, codes)
  if(!any(matches)) return(FALSE)
  
  # Get corresponding dates
  matched_dates <- as.Date(dates[matches])
  baseline <- as.Date(assessment_date)
  
  # Check if any occurred before or at baseline
  any(matched_dates <= baseline, na.rm = TRUE)
}

#===============================================================================
# Identify competing liver disease at baseline, apply inclusion and exclusion criteria
#===============================================================================

data_processed <- data_clean_with_gt %>%
  
  dplyr::filter(
    # =========================================================================
    # Exclude missing variables that need to be used in the model
    # =========================================================================
    !is.na(age), !is.na(sex), !is.na(bmi), !is.na(has_t2dm),
    !is.na(hba1c), !is.na(hdl), !is.na(triglycerides), !is.na(has_hypertension),
    !is.na(alcohol_grams_week), 
    smoking %in% c("Never","Previous","Current"),
    diabetes_dx %in% c("Yes","No"),
    !is.na(ggt), !is.na(hip_circ), !is.na(waist_circ), 
    !is.na(PNPLA3_rs738409_G), !is.na(TM6SF2_rs58542926_T), !is.na(HSD17B13_rs9992651_A),
  )  %>%
  rowwise() %>%
  mutate(
    
    # =========================================================================
    # COMPETING LIVER DISEASES (must occur BEFORE or AT assessment date)
    # =========================================================================
    
    # 1. Chronic Viral Hepatitis
    has_viral_hepatitis = has_diagnosis_before_baseline(
      icd10_dates_all, icd10_codes, 
      "B16\\.|B17\\.1|B18\\.|B19\\.", 
      assessment_date
    ) | has_diagnosis_before_baseline(
      icd9_dates_all, icd9_codes, 
      "070\\.2|070\\.3|070\\.4|070\\.5", 
      assessment_date
    ),
    
    # 2. Autoimmune / Cholestatic Liver Disease
    has_autoimmune_liver = has_diagnosis_before_baseline(
      icd10_dates_all, icd10_codes, 
      "K75\\.4|K74\\.3|K83\\.0", 
      assessment_date
    ) | has_diagnosis_before_baseline(
      icd9_dates_all, icd9_codes, 
      "571\\.42|571\\.6|576\\.1", 
      assessment_date
    ),
    
    # 3. Metabolic Inherited Liver Disorders
    has_metabolic_liver = has_diagnosis_before_baseline(
      icd10_dates_all, icd10_codes, 
      "E83\\.1|E83\\.0|E88\\.0", 
      assessment_date
    ) | has_diagnosis_before_baseline(
      icd9_dates_all, icd9_codes, 
      "275\\.0|275\\.1|273\\.4", 
      assessment_date
    ),
    
    # 4. Toxic Liver Disease
    has_toxic_liver = has_diagnosis_before_baseline(
      icd10_dates_all, icd10_codes, 
      "K71\\.", 
      assessment_date
    ) | has_diagnosis_before_baseline(
      icd9_dates_all, icd9_codes, 
      "573\\.3", 
      assessment_date
    ),
    
    # 5. Other (Budd-Chiari)
    has_budd_chiari = has_diagnosis_before_baseline(
      icd10_dates_all, icd10_codes, 
      "I82\\.0", 
      assessment_date
    ) | has_diagnosis_before_baseline(
      icd9_dates_all, icd9_codes, 
      "453\\.0", 
      assessment_date
    ),
    
    # 6a. Cirrhosis & Fibrosis at baseline
    has_cirrhosis = has_diagnosis_before_baseline(
      icd10_dates_all, icd10_codes, 
      "K70\\.3|K74\\.0|K74\\.1|K74\\.2|K74\\.6|K76\\.6", 
      assessment_date
    ) | has_diagnosis_before_baseline(
      icd9_dates_all, icd9_codes, 
      "571\\.2|571\\.5|572\\.3", 
      assessment_date
    ),
    
    # 6b. Hepatic Decompensation at baseline
    has_decompensation = has_diagnosis_before_baseline(
      icd10_dates_all, icd10_codes, 
      "I85\\.0|I85\\.9|I86\\.4|K70\\.4|K72\\.|K76\\.7|R18", 
      assessment_date
    ) | has_diagnosis_before_baseline(
      icd9_dates_all, icd9_codes, 
      "456\\.0|456\\.1|456\\.2|570|572\\.2|572\\.4|789\\.5", 
      assessment_date
    ),
    
    # 6c. Hepatocellular Carcinoma at baseline
    has_hcc = has_diagnosis_before_baseline(
      icd10_dates_all, icd10_codes, 
      "C22\\.0", 
      assessment_date
    ) | has_diagnosis_before_baseline(
      icd9_dates_all, icd9_codes, 
      "155\\.0", 
      assessment_date
    ) | has_diagnosis_before_baseline(
      cancer_dates_all, cancer_icd10_all, 
      "C22\\.0", 
      assessment_date
    ),
    
    # 6d. Liver Transplant at baseline
    has_liver_transplant = has_diagnosis_before_baseline(
      icd10_dates_all, icd10_codes, 
      "Z94\\.4", 
      assessment_date
    ) | has_diagnosis_before_baseline(
      icd9_dates_all, icd9_codes, 
      "V42\\.7|996\\.81", 
      assessment_date
    ) | (!is.na(opcs_procedure) && grepl("J01", opcs_procedure)),
    
    # Any competing liver disease
    has_competing_liver_disease = has_viral_hepatitis | has_autoimmune_liver | 
      has_metabolic_liver | has_toxic_liver | 
      has_budd_chiari | has_cirrhosis | 
      has_decompensation | has_hcc | has_liver_transplant
    
  ) %>%
  ungroup() %>%
  mutate(
    
    # =========================================================================
    # INCLUSION CRITERIA
    # =========================================================================
    
    meets_inclusion = has_t2dm | has_complicated_obesity | 
      has_metabolic_syndrome | has_excess_alcohol,
    
    # =========================================================================
    # EXCLUSION CRITERIA
    # =========================================================================
    
    # Missing FIB-4 components
    missing_fib4 = is.na(age) | is.na(ast) | is.na(alt) | is.na(platelets),
    
    # Apply exclusions
    exclude_missing_fib4 = missing_fib4,
    exclude_high_fib4 = fib4_high,
    exclude_competing_disease = has_competing_liver_disease,
    
    # Final inclusion/exclusion flag
    excluded = exclude_missing_fib4 | exclude_high_fib4 | exclude_competing_disease,
    included_in_cohort = meets_inclusion & !excluded
  )

#===============================================================================
# Generate Summary Statistics Data Frame
#===============================================================================

# Calculate counts for flow diagram
total_n <- nrow(data_processed)
meets_inclusion_n <- sum(data_processed$meets_inclusion, na.rm = TRUE)
excluded_among_eligible_n <- sum(data_processed$meets_inclusion & data_processed$excluded, na.rm = TRUE)
final_cohort_n <- sum(data_processed$included_in_cohort, na.rm = TRUE)

# Inclusion criteria counts
inclusion_summary <- data.frame(
  Criterion = c("Type 2 Diabetes", 
                "Medically Complicated Obesity", 
                "Metabolic Syndrome", 
                "Excess Alcohol Consumption",
                "Any Inclusion Criterion Met"),
  N = c(
    sum(data_processed$has_t2dm, na.rm = TRUE),
    sum(data_processed$has_complicated_obesity, na.rm = TRUE),
    sum(data_processed$has_metabolic_syndrome, na.rm = TRUE),
    sum(data_processed$has_excess_alcohol, na.rm = TRUE),
    meets_inclusion_n
  ),
  Percentage_of_Total = c(
    round(100 * sum(data_processed$has_t2dm, na.rm = TRUE) / total_n, 2),
    round(100 * sum(data_processed$has_complicated_obesity, na.rm = TRUE) / total_n, 2),
    round(100 * sum(data_processed$has_metabolic_syndrome, na.rm = TRUE) / total_n, 2),
    round(100 * sum(data_processed$has_excess_alcohol, na.rm = TRUE) / total_n, 2),
    round(100 * meets_inclusion_n / total_n, 2)
  )
)

# Exclusion criteria counts (among those meeting inclusion)
exclusion_summary <- data.frame(
  Criterion = c("Missing FIB-4 Components",
                "FIB-4 > 2.67",
                "Competing Liver Disease",
                "  - Chronic Viral Hepatitis",
                "  - Autoimmune/Cholestatic Liver Disease",
                "  - Metabolic Inherited Liver Disorders",
                "  - Toxic Liver Disease",
                "  - Budd-Chiari Syndrome",
                "  - Cirrhosis/Fibrosis",
                "  - Hepatic Decompensation",
                "  - Hepatocellular Carcinoma",
                "  - Liver Transplant",
                "Any Exclusion (Among Eligible)"),
  N = c(
    sum(data_processed$meets_inclusion & data_processed$exclude_missing_fib4, na.rm = TRUE),
    sum(data_processed$meets_inclusion & data_processed$exclude_high_fib4, na.rm = TRUE),
    sum(data_processed$meets_inclusion & data_processed$exclude_competing_disease, na.rm = TRUE),
    sum(data_processed$meets_inclusion & data_processed$has_viral_hepatitis, na.rm = TRUE),
    sum(data_processed$meets_inclusion & data_processed$has_autoimmune_liver, na.rm = TRUE),
    sum(data_processed$meets_inclusion & data_processed$has_metabolic_liver, na.rm = TRUE),
    sum(data_processed$meets_inclusion & data_processed$has_toxic_liver, na.rm = TRUE),
    sum(data_processed$meets_inclusion & data_processed$has_budd_chiari, na.rm = TRUE),
    sum(data_processed$meets_inclusion & data_processed$has_cirrhosis, na.rm = TRUE),
    sum(data_processed$meets_inclusion & data_processed$has_decompensation, na.rm = TRUE),
    sum(data_processed$meets_inclusion & data_processed$has_hcc, na.rm = TRUE),
    sum(data_processed$meets_inclusion & data_processed$has_liver_transplant, na.rm = TRUE),
    excluded_among_eligible_n
  ),
  Percentage_of_Eligible = c(
    round(100 * sum(data_processed$meets_inclusion & data_processed$exclude_missing_fib4, na.rm = TRUE) / meets_inclusion_n, 2),
    round(100 * sum(data_processed$meets_inclusion & data_processed$exclude_high_fib4, na.rm = TRUE) / meets_inclusion_n, 2),
    round(100 * sum(data_processed$meets_inclusion & data_processed$exclude_competing_disease, na.rm = TRUE) / meets_inclusion_n, 2),
    round(100 * sum(data_processed$meets_inclusion & data_processed$has_viral_hepatitis, na.rm = TRUE) / meets_inclusion_n, 2),
    round(100 * sum(data_processed$meets_inclusion & data_processed$has_autoimmune_liver, na.rm = TRUE) / meets_inclusion_n, 2),
    round(100 * sum(data_processed$meets_inclusion & data_processed$has_metabolic_liver, na.rm = TRUE) / meets_inclusion_n, 2),
    round(100 * sum(data_processed$meets_inclusion & data_processed$has_toxic_liver, na.rm = TRUE) / meets_inclusion_n, 2),
    round(100 * sum(data_processed$meets_inclusion & data_processed$has_budd_chiari, na.rm = TRUE) / meets_inclusion_n, 2),
    round(100 * sum(data_processed$meets_inclusion & data_processed$has_cirrhosis, na.rm = TRUE) / meets_inclusion_n, 2),
    round(100 * sum(data_processed$meets_inclusion & data_processed$has_decompensation, na.rm = TRUE) / meets_inclusion_n, 2),
    round(100 * sum(data_processed$meets_inclusion & data_processed$has_hcc, na.rm = TRUE) / meets_inclusion_n, 2),
    round(100 * sum(data_processed$meets_inclusion & data_processed$has_liver_transplant, na.rm = TRUE) / meets_inclusion_n, 2),
    round(100 * excluded_among_eligible_n / meets_inclusion_n, 2)
  )
)

# Final cohort summary (FLOW DIAGRAM)
final_summary <- data.frame(
  Step = c("1. Total Participants After Data Cleaning",
           "2. Meeting ≥1 Inclusion Criterion",
           "3. Excluded (Missing FIB-4, FIB-4≥2.67, or Competing Disease)",
           "4. Final Analytical Cohort"),
  N = c(
    total_n,
    meets_inclusion_n,
    excluded_among_eligible_n,
    final_cohort_n
  ),
  Percentage_of_Previous_Step = c(
    100,
    round(100 * meets_inclusion_n / total_n, 2),
    round(100 * excluded_among_eligible_n / meets_inclusion_n, 2),
    round(100 * final_cohort_n / meets_inclusion_n, 2)
  ),
  Cumulative_Percentage = c(
    100,
    round(100 * meets_inclusion_n / total_n, 2),
    round(100 * (meets_inclusion_n - excluded_among_eligible_n) / total_n, 2),
    round(100 * final_cohort_n / total_n, 2)
  )
)

# Print tables
cat("\n=== COHORT SELECTION SUMMARY TABLES ===\n\n")

cat("STEP 1: INCLUSION CRITERIA\n")
print(inclusion_summary, row.names = FALSE)

cat("\n\nSTEP 2: EXCLUSION CRITERIA (Among those meeting inclusion)\n")
print(exclusion_summary, row.names = FALSE)

cat("\n\nCOHORT FLOW DIAGRAM:\n")
print(final_summary, row.names = FALSE)

# Verification
cat("\n\nVERIFICATION:\n")
cat("Meeting inclusion:", meets_inclusion_n, "\n")
cat("Excluded among eligible:", excluded_among_eligible_n, "\n")
cat("Final cohort:", final_cohort_n, "\n")
cat("Sum check:", excluded_among_eligible_n + final_cohort_n, 
    "(should equal", meets_inclusion_n, ")\n")

# Create final cohort dataset
final_cohort <- data_processed %>%
  filter(included_in_cohort)

cat("  Final N =", nrow(final_cohort), "\n\n")

#===============================================================================
# Define Outcomes: MALOs (Major Adverse Liver Outcomes)
# Include 6-month washout period from baseline assessment
#===============================================================================

# Main helper function to get first date after washout of occurrence
get_first_date_after_washout <- function(date_string, icd_string, pattern, assessment_date, washout_days = 183) {
  # Handle NA inputs
  if(is.na(date_string) || is.na(icd_string) || date_string == "" || icd_string == "") {
    return(as.Date(NA_character_))
  }
  
  dates <- strsplit(as.character(date_string), ",")[[1]]
  codes <- strsplit(as.character(icd_string), ",")[[1]]
  
  if(length(dates) == 0 || length(codes) == 0) return(as.Date(NA_character_))
  
  # Find matching codes
  matches <- grepl(pattern, codes)
  if(!any(matches)) return(as.Date(NA_character_))
  
  # Get corresponding dates
  matched_dates <- suppressWarnings(as.Date(dates[matches]))
  washout_end <- as.Date(assessment_date) + washout_days
  
  # Filter to valid dates after washout
  valid_dates <- matched_dates[!is.na(matched_dates) & matched_dates > washout_end]
  
  if(length(valid_dates) == 0) return(as.Date(NA_character_))
  
  return(min(valid_dates))
}

# HCC helper to check both cancer registry and hospital records
get_hcc_date <- function(icd10_dates_all, icd10_codes, cancer_dates_all, cancer_icd10_all, assessment_date) {
  # Initialize as Date NA
  hcc_hospital <- as.Date(NA_character_)
  hcc_cancer <- as.Date(NA_character_)
  
  # Try hospital records
  if(!is.na(icd10_dates_all) && !is.na(icd10_codes)) {
    hcc_hospital <- get_first_date_after_washout(icd10_dates_all, icd10_codes, 
                                                 "C22\\.0", assessment_date)
  }
  
  # Try cancer registry
  if(!is.na(cancer_dates_all) && !is.na(cancer_icd10_all)) {
    hcc_cancer <- get_first_date_after_washout(cancer_dates_all, cancer_icd10_all, 
                                               "C22\\.0", assessment_date)
  }
  
  # Combine results
  all_dates <- c(hcc_hospital, hcc_cancer)
  all_dates <- all_dates[!is.na(all_dates)]
  
  if(length(all_dates) == 0) {
    return(as.Date(NA_character_))
  }
  
  return(min(all_dates))
}

# Ascites helper - need to check for empty strings
get_ascites_date <- function(icd10_dates_all, icd10_codes, assessment_date) {
  if(is.na(icd10_dates_all) || is.na(icd10_codes) || 
     icd10_dates_all == "" || icd10_codes == "") {
    return(as.Date(NA_character_))
  }
  
  dates <- strsplit(as.character(icd10_dates_all), ",")[[1]]
  codes <- strsplit(as.character(icd10_codes), ",")[[1]]
  
  # Find R18 codes
  r18_idx <- which(grepl("R18", codes))
  if(length(r18_idx) == 0) return(as.Date(NA_character_))
  
  r18_dates <- suppressWarnings(as.Date(dates[r18_idx]))
  r18_dates <- r18_dates[!is.na(r18_dates)]
  
  # Check if liver anchor exists
  liver_anchor_idx <- which(grepl("K70\\.|K76\\.0", codes))
  if(length(liver_anchor_idx) == 0) return(as.Date(NA_character_))
  
  liver_anchor_dates <- suppressWarnings(as.Date(dates[liver_anchor_idx]))
  liver_anchor_dates <- liver_anchor_dates[!is.na(liver_anchor_dates)]
  
  # Find first R18 after washout
  washout_end <- as.Date(assessment_date) + 183
  valid_r18 <- r18_dates[r18_dates > washout_end]
  if(length(valid_r18) == 0) return(as.Date(NA_character_))
  
  # For each valid R18, check if there's a liver anchor on or before it
  has_anchor <- sapply(valid_r18, function(r18_date) {
    any(liver_anchor_dates <= r18_date)
  })
  
  valid_with_anchor <- valid_r18[has_anchor]
  
  if(length(valid_with_anchor) == 0) {
    return(as.Date(NA_character_))
  }
  
  return(min(valid_with_anchor))
}

# Liver death helper
get_liver_death_date <- function(date_of_death, death_cause_primary, death_cause_secondary, assessment_date) {
  # Check for NA or empty first
  if(length(date_of_death) == 0) return(as.Date(NA_character_))
  if(is.na(date_of_death)) return(as.Date(NA_character_))
  if(as.character(date_of_death) == "") return(as.Date(NA_character_))
  
  death_date <- as.Date(date_of_death)
  if(is.na(death_date)) return(as.Date(NA_character_))
  
  # Only count if after washout
  washout_end <- as.Date(assessment_date) + 183
  if(death_date <= washout_end) return(as.Date(NA_character_))
  
  # Check if primary or secondary cause is liver-related
  primary <- ifelse(is.na(death_cause_primary), "", as.character(death_cause_primary))
  secondary <- ifelse(is.na(death_cause_secondary), "", as.character(death_cause_secondary))
  
  liver_pattern <- "K70\\.|K74\\.|K76\\.|I85\\.|I86\\.4|K72\\.|C22\\.|R18"
  
  is_liver_death <- grepl(liver_pattern, primary) | grepl(liver_pattern, secondary)
  
  if(is_liver_death) {
    return(death_date)
  }
  
  return(as.Date(NA_character_))
}

#===============================================================================
# Make the FINAL COHORT WITH OUTCOMES
#===============================================================================

final_cohort_with_outcomes <- final_cohort %>%

  rowwise() %>%
  mutate(
    
    # Washout end date (6 months = 183 days after assessment)
    washout_end = assessment_date + 183,
    
    # =========================================================================
    # INCIDENT CIRRHOSIS
    # =========================================================================
    date_cirrhosis = get_first_date_after_washout(
      icd10_dates_all, icd10_codes, 
      "K70\\.2|K70\\.3|K74\\.0|K74\\.1|K74\\.2|K74\\.6|K76\\.6|I85\\.9|I86\\.4", 
      assessment_date
    ),
    
    # =========================================================================
    # DECOMPENSATION
    # =========================================================================
    
    # Ascites (R18) - only count if liver anchor code exists
    date_ascites = get_ascites_date(icd10_dates_all, icd10_codes, assessment_date),
    
    # Hepatorenal syndrome
    date_hrs = get_first_date_after_washout(icd10_dates_all, icd10_codes, 
                                            "K76\\.7", assessment_date),
    
    # Hepatic encephalopathy
    date_he = get_first_date_after_washout(icd10_dates_all, icd10_codes, 
                                           "K70\\.4|K72\\.0|K72\\.1|K72\\.9", 
                                           assessment_date),
    
    # Acute variceal bleeding
    date_avb = get_first_date_after_washout(icd10_dates_all, icd10_codes, 
                                            "I85\\.0", assessment_date),
    
    # Alcoholic hepatitis / ALD
    date_ald = get_first_date_after_washout(icd10_dates_all, icd10_codes,
                                            "K70\\.1|K70\\.9", assessment_date),
    
    # Any decompensation (earliest of ascites, HRS, HE, AVB, ALD)
    date_decompensation = {
      dates <- c(date_ascites, date_hrs, date_he, date_avb, date_ald)
      dates <- dates[!is.na(dates)]
      if(length(dates) == 0) NA else min(dates)
    },
    
    # =========================================================================
    # HEPATOCELLULAR CARCINOMA (HCC)
    # =========================================================================
    date_hcc = get_hcc_date(icd10_dates_all, icd10_codes, 
                            cancer_dates_all, cancer_icd10_all, 
                            assessment_date),
    
    # =========================================================================
    # LIVER TRANSPLANT
    # =========================================================================
    date_liver_transplant = get_first_date_after_washout(
      icd10_dates_all, icd10_codes, 
      "Z94\\.4", assessment_date
    ),
    
    # =========================================================================
    # LIVER-RELATED DEATH
    # =========================================================================
    date_liver_death = get_liver_death_date(
      date_of_death, death_cause_primary_i0, 
      death_cause_secondary_all, assessment_date
    ),
    
    # =========================================================================
    # COMPOSITE MALO OUTCOME
    # =========================================================================
    date_malo = {
      dates <- c(date_cirrhosis, date_decompensation, date_hcc, 
                 date_liver_transplant, date_liver_death)
      dates <- dates[!is.na(dates)]
      if(length(dates) == 0) NA else min(dates)
    },
    
    # =========================================================================
    # NON-LIVER DEATH (COMPETING EVENT)
    # =========================================================================
    date_non_liver_death = {
      if(is.na(date_of_death)) {
        NA
      } else {
        death_date <- as.Date(date_of_death)
        
        # Only count if after washout and NOT a liver death
        if(death_date <= (assessment_date + 183)) {
          NA
        } else if(!is.na(date_liver_death)) {
          NA
        } else {
          death_date
        }
      }
    },
    
    # =========================================================================
    # EVENT FLAGS
    # =========================================================================
    event_cirrhosis = !is.na(date_cirrhosis),
    event_decompensation = !is.na(date_decompensation),
    event_hcc = !is.na(date_hcc),
    event_liver_transplant = !is.na(date_liver_transplant),
    event_liver_death = !is.na(date_liver_death),
    event_malo = !is.na(date_malo),
    event_non_liver_death = !is.na(date_non_liver_death),
    
    # =========================================================================
    # TIME TO EVENT (in years)
    # =========================================================================
    
    # Censoring date (last known date or end of follow-up)
    date_censor = {
      # Use latest available date (death or last contact)
      if(!is.na(date_of_death)) {
        as.Date(date_of_death)
      } else {
        # Use date when UK Biobank released data
        as.Date("2026-02-10")
      }
    },
    
    # Time to MALO or censoring (in years from assessment)
    time_to_event = {
      if(event_malo) {
        as.numeric(difftime(date_malo, assessment_date, units = "days")) / 365.25
      } else if(event_non_liver_death) {
        as.numeric(difftime(date_non_liver_death, assessment_date, units = "days")) / 365.25
      } else {
        as.numeric(difftime(date_censor, assessment_date, units = "days")) / 365.25
      }
    }
    
  ) %>%
  ungroup()

# =========================================================================
# OUTCOME SUMMARY
# =========================================================================

# MALO Events Summary
malo_summary <- data.frame(
  Outcome = c(
    "Incident Cirrhosis",
    "Decompensation (Any)",
    "  - Ascites",
    "  - Hepatorenal Syndrome (HRS)",
    "  - Hepatic Encephalopathy (HE)",
    "  - Acute Variceal Bleeding (AVB)",
    "  - Alcoholic Liver Disease (ALD)",
    "Hepatocellular Carcinoma (HCC)",
    "Liver Transplant",
    "Liver-Related Death",
    "Composite MALO"
  ),
  N_Events = c(
    sum(final_cohort_with_outcomes$event_cirrhosis, na.rm = TRUE),
    sum(final_cohort_with_outcomes$event_decompensation, na.rm = TRUE),
    sum(!is.na(final_cohort_with_outcomes$date_ascites)),
    sum(!is.na(final_cohort_with_outcomes$date_hrs)),
    sum(!is.na(final_cohort_with_outcomes$date_he)),
    sum(!is.na(final_cohort_with_outcomes$date_avb)),
    sum(!is.na(final_cohort_with_outcomes$date_ald)),
    sum(final_cohort_with_outcomes$event_hcc, na.rm = TRUE),
    sum(final_cohort_with_outcomes$event_liver_transplant, na.rm = TRUE),
    sum(final_cohort_with_outcomes$event_liver_death, na.rm = TRUE),
    sum(final_cohort_with_outcomes$event_malo, na.rm = TRUE)
  ),
  Percentage = c(
    round(100 * sum(final_cohort_with_outcomes$event_cirrhosis, na.rm = TRUE) / nrow(final_cohort_with_outcomes), 2),
    round(100 * sum(final_cohort_with_outcomes$event_decompensation, na.rm = TRUE) / nrow(final_cohort_with_outcomes), 2),
    round(100 * sum(!is.na(final_cohort_with_outcomes$date_ascites)) / nrow(final_cohort_with_outcomes), 2),
    round(100 * sum(!is.na(final_cohort_with_outcomes$date_hrs)) / nrow(final_cohort_with_outcomes), 2),
    round(100 * sum(!is.na(final_cohort_with_outcomes$date_he)) / nrow(final_cohort_with_outcomes), 2),
    round(100 * sum(!is.na(final_cohort_with_outcomes$date_avb)) / nrow(final_cohort_with_outcomes), 2),
    round(100 * sum(!is.na(final_cohort_with_outcomes$date_ald)) / nrow(final_cohort_with_outcomes), 2),
    round(100 * sum(final_cohort_with_outcomes$event_hcc, na.rm = TRUE) / nrow(final_cohort_with_outcomes), 2),
    round(100 * sum(final_cohort_with_outcomes$event_liver_transplant, na.rm = TRUE) / nrow(final_cohort_with_outcomes), 2),
    round(100 * sum(final_cohort_with_outcomes$event_liver_death, na.rm = TRUE) / nrow(final_cohort_with_outcomes), 2),
    round(100 * sum(final_cohort_with_outcomes$event_malo, na.rm = TRUE) / nrow(final_cohort_with_outcomes), 2)
  )
)

# Competing Events Summary
competing_summary <- data.frame(
  Event = "Non-Liver Death",
  N_Events = sum(final_cohort_with_outcomes$event_non_liver_death, na.rm = TRUE),
  Percentage = round(100 * sum(final_cohort_with_outcomes$event_non_liver_death, na.rm = TRUE) / nrow(final_cohort_with_outcomes), 2)
)

# Follow-up Summary
followup_summary <- data.frame(
  Metric = c("Total Cohort", "Median Follow-up (years)", "Mean Follow-up (years)", "Max Follow-up (years)"),
  Value = c(
    nrow(final_cohort_with_outcomes),
    round(median(final_cohort_with_outcomes$time_to_event, na.rm = TRUE), 2),
    round(mean(final_cohort_with_outcomes$time_to_event, na.rm = TRUE), 2),
    round(max(final_cohort_with_outcomes$time_to_event, na.rm = TRUE), 2)
  )
)

# Print summaries
cat("\n=== OUTCOME SUMMARY ===\n\n")
cat("MALO EVENTS:\n")
print(malo_summary, row.names = FALSE)

cat("\n\nCOMPETING EVENTS:\n")
print(competing_summary, row.names = FALSE)

cat("\n\nFOLLOW-UP:\n")
print(followup_summary, row.names = FALSE)

#=============================================
# SAVE ALL FILES
#=============================================
 
summary_list <- list("inclusion" = inclusion_summary,
                     "exclusion" = exclusion_summary,
                     "final_cohort" = final_summary,
                     "MALO" = malo_summary,
                     "competing_event" = competing_summary,
                     "follow_up" = followup_summary)
write_xlsx(summary_list, file.path(data_dir, "BFA_cohort_summaries.xlsx"))
write_feather(final_cohort_with_outcomes, file.path(data_dir, "BFA_principal_data.feather"))

cat("=============== DONE ===============\n\n")
