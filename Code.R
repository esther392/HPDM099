# In command prompt type: ssh slade

# Setup
library(tidyverse)
library(aurum)
library(lubridate)
library(dplyr)
library(tidyr)
library(gtsummary)
library(table1)
library(ggplot2)
library(data.table)
library(Hmisc)
library(patchwork)

rm(list=ls())

#:--------------------------------------
# load t2d_1stinstance

# Select columns needed and save as new table
t2d_selected <- t2d_1stinstance %>%
  select(patid, gender, dob, dstartdate, dstartdate_age, dstartdate_dm_dur_all,
         preweight, death_date, diabetes_type, ethnicity_5cat, predbp, presbp, 
         prebmi, preegfr, prealt, prehdl, qrisk2_10yr_score, predrug_primary_hhf, 
         predrug_stroke, predrug_latest_stroke, predrug_tia, predrug_latest_tia, 
         predrug_latest_primary_hhf, predrug_primary_incident_mi, predrug_latest_primary_incident_mi,
         predrug_primary_incident_stroke, predrug_latest_primary_incident_stroke,
         predrug_myocardialinfarction, predrug_latest_myocardialinfarction, 
         predrug_revasc, predrug_latest_revasc, predrug_angina, predrug_latest_angina,
         predrug_amputation, predrug_efi_peripheral_vascular_disease, prefastingglucose, 
         prefastingglucosedate, predrug_efi_thyroid_disease, Empagliflozin, Acarbose, 
         DPP4, Glinide, GIPGLP1, GLP1, MFN, SGLT2, SU, TZD, INS, prehba1c, 
         predrug_haem_cancer, predrug_latest_haem_cancer, predrug_solid_cancer, 
         predrug_latest_solid_cancer, alcohol_cat, drug_class, drug_class_start, 
         dstopdate_class, drug_order, drugline, drug_instance, drug_substance, 
         dstopdate_substance, imd_decile, pretotalcholesterol, pretriglyceride,
         prealbumin_blood, precreatinine_blood,
         predrug_latest_oralsteroids, predrug_medspecific_gi, predrug_latest_medspecific_gi, 
         predrug_latest_prodspecific_gi, predrug_unspecific_gi, predrug_latest_unspecific_gi,
         predrug_af, predrug_latest_af, predrug_hypertension, predrug_latest_hypertension,
         predrug_ihd, predrug_latest_ihd, predrug_pad, predrug_latest_pad, predrug_heartfailure, 
         predrug_latest_heartfailure) %>%
  mutate(dstartdate = as.Date(dstartdate)) %>%
  analysis$cached("test") # name: ys_test; path: cprd_jun24dm_analysis.ys_test

#View(t2d_selected %>% head(100) %>% collect())
#2893016

# Made a new column for patients with drug_order == 2 & afterwards into one individual
t2d_selected2 <- t2d_selected %>%
  mutate(
    dstartdate = as.Date(dstartdate),
    dstopdate_substance = as.Date(dstopdate_substance),
    duration_days = as.numeric(dstopdate_substance - dstartdate), # Treatment duration
    new_patid = paste0(patid, "_", drug_order) # Create unique ID for each drug orders per patient
  ) %>%
  relocate(new_patid, .after = patid) %>%
  relocate(duration_days, .after = new_patid)


#:--------------------------------------
# TARGET POPULATION
# Select patients who could potentially benefit from SGLT2i

# 1. Keep patients who started second-line treatment in and after 2017
# Identify patients who start 2nd-line BEFORE 2017
exclude_patids <- t2d_selected2 %>%
  filter(drug_order == 2,
         dstartdate < as.Date("2017-01-01")) %>%
  distinct(new_patid) %>%
  collect()

# Exclude those patients from the dataset
t2d_selected3 <- t2d_selected2 %>%
  filter(!new_patid %in% exclude_patids$new_patid)

# Get the numbers of patients excluding 1st-line treatment
t2d_selected3 %>% filter(drug_order != 1) %>% count()
#1333576

## MIGHT NEED SOME TIME TO RUN
# Bring relevant columns into local memory for faster processing
t2d_selected3_local <- t2d_selected3 %>% collect()

## MIGHT NEED SOME TIME TO RUN
# Create string summaries of all drug start/stop dates per patient
dstartstop_summary <- t2d_selected3_local %>%
  group_by(patid) %>%
  summarise(
    all_dstartdates = paste(sort(as.character(dstartdate)), collapse = ","),
    all_dstopdates = paste(sort(as.character(dstopdate_class)), collapse = ","),
    .groups = "drop"
  )

# Merge treatment timeline back to main dataset
t2d_startstop <- t2d_selected3_local %>%
  left_join(dstartstop_summary, by = "patid")

# Remove all 1st-line drug
t2d_no_drg1 <- t2d_startstop %>%
  filter(drug_order != 1)
#1333576

# 2. Exclude Patients who take Insulin, Acarbose, Glinide & Metformin
# Identify new_patids with unwanted drugs
exclude_patids2 <- t2d_no_drg1 %>%
  filter(drug_class %in% c("Acarbose", "Glinide", "MFN", "INS")) %>%
  distinct(new_patid) %>%
  collect()

# Filter out those patients from full dataset
t2d_selected4 <- t2d_no_drg1 %>%
  filter(!new_patid %in% exclude_patids2$new_patid)
#1048823

## MIGHT NEED SOME TIME TO RUN
# 3. If multiple records per new_patid, keep only the one with longest treatment duration
t2d_selected5 <- t2d_selected4 %>%
  group_by(new_patid) %>%
  slice_max(order_by = duration_days, n = 1, with_ties = FALSE) %>%
  ungroup()
#943774

# 4. Keep patients with eGFR >= 30 (not contraindication to SGLT2i)
t2d_egfr_filtered <- t2d_selected5 %>%
  filter(preegfr >= 30)
#913083

# 5. Define patients with cardiovascular risk (per NICE guidelines)
# Add a binary flag for patients meeting CVD risk criteria:
## Includes established CVD, heart failure, or QRISK2 score ≥ 10
t2d_cvd_flagged <- t2d_egfr_filtered %>%
  mutate(
    cvd_risk_flag = (
      # A. Established CVD
      (predrug_angina == 1 & !is.na(predrug_latest_angina) & predrug_latest_angina < dstartdate) |
        (predrug_af == 1 & !is.na(predrug_latest_af) & predrug_latest_af < dstartdate) |
        (predrug_hypertension == 1 & !is.na(predrug_latest_hypertension) & predrug_latest_hypertension < dstartdate) |
        (predrug_revasc == 1 & !is.na(predrug_latest_revasc) & predrug_latest_revasc < dstartdate) |
        (predrug_ihd == 1 & !is.na(predrug_latest_ihd) & predrug_latest_ihd < dstartdate) |
        (predrug_myocardialinfarction == 1 & !is.na(predrug_latest_myocardialinfarction) & predrug_latest_myocardialinfarction < dstartdate) |
        (predrug_pad == 1 & !is.na(predrug_latest_pad) & predrug_latest_pad < dstartdate) |
        (predrug_stroke == 1 & !is.na(predrug_latest_stroke) & predrug_latest_stroke < dstartdate) |
        (predrug_tia == 1 & !is.na(predrug_latest_tia) & predrug_latest_tia < dstartdate) |
        
        # B. Heart failure
        (predrug_heartfailure == 1 & !is.na(predrug_latest_heartfailure) & predrug_latest_heartfailure < dstartdate) |
        
        # C. High QRISK2 score
        (!is.na(qrisk2_10yr_score) & qrisk2_10yr_score >= 10)
    )
  )

# Keep only patients with at least one CVD risk factor
t2d_cvd_risk <- t2d_cvd_flagged %>%
  filter(cvd_risk_flag == TRUE)
#819784


#:--------------------------------------
# TRIAL ELIGIBLE
# Code for Exclusion criteria
# Exclude patients with age < 18 and BMI > 45
t2d_exclusion1 <- t2d_cvd_risk %>%
  filter(!(dstartdate_age < 18 & prebmi > 45))
#819780

# Exclude patients who DO NOT have eligible cardiovascular conditions (Accordin g to Section C)
t2d_no_cvd <- t2d_exclusion1 %>%
  filter(
    !(
      # A. Myocardial infarction more than 2 months prior to drug start
      (predrug_myocardialinfarction == 1 &
         as.numeric(dstartdate - predrug_latest_myocardialinfarction) > 60) |
        
        # B. Coronary revascularization > 2 months ago
        (predrug_revasc == 1 &
           as.numeric(dstartdate - predrug_latest_revasc) > 60) |
        
        # C. Recent hospital discharge for major CV events (within past 12 months)
        (predrug_primary_hhf == 1 &
           as.numeric(dstartdate - predrug_latest_primary_hhf) <= 365) |
        (predrug_primary_incident_mi == 1 &
           as.numeric(dstartdate - predrug_latest_primary_incident_mi) <= 365) |
        (predrug_primary_incident_stroke == 1 &
           as.numeric(dstartdate - predrug_latest_primary_incident_stroke) <= 365) |
        
        # D. Unstable angina > 2 months prior
        (predrug_angina == 1 &
           as.numeric(dstartdate - predrug_latest_angina) > 60) |
        
        # E. Stroke > 2 months prior
        (predrug_stroke == 1 &
           as.numeric(dstartdate - predrug_latest_stroke) > 60) |
        
        # F. Limb amputation or peripheral artery disease
        predrug_amputation == 1 |
        predrug_efi_peripheral_vascular_disease == 1
    )
  )
#577426


# Exclude patients that:
## (1) Had a previous drug start or stop within 12 weeks before the current treatment start date &
## (2) Had an HbA1c percent value NOT between 7.0 and 10.0 at the current treatment start
# Also remove patient who only have 1 drugline (only one date in all_dstartdates/all_dstopdates)

# Convert dataset to data.table
dt <- as.data.table(t2d_no_cvd)

# Count how many treatment start/stop dates each patient had (based on comma-separated strings)
dt[, n_start_dates := lengths(str_split(all_dstartdates, ","))]
dt[, n_stop_dates  := lengths(str_split(all_dstopdates, ","))]

# Keep only patients who had more than one drug line (i.e. more than one start and stop date)
dt <- dt[n_start_dates > 1 & n_stop_dates > 1]

# Convert comma-separated start/stop date strings into proper Date vectors
dt_unique <- unique(dt[, .(patid, all_dstartdates, all_dstopdates)])
## MIGHT NEED SOME TIME TO RUN
dt_unique[, all_start_dates := lapply(all_dstartdates, function(x) as.Date(str_split(x, ",")[[1]]))]
dt_unique[, all_stop_dates  := lapply(all_dstopdates,  function(x) as.Date(str_split(x, ",")[[1]]))]

# Flag patients who started a new drug within 12 weeks of a prior start/stop date
within_12weeks <- function(date_vec, ref_date) {
  diffs <- as.numeric(ref_date - date_vec)
  any(diffs > 0 & diffs <= 84) # 12 weeks = 84 days
}

# Apply flag logic for each patient
dt_unique[, exclude_due_to_timing := {
  sdates <- all_start_dates[[1]]  # all start dates for this patient
  stdops <- all_stop_dates[[1]] # all stop dates for this patient
  
  any(sapply(2:length(sdates), function(i) {
    cur_start <- sdates[i]   # current start date
    prev_starts <- sdates[1:(i - 1)] # all prior start dates
    prev_stops  <- stdops[1:(i - 1)] # all prior stop dates
    
    # Exclude if current start is within 12 weeks of any prior start or stop
    within_12weeks(prev_starts, cur_start) || within_12weeks(prev_stops, cur_start)
  }))
}, by = patid]

# Identify patients to exclude due to recent prior treatment
exclude_timing <- dt_unique[exclude_due_to_timing == TRUE, patid]

# Convert HbA1c mmol/mol to percent
dt[, hba1c_percent := (prehba1c / 10.929) + 2.15]

# Final filter — exclude rows for patients who meet both:
# (a) Were flagged due to timing overlap AND
# (b) Had an HbA1c value outside 7.0–10.0%
result_dt <- dt[!(patid %in% exclude_timing & (hba1c_percent < 7.0 | hba1c_percent > 10.0))]
#420440


# 1.Exclude Uncontrolled hyperglycemia with glucose > 240 mg/dL
# Convert to mg/dL, filter > 240, and remove missing dates
t2d_hyperglycemia_exclude <- result_dt %>%
  mutate(
    prefastingglucose_mgdl = prefastingglucose * 18.0182,
    prefastingglucosedate = as.Date(prefastingglucosedate)
  ) %>%
  filter(
    prefastingglucose_mgdl > 240, 
    !is.na(prefastingglucosedate) # Exclude missing dates
  ) %>% 
  distinct(patid, prefastingglucosedate) %>% # Keep one row per patient per date
  group_by(patid) %>%
  filter(n() >= 2) %>% # Keep only patients with 2+ high glucose readings on different dates
  ungroup()

# Get list of patients to exclude
exclude_ids <- t2d_hyperglycemia_exclude %>%
  distinct(patid) %>%
  collect()  # Pull into memory

# Filter them out from the main dataset
t2d_exclusion3 <- result_dt %>%
  filter(!patid %in% exclude_ids$patid)
#419553

# 2. Exclude patients with elevated ALT (>100 IU/L) indicating possible liver dysfunction
# U/L * 1 = IU/L
t2d_exclude_alt <- t2d_exclusion3 %>%
  filter(prealt <= 100 | is.na(prealt))
#414919

# 3. Exclude patients with bariatric surgery within 2 years before drug start
bariatric_data <- bariatric_data %>% 
  analysis$cached("bariatric_data")

# Get most recent bariatric surgery end date per patient
bariatric_latest <- bariatric_data %>%
  group_by(patid) %>%
  summarise(latest_epiend = max(epiend, na.rm = TRUE)) %>%
  ungroup()

t2d_exclude_bariatric <- t2d_exclude_alt %>%
  left_join(bariatric_latest %>%
              rename(latest_bariatric = latest_epiend) %>%
              select(patid, latest_bariatric) %>%
              collect() ,
            by = "patid") %>%
  # Calculate time difference in days
  mutate(
    days_since_bariatric = as.numeric(dstartdate - latest_bariatric)
    ) %>%
  # Keep if no bariatric data (NA in epiend) or surgery ended more than 2 years (730 days) before dstartdate
  filter(
    is.na(days_since_bariatric) | days_since_bariatric > 730
    ) 
#412932

# 4. Exclude patients with history of blood dyscrasias
dyscrasias_data <- dyscrasias_data %>%
  analysis$cached("dycrasia_data") %>%
  collect()

# Get unique patients with dyscrasias observations
patients_with_dyscrasias <- dyscrasias_data %>%
  distinct(patid)

# Exclude these patients
t2d_exclude_dyscrasias <- t2d_exclude_bariatric %>%
  anti_join(patients_with_dyscrasias, by = "patid")
#403555

# 5.Exclude patients with history of cancer within last 5 years
t2d_exclude_cancer <- t2d_exclude_dyscrasias %>%
  filter(
    is.na(predrug_latest_haem_cancer) | dstartdate - predrug_latest_haem_cancer >= 5 * 365, 
    is.na(predrug_latest_solid_cancer) | dstartdate - predrug_latest_solid_cancer >= 5 * 365
  )
#381334

# 6. Exclude patients:
## take high-dose semaglutide --> antiobesity treatment
## their following treatment is within 3 months
# Split new_patid into patid and drug_order
t2d_split2 <- t2d_exclude_cancer %>%
  mutate(
    patid_extracted = sub("_.*", "", new_patid),
    drug_order_extracted = as.integer(sub(".*_", "", new_patid))
  )

# Find rows where drug substance is high-dose semaglutide
high_dose_rows <- t2d_split2 %>%
  filter(drug_substance == "High-dose semaglutide") %>%
  select(
    patid_extracted,
    drug_order_extracted,
    high_dose_start = dstartdate,
    new_patid
  )

# Find the very next drug order row for each high-dose semaglutide row
next_rows <- high_dose_rows %>%
  mutate(next_order = drug_order_extracted + 1) %>% # Next drug order after the high-dose semaglutide
  left_join(
    t2d_split2,
    by = c("patid_extracted" = "patid_extracted", "next_order" = "drug_order_extracted"),
    suffix = c("_highdose", "_next")
  ) %>%
  filter(!is.na(dstartdate)) %>%
  mutate(
    days_between = as.numeric(dstartdate - high_dose_start), # Days difference between next treatment start and semaglutide start
    within_3_months = days_between >= 0 & days_between <= 84 # Flag if next treatment starts within 84 days (3 months)
  ) %>%
  filter(within_3_months) # Keep only next treatments that started within 3 months

# Combine rows to exclude: both the high-dose semaglutide rows AND their immediate next rows starting within 3 months
rows_to_exclude <- bind_rows(
  high_dose_rows %>% select(new_patid), # high-dose semaglutide rows
  next_rows %>% select(new_patid = new_patid_next) # immediate next rows (note: renamed new_patid_next)
) %>% distinct() # Remove any duplicates

# Exclude both the high-dose and the immediate follow-up rows
t2d_exclude_semaglutide <- t2d_split2 %>%
  filter(!new_patid %in% rows_to_exclude$new_patid)
#381317

# 7. Exclude if oral steroid use was within 3 months before drug start
t2d_exclude_steroids <- t2d_exclude_semaglutide %>%
  filter(
    is.na(predrug_latest_oralsteroids) | # Keep if no oral steroid prescription
      ((dstartdate - predrug_latest_oralsteroids) > 90) | # Keep if steroid use was more than 3 months ago
      (dstartdate < predrug_latest_oralsteroids) # Keep if steroid use is in the future
  )
#367022

# 8. Exclude patients classified as "Heavy" alcohol users
t2d_exclude_alcohol <- t2d_exclude_steroids %>%
  filter(alcohol_cat != "Heavy")
#339393

# 9. Exclude genito-urinary infection within 2 weeks prior to drug start
t2d_exclude_gi <- t2d_exclude_alcohol %>%
  mutate(
    gi_exclude = case_when(
      predrug_medspecific_gi == 1 & !is.na(predrug_latest_medspecific_gi) &
        (dstartdate - predrug_latest_medspecific_gi <= 14 &
           dstartdate - predrug_latest_medspecific_gi >= 0) ~ TRUE,
      
      !is.na(predrug_latest_prodspecific_gi) &
        (dstartdate - predrug_latest_prodspecific_gi <= 14 &
           dstartdate - predrug_latest_prodspecific_gi >= 0) ~ TRUE,
      
      predrug_unspecific_gi == 1 & !is.na(predrug_latest_unspecific_gi) &
        (dstartdate - predrug_latest_unspecific_gi <= 14 &
           dstartdate - predrug_latest_unspecific_gi >= 0) ~ TRUE,
      
      TRUE ~ FALSE
    )
  ) %>%
  filter(!gi_exclude) %>%
  select(-gi_exclude)
#337156

# 10. Exclude patients with stroke in the 2 months before treatment
t2d_exclude_stroke <- t2d_exclude_gi %>%
  filter(
    predrug_stroke == 0 |
      is.na(predrug_latest_stroke) |
      (dstartdate - predrug_latest_stroke) > 60
  )
#336234

# 11. Exclude TIA (transient ischemic attack) in the past 2 months
t2d_exclude_tia <- t2d_exclude_stroke %>%
  filter(
    predrug_tia == 0 |
      is.na(predrug_latest_tia) |
      (dstartdate - predrug_latest_tia) > 60
  )
#335779

# 12. Exclude if blood pressure is > 160 systolic OR > 100 diastolic
t2d_exclude_bp <- t2d_exclude_tia %>%
  filter(
    !(presbp > 160 | predbp > 100) 
  )
#317053


#:--------------------------------------
# COMPARISON TABLE: Generate Table 1 for Target, Included, and Excluded groups

# Derive df_target
df_target <- t2d_cvd_risk %>% 
  filter(drug_order != 1) %>% # Exclude 1st-line treatment, focus on 2nd-line or later
  mutate(new_patid = as.character(new_patid)) %>% # Ensure consistent ID type for joins
  collect()

# Define df_include (patients who meet trial eligibility criteria)
df_include <- t2d_exclude_bp %>% 
  mutate(new_patid = as.character(new_patid)) %>% 
  collect()

# Define df_exclude as those in df_target but not in df_include
df_exclude <- anti_join(df_target, df_include, by = "new_patid") # Excluded = not trial-eligible

# Define the preprocessing function for Table 1
prep_table1_data <- function(df) {
  df <- df %>%
    mutate(
      # Convert gender to factor with labels
      gender = factor(gender, levels = c(1 ,2), labels = c("Male", "Female")),
      
      # Recode 5-category ethnicity
      ethnicity = case_when(
        ethnicity_5cat == 0 ~ "White",
        ethnicity_5cat == 1 ~ "South Asian",
        ethnicity_5cat == 2 ~ "Black",
        ethnicity_5cat %in% c(3, 4) ~ "Other/Mixed",
        TRUE ~ "Unknown"
      ),
      ethnicity = factor(ethnicity, levels = c("White", "South Asian", "Black", "Other/Mixed", "Unknown")),
      
      # Assign diabetes duration
      t2d_duration = dstartdate_dm_dur_all,
      
      # Categorise treatment class initiated
      treatment_initiated = case_when(
        drug_class == "GIPGLP1" ~ "GLP1",
        drug_class %in% c("DPP4", "GLP1", "SGLT2", "SU", "TZD") ~ drug_class,
        TRUE ~ "Other/Missing"
      ),
      treatment_initiated = factor(treatment_initiated,
                                   levels = c("DPP4", "GLP1", "SGLT2", "SU", "TZD", "Other/Missing")),
      
      # Group IMD deciles into bins
      imd_decile_group = case_when(
        imd_decile %in% 1:2 ~ "1 or 2",
        imd_decile %in% 3:4 ~ "3 or 4",
        imd_decile %in% 5:6 ~ "5 or 6",
        imd_decile %in% 7:8 ~ "7 or 8",
        imd_decile %in% 9:10 ~ "9 or 10",
        TRUE ~ "Missing"
      ),
      imd_decile_group = factor(imd_decile_group,
                                levels = c("1 or 2", "3 or 4", "5 or 6", "7 or 8", "9 or 10", "Missing")),
      
      # Group drug_order into "4+" category
      drug_order_grouped = case_when(
        as.numeric(drug_order) <= 4 ~ as.character(drug_order),
        as.numeric(drug_order) > 4 ~ "4+",
        TRUE ~ NA_character_
      ),
      drug_order = factor(drug_order_grouped, levels = c("2", "3", "4", "4+")),
      
      # Categorise QRISK2 score
      qrisk2_category = case_when(
        qrisk2_10yr_score < 10 ~ "<10% (Low)",
        qrisk2_10yr_score >= 10 & qrisk2_10yr_score <= 20 ~ "10–20% (Moderate)",
        qrisk2_10yr_score > 20 ~ ">20% (High)",
        TRUE ~ NA_character_
      )
    )
  
  # Convert all predrug cardiovascular conditions to binary factors
  # 0 = no history, 1 = history present
  cvd_vars <- df %>% select(starts_with("predrug_")) %>% names()
  
  df <- df %>%
    mutate(across(all_of(cvd_vars), ~ factor(ifelse(is.na(.) | . == 0, 0, 1), levels = c(0, 1))))
  
  return(df)
}

## MIGHT NEED SOME TIME TO RUN
# Preprocess each group using the function above
df_target2  <- prep_table1_data(df_target)  %>% mutate(group = "Target")
df_include2 <- prep_table1_data(df_include) %>% mutate(group = "Included")
df_exclude2 <- prep_table1_data(df_exclude) %>% mutate(group = "Excluded")

# Combine all groups into one dataframe for Table 1
df_combined <- bind_rows(df_target2, df_include2, df_exclude2)

# Convert group to factor to set display order in Table 1
df_combined$group <- factor(df_combined$group, levels = c("Target", "Included", "Excluded"))

# Label variables for display in Table 1 (table1 package)
label(df_combined$dstartdate_age) <- "Age"
label(df_combined$gender) <- "Gender"
label(df_combined$ethnicity) <- "Ethnicity"
label(df_combined$t2d_duration) <- "Time since diagnosis of type 2 diabetes (years)"
label(df_combined$prebmi) <- "Body mass index (kg/m²)"
label(df_combined$preweight) <- "Weight (kg)"
label(df_combined$preegfr) <- "Estimated glomerular filtration rate (mL/min/1.73m²)"
label(df_combined$prehba1c) <- "Glycated hemoglobin (mmol/mol)"
label(df_combined$prehdl) <- "HDL cholesterol (mmol/L)"
label(df_combined$prealt) <- "Alanine aminotransferase (IU/L)"
label(df_combined$presbp) <- "Systolic blood pressure (mm Hg)"
label(df_combined$predbp) <- "Diastolic blood pressure (mm Hg)"
label(df_combined$pretotalcholesterol) <- "Total cholesterol (mmol/L)"
label(df_combined$pretriglyceride) <- "Triglycerides (mmol/L)"
label(df_combined$qrisk2_category) <- "QRISK2 category"
label(df_combined$treatment_initiated) <- "Treatment initiated"
label(df_combined$imd_decile_group) <- "IMD decile group"
label(df_combined$drug_order) <- "Drug order"

# Label CVD history variables
label(df_combined$predrug_angina) <- "Angina"
label(df_combined$predrug_af) <- "Atrial fibrillation"
label(df_combined$predrug_hypertension) <- "Hypertension"
label(df_combined$predrug_revasc) <- "Revascularisation"
label(df_combined$predrug_ihd) <- "Ischemic heart disease"
label(df_combined$predrug_myocardialinfarction) <- "Myocardial infarction"
label(df_combined$predrug_pad) <- "Peripheral arterial disease"
label(df_combined$predrug_stroke) <- "Stroke"
label(df_combined$predrug_tia) <- "Transient ischemic attack"
label(df_combined$predrug_heartfailure) <- "Heart failure"

# Define the variables to include in Table 1
variables_table1 <- c(
  "gender", "ethnicity", "treatment_initiated", "imd_decile_group", "drug_order",
  "dstartdate_age", "t2d_duration", "prebmi", "preweight", "preegfr", "prehba1c",
  "prehdl", "prealt", "presbp", "predbp", "pretotalcholesterol", "pretriglyceride",
  "predrug_angina", "predrug_af", "predrug_hypertension", "predrug_revasc",
  "predrug_ihd", "predrug_myocardialinfarction", "predrug_pad", "predrug_stroke",
  "predrug_tia", "predrug_heartfailure", "qrisk2_category"
)

# Create formula for the table1 function to compare by group
formula_data_description <- as.formula(
  paste("~", paste(variables_table1, collapse = " + "), "| group")
)

## MIGHT NEED SOME TIME TO RUN
# Generate the table using the 'table1' package
table1_output <- table1(
  formula_data_description,
  data = df_combined,
  droplevels = TRUE,
  render.continuous = "Mean (SD)",  # Show mean and standard deviation for continuous vars
  render.missing = NULL, # Do not show missing value counts separately
  overall = FALSE # Do not include an "Overall" column
)

# Print the table
print(table1_output)




#:------------------------------------------------------------------------------
# Visualization graphs

# Categorical Bar Chart
# Define a list of categorical variables to plot
categorical_vars <- c(
  "gender", "ethnicity", "treatment_initiated", "imd_decile_group", 
  "drug_order", "qrisk2_category"
)

# Define a function to create side-by-side bar charts for each categorical variable
plot_categorical <- function(var) {
  var_label <- Hmisc::label(df_combined[[var]]) # Retrieve human-readable label
  
  df_combined %>%
    count(group, !!sym(var)) %>% # Count occurrences by group and variable level
    ggplot(aes(x = !!sym(var), y = n, fill = group)) +
    geom_col(position = "dodge") +  # Use side-by-side bars
    labs(
      title = var_label,  # Use descriptive title
      x = NULL, # No x-axis label
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1), # Tilt x labels
      legend.title = element_blank()
    )
}

# Apply the plotting function to each variable and arrange results in a 2-column grid
cat_plots <- map(categorical_vars, plot_categorical)
wrap_plots(cat_plots, ncol = 2)



# Selected Continuous Boxplot
# List of selected continuous variables
cont_vars <- c("dstartdate_age", "prehba1c", "preegfr", "presbp")

# Map variable labels for display (fall back to name if no label)
label_lookup <- map_chr(cont_vars, ~ Hmisc::label(df_combined[[.x]]) %||% .x) |> set_names(cont_vars)

## MIGHT NEED SOME TIME TO RUN
# Prepare the dataset for plotting
df_cont <- df_combined %>%
  select(group, all_of(cont_vars)) %>%
  mutate(across(all_of(cont_vars), ~ unclass(.))) %>%  # Remove labels to allow pivot
  pivot_longer(-group, names_to = "variable", values_to = "value") %>%
  filter(is.finite(value)) %>%  # Keep only numeric and non-missing values
  mutate(label = label_lookup[variable]) # Attach label for each variable

# Create boxplots, faceted by variable label
ggplot(df_cont, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.alpha = 0.1) + # Lower outlier opacity for readability
  facet_wrap(~label, scales = "free_y") + # Allow separate y-scales
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = "none")




# Binary CVD Variables Stacked Bars
# All binary cardiovascular variables
cvd_vars <- c(
  "predrug_angina", "predrug_af", "predrug_hypertension",
  "predrug_revasc", "predrug_ihd", "predrug_myocardialinfarction",
  "predrug_pad", "predrug_stroke", "predrug_tia", "predrug_heartfailure"
)

# Retrieve variable labels for titles
cvd_labels <- sapply(cvd_vars, function(var) Hmisc::label(df_combined[[var]]))

# Convert labelled variables to numeric (or logical if preferred)
df_combined2 <- df_combined %>%
  mutate(across(all_of(cvd_vars), as.numeric)) # Avoids pivot errors

## MIGHT NEED SOME TIME TO RUN
# Reshape and compute prevalence
df_cvd_long <- df_combined2 %>%
  select(group, all_of(cvd_vars)) %>%
  pivot_longer(cols = -group, names_to = "cvd", values_to = "present") %>%
  group_by(group, cvd, present) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group, cvd) %>%
  mutate(
    pct = n / sum(n),
    cvd_label = cvd_labels[cvd]  # Add label back here
  ) %>%
  filter(present == 1) # Keep only those with the condition

# Plot prevalence of all CVD history conditions
ggplot(df_cvd_long, aes(x = reorder(cvd_label, -pct), y = pct, fill = group)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Prevalence of CVD history by group",
    x = NULL,
    y = "Prevalence (%)"
  ) +
  coord_flip() +
  theme_minimal() +
  theme(
    legend.title = element_blank()
  )






### MIGHT NEED LONGER TIME TO RUN
# Full Continuous Box plot
continuous_vars <- c(
  "dstartdate_age", "t2d_duration", "prebmi", "preweight", "preegfr",
  "prehba1c", "prehdl", "prealt", "presbp", "predbp",
  "pretotalcholesterol", "pretriglyceride"
)

# Plotting function for each continuous variable
plot_continuous <- function(var) {
  df_plot <- df_combined %>%
    filter(is.finite(!!sym(var)))  # Remove NA, NaN, Inf values
  
  var_label <- Hmisc::label(df_combined[[var]]) # Get label
  
  ggplot(df_plot, aes(x = group, y = !!sym(var), fill = group)) +
    geom_boxplot(outlier.alpha = 0.1) +
    labs(
      title = var_label,
      y = NULL,
      x = NULL
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
}

# Generate all plots and arrange in 4-column layout
cont_plots <- map(continuous_vars, plot_continuous)
wrap_plots(cont_plots, ncol = 4)  


