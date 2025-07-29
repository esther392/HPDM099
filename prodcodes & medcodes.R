# libraries
library(tidyverse)
library(dplyr)
library(aurum)
library(readr)
library(readxl)
library(rlang)

rm(list=ls())

# load cprd

# load own analysis
analysis = cprd$analysis("ys")

## Bariatric Surgery (OPCS)
## Primary procedure: G281, G282, G283, G284, G285, G301, G302, G303, G304, G312, G321, G331, G716. 
## Revision procedure: G305, G315, G316, G322, G323, G324, G325, G332, G387, G717 
## Gastric balloons and bubbles: G481, G485. 
## OPCS codes from:
## https://digital.nhs.uk/data-and-information/publications/statistical/national-obesity-audit/weight-management-services-hes-22-23-final-and-csds-dq-q1-23-24/appendices

# OPCS codes for bariatric surgery
bariatric_opcs <- c(
  "G281", "G282", "G283", "G284", "G285", "G301", "G302", "G303", "G304", 
  "G312", "G321", "G331", "G716", "G305", "G315", "G316", "G322", "G323", 
  "G324", "G325", "G332", "G387", "G717", "G481", "G485"
)

hes_bariatric_procedures <- cprd$tables$hesProceduresEpi %>%
  filter(OPCS %in% bariatric_opcs) %>%
  select(patid, OPCS, epistart, epiend) %>%
  analysis$cached("bariatric_data")
# 18669
View(hes_bariatric_procedures %>% head(100) %>% collect())





# Anti-obesity (Prodcodes)
## Orlistat (Xenical)
## Alli
## Semaglutide (Wegovy / Ozempic)
## Naltrexone / Bupropion (Mysimba)

# Import Aurum prodcode lookup table
aurum_prodcodes <- read_delim("D:/Documents/Exeter/Dissertation/CPRD/202406_EMISProductDictionary.txt", 
                              col_types = cols(.default=col_character()))

# Find by drug substances name (case-insensitive)
antiobesity_substance_match <- aurum_prodcodes %>%
  filter(grepl("orlistat", DrugSubstanceName, ignore.case = TRUE) |
           grepl("semaglutide", DrugSubstanceName, ignore.case = TRUE) |
           grepl("naltrexone", DrugSubstanceName, ignore.case = TRUE) |
           grepl("bupropion", DrugSubstanceName, ignore.case = TRUE))
# 44

#Find by brand names
antiobesity_name_match <- aurum_prodcodes %>%
  filter(grepl("xenical", `Term from EMIS`, ignore.case = TRUE) |
           grepl("alli", `Term from EMIS`, ignore.case = TRUE) | # OTC, wont able to identify
           grepl("wegovy", `Term from EMIS`, ignore.case = TRUE) |
           grepl("ozempic", `Term from EMIS`, ignore.case = TRUE) |
           grepl("mysimba", `Term from EMIS`, ignore.case = TRUE))
# 315

# Join both datasets
antiobesity_combined <- antiobesity_name_match %>%
  union_all(antiobesity_substance_match)%>%
  rename(prodcodeid = ProdCodeId)
#359


# Pull the relevant drugIssue rows into local memory first
drugissue_local <- cprd$tables$drugIssue %>%
  filter(prodcodeid %in% antiobesity_prodcodes$prodcodeid) %>%
  select(patid, issuedate, drugrecid, prodcodeid) %>%
  analysis$cached("drugissue_local")


# Then join locally
antiobesity_drugissue <- drugissue_local %>%
  inner_join(antiobesity_combined, by = "prodcodeid", copy = TRUE) %>%
  select(patid, issuedate, prodcodeid, DrugSubstanceName, `Term from EMIS`) %>%
  analysis$cached("antiobesity_data")
#5355193
View(antiobesity_drugissue %>% head(100) %>% collect())





# Pregnant (Medcodes)
# Vector of Read code stems for pregnancy
# https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fpds.5584&file=pds5584-sup-0001-supinfo.docx
pregnancy_stems <- c(
  "L", "27", "444", "4H", "584", "62", "63", "68b", "7E08", "7E0J", "7E24", "7F", 
  "939", "95", "9OqC", "9OqJ", "Q", "ZV22", "ZV23", "ZV24", "ZV27", "ZV28", "ZV29", "ZV3"
)

# Create regex pattern that matches any of the stems at the beginning of a string
prefix_regex <- paste0("^(", paste(pregnancy_stems, collapse = "|"), ")")

# Collect the data into memory
medDict_df <- collect(cprd$tables$medDict)

# GFilter using str_detect (case-sensitive by default in stringr)
pregnancy_medcodes_from_stems <- medDict_df %>%
  filter(
    str_detect(cleansedreadcode, prefix_regex) |
      str_detect(originalreadcode, prefix_regex)
  ) %>%
  select(medcodeid, cleansedreadcode, originalreadcode)
#4749


# Read xlsx file containing pregnancy-related Read codes
# https://pmc.ncbi.nlm.nih.gov/articles/instance/6618019/bin/PDS-28-923-s003.xlsx
df <- read_xlsx("D:Downloads/PDS-28-923-s003.xlsx", sheet = "Data")

# Filter to keep only pregnancy codes for mothers, excluding generic and postnatal codes
filtered_readcodes <- df %>%
  filter(
    mother == 1,
    is.na(lmp) | lmp != 1,
    is.na(postnatal_8wk) | postnatal_8wk != 1,
    is.na(postnatal_other) | postnatal_other != 1,
    is.na(preg_related) | preg_related != 1
  ) %>%
  mutate(`read code` = str_trim(`read code`)) # Trim whitespace from read codes

# Trim whitespace from original read codes in medDict and keep relevant columns
medDict_df_trim <- medDict_df %>%
  mutate(originalreadcode = str_trim(originalreadcode)) %>%
  select(medcodeid, cleansedreadcode, originalreadcode)
#253840 from medDcit table

# Find filtered read codes matching original read codes in medDict
matches_orig <- filtered_readcodes %>%
  semi_join(medDict_df_trim, by = c("read code" = "originalreadcode"))

# Find filtered read codes matching cleansed read codes in medDict
matches_cleansed <- filtered_readcodes %>%
  semi_join(medDict_df_trim, by = c("read code" = "cleansedreadcode"))

# Combine matches from both joins and remove duplicates
filtered_matches <- bind_rows(matches_orig, matches_cleansed) %>%
  distinct()
#2482

# Join to medDict to get medcodeid, handling missing values by matching cleansedreadcode if needed
filtered_with_medcodeid <- filtered_matches %>%
  left_join(medDict_df_trim, by = c("read code" = "originalreadcode")) %>%
  mutate(medcodeid = coalesce(medcodeid, medDict_df_trim$medcodeid[match(`read code`, medDict_df_trim$cleansedreadcode)]))

# Combine medcodes found by stems and by filtered read codes, keeping unique medcodeid only
combined_df <- bind_rows(
  pregnancy_medcodes_from_stems,
  filtered_with_medcodeid
) %>%
  distinct(medcodeid, .keep_all = TRUE)  # Keep only unique medcodeid rows
#4943

# Extract vector of medcodeid for pregnancy-related codes
pregnancy_medcodes_vec <- combined_df$medcodeid

# Filter observations to keep only those with pregnancy-related medcodes
observation_medcodes <- cprd$tables$observation %>%
  filter(medcodeid %in% pregnancy_medcodes_vec) %>%
  select(patid, obsdate, medcodeid) %>%
  analysis$cached("observation_medcodes_pregnant")
#6865413

# Retrieve maternity events from Hospital Episode Statistics (HES)
pregnancy_maternity <- cprd$tables$hesMaternity %>%
  select(patid, epistart, epiend) %>%
  analysis$cached("pregnancy_maternity")
#724162

# Add missing columns to observation_medcodes
obs_events <- observation_medcodes %>%
  mutate(
    epistart = as.Date(NA),
    epiend = as.Date(NA)
  ) %>%
  select(patid, obsdate, epistart, epiend, medcodeid)

# Add missing columns to pregnancy_maternity
maternity_events <- pregnancy_maternity %>%
  mutate(
    obsdate = as.Date(NA),
    medcodeid = NA_integer_
  ) %>%
  select(patid, obsdate, epistart, epiend, medcodeid) 

# Combine observation and maternity events into one dataset
pregnancy_events <- obs_events %>%
  union_all(maternity_events) %>%
  analysis$cached("pregnancy_data")
#7589575

View(pregnancy_events %>% head(100) %>% collect())






# Blood dyscrasia (Medcodes)
# Find other name
# Import Aurum medcode lookup table
aurum_medcodes <- read_delim("D:/Documents/Exeter/Dissertation/CPRD/202406_EMISMedicalDictionary.txt", 
                             col_types = cols(.default=col_character()))

# Find medcodes potentially indicating blood dyscrasia conditions (case-insensitive)
blood_dyscrasia_terms <- c(
  "aplastic anemia",
  "leukemia",
  "lymphoma",
  "multiple myeloma",
  "myelodysplastic",
  "thrombocytopenia",
  "neutropenia",
  "pancytopenia",
  "hemolytic anemia",
  "polycythemia vera",
  "waldenstrom",
  "macroglobulinemia",
  "blood dyscrasia", "dyscrasia", "dyscrasias" # general terms
)

# Use a single grepl() with collapse = "|"
blood_dyscrasia_substance_match <- aurum_medcodes %>%
  filter(grepl(paste(blood_dyscrasia_terms, collapse = "|"), Term, ignore.case = TRUE))


# Rename MedCodeId for consistent join
blood_dyscrasia_substance_match <- blood_dyscrasia_substance_match %>%
  rename(medcodeid = MedCodeId)

# extract vector of medcodeids for blood dyscrasia
blood_medcode_ids <- blood_dyscrasia_substance_match$medcodeid

# Observation table
observation_medcodes <- cprd$tables$observation %>%
  filter(medcodeid %in% blood_medcode_ids) %>%
  select(patid, obsdate, medcodeid) %>%
  mutate(medcodeid = as.character(medcodeid)) %>%
  analysis$cached("dycrasia_medcodes")

# Perform inner join
dyscrasia_observations <- observation_medcodes %>%
  inner_join(blood_dyscrasia_substance_match, by = "medcodeid", copy = TRUE) %>%
  analysis$cached("dycrasia_data")
#253011
View(dyscrasia_observations %>% head(100) %>% collect())







