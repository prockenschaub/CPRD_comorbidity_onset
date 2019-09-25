###########################################################################
# Author:   Patrick Rockenschaub
# Project:  Preserve Antibiotics through Safe Stewardship (PASS)
#           Primary Care work package 1
#           Comorbidity onset analysis
#
# File:     01_extract_data.R
# Date:     11/06/2019
# Task:     Extract all necessary records from the database and store them
#           in local R files for reproducability
#
###########################################################################


# Path from project directory to this file
# NOTE: must be set in each programme separately
subfolder <- "05_comorb_onset"

# Initialise the workspace
source(file.path(subfolder, "00_init.R"))
source(file.path(subfolder, "00_functions.R"))


# Select practices --------------------------------------------------------
## @knitr practices
#+ practices, include = FALSE

# Obtain information on when practices provided data
def_practice <- 
  practice_db() %>% 
  filter(linked == 1L) %>% 
  compute(name = "practices")

practices <- collect_dt(def_practice, convert = TRUE)



# Select patient population -----------------------------------------------
## @knitr patients
#+ patients, include = FALSE

# Define start and end dates
def_patients <-
  study_population_db(link = TRUE) %>% 
  in_date(study_start, study_end) %>% 
  semi_join(def_practice, by = "pracid") %>% 
  left_join(imd_db(), by = c("patid", "pracid"))

# Apply exclusion criteria
def_patients <-
  def_patients %>% 
  filter(
    !is.na(birth_date), 
    !is.na(female),
    !is.na(imd)
  ) %>% 
  select(patid, pracid, female, imd, birth_date, 
         death_date, enter_date, leave_date)

def_patients <- compute(def_patients, name = "patients")

patients <- collect_dt(def_patients, convert = TRUE)




# Select all systemic antibiotic prescribing ------------------------------
## @knitr abx
#+ abx, include = FALSE

# Define the included and excluded BNF chapters
bnf_systemic <- "0501"
bnf_excl <- c("050109", "050110")

# Get all prescription data for the above defined study period  
def_abx <-
  abx_bnf_db(bnf_systemic, bnf_excl) %>% 
  semi_join(def_patients, by = "patid") %>% 
  in_date(study_start, study_end) %>% 
  select(patid, prescdate = eventdate, prodcode, consid, 
         qty, ndd, issueseq) %>% 
  arrange(patid, eventdate, issueseq)

abx <- collect_dt(def_abx, convert = TRUE)


# Some drugs are classified in 5.1., but are crossover-products
# Exclude those
abx_info <- antibiotics() %>% collect_dt()

atc_a02b <- abx_info[atcchapter == "A02B"]
abx %<>% .[!atc_a02b, on = "prodcode"]
remove(atc_a02b)


abx %<>% .[!(abx_info[group %in% c("Antifungal", "Antileprotic", 
                                   "Antituberculosis", "No antibiotic")]),
           on = "prodcode"]


# No unique identifier of prescription exists in the database, so create 
# one for this project.
#
# NOTE: This identifier might change with each new database query, so the 
#       identifier is only valid within the same session (and not between
#       sessions)
abx <- abx[order(patid, prescdate, prodcode, qty, issueseq)]
abx[, "abx_id" := .I]





# Obtain comorbidities ----------------------------------------------------
## @knitr comorbs
#+ comorbs, include = FALSE

# Use all other comorbidity code
codes_qof <- 
  codes_qof_db() %>% 
  filter(!(subclass %in% c("ckd_stage_1", "ckd_stage_2"))) %>% 
  collect_dt()


# Select all records with one of those codes in the study population
# NOTE: exclude COPD (by definition) and CKD stage 1&2 (unreliable coding)
comorb <- 
  qof_db() %>% 
  records_db(codes_qof) %>% 
  semi_join(def_patients, by = "patid") %>% 
  filter(!is.na(eventdate)) %>% 
  select(-num_rows) %>% 
  collect_dt(convert = TRUE)


# Add the disease labels
comorb[codes_qof, on = "medcode", comorb := comorbidity]

# Keep the earliest for each comorbidity
comorb %<>% 
  split(., f = .$comorb) %>% 
  map_df(first_record)
  


# Get all consultations in the study period -------------------------------
## @knitr comorbs
#+ comorbs, include = FALSE

def_consultations <-
  consult_db() %>% 
  in_date(study_start , study_end) %>% 
  semi_join(def_patients, by = "patid") %>% 
  select(patid, eventdate)

cons <- collect_dt(def_consultations, convert = TRUE)
cons %<>% unique()



# Save the datasets -------------------------------------------------------

mget(c("practices", "patients", "abx", "comorb", "cons")) %>% 
  walk2(., names(.), save_derived, compress = "gz")


