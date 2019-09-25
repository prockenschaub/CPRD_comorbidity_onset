###########################################################################
# Author:   Patrick Rockenschaub
# Project:  Preserve Antibiotics through Safe Stewardship (PASS)
#           Primary Care work package 1
#           Comorbidity onset analysis
#
# File:     00_init.R
# Date:     16/05/2019
# Task:     Initialise the workspace for a clean new analysis of antibiotic
#           prescribing related to the onset/diagnosis of chronic disease
#
###########################################################################


# Load the general environment (only if not called before)
if(exists(".conn")){
   disconnect_db()
}

# Load the base settings and functions
suppressMessages({
  source("00_init.R")
  source("00_basic_tables.R")
  source("00_code_lists.R")
})

# Load local functions
subfolder <- "05_comorb_onset"

# Set path to store the derived datasets
der_dir <- "01_derived"
res_dir <- "99_results"

study_start <- ymd("2008-01-01")
study_end <- ymd("2015-12-31")