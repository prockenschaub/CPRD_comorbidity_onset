###########################################################################
# Author:   Patrick Rockenschaub
# Project:  Preserve Antibiotics through Safe Stewardship (PASS)
#           Primary Care work package 1
#           Comorbidity onset analysis
#
# File:     02b_match.R
# Date:     17/06/2019
# Task:     Match a control group of non-comorbid patients
#
###########################################################################


subfolder <- "05_comorb_onset"

# Initialise the workspace
source(file.path(subfolder, "00_init.R"))
source(file.path(subfolder, "00_functions.R"))

library(MatchIt)
library(parallel)

# Load the data
data <- read_rds(file.path(subfolder, "01_derived", "data.rds"))
data <- unique(data[, .(patid, comorb, age, female, imd)])
data[, has_comorb := comorb != "none"]


# Set up the parallel computing environment and run matching 
# within the substrata
cl <- makeCluster(detectCores())


em <- data %>% 
  split(str_c("f", data$female, "i", data$imd)) %>% 
  parLapply(cl = cl, matchit, formula = has_comorb ~ age, 
            method = "nearest", ratio = 1)

# Look for balance in each strata

strat.data <- function(object, subset){
  # This method is heavily adapted from MatchIt's match.data to allow
  # the extraction of data from the stratified matching 
  
  data <- subset
  treat <- object$treat
  wt <- object$weights
  vars <- names(data)
  
  dta <- data.frame(cbind(data, object$weights))
  names(dta) <- c(names(data), weights)
  data <- dta
  
  return(as.data.table(data[wt > 0, ]))
}

match_dts <- data %>% 
  split(str_c("f", data$female, "i", data$imd)) %>% 
  map2(em, ~strat.data(.y, .x))

match_dts %>% map(~ .[, mean(age), by = has_comorb])
data_matched <- rbindlist(match_dts)
match_idx <- data_matched[, .(patid)]

save_derived(match_idx, "match_idx", compress = "gz")

