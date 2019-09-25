###########################################################################
# Author:   Patrick Rockenschaub
# Project:  Preserve Antibiotics through Safe Stewardship (PASS)
#           Primary Care work package 1
#           Analysis of comorbidity diagnosis and prescribing
#
# File:     02a_derive_tables.R
# Date:     11/06/2019
# Task:     Use the extracted data and derive the necessary data tables
#
###########################################################################



subfolder <- "05_comorb_onset"

# Initialise the workspace
source(file.path(subfolder, "00_init.R"))
source(file.path(subfolder, "00_functions.R"))

library(forcats)
library(ggplot2)

# Load the data
char(practices, patients, abx, comorb, cons) %>% 
  walk(load_derived)


type <- "all"
include_index <- FALSE


attrition <- list()



# Limit to patients aged 18 at some point ---------------------------------

patients %<>% .[time_length(birth_date %--% leave_date, u = "years") >= 18]
attrition$`1_eligible` <- nrow(patients)


# Define the incidence date for each patient/condition --------------------

comorb <- comorb[(patients[, .(patid)]), on = "patid", nomatch = 0]
setnames(comorb, "eventdate", "inci_date")
setorder(comorb, patid, inci_date)

if(type == "first"){
  # If the `first` is flag, only consider the first comorbidity for each 
  # patient
  comorb %<>% .[, .SD[1], by = patid]
}

attrition$`2a_comorb` <- length(unique(comorb$patid))
attrition$`2b_non`    <- nrow(patients[!comorb, on = "patid"])

# Add enter and leave date to be able to decide inclusion
df <- patients[, .(patid, pracid, birth_date, enter_date, leave_date)]
df %<>% merge(comorb, by = "patid", all.x = TRUE)
df[, c("consid", "staffid", "medcode") := NULL]

attrition$`3_excl_before_co` <- 
  nrow(df[!is.na(inci_date) & 
           (inci_date < enter_date | 
            time_length(birth_date %--% inci_date, u = "years") < 18)])
attrition$`4_in_obs_co` <- attrition$`2a_comorb` - attrition$`3_excl_before_co`

# Every patient without a diagnosis is labelled as none and given a 
# pseudo-index date based on a randomly drown consultation date
df[is.na(comorb), comorb := "none"]

set.seed(2006)
rnd_cons <- cons[, .(cons_date = sample(eventdate, size = 1)), by = patid]
rnd_cons[, inci_date := as_date(NA)]

df[rnd_cons, on = .(patid, inci_date), inci_date := cons_date]
df %<>% .[!is.na(inci_date)] # remove patients without a consultation

# Finally, set the age at (pseudo-)onset for every patient
df[, age := time_length(birth_date %--% inci_date, u = "year")]

attrition$`3_excl_before_no` <- 
  nrow(df[comorb == "none" & 
            (inci_date < enter_date | 
               time_length(birth_date %--% inci_date, u = "years") < 18)])
attrition$`4_in_obs_no` <- attrition$`2b_non` - attrition$`3_excl_before_no`





# Limit incidence dates to adults with one year prior and after -----------
fu <- df[age > 18 & 
         inci_date %between% list(enter_date %m+% days(360),  # 12 months a 30 days
                                  leave_date %m-% days(360))] # 12 months a 30 days
fu[, inci_season := cut(month(inci_date), 
                        c(0, 3, 6, 9, 12), 
                        c("winter", "spring", "summer", "fall"))]

attrition$`6_fu_co` <- nrow(fu[comorb != "none"])
attrition$`5_excl_fu_co` <- attrition$`4_in_obs_co` - attrition$`6_fu_co`
attrition$`6_fu_no` <- nrow(fu[comorb == "none"])
attrition$`5_excl_fu_no` <- attrition$`4_in_obs_no` - attrition$`6_fu_no`


# Plot the distributions of incidence dates
ggplot(fu, aes(inci_date)) + 
  geom_histogram() + 
  facet_wrap(~ comorb, scales = "free_y") + 
  theme_minimal()

# For each prescription, count the days to or from disease onset
abx_fu <- abx[, .(patid, prescdate)]
abx_fu %<>% .[(fu[, .(patid, inci_date)]), on = "patid", nomatch = 0]
abx_fu[, time_to := time_length(inci_date %--% prescdate, u = "days")]

if(!include_index)
  abx_fu %<>% .[time_to != 0] # Exclude prescriptions at the index date


abx_fu[, rnded := floor((time_to - as.integer(!include_index) * 
                         as.integer(time_to > 0)) / 30)]

abx_fu %<>% .[rnded %between% list(-13, 11)]
setnames(abx_fu, "rnded", "month")

rel_abx <- abx_fu[, .N, by = .(patid, month)]



# Combine the follw-up time and the abx -----------------------------------

data <- fu[, .(month = (-13):11),, by = .(patid, comorb, inci_date, inci_season)]
data[, season := cut(month(inci_date %m+% months(month - as.integer(month > 0))), 
                      c(0, 3, 6, 9, 12), 
                      c("winter", "spring", "summer", "fall"))]
data[rel_abx, on = .(patid, month), "N" := N]
data[is.na(N), N := 0L]


# Set demographic variables
data[df, on = "patid", age := age]
data[patients, on = "patid", female := female]
data[patients, on = "patid", imd := imd]

data[, comorb := factor(comorb, sort(unique(comorb)))]
data[, comorb := fct_relevel(comorb, "none", after = Inf)]



# Save to disk ------------------------------------------------------------

name <- "data"

if(type != "first")
  name <- str_c(name, "_all")

if(!include_index)
  name <- str_c(name, "_noindex")

write_rds(data, file.path(subfolder, der_dir, str_c(name, ".rds")), compress = "gz")


