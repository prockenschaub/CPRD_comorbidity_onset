###########################################################################
# Author:   Patrick Rockenschaub
# Project:  Preserve Antibiotics through Safe Stewardship (PASS)
#           Primary Care work package 1
#           Comorbidity onset analysis
#
# File:     03_univariable.R
# Date:     11/06/2019
# Task:     Create a forest plot of univariable associations with all
#           included variables
#
###########################################################################


subfolder <- "05_comorb_onset"

# Initialise the workspace
source(file.path(subfolder, "00_init.R"))
source(file.path(subfolder, "00_functions.R"))

library(ggplot2)
library(gtable)
library(gridExtra)

library(forcats)
library(broom)

# Load the data
data <- read_rds(file.path(subfolder, "01_derived", "data.rds"))
data %<>% .[month > -13]

match_idx <- read_rds(file.path(subfolder, "01_derived", "match_idx.rds"))
data[, comorb := fct_expand(comorb, "control")]
data[match_idx, on = "patid", 
     comorb := factor(ifelse(comorb == "none", "control", as.character(comorb)),
                       levels(comorb))]


# Aggregate the data to one 24 month long observation per pat -------------

data %<>% .[, .(n_abx = sum(N), n_month = .N), 
            by = .(patid, age, female, imd, inci_season, comorb)]
data[, age_cat := cut(age, c(0, 40, 60, 80, Inf), 
                           c("18-39", "40-59", "60-79", "80+"), 
                      right = FALSE)]
data[, sex := factor(female, 0:1, c("Men", "Women"))]
data[, imd_cat := factor(imd)]
levels(data$comorb) <- c("Asthma", "CHD", "CKD", "COPD", "Diabetes", 
                         "Heart failure", "PAD", "Stroke", 
                         "No comorbidity", "Control")

# Set reference category to the largest group of each level
lbls_age <- levels(data$age_cat)
lbls_sex <- levels(data$sex)
lbls_imd <- levels(data$imd_cat)
lbls_comorb <- levels(data$comorb)
lbls <- list(lbls_age, lbls_sex, lbls_imd, lbls_comorb)

data[, age_cat := fct_relevel(age_cat, "60-79")]
data[, imd_cat := fct_relevel(imd_cat, "3")]
data[, comorb := fct_relevel(comorb, "No comorbidity")]



# Tabulate counts ---------------------------------------------------------

mean_sd <- function(x){
  # Calculate the mean and standard deviation of a vector
  #
  # Parameters
  # ----------
  #  x : numeric vector
  #    data to summarise
  #
  # Result
  # ------
  #  : character vector
  #    mean and standard deviation as character
  
  mu <- mean(x)
  sig <- sd(x)
  
  str_c(prty(mu, 1), " (", prty(sig, 1), ")")
}

data[, .(age = mean_sd(age))]
data[, .(.N, .N / nrow(data)), by = sex]
data[, .(imd = mean_sd(imd))]


tab_com <- data[, .(.N, perc = perc(.N, nrow(data))), by = comorb]
tab_age <- data[, .(age = mean_sd(age)), by = comorb]
tab_sex <- data[, .N, by = .(comorb, female)] %>% 
              .[, .(female, sex = perc(N, sum(N))), by = comorb] %>% 
              .[female == 1, .(comorb, sex)]
tab_imd <- data[, .(imd = mean_sd(imd)), by = comorb]

# Combine into a single table
tbl_counts <- list(tab_com, tab_age, tab_sex, tab_imd) %>% 
  reduce(merge, by = "comorb", sort = FALSE)
setorder(tbl_counts, comorb)



# Run and present the estimates -------------------------------------------

univar_glm <- function(dt, var){
  # Run a generalised linear quasi-poisson model with one
  # explanatory variable
  #
  # Parameters
  # ----------
  #  dt : data.table 
  #    data to train on
  #  var : character
  #    name of the explanatory variable/column as vector of length 1
  #
  # Return
  # ------
  #  : glm
  #    fitted model
  
  glm(as.formula(str_c("n_abx ~ ", var)), data = dt, 
      family = "quasipoisson", offset = log(n_month))
}


coef_ci <- function(mod, var, lbls, per = 1){
  # Extract the coefficients from a univariable (quasi-)poisson generalised linear model
  # fit via `glm()`
  #
  # Parameters
  # ----------
  #  mod : glm
  #    fitted model object
  #  var : character
  #    name of the explanatory variable/column as vector of length 1
  #  lbls : character
  #    levels of the variable in the desired order
  #  per : int
  #    factor to multiply the estimates by (e.g. to get rates per 1000)
  #
  # Return
  # ------
  #  : data.table
  #    point estimates and 95% confidence boundaries
  
  # Build the design matrix
  pred_dt <- unique(data[, .SD, .SDcols = var])
  C <- model.matrix(data = pred_dt, as.formula(str_c("~", var)))
  
  # Bring into the desired order
  pred_dt[, (var) := fct_relevel(get(var), lbls)]
  C <- C[order(pred_dt[[var]]), ]
  
  # Calculate means and standard errors  
  preds <- data.table(
    term = lbls,
    eta = as.vector(C %*% coefficients(mod)),
    se = sqrt(diag(C %*% vcov(mod) %*% t(C)))
  )                
  
  preds[, est := per * exp(eta)]
  preds[, lower := per * exp(eta + qnorm(0.025) * se)]
  preds[, upper := per * exp(eta + qnorm(0.975) * se)]
  preds[, ci := str_c(prty(est, 1), " (", 
                      prty(lower, 1), "-", 
                      prty(upper, 1), ")")]
  
  preds[]
}



sum(data$n_abx)

# Run a univariable analysis for all patient variables
vars <- c("age_cat", "sex", "imd_cat", "comorb")
univar <-  c("1", vars) %>% 
  map(univar_glm, dt = data)

interc <- tidy(univar[[1]])
total <- data.table(
  est = 1000 * exp(interc[1, "estimate"]),
  lower = 1000 * exp(interc[1, "estimate"] - 1.96 * interc[1, "std.error"]),
  upper = 1000 * exp(interc[1, "estimate"] + 1.96 * interc[1, "std.error"])                
)

rates <- pmap_df(list(univar[-1], vars, lbls), coef_ci, per = 1000)
rates[, term := factor(term, levels = rev(unlist(lbls)))]



# Calculate the adjusted values -------------------------------------------

adj_mod <- data[comorb != "No comorbidity"] %>% 
  split(.$comorb) %>% .[-1] %>% # -1 to remove the empty factor 
  map(univar_glm, var = "age_cat + sex + imd_cat")

pred_dt <- data.table(
  age_cat = sort(data$age_cat)[1],
  sex = sort(data$sex)[1],
  imd_cat = sort(data$imd_cat)[1],
  n_month = 1
)

make_prediction <- function(mod, newdata){
  # Make a prediction of the average prescribing
  #
  # Parameters
  # ----------
  #  mod : glm
  #    a fit generalised linear model
  #  newdata : data.table
  #    all covariates for the predictions
  #
  # Return
  # ------
  #  : data.table
  #    predicted means for all given covariate combinations
  
  preds <- predict(mod, pred_dt, se.fit = TRUE)
  preds["residual.scale"] <- NULL
  setDT(preds)
  
  preds
}

preds <- map_df(adj_mod, make_prediction, pred_dt, .id = "comorb") %>% copy()

preds[, est := 1000 * exp(fit)]
preds[, lower := 1000 * exp(fit + qnorm(0.025) * se.fit)]
preds[, upper := 1000 * exp(fit + qnorm(0.975) * se.fit)]
preds[, ci := str_c(prty(est, 2), " (", prty(lower, 2), "-", prty(upper, 2), ")")]

preds[, c("fit", "se.fit") := NULL]


