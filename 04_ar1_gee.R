###########################################################################
# Author:   Patrick Rockenschaub
# Project:  Preserve Antibiotics through Safe Stewardship (PASS)
#           Primary Care work package 1
#           Analysis of comorbidity diagnosis and prescribing
#
# File:     04_ar1_gee.R
# Date:     11/06/2019
# Task:     Caclulate the average monthly prescribing using generalised 
#           estimating equations with an AR1 correlation structure
#
###########################################################################



subfolder <- "05_comorb_onset"

# Initialise the workspace
source(file.path(subfolder, "00_init.R"))
source(file.path(subfolder, "00_functions.R"))

library(ggplot2)
library(geepack)
library(broom)
library(forcats)

sensitivity <- FALSE
exclude_index <- TRUE

# Load the data
if(sensitivity & exclude_index){
  data <- read_rds(file.path(subfolder, "01_derived", "data_all_noindex.rds"))
} else if(sensitivity){
  data <- read_rds(file.path(subfolder, "01_derived", "data_all.rds"))
} else if(exclude_index){
  data <- read_rds(file.path(subfolder, "01_derived", "data_noindex.rds"))
} else {
  data <- read_rds(file.path(subfolder, "01_derived", "data.rds"))
}

data %<>% .[month > -13]

match_idx <- read_rds(file.path(subfolder, "01_derived", "match_idx.rds"))
data[, comorb := fct_expand(comorb, "control")]
data[match_idx, on = "patid", 
     comorb := factor(ifelse(comorb == "none", "control", as.character(comorb)),
                      levels(comorb))]


# Set factors
data[, age_cat := cut(age, c(0, 40, 60, 80, Inf), 
                           c("18-39", "40-59", "60-79", "80+"), 
                      right = FALSE)]
data[, imd_cat := factor(imd)]

data[, age_cat := fct_relevel(age_cat, "60-79")]
data[, imd_cat := fct_relevel(imd_cat, "3")]

# Set month and quarter indicators
data[, month := factor(month + 13L, 1:24)]
data[, quarter := factor(ceiling(as.integer(month) / 3), 1:8)]


# Run the models ----------------------------------------------------------

set_month_contrasts <- function(dt, formula, n_quarter = 8){
  # Create a model matrix with contrasts for months oscillating around 
  # their quarterly averages. 
  #
  # If both months and quarters are naively included in the creation of the
  # model matrix, they are linear combinations of each other and the model
  # matrix is rank deficient. This function changes the contrasts manually
  # so that the sum of the months in each quarter sum to the quarterly
  # average, thus making the model identifiable again.
  #
  # Parameters
  # ----------
  #  dt : data.frame
  #    source data used to create the model matrix
  #  formula : formula
  #    right-sided formula object specifying the structure of the model;
  #    must include month and quarter
  #  n_quarter : int
  #    the number of quarters present in the data
  #
  # Result
  # ------
  #  : matrix
  #    model matrix with correct contrasts
  
  m <- model.matrix(formula, data = dt)
  
  # Add the first month as a separate column
  m <- cbind(m, 0)
  colnames(m) <- c(colnames(m)[-ncol(m)], "month1")
  m[rowSums(m[, grep("month", colnames(m)[-ncol(m)])]) == 0, "month1"] <- 1
  
  # Set contrasts to sum to zero for each quarter
  for(i in 1:n_quarter){
    base <- 3 * (i-1)
    
    month <- str_c("month", base + 1:3)
    m[m[, month[3]] == 1, month[1:2]] <- c(-1, -1)
  }
  
  # Delete the superfluous months in each quarter and order the columns
  m <- m[, !grepl(str_c("month(", str_c((1:n_quarter) * 3, collapse = "|"), ")"), colnames(m))]
  idx_months <- grep("month", colnames(m))
  col_months <- colnames(m)[idx_months]
  ord_months <- order(as.numeric(str_extract(col_months, "[0-9].*")))
  
  cbind(m[, -idx_months], m[, col_months[ord_months]])
}


fit_gee_ar1 <- function(dt, formula){
  # Fit a single generalised estimating equation to estimate monthly
  # and quarterly averages
  #
  # Parameters
  # ----------
  #  dt : data.frame 
  #    source data for the model
  #  formula : formula
  #    the model definition as a formula; must include month and quarter
  #
  # Return
  # ------
  #  : geeglm
  #    the fit model
  
  # Set the correct contrasts
  m <- set_month_contrasts(dt, formula)[, -1]

  geeglm(dt$outcome ~ ., data = as.data.frame(m), corstr = "ar1",
         family = "poisson", id = dt$patid, waves = dt$month)
}

data %<>% .[comorb != "none"]
data[, comorb := fct_drop(comorb)]
data$outcome <- data$N
form <- ~ age_cat + female + imd_cat + season + month + quarter

ar1 <- data %>% 
  split(.$comorb) %>%    # Run the models separately for each comorbidity
  map(fit_gee_ar1, formula = form)



# Make predictions for the monthly means ----------------------------------

pred_dt <- data.table(
  age_cat = data[age_cat == "60-79"]$age_cat[1],
  female = 0,
  imd_cat = data[imd_cat == "3"]$imd_cat[1],
  season = data[season == "winter"]$season[1],
  month = factor(1:24),
  quarter = factor(rep(1:8, each = 3))
)


make_prediction <- function(dt, mod, formula){
  # Make a prediction of the average prescribing
  #
  # Parameters
  # ----------
  #  dt : data.table 
  #    artificial or new data used to make predictions
  #  mod : geeglm
  #    a fit generalised estimating equation
  #  formula : formula
  #    the formula used in the model
  #
  # Return
  # ------
  #  : data.table
  #    the input plus additional columns for the monthly mean and standard 
  #    error and the quarterly mean
  
  C <- set_month_contrasts(dt, formula)
  dt$mu_month <- C %*% coefficients(mod)
  dt$se_month <- diag(sqrt(C %*% mod$geese$vbeta %*% t(C)))
  
  C[, grep("month", colnames(C))] <- 0
  dt$mu_quarter <- C %*% coefficients(mod)
  dt$se_quarter <- diag(sqrt(C %*% mod$geese$vbeta %*% t(C)))
  
  dt
}


preds <- ar1 %>% 
  map_df(make_prediction, dt = pred_dt, formula = form, .id = "comorb") %>% 
  copy() # Fix shallow copy warning

preds[, comorb := factor(comorb, levels = levels(data$comorb))]

# Get rates, ratios and p-values for the quarters --------------------------

# Extract the coefficients for each comorbidty
coe <- ar1 %>% map_df(tidy, .id = "comorb")
setDT(coe)

# Get the absolute rate at baseline
rate <- preds[age_cat == "60-79" & month == "1"]
rate[, est := 1000 * exp(mu_quarter)]
rate[, lower := 1000 * exp(mu_quarter + qnorm(0.025) * se_quarter)]
rate[, upper := 1000 * exp(mu_quarter + qnorm(0.975) * se_quarter)]
rate %<>% .[, .(comorb, est, lower, upper)]
rate[, ci := str_c(prty(est, 2), " (", prty(lower, 2), "-", prty(upper, 2), ")")]

# Get rate ratios for the quarters
rr <- coe[str_detect(term, "quarter")]
rr[, quarter := as.integer(str_extract(term, "[2-8]"))]
rr[, quarter := factor(quarter, 1:8)]

# Calculate a point estimate and confidence interval
rr[, est := exp(estimate)]
rr[, lower := exp(estimate + qnorm(0.025) * std.error)]
rr[, upper := exp(estimate + qnorm(0.975) * std.error)]

# Class the p-value into discrete categories for plotting
rr[, p_cat := cut(p.value, c(-Inf, 0.001, 0.01, 0.05, Inf),
                           c("<0.001", "0.001-0.01", "0.01-0.05", ">0.05"))]

# Add the p-value to the predictions
preds[rr, on = .(comorb, quarter), p_cat := p_cat]
preds[, p_cat := fct_explicit_na(p_cat, na_level = "Reference quarter")]

# Bring the four quarters around diagnosis into table format
rr %<>% .[, .(comorb, quarter, est, lower, upper)]
rr[, ci := str_c(prty(est, 2), " (", prty(lower, 2), "-", prty(upper, 2), ")")]
rr[, c("est", "lower", "upper") := NULL]
rr %<>% dcast(comorb ~ quarter, value.var = "ci")


# Plot the monthly means and 95%-CIs --------------------------------------

per_1000 <- function(x){
  # Format to disply the y-axis per 1,000 patients
  x * 1000
}


limits <- expand_limits(
  # Set y-axis limits separately for each comorbidity facet
  comorb = sort(unique(preds$comorb)),
  y = c(175, 105, 105, 255, 105, 205, 105, 105, 105) / 1000
)


# Black and white version
p <- ggplot(preds, aes(x = as.integer(month) - 12.5, y = exp(mu_month), 
                      ymin = exp(mu_month + qnorm(0.025) * se_month), 
                      ymax = exp(mu_month + qnorm(0.975) * se_month))) + 
  geom_vline(xintercept = 0, size = 2, colour = "grey92") + 
  geom_pointrange(colour = "grey40", fatten = 0.5) + 
  geom_point(colour = "grey40", size = 1) + 
  geom_line(aes(y = exp(mu_quarter), group = quarter), colour = "white", size = 1.4) +
  geom_line(aes(y = exp(mu_quarter), group = quarter), colour = "black", size = 1) + 
  limits + 
  facet_wrap(~ comorb, ncol = 2, scales = "free_y", labeller = comorb_labeller) + 
  scale_x_continuous(breaks = c(-12, -9, -6, -3, 0, 3, 6, 9, 12)) + 
  scale_y_continuous(limits = c(0, NA), breaks = scales::pretty_breaks(4),
                     labels = scales::trans_format("per_1000", scales::number)) +
  labs(x = "\n\nMonths to diagnosis", y = "Antibiotics per 1,000 patients\n\n") + 
  coord_cartesian(xlim = c(-12, 12), expand = FALSE) + 
  theme_minimal() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.5, "lines"),
    strip.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

# Single colour version
p_dark <- p + 
  geom_line(aes(y = exp(mu_quarter), group = quarter), colour = "#046C94", size = 1)

p_bright <- p + 
  geom_line(aes(y = exp(mu_quarter), group = quarter), colour = "#D69C4E", size = 1)

# Colourful version according to p-value
p_multi <- p + 
  geom_line(aes(y = exp(mu_quarter), group = quarter, colour = p_cat), size = 1) + 
  scale_colour_manual(values = c("#e66101", "#fdb863", "#b2abd2", "#5e3c99", "#111111"),
                      drop = FALSE) + 
  guides(colour = guide_legend("p-value:", ncol = 2, nrow = 4))
p_multi <- gridExtra::grid.arrange(shift_legend(p_multi))

p


name <- "fig2-rates" # Excluding baseline prescribing is main analysis

if(sensitivity){
  name <- str_c("supp-", name, "-all")
} else if(!exclude_index){
  name <- str_c("supp-", name, "-index")
}

ggsave(file.path(subfolder, res_dir, str_c(name, ".jpg")), p,
       dpi = 600, unit = "in", width = 6.1, height = 9)  
ggsave(file.path(subfolder, res_dir, str_c(name, "_dark.jpg")), p_dark,
       dpi = 600, unit = "in", width = 6.1, height = 9)
ggsave(file.path(subfolder, res_dir, str_c(name, "_bright.jpg")), p_bright,
       dpi = 600, unit = "in", width = 6.1, height = 9)
ggsave(file.path(subfolder, res_dir, str_c(name, "_multi.jpg")), p_multi,
       dpi = 600, unit = "in", width = 6.1, height = 9)
ggsave(file.path(subfolder, res_dir, str_c(name, ".svg")), p,
       unit = "in", width = 6.1, height = 9)  












