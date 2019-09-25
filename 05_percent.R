###########################################################################
# Author:   Patrick Rockenschaub
# Project:  Preserve Antibiotics through Safe Stewardship (PASS)
#           Primary Care work package 1
#           Analysis of comorbidity diagnosis and prescribing
#
# File:     04_percent.R
# Date:     12/06/2019
# Task:     Caclulate and plot the proportion of patients in each month 
#           that received X prescriptions in the last 2 months
#
###########################################################################



subfolder <- "05_comorb_onset"

# Initialise the workspace
source(file.path(subfolder, "00_init.R"))
source(file.path(subfolder, "00_functions.R"))

library(ggplot2)
library(forcats)

# Load the data
data <- read_rds(file.path(subfolder, "01_derived", "data_noindex.rds"))

match_idx <- read_rds(file.path(subfolder, "01_derived", "match_idx.rds"))
data[, comorb := fct_expand(comorb, "control")]
data[match_idx, on = "patid", 
     comorb := factor(ifelse(comorb == "none", "control", as.character(comorb)),
                      levels(comorb))]

data <- data[comorb != "none"]
data %<>% .[month > -13]

# Aggregate the cumulative number of antibiotics in each month ------------

#data[, cum := N + shift(N), by = .(patid, comorb)]
data[, cum := N, by = .(patid, comorb)]
data[cum >= 3, cum := 3]
data[, cum := factor(cum, rev(0:3))]

hist <- data[, .N, by = .(comorb, month, cum)] %>% 
           .[, .(cum, N, perc = N / sum(N)), by = .(comorb, month)]


hist[cum != "0", .(N = sum(N), perc = sum(perc)), by = .(comorb, month)]




# Plot the monthly percent and 95%-CIs --------------------------------------

limits <- expand_limits(
  # Set y-axis limits separately for each comorbidity facet
  comorb = sort(unique(data$comorb)),
  #y = c(0.4, 0.2, 0.2, 0.4, 0.2, 0.4, 0.2, 0.2, 0.2)
  y = c(0.2, 0.1, 0.1, 0.2, 0.1, 0.2, 0.1, 0.1, 0.1)
)


p <- ggplot(hist[cum != "0"], aes(x = month + 0.5, y = perc, fill = cum)) + 
  geom_vline(xintercept = 0, size = 1, colour = "grey") + 
  geom_col(width = 0.8) + 
  scale_x_continuous(breaks = c(-12, -9, -6, -3, -0.5, 2, 5, 8, 11) + 0.5,
                     labels = (-4):4 * 3) + 
  scale_y_continuous(breaks = scales::pretty_breaks(4),
                     labels = scales::percent_format(accuracy = 1)) + 
  scale_fill_manual(values = c("#e66101", "#fdb863", "#5e3c99"),
                    labels = c("Three or more antibiotics", "Two antibiotics", "One antibiotic")) + 
  limits + 
  guides(fill = guide_legend(title.position = "top", reverse = TRUE)) + 
  facet_wrap(~ comorb, ncol = 2, scales = "free", labeller = comorb_labeller) + 
  labs(x = "\n\nMonths to diagnosis", y = "Percent of patients\n\n", 
       fill = "# of antibiotics in previous 2 months:") + 
  theme_minimal() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

p <- gridExtra::grid.arrange(shift_legend(p))
p

ggsave(file.path(subfolder, res_dir, "fig3-percent.tiff"), p,
       dpi = 600, unit = "in", width = 6.1, height = 9)  
ggsave(file.path(subfolder, res_dir, "fig3-percent.svg"), p,
       unit = "in", width = 6.1, height = 9)  
