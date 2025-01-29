##### THIS SCRIPT COMPILES ALL REQUIRED FUNCTIONS FOR THE EFFECTIVE N PROJ #####

### Packages required ###
pkgs <- c(
  "survival",
  "tinytex",
  "tidyverse",
  "dplyr",
  "ggplot2",
  "simsurv",
  "matrixStats",
  "ggsurvfit",
  "svglite"
)

# Check if the package is installed, if yes load, if no: install + load 
vapply(pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  require(pkg, character.only = TRUE, quietly = TRUE)
}, FUN.VALUE = logical(length = 1L))

rm(pkgs)


source("code/modified.R")
source("code/calculate_ess.R")
source("code/survfit_n.R")
source("code/plot_km_eff.R")
source("code/plot_km_eff2.R")
source("code/plot_effective_n.R")
source("code/sf_to_df.R")
source("code/trial_report.R")