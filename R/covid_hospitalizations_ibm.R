########################################################
### Title: COVID hospitalizations from IBM output ######
########################################################

setwd("C:/Users/KL1215/Box Sync/COVID overflows/IBM/ibm_formatting")

# import excel sheet: Let's try good compliance first
excel_sheet <- "GoodCompliance_R0=2.8.xlsx"

# uses the following packages and functions
library(here)
source(here::here("scripts","read_ibm_excel.R"))
source(here::here("scripts", "packages.R"))
source(here::here("scripts", "rolling_sums.R"))
source(here::here("scripts", "symptom_inc_calc.R"))
source(here::here("scripts", "prevalence_plotter.R"))
source(here::here("scripts", "cumulative_death.R"))

# read excel sheet
ibm_data <- read_ibm_excel(excel_sheet)
