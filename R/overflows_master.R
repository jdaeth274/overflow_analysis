###############################################################################
## Master Overflows analysis script ###########################################
###############################################################################

require(vroom)
require(dplyr)
require(plyr)
require(ggplot2)
require(ggpubr)
require(snow)
require(tidyverse)
require(lubridate)
require(KFAS)
require(forecast)
require(pryr)


## First step to load up HES dataset ##

input_args <- commandArgs(trailingOnly = TRUE)

hes_data <- vroom(input_args[1], delim = ",")


cohort_allocation <- cohort_set_up(num_cores = 12,
                                   data_loc = "D:/Overflows/HES_APC_CC_0912_Neoplasmsv2.csv")
transitions_data <- hes_transitions(cohort_allocation)

survival_res <- survival_analysis_set_up(transitions_data, single_ICD = TRUE)

failure_func_df <- survival_res[[1]]
emergency_ga <- survival_res[[2]]
emergency_cc <- survival_res[[3]]
elective_ga <- survival_res[[4]]
elective_cc <- survival_res[[5]]


test_1 <- cohort_allocation[cohort_allocation$hesid == "09A941FAB40A2072CEBE9064E9241B5E",]


neo_vals <- vroom("D:/Overflows/HES_APC_CC_0912_Neoplasmsv2.csv", delim = ",")
###############################################################################
## Function list here: 

source(file = "D:/Overflows/cohort_identification.R")
source(file = "D:/Overflows/survival analysis.R")
source(file = "D:/Overflows/transitions_coding.R")
