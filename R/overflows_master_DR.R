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

devtools::install_github("erickawaguchi/fastcmprsk",ref = "developer",
                         lib="C:/R-3.6.2/library")


require(fastcmprsk)


setwd("E:/HES/COVID/")

source(file = "D:/Dropbox/COVID19/Overflow/Rscript/R/cohort_identification.R")
source(file = "D:/Dropbox/COVID19/Overflow/Rscript/R/transitions_coding.R")
source(file = "D:/Dropbox/COVID19/Overflow/Rscript/R/survival analysis.R")
source(file = "D:/Dropbox/COVID19/Overflow/Rscript/R/time_series_creator.R")
source(file = "D:/Dropbox/COVID19/Overflow/Rscript/R/time_series_forecast.R")

neos_one <- vroom(file = "./HES_APC_CC_0912_Neoplasms1.csv", delim = ",")
subset_dat <- read.csv("./HES_APC_CC_0913_TEMP02_s1.csv", stringsAsFactors = FALSE)

## First step to load up HES dataset ##

input_args <- commandArgs(trailingOnly = TRUE)

hes_data <- vroom(input_args[1], delim = ",")


cohort_allocation <- cohort_set_up(num_cores = 14,
                                   data_loc = "E:/HES/COVID/HES_APC_CC_0913_TEMP02_s1.csv")
                                    
transitions_data <- hes_transitions(cohort_allocation)
neo_trans_data <- vroom(file = "E:/HES/COVID/neo_trans_data_07_06_2020.csv",
                        delim = ",")

## Just for neoplasms ##
vroom_write(transitions_data,
            path = "E:/HES/COVID/HES_APC_CC_0913_transitions_all_ICD.csv",
            delim = ",",
            num_threads = 6,
            progress = TRUE)

source(file = "D:/Dropbox/COVID19/Overflow/Rscript/R/survival analysis.R")
survival_res <- survival_analysis_set_up(transitions_data, single_ICD = FALSE, base_dir = "D:/Dropbox/COVID19/Overflow/",
                                         core_18 = TRUE)


survival_res_neos <- survival_analysis_set_up(neo_trans_data, single_ICD = TRUE, base_dir = "D:/Dropbox/COVID19/Overflow/",
                                              core_18 = TRUE, single_icd = 2)


failure_func_df <- survival_res[[1]]
emergency_ga <- survival_res[[2]]
emergency_cc <- survival_res[[3]]
elective_ga <- survival_res[[4]]
elective_cc <- survival_res[[5]]

## time_series ##
source(file = "D:/Dropbox/COVID19/Overflow/Rscript/R/time_series_creator.R")
time_series_data <- time_series_creator(hes_data = transitions_data, num_cores = 12, forecast_date = "2012-03-01")

times_series_forecasts <- running_forecasts(total_cohort_data = time_series_data, train_date = as.Date("2012-03-01"),
                                       forecast_period = 52,  base_dir = "E:/HES/COVID/")


emergency_admis <- times_series_forecasts[[1]]
emergency_frail <- times_series_forecasts[[2]]
elective_admis <- times_series_forecasts[[3]]
elective_median <- times_series_forecasts[[4]]
elective_mean <- times_series_forecasts[[5]]
elective_frail <- times_series_forecasts[[6]]


#test_1 <- cohort_allocation[cohort_allocation$hesid == "09A941FAB40A2072CEBE9064E9241B5E",]


neo_vals <- vroom("E:/HES/COVID/HES_APC_CC_0913_TEMP02.csv", delim = ",")
###############################################################################
## Function list here: 


neos_one <- neo_trans_data[neo_trans_data$cohort == 1,c("ga_transitions","agegrp_v3", "WaitingTime", "GA_LoS")]


fail_type <- "2-ga-3"

current_fail <- as.integer(str_split_fixed(fail_type, "-",3)[1])
current_type <- as.character(str_split_fixed(fail_type,"-",3)[2])
current_age <- as.character(str_split_fixed(fail_type,"-",3)[3])

cohort_data <- neos_one[neos_one$agegrp_v3 == as.integer(current_age),]
plyr::count(cohort_data$ga_transitions)

ga_function <- function(cohort_data){
print(current_type)
print(current_fail)
print(current_age)
if(current_type == "ga"){
  agegrp_3_2_time <- cohort_data$GA_LoS
  agegrp_3_2_trans <- cohort_data$ga_transitions
  agegrp_3_2_WT <- cohort_data$WaitingTime
  
  na_tims <- which(is.na(agegrp_3_2_time))
  if(length(na_tims) > 0){
    agegrp_3_2_time <- agegrp_3_2_time[-na_tims]
    agegrp_3_2_trans <- agegrp_3_2_trans[-na_tims]
    agegrp_3_2_WT <- agegrp_3_2_WT[-na_tims]
  }
  old_na <- which(is.na(agegrp_3_2_trans))
  if(length(old_na) > 0){
    agegrp_3_2_time <- agegrp_3_2_time[-old_na]
    agegrp_3_2_trans <- agegrp_3_2_trans[-old_na]
    agegrp_3_2_WT <- agegrp_3_2_WT[-old_na]
  }
  if(current_fail != 1){
    old_2ers <- which(agegrp_3_2_trans == current_fail)
    old_1ers <- which(agegrp_3_2_trans == 1)
    
    agegrp_3_2_trans[old_2ers] <- 1
    agegrp_3_2_trans[old_1ers] <- 2
    
  }
  
  
  CI.agegrp1_t1 <- fastCrr(Crisk(agegrp_3_2_time, agegrp_3_2_trans, failcode = 1) ~ agegrp_3_2_WT,
                           variance = TRUE, returnDataFrame = TRUE)
  
  # cif1_pred <- predict(CI.agegrp1_t1, newdata = 0, tL = 1)
  # cif1_pred <- cbind(cif1_pred$ftime, cif1_pred$CIF)
  cif1_pred <- fastcmprsk:::predict.fcrr(CI.agegrp1_t1, newdata = 0, tL = 1, type = "interval")
  cif1_pred <- cbind(cif1_pred$ftime, cif1_pred$CIF, cif1_pred$lower, cif1_pred$upper )
  
}

}
ga_function(cohort_data = cohort_data)




