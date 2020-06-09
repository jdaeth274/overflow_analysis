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

setwd("E:/HES/COVID/")

source(file = "C:/Users/dr314/Dropbox/COVID19/Overflow/Rscript/R/cohort_identification.R")
source(file = "C:/Users/dr314/Dropbox/COVID19/Overflow/Rscript/R/transitions_coding.R")
source(file = "C:/Users/dr314/Dropbox/COVID19/Overflow/Rscript/R/survival analysis.R")


tesdt_df <- data.frame(matrix(ncol = 100, nrow = 1000, data = 100.9))

for(k in 1:nrow(tesdt_df)){
  current_val <- tesdt_df[k,1]
  print(current_val)
  
}



neos_one <- vroom(file = "./HES_APC_CC_0912_Neoplasms1.csv", delim = ",")



## First step to load up HES dataset ##

input_args <- commandArgs(trailingOnly = TRUE)

hes_data <- vroom(input_args[1], delim = ",")


cohort_allocation <- cohort_set_up(num_cores = 14,
                                   data_loc = "E:/HES/COVID/HES_APC_CC_0912_Neoplasms1.csv")
                                    
transitions_data <- hes_transitions(cohort_allocation)
system.time(save(transitions_data, file = "E:/HES/COVID/HES_APC_CC_09_13_Temp02_transitions.RData"))


## Just for neoplasms ##
neos_transitions <- transitions_data[transitions_data$ICD == 2,]



survival_res <- survival_analysis_set_up(transitions_data, single_ICD = TRUE, base_dir = "C:/Users/dr314/Dropbox/COVID19/Overflow/")

failure_func_df <- survival_res[[1]]
emergency_ga <- survival_res[[2]]
emergency_cc <- survival_res[[3]]
elective_ga <- survival_res[[4]]
elective_cc <- survival_res[[5]]


#test_1 <- cohort_allocation[cohort_allocation$hesid == "09A941FAB40A2072CEBE9064E9241B5E",]


neo_vals <- vroom("E:/HES/COVID/HES_APC_CC_0913_TEMP02.csv", delim = ",")
###############################################################################
## Function list here: 


test_data <- vroom("E:/HES/COVID/HES_APC_CC_0913_TEMP02.csv",
                   delim = ",")
test_rows <- plyr::count(test_hesid)
test_data_narrow <- test_data[test_data$hesid == test_hesid,]
test_data_narrow$index <- seq(1, nrow(test_data_narrow))
test_data_narrow$cohort <- 3

electives <- which(test_data_narrow$admimeth_C == 1)

test_data_narrow$cohort[electives] <- 1
parallel_cols <- which(colnames(test_data_narrow) %in% c("hesid","index","diag_01","cohort","rttstart",
                                             "admidate_MDY"))
parallel_data <- test_data_narrow[,parallel_cols]
paralell_admi_nas <- which(is.na(parallel_data$admidate_MDY))
if(length(paralell_admi_nas) > 0){
  parallel_data <- parallel_data[-paralell_admi_nas,]
}
elective_nas <- parallel_data[parallel_data$cohort == 1,]
elective_na_rows <- which(is.na(elective_nas$rttstart))
remove_rtt_missing <- elective_nas$index[elective_na_rows]
if(length(remove_rtt_missing)>0)
  parallel_data <- parallel_data[-remove_rtt_missing,]
parallel_data$rttstart <- as.Date(parallel_data$rttstart, format = "%d%b%Y")
parallel_data$admidate_MDY <- as.Date(parallel_data$admidate_MDY, format = "%d%b%Y")    


cohort_allocator(hes_data = parallel_data, hes_ids = test_rows)





