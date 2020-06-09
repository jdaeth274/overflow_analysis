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

source(file = "D:/Overflows/cohort_identification.R")
source(file = "D:/Overflows/survival analysis.R")
source(file = "D:/Overflows/transitions_coding.R")
source(file = "D:/Overflows/time_series_creator.R")
source(file = "D:/Overflows/time_series_forecast.R")



## First step to load up HES dataset ##

input_args <- commandArgs(trailingOnly = TRUE)

hes_data <- vroom(input_args[1], delim = ",")


cohort_allocation <- cohort_set_up(num_cores = 12,
                                   data_loc = "D:/Overflows/HES_APC_CC_0912_Neoplasmsv2.csv")
transitions_data <- hes_transitions(cohort_allocation)

transitions_data <- vroom("D:/Overflows/neo_trans.csv", delim = ",", num_threads = 8)


survival_res <- survival_analysis_set_up(transitions_data, single_ICD = TRUE, base_dir = "D:/Overflows",
                                         single_icd = "2")

failure_func_df <- survival_res[[1]]
emergency_ga <- survival_res[[2]]
emergency_cc <- survival_res[[3]]
elective_ga <- survival_res[[4]]
elective_cc <- survival_res[[5]]


write.csv(failure_func_df,
          file = "D:/Overflows/cohort_1_to_2_failure_function.csv",
          row.names = FALSE,quote = FALSE)

write.csv(emergency_ga,
          file = "D:/Overflows/cohort_3_GA_transitions.csv",
          row.names = FALSE,quote = FALSE)

write.csv(emergency_cc,
          file = "D:/Overflows/cohort_3_CC_transitions.csv",
          row.names = FALSE,quote = FALSE)

write.csv(elective_ga,
          file = "D:/Overflows/cohort_1_GA_transitions_coeff.csv",
          row.names = FALSE,quote = FALSE)

write.csv(elective_cc,
          file = "D:/Overflows/cohort_1_CC_transitions_coeff.csv",
          row.names = FALSE,quote = FALSE)


transitions_data$ICD <- "2"
source(file = "D:/Overflows/time_series_creator.R")
time_series <- time_series_creator(transitions_data, num_cores = 12)
source(file = "D:/Overflows/time_series_forecast.R")
times_series_data <- running_forecasts(total_cohort_data = time_series, train_date = as.Date("2012-03-01"),
                                       forecast_period = 52, single_icd = 2, base_dir = "D:/Overflows/")
emergency_admis <- times_series_data[[1]]
emergency_frail <- times_series_data[[2]]
elective_admis <- times_series_data[[3]]
elective_median <- times_series_data[[4]]
elective_mean <- times_series_data[[5]]
elective_frail <- times_series_data[[6]]

write.csv(emergency_admis, 
          file = "D:/Overflows/neoplasms_cohort_1_Admissions_forecast.csv",
          row.names = FALSE, quote = FALSE)
write.csv(emergency_frail, 
          file = "D:/Overflows/neoplasms_cohort_1_Frail_forecast.csv",
          row.names = FALSE, quote = FALSE)
write.csv(elective_admis, 
          file = "D:/Overflows/neoplasms_cohort_3_new_pool_forecast.csv",
          row.names = FALSE, quote = FALSE)
write.csv(elective_median, 
          file = "D:/Overflows/neoplasms_cohort_3_median_wait_forecast.csv",
          row.names = FALSE, quote = FALSE)
write.csv(elective_mean, 
          file = "D:/Overflows/neoplasms_cohort_3_mean_wait_forecast.csv",
          row.names = FALSE, quote = FALSE)
write.csv(elective_frail, 
          file = "D:/Overflows/neoplasms_cohort_3_frailty_forecast.csv",
          row.names = FALSE, quote = FALSE)


test_1 <- cohort_allocation[cohort_allocation$hesid == "09A941FAB40A2072CEBE9064E9241B5E",]


neo_vals <- vroom("D:/Overflows/HES_APC_CC_0912_Neoplasmsv2.csv", delim = ",")
###############################################################################
## Function list here: 


test_data <- vroom("D:/Overflows/HES_APC_CC_0912_Neoplasmsv2.csv",
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


if(7 %in% agegrp1_pred[,1]){
  seven_t1_cc <- agegrp1_pred[agegrp1_pred[,1] == 7, 2]
}else if(nrow(agegrp1_pred[agegrp1_pred[,1] <7,,drop = FALSE ]) > 0){
  seven_t1_cc <- agegrp1_pred[agegrp1_pred[,1] < 7, 2]
  seven_t1_cc <- max(seven_t1_cc)
}else{
  seven_t1_cc <- 0
}
seven_t1_cc



## testing out fastcrr ##


agegrp1 <- transitions_data[transitions_data$cohort == 1 &
                              transitions_data$agegrp_v3 == 1,]
agegrp1 <- agegrp1[-which(is.na(agegrp1$GA_LoS)),]
agegrp1[which(is.na(agegrp1$ga_transitions)),"ga_transitions"] <- 0
require(tictoc)
require(fastcmprsk)
tic("fast crr - use for transition 1")
  ci_2 <- fastcmprsk::fastCrr(Crisk(agegrp1$GA_LoS, agegrp1$ga_transitions,failcode = 1) ~ agegrp1$WaitingTime,
                            variance = TRUE, var.control = varianceControl(useMultipleCores = TRUE), returnDataFrame = TRUE)
toc()


### getting predict.fcrr from github ###
source("D:/Overflows/pedict_fcrr.R")

tic("Pred step")
ci_pred <- predict(ci_2, newdata = 0, tL = 1)
toc()
fcrrresults <- do.call(rbind, Map(data.frame, A=ci_pred$ftime, B=ci_pred$CIF))



### compare with old crr function 
agegrp1 <- transitions_data[transitions_data$cohort == 1 &
                              transitions_data$agegrp_v3 == 1,]

tic("old CRR - for comparison")
CI.agegrp1_t1 <- crr(ftime = agegrp1$GA_LoS, fstatus = agegrp1$ga_transitions,
                     cov1 = agegrp1$WaitingTime, failcode = 1)
toc()

t1_pred <- predict.crr(CI.agegrp1_t1, cov1 = 0)
crrresults <- do.call(rbind, Map(data.frame, A=t1_pred[,1], B=t1_pred[,2]))

tic("old CRR - use for transition 2")
CI.agegrp1_t2 <- crr(ftime = agegrp1$GA_LoS, fstatus = agegrp1$ga_transitions,
                     cov1 = agegrp1$WaitingTime, failcode = 2)
toc()

t2_pred <- predict.crr(CI.agegrp1_t2, cov1 = 0)
crrresults_t2 <- do.call(rbind, Map(data.frame, A=t2_pred[,1], B=t2_pred[,2]))

### use fastcrr for the largest transition group; use crr for everything else
agegrp1 <- transitions_data[transitions_data$cohort == 1 &
                              transitions_data$agegrp_v3 == 1,]
agegrp2 <- transitions_data[transitions_data$cohort == 1 &
                              transitions_data$agegrp_v3 == 2,]
agegrp3 <- transitions_data[transitions_data$cohort == 1 &
                              transitions_data$agegrp_v3 == 3,]
plyr::count(agegrp2$ga_transitions)

### crr with cc transition 2
tic("GA crr")
CI.agegrp2_t1 <- crr(ftime = agegrp2$GA_LoS, fstatus = agegrp2$ga_transitions,
                        cov1 = agegrp2$WaitingTime, failcode = 1)
toc()
t1_pred_cc <- predict.crr(CI.agegrp1_t1_cc, cov1 = 0)
crrresults_t1_cc <- do.call(rbind, Map(data.frame, A=t1_pred_cc[,1], B=t1_pred_cc[,2]))

### fcrr with cc transition 2
agegrp1 <- agegrp1[-which(is.na(agegrp1$cc_LoS)),]
agegrp1[which(is.na(agegrp1$cc_transitions)),"cc_transitions"] <- 0

ci_2_cc <- fastcmprsk::fastCrr(Crisk(agegrp1$cc_LoS, agegrp1$cc_transitions,failcode = 2) ~ agegrp1$WaitingTime,
                            variance = TRUE, var.control = varianceControl(useMultipleCores = TRUE), returnDataFrame = TRUE)

ci_pred_cc <- predict(ci_2_cc, newdata = 0, tL = 1)

fcrrresults_cc <- do.call(rbind, Map(data.frame, A=ci_pred$ftime, B=ci_pred$CIF))