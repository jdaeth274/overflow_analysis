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


require(fastcmprsk, lib.loc = "C:/R-3.6.2/library")


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
            path = "E:/HES/COVID/HES_APC_CC_0913_transitions_all_ICD_S1.csv",
            delim = ",")

system.time(transitions_data <- vroom(file = "E:/HES/COVID/HES_APC_CC_0913_transitions_all_ICD.csv", delim = ","))

source(file = "D:/Dropbox/COVID19/Overflow/Rscript/R/survival analysis.R")
survival_res <- survival_analysis_set_up(transitions_data, single_ICD = FALSE, base_dir = "D:/Dropbox/COVID19/Overflow/",
                                         core_18 = TRUE, elective_run = TRUE, emergency_run = FALSE, crr_try = FALSE)


nsurvival_res_neos <- survival_analysis_set_up(neo_trans_data, single_ICD = TRUE, base_dir = "D:/Dropbox/COVID19/Overflow/",
                                              core_18 = TRUE, single_icd = 2)


failure_func_df <- survival_res[[1]]
emergency_ga <- survival_res[[2]]
emergency_cc <- survival_res[[3]]
elective_ga <- survival_res[[4]]
elective_cc <- survival_res[[5]]


write.csv(emergency_ga,
          file = "D:/Dropbox/COVID19/Overflow/cohort_3_ga_transitions_all_ICD.csv",
          row.names = FALSE,
          quote = FALSE)

write.csv(emergency_cc,
          file = "D:/Dropbox/COVID19/Overflow/cohort_3_cc_transitions_all_ICD.csv",
          row.names = FALSE,
          quote = FALSE)

## make elective na mean & medians the day 7 vals 
elective_ga_na <- which(is.na(elective_ga$mean_7))
elective_ga_na_med <- which(is.na(elective_ga$median_7))
elective_ga[c(elective_ga_na, elective_ga_na_med),] <- elective_ga[c(elective_ga_na, elective_ga_na_med), "day_7"]

write.csv(elective_ga,
          file = "D:/Dropbox/COVID19/Overflow/Wolfram/cohort_1_ga_all_ICDS.csv",
          row.names = FALSE,
          quote = FALSE)

elective_ga_na <- which(is.na(elective_cc_2$mean_7))
elective_ga_na_med <- which(is.na(elective_cc_2$median_7))
elective_cc_2[elective_ga_na,c("mean_7","median_7")] <- elective_cc[c(elective_ga_na,
                                                                    elective_ga_na_med), "day_7"]
write.csv(elective_cc_2,
          file = "D:/Dropbox/COVID19/Overflow/Wolfram/cohort_1_cc_all_ICDS.csv",
          row.names = FALSE,
          quote = FALSE)


## time_series ##
source(file = "D:/Dropbox/COVID19/Overflow/Rscript/R/time_series_creator.R")
time_series_data <- time_series_creator(hes_data = transitions_data, num_cores = 12, forecast_date = "2012-03-01")

waiting_pool <- time_series_data[[3]]
write.csv(waiting_pool,
          file = "D:/Dropbox/COVID19/Overflow/Wolfram/waiting_pool.csv",
          row.names = FALSE,
          quote = FALSE)

source(file = "D:/Dropbox/COVID19/Overflow/Rscript/R/time_series_forecast.R")
times_series_forecasts <- running_forecasts(total_cohort_data = time_series_data, train_date = as.Date("2012-03-01"),
                                       forecast_period = 52,  base_dir = "D:/Dropbox/COVID19/Overflow/Results/TimeSeries/")



emergency_admis <- times_series_forecasts[[1]]
emergency_frail <- times_series_forecasts[[2]]
emergency_cc <- times_series_forecasts[[3]]

## cohort 3, icd 2, age 3
emerg <- time_series_data[[1]]
cc_actual <- emerg[emerg$ICD == 9 &
                     emerg$agegrp_v3 == 3,]
cc_actual$date <- as.Date(cc_actual$date)
pred_cc <- emergency_frail[emergency_frail$icd_name == "9" &
                          emergency_frail$age == "3",]
cc_actual <- cc_actual[cc_actual$date >= pred_cc$date[1] ,]
cc_actual <- cc_actual[order(cc_actual$date),]

plot(cc_actual$prop_Frail)
lines(pred_cc$median)
plot(pred_cc$median)
lines(cc_actual$prop_cc)


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





###############################################################################
## quick test of the fastcrr vs crr on the neoplasms cc data ##################
###############################################################################


plyr::count(neo_trans_data$cc_transitions)
neo_trans_data_1 <- neo_trans_data[neo_trans_data$cohort == 1,]
icd_12_cohort_1 <- transitions_data[transitions_data$cohort == 1 &
                                      transitions_data$ICD == 12,]

cc_failcode_3 <- crr(ftime = neo_trans_data_1$cc_LoS, fstatus = neo_trans_data_1$cc_transitions,
                     cov1 = neo_trans_data_1$WaitingTime, failcode = 3)
cc_failcode_3$coef

age_1_icd_12 <- icd_12_cohort_1[icd_12_cohort_1$agegrp_v3 == 1,]

agegrp_3_2_time <- age_1_icd_12$GA_LoS
agegrp_3_2_trans <- age_1_icd_12$ga_transitions
agegrp_3_2_WT <- age_1_icd_12$WaitingTime
current_fail <- 1

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
  agegrp_3_2_trans[old_1ers] <- current_fail
  
}

cc_failcode_3_fast <- fastCrr(Crisk(agegrp_3_2_time, agegrp_3_2_trans, failcode = 1) ~ agegrp_3_2_WT,
                              variance = TRUE, returnDataFrame = TRUE)



###############################################################################
## Write over a icd 12 age grp 1 transitions 2 for ga code #####################
###############################################################################

typpei <- "2-ga-2"
icd_12_dat <- transitions_data[transitions_data$ICD == 12, c("ga_transitions","agegrp_v3","WaitingTime","GA_LoS")]

icd_12_crr <- crr_cluster_run_18(fail_type = typpei, icd_12_dat)



CI.agegrp1_t1 <- icd_12_crr[[1]]

tic("Predicting trans 1")
cif1_pred <- icd_12_crr[[2]]
toc()

cif1 <- icd_12_crr[[2]]




if(7 %in% cif1[,1]){
  seven_t1 <- cif1[cif1[,1] == 7, 2]
}else if(nrow(cif1[cif1[,1] <7, ,drop = FALSE]) > 0){
  seven_t1 <- cif1[cif1[,1] < 7, 2]
  seven_t1 <- max(seven_t1)
}else{
  seven_t1 <- 0
}

max_time <- max(cif1[,1])
if(max_time <= 7){
  multi_df <- seven_t1
}else{
  
  seven_mults <- seq(7, (max_time + 6), by = 7)
  start_sevens <- seven_t1
  
  diff_vec <- start_sevens
  actual_vec <- start_sevens
  
  for(j in 2:length(seven_mults)){
    current_mult_vals <- cif1[cif1[,1] > seven_mults[j-1] &
                                cif1[,1] <= seven_mults[j],,drop = FALSE]
    if(nrow(current_mult_vals) == 0){
      diff_vec <- append(diff_vec, 0)
      actual_vec <- append(actual_vec, actual_vec[length(actual_vec)])
    }else{
      diff_val <- max(current_mult_vals[,2]) - actual_vec[length(actual_vec)]
      diff_vec <- append(diff_vec, diff_val)
      actual_vec <- append(actual_vec, max(current_mult_vals[,2]))
    }
    
  }
}

elective_ga$coeff[65] <- CI.agegrp1_t1$coef
elective_ga$variance[65] <- CI.agegrp1_t1$var[1,1]
elective_ga$day_7[65] <- seven_t1
elective_ga$mean_7[65] <- mean(diff_vec, na.rm = TRUE)
elective_ga$median_7[65] <- median(diff_vec, na.rm = TRUE)
elective_ga$age[65] <- 1
elective_ga$ICD[65] <- 12
elective_ga$transition[65] <- 2






############################################ to be inserted into one of the source files
# calculates patients currently in hospital at start of pandemic
transitions_data_whole$disdate_MDY <- as.Date(transitions_data_whole$disdate_MDY, format = "%d%b%Y")

transitions_data_whole$currentpt <- ifelse(transitions_data_whole$admidate_MDY < "2012-03-05"
                                           & transitions_data_whole$disdate_MDY >= "2012-03-05", 1, 0)

y0 <- subset(transitions_data_whole, currentpt == 1)
y0_collapsed <- do.call(data.frame, aggregate(y0$One, list(y0$ICD, y0$agegrp_v3, y0$admimeth_C, y0$cc), sum))

colnames(y0_collapsed) <- c("ICD","agegrp_v3","admimeth","cc","y0")

y0_collapsed$disease <- "ICD"

y0_collapsed$age <- paste("_AGE",y0_collapsed$agegrp_v3)

y0_collapsed$p <- paste0(y0_collapsed$disease, y0_collapsed$ICD,y0_collapsed$age)
y0_collapsed$S <- ifelse(y0_collapsed$cc == 1, "C", "G")
y0_collapsed$a <- ifelse(y0_collapsed$admimeth == 1, "N", "E")

export <- select(y0_collapsed, p,S,a,y0)

write.csv(export,
          file = "D:/Dropbox/COVID19/Overflow/Wolfram/y0.csv", 
          row.names = FALSE,
          quote = FALSE)


#####waiting patients

transitions_data_whole$disdate_MDY <- as.Date(transitions_data_whole$disdate_MDY, format = "%d%b%Y")

transitions_data_whole$currentpt <- ifelse(transitions_data_whole$rttstart < "2012-03-05"
                                           & transitions_data_whole$admidate_MDY >= "2012-03-05" 
                                           & transitions_data_whole$admimeth_C==1, 1, 0)

x0 <- subset(transitions_data_whole, currentpt == 1)
x0_collapsed <- do.call(data.frame, aggregate(x0$One, list(x0$ICD, x0$agegrp_v3, x0$admimeth_C), sum))

colnames(y0_collapsed) <- c("ICD","agegrp_v3","admimeth","x0")

x0_collapsed$disease <- "ICD"

x0_collapsed$age <- paste("_AGE",x0_collapsed$agegrp_v3)

x0_collapsed$p <- paste0(x0_collapsed$disease, x0_collapsed$ICD,x0_collapsed$age)
#x0_collapsed$S <- ifelse(y0_collapsed$cc == 1, "C", "G")
x0_collapsed$a <- ifelse(x0_collapsed$admimeth == 1, "N", "E")

export <- select(x0_collapsed, p,S,a,x0)

write.csv(export,
          file = "D:/Dropbox/COVID19/Overflow/Wolfram/x0.csv", 
          row.names = FALSE,
          quote = FALSE)












