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
                                         core_18 = TRUE, elective_run = TRUE, emergency_run = TRUE, crr_try = TRUE)
survival_res_emerg <- survival_analysis_set_up(transitions_data, single_ICD = FALSE, base_dir = "D:/Dropbox/COVID19/Overflow/",
                                         core_18 = TRUE, elective_run = FALSE, emergency_run = TRUE, crr_try = TRUE)
failure_func_df <- survival_res[[1]]
emergency_ga <- survival_res[[2]]
emergency_cc <- survival_res[[3]]
elective_ga <- survival_res[[4]]
elective_cc <- survival_res[[5]]

save(survival_res, file = "D:/Dropbox/COVID19/Overflow/SA_res_2020_7_7.Rdata")
load(file = "D:/Dropbox/COVID19/Overflow/SA_res_2020_7_7.Rdata", verbose = TRUE)

write.csv(failure_func_df,
          file = "D:/Dropbox/COVID19/Overflow/cohort_1_to_3_5_7_2020.csv",
          row.names = FALSE,
          quote = FALSE)


write.csv(emergency_ga,
          file = "D:/Dropbox/COVID19/Overflow/cohort_3_ga_transitions_all_ICD_5_7_2020.csv",
          row.names = FALSE,
          quote = FALSE)

write.csv(emergency_cc,
          file = "D:/Dropbox/COVID19/Overflow/cohort_3_cc_transitions_all_ICD_5_7_2020.csv",
          row.names = FALSE,
          quote = FALSE)

## make elective na mean & medians the day 7 vals 
elective_ga_na <- which(is.na(elective_ga$mean_7))
elective_ga_na_med <- which(is.na(elective_ga$median_7))
elective_ga[c(elective_ga_na, elective_ga_na_med),] <- elective_ga[c(elective_ga_na, elective_ga_na_med), "day_7"]

write.csv(elective_ga,
          file = "D:/Dropbox/COVID19/Overflow/Wolfram/cohort_1_ga_all_ICDS_5_7_2020.csv",
          row.names = FALSE,
          quote = FALSE)

elective_ga_na <- which(is.na(elective_cc_2$mean_7))
elective_ga_na_med <- which(is.na(elective_cc_2$median_7))
elective_cc_2[elective_ga_na,c("mean_7","median_7")] <- elective_cc[c(elective_ga_na,
                                                                    elective_ga_na_med), "day_7"]
write.csv(elective_cc,
          file = "D:/Dropbox/COVID19/Overflow/Wolfram/cohort_1_cc_all_ICDS_5_7_2020.csv",
          row.names = FALSE,
          quote = FALSE)


icd_7 <- transitions_data[transitions_data$ICD == 7,]
source(file = "D:/Dropbox/COVID19/Overflow/Rscript/R/survival analysis.R")
test_ICD_3_sa <- survival_analysis_set_up(icd_7,base_dir = "D:/Dropbox/COVID19/Overflow/",
                                         core_18 = TRUE, elective_run = TRUE, emergency_run = FALSE, crr_try = TRUE,
                                         single_ICD = TRUE, single_icd = 7)


## time_series ##

source(file = "D:/Dropbox/COVID19/Overflow/Rscript/R/time_series_creator.R")
time_series_data <- time_series_creator(hes_data = transitions_data, num_cores = 12, forecast_date = "2012-03-01",
                                        emergency_run = TRUE, elective_ts = TRUE, forecast_cutoff = "2012-03-05")

waiting_pool <- time_series_data[[3]]
write.csv(waiting_pool,
          file = "D:/Dropbox/COVID19/Overflow/Wolfram/waiting_pool.csv",
          row.names = FALSE,
          quote = FALSE)

source(file = "D:/Dropbox/COVID19/Overflow/Rscript/R/time_series_forecast.R")
times_series_forecasts <- running_forecasts(total_cohort_data = time_series_data, train_date = as.Date("2012-03-01"),
                                       forecast_period = 52,  base_dir = "D:/Dropbox/COVID19/Overflow/Results/TimeSeries/",
                                       run_admis = TRUE)



save(times_series_forecasts, file = "D:/Dropbox/COVID19/Overflow/time_series_data.Rdata")


## Sum up file 

source(file = "D:/Dropbox/COVID19/Overflow/Rscript/R/final_output_creator.R")

pi_y_df <- sum_up_function(SA_res = survival_res, time_series_data_res = time_series_data,
                time_series_forecasts = times_series_forecasts, forecast_length = 52)


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


############################################ to be inserted into one of the source files
# calculates patients currently in hospital at start of pandemic
transitions_data_whole$disdate_MDY <- as.Date(transitions_data_whole$disdate_MDY, format = "%d%b%Y")

transitions_data$currentpt <- ifelse(transitions_data_whole$admidate_MDY < "2012-03-05"
                                           & transitions_data_whole$disdate_MDY >= "2012-03-05", 1, 0)

y0 <- subset(transitions_data, currentpt == 1)
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

transitions_data$currentpt <- ifelse(transitions_data$rttstart < "2012-03-05"
                                           & transitions_data$admidate_MDY >= "2012-03-05" 
                                           & transitions_data$admimeth_C==1, 1, 0)

x0 <- transitions_data[transitions_data$currentpt == 1,]
x0_collapsed <- do.call(data.frame, aggregate(x0$One, list(x0$ICD, x0$agegrp_v3), sum))

colnames(x0_collapsed) <- c("ICD","agegrp_v3","x0")

x0_collapsed$disease <- "ICD"

x0_collapsed$age <- paste("_AGE",x0_collapsed$agegrp_v3)

x0_collapsed$p <- paste0(x0_collapsed$disease, x0_collapsed$ICD, x0_collapsed$age)
#x0_collapsed$S <- ifelse(y0_collapsed$cc == 1, "C", "G")
export <- select(x0_collapsed, p,x0)

write.csv(export,
          file = "D:/Dropbox/COVID19/Overflow/Wolfram/x0.csv", 
          row.names = FALSE,
          quote = FALSE)




icd_7_dat <- transitions_data[transitions_data$ICD == 7 &
                                transitions_data$cohort == 1,]

agegrp_3_2_time <- icd_7_dat$GA_LoS
agegrp_3_2_trans <- icd_7_dat$ga_transitions
agegrp_3_2_WT <- icd_7_dat$WaitingTime
current_fail <- 2


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
      
      if(!(1 %in% unique(agegrp_3_2_trans)))
        stop("1 not in agegrp 3_2 trans for some reason")
      
    }
    
    CI.agegrp1_t1 <- NULL
    print(unique(agegrp_3_2_trans))
    try(CI.agegrp1_t1 <- fastCrr(Crisk(agegrp_3_2_time, agegrp_3_2_trans, failcode = 1) ~ agegrp_3_2_WT,
                                 variance = TRUE, returnDataFrame = TRUE))
    if(length(CI.agegrp1_t1) != 0){
      
      ## Sometimes the bootstrapping in the predict.fcrr will mean with low event numbers it can't 
      ## estimate them due to using Crisk again, we'll try turning this off if we get an error
      cif22_pred <- NULL
      no_var <- FALSE
      try(cif22_pred <- fastcmprsk:::predict.fcrr(CI.agegrp1_t1, newdata = 0, tL = 1, type = "interval"), silent = TRUE)
      if(length(cif1_pred) == 0){
        try(cif22_pred <- fastcmprsk:::predict.fcrr(CI.agegrp1_t1, newdata = 0, tL = 1, type = "interval",
                                                   getBootstrapVariance = FALSE), silent = TRUE)




  }
}
    
cif1_pred    
cif2_pred    

test_crr_1 <- crr(ftime = icd_3_age_3$cc_LoS, fstatus = icd_3_age_3$cc_transitions, cov1 = icd_3_age_3$WaitingTime, failcode = 1)
test_2_crr <- crr(ftime = icd_3_age_3$cc_LoS, fstatus = icd_3_age_3$cc_transitions, cov1 = icd_3_age_3$WaitingTime, failcode = 2)

predict_1_crr <- predict.crr(test_crr_1, cov1 = 0)
predict_2_crr <- predict.crr(test_2_crr, cov1 = 0)
predict(test_crr_1, 0)
predict(test_2_crr, 0)

plot.predict.crr(predict_1_crr)

icd_3_cohort_1 <- icd_3[icd_3$cohort == 1,]
icd_3_age_3 <- icd_3_cohort_1
icd_3_cohort_1$agegrp_v3 <- icd_3_cohort_1$agegrp_v3 * 10 

CI.byagegrp <- cuminc(ftime = icd_3_cohort_1$cc_LoS, fstatus = icd_3_cohort_1$cc_transitions, group = icd_3_cohort_1$agegrp_v3)
cuminc_only <- cuminc_list_to_df(CI.byagegrp,ICD = 3)
cuminc_only[[2]]
    