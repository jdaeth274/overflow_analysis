###############################################################################
## Transitions labelling in HES ###############################################
###############################################################################



###############################################################################
## Writing out the transitions function #######################################
###############################################################################

hes_transitions <- function(hes_dataset){
  
  if(tibble::is_tibble(hes_dataset)){
    hes_dataset <- as.data.frame(hes_dataset)
  }
  
  ## cohort 3 WT 
  if(!("WaitingTime" %in% colnames(hes_dataset))){
  
    print("Getting the WT data")
    hes_dataset$WaitingTime <- as.numeric(admidate_MDY - rttstart)
    elective_df <- hes_dataset[hes_dataset$cohort != 3,]
    one_year_under <- which(as.integer(elective_df$WaitingTime) <= 365)
    elective_df$WaitingTime[one_year_under] <- elective_df$WaitingTime[one_year_under]
    hes_dataset[hes_dataset$cohort == 1, "WaitingTime"] <- elective_df$WaitingTime
    
  
  }
  
  ## cohort 12 waiting time 
  print("Getting the WT data cohort 2")
  cohort2_indies <- which(hes_dataset$cohort == 2)
  wait_time <- hes_dataset$admidate_MDY[cohort2_indies] - hes_dataset$rttstart[cohort2_indies]
  hes_dataset$WaitingTime[cohort2_indies] <- wait_time
  
  if(!("WT" %in% colnames(hes_dataset))){
    hes_dataset$WT <- hes_dataset$WaitingTime
  }else{
    hes_dataset$WT_old <- hes_dataset$WT
    hes_dataset$WT <- hes_dataset$WaitingTime
  }
  
  
  
  
  
  ## Elective to Emergency 
  print("Electives to emergencies")
  hes_dataset$Elective2Emergency <- 0
  hes_dataset$Elective2Emergency[which(hes_dataset$cohort == 2)] <- 1
  
  ## ga to cc ##
  print("GA > CC")
  hes_dataset$ga_cc <- 0
  hes_dataset$ga_cc[which(hes_dataset$cc == 1)] <- 1
  odd_dismeth <- which(is.na(hes_dataset$dismeth))
  hes_dataset$ga_cc[odd_dismeth] <- NA
  
  ## CC to G&A ##
  print("CC > GA")
  hes_dataset$cc_ga <- 0
  hes_dataset[(hes_dataset$ccdisdest %in% c(1,2,3)) & (hes_dataset$cc == 1) & (hes_dataset$cc_dis_flg != 1),"cc_ga"] <- 1
  hes_dataset[!(is.na(hes_dataset$ccdisdest)) & (hes_dataset$ccdisdest > 3) & (hes_dataset$cc == 1), "cc_ga"] <- 0
  hes_dataset[hes_dataset$cc == 0, "cc_ga"] <- NA
  
  ## transitions to the grave ##
  print("morir")  
  hes_dataset[!is.na(hes_dataset$dismeth) & hes_dataset$dismeth == 5, "dismeth"] <- NA
  
  hes_dataset$death <- 0
  hes_dataset$death[which(hes_dataset$dismeth == 4 | hes_dataset$ccdisdest == 6)] <- 1
  hes_dataset$death[which(is.na(hes_dataset$dismeth))] <- NA
  
  ## cc deaths ##
  print("CC deaths")
  hes_dataset$cc_death <- 0
  hes_dataset[(hes_dataset$dismeth == 4) & (hes_dataset$cc == 1) & (hes_dataset$cc_dis_flg == 1),"cc_death"] <- 1
  hes_dataset[!is.na(hes_dataset$ccdisdest) & (hes_dataset$ccdisdest == 6) & (hes_dataset$cc == 1) & (hes_dataset$cc_dis_flg == 1),"cc_death"] <- 1
  hes_dataset[hes_dataset$cc == 0 | is.na(hes_dataset$death),"cc_death" ] <- NA
  
  ## ga deaths ##
  print("GA deaths")
  hes_dataset$ga_death <- 0
  hes_dataset[!is.na(hes_dataset$dismeth) &hes_dataset$dismeth == 4 & hes_dataset$cc == 0, "ga_death"] <- 1
  hes_dataset$ga_death[which(is.na(hes_dataset$death))] <- NA
  
  ## ga to cc then die in ga ##
  print("GA>CC>GA deaths")
  hes_dataset$gaccga_death <- NA
  hes_dataset[!is.na(hes_dataset$dismeth) & !is.na(hes_dataset$ccdisdest) & (hes_dataset$dismeth == 4 | hes_dataset$ccdisdest == 6) & hes_dataset$cc == 1 & hes_dataset$cc_dis_flg != 1, "gaccga_death"] <- 1
  
  ## ga recover ##
  
  hes_dataset$ga_recover <- 0
  hes_dataset[!is.na(hes_dataset$dismeth) & (hes_dataset$dismeth %in% c(1,2,3)) & hes_dataset$ga_cc == 0, "ga_recover"] <- 1
  hes_dataset$ga_recover[which(is.na(hes_dataset$dismeth))] <- NA

  
  ## ga to cc back to gg then recover ##
  
  hes_dataset$gaccga_recover <- NA
  hes_dataset[!is.na(hes_dataset$dismeth) & !is.na(hes_dataset$cc_ga) & (hes_dataset$dismeth %in% c(1,2,3)) & hes_dataset$cc_ga == 1,"gaccga_recover"] <- 1
  
  ## CC to recover ##
  
  hes_dataset$cc_recover <- 0
  hes_dataset[!is.na(hes_dataset$ccdisdest) &(hes_dataset$ccdisdest %in% c(4,5)) & (hes_dataset$cc_dis_flg == 1 | hes_dataset$cc_ga == 0), "cc_recover"] <- 1
  hes_dataset$cc_recover[which(hes_dataset$cc == 0)] <- NA
  
  ## competing risk variables GA TRANSITIONS ##
  
  hes_dataset$ga_transitions <- NA
  hes_dataset$ga_transitions[which(hes_dataset$ga_recover == 1)] <- 1
  hes_dataset$ga_transitions[which(hes_dataset$ga_cc == 1)] <- 2
  hes_dataset$ga_transitions[which(hes_dataset$ga_death == 1)] <- 3
  
  ## cc transitions 
  
  hes_dataset$cc_transitions <- NA
  hes_dataset$cc_transitions[which(hes_dataset$cc_recover == 1)] <- 1
  hes_dataset$cc_transitions[which(hes_dataset$cc_ga == 1)] <- 2
  hes_dataset$cc_transitions[which(hes_dataset$cc_death == 1)] <- 3
  
  ## straight to cc 
  
  hes_dataset$straighttocc <- NA
  hes_dataset[hes_dataset$cc_start_flg == 1 & hes_dataset$cc == 1, "straighttocc"] <- 1
  hes_dataset[hes_dataset$cc_start_flg == 0 | hes_dataset$cc == 0, "straighttocc"] <- 0
  
  ## update LOS
  hes_dataset$GA_LoS_old <- hes_dataset$GA_LoS
  hes_dataset$GA_LoS <- ga_length(hes_dataset)
  hes_dataset <- hes_dataset[hes_dataset$GA_LoS >= 0,]
  
  return(hes_dataset)
  
}

ga_length <- function(hes_data){
  
  ## take this as the epistart and epiend, due to patients with multiple episodes in one admission
  
  hes_data <- hes_data[,c("GA_LoS","cc_LoS","epiend_str","epistart_str","cc","ccstartdate_MDY","ccdisdate_MDY")]
  
  hes_data$GA_LoS_old <- hes_data$GA_LoS
  hes_data$cc_LoS_old <- hes_data$cc_LoS
  
  ## 
  
  hes_data$GA_LoS <- ifelse(hes_data$cc == 1, as.Date(hes_data$ccstartdate_MDY,format = "%d%b%Y") - hes_data$epistart_str, hes_data$epiend_str - hes_data$epistart_str)
  
  return(hes_data$GA_LoS)  
  
}










