###############################################################################
## Transitions labelling in HES ###############################################
###############################################################################

hes_dataset <- read.csv("E:/HES/COVID/HES_APC_CC_0913_TEMP02_s1.csv",
                        stringsAsFactors = FALSE)


###############################################################################
## Writing out the transitions function #######################################
###############################################################################

hes_transitions <- function(hes_dataset){

  ## ga to cc ##
  
  hes_dataset$ga_cc <- 0
  hes_dataset$ga_cc[which(hes_dataset$cc == 1)] <- 1
  odd_dismeth <- which(is.na(hes_dataset$dismeth))
  hes_dataset$ga_cc[odd_dismeth] <- NA
  
  ## CC to G&A ##
  
  hes_dataset$cc_ga <- 0
  hes_dataset[(hes_dataset$ccdisdest %in% c(1,2,3)) & (hes_dataset$cc == 1) & (hes_dataset$cc_dis_flg != 1),"cc_ga"] <- 1
  hes_dataset[!(is.na(hes_dataset$ccdisdest)) & (hes_dataset$ccdisdest > 3) & (hes_dataset$cc == 1), "cc_ga"] <- 0
  hes_dataset[hes_dataset$cc == 0, "cc_ga"] <- NA
  
  ## transitions to the grave ##
  
  hes_dataset[!is.na(hes_dataset$dismeth) & hes_dataset$dismeth == 5, "dismeth"] <- NA
  
  hes_dataset$death <- 0
  hes_dataset$death[which(hes_dataset$dismeth == 4 | hes_dataset$ccdisdest == 6)] <- 1
  hes_dataset$death[which(is.na(hes_dataset$dismeth))] <- NA
  
  ## cc deaths ##
  
  hes_dataset$cc_death <- 0
  hes_dataset[(hes_dataset$dismeth == 4) & (hes_dataset$cc == 1) & (hes_dataset$cc_dis_flg == 1),"cc_death"] <- 1
  hes_dataset[!is.na(hes_dataset$ccdisdest) & (hes_dataset$ccdisdest == 6) & (hes_dataset$cc == 1) & (hes_dataset$cc_dis_flg == 1),"cc_death"] <- 1
  hes_dataset[hes_dataset$cc == 0 | is.na(hes_dataset$death),"cc_death" ] <- NA
  
  ## ga deaths ##
  
  hes_dataset$ga_death <- 0
  hes_dataset[!is.na(hes_dataset$dismeth) &hes_dataset$dismeth == 4 & hes_dataset$cc == 0, "ga_death"] <- 1
  hes_dataset$ga_death[which(is.na(hes_dataset$death))] <- NA
  
  ## ga to cc then die in ga ##
  
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
  
  return(hes_dataset)
  
}

hes_data_avec_trans <- hes_transitions(hes_dataset)










