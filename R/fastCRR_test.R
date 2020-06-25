require(cmprsk)
require(devtools)

install_github("erickawaguchi/fastcmprsk", ref = "developer")
require(fastcmprsk)

## Data here is cc transitions for neoplasm patients ##
neo_trans_data_1 <- read.csv()

## cmprsk crr run ##
cc_failcode_3 <- crr(ftime = neo_trans_data_1$cc_LoS, fstatus = neo_trans_data_1$cc_transitions,
                     cov1 = neo_trans_data_1$WaitingTime, failcode = 3)
cc_failcode_3$coef

## Setting up run for fastcmprsk, removing NAs from ftime (cc_LoS) & fstatus (cc_transitions)

agegrp_3_2_time <- neo_trans_data_1$cc_LoS
agegrp_3_2_trans <- neo_trans_data_1$cc_transitions
agegrp_3_2_WT <- neo_trans_data_1$WaitingTime
current_fail <- 3

## remove values with no Length of stay in CC
na_tims <- which(is.na(agegrp_3_2_time))
if(length(na_tims) > 0){
  agegrp_3_2_time <- agegrp_3_2_time[-na_tims]
  agegrp_3_2_trans <- agegrp_3_2_trans[-na_tims]
  agegrp_3_2_WT <- agegrp_3_2_WT[-na_tims]
}

## Step above should catch all patients who didn't go to cc 
## but if below catches those that might somehow have been miscoded

old_na <- which(is.na(agegrp_3_2_trans))
if(length(old_na) > 0){
  agegrp_3_2_time <- agegrp_3_2_time[-old_na]
  agegrp_3_2_trans <- agegrp_3_2_trans[-old_na]
  agegrp_3_2_WT <- agegrp_3_2_WT[-old_na]
}

## Initially when using fastCrr, would get odd ouput if the failcode input for Crisk
## wasn't 1, this step takes current fail code and switches if with 1, if not 1 alread

if(current_fail != 1){
  old_2ers <- which(agegrp_3_2_trans == current_fail)
  old_1ers <- which(agegrp_3_2_trans == 1)
  
  agegrp_3_2_trans[old_2ers] <- 1
  agegrp_3_2_trans[old_1ers] <- current_fail
  
}

## Running fast CRR

cc_failcode_3_fast <- fastCrr(Crisk(agegrp_3_2_time, agegrp_3_2_trans, failcode = 1) ~ agegrp_3_2_WT,
                              variance = TRUE, returnDataFrame = TRUE)
