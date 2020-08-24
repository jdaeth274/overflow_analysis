##########################################################
###### LOS frequency tables ##############################
##########################################################

##############################################################################################################
# Required packages
require(data.table)
require(dplyr)
##############################################################################################################
## Load HES data
system.time(transitions_data_whole <- fread("D:/Overflows/data/HES_APC_CC_0913_transitions_all_ICD.csv"))

# Take only 1 year of data (11-12 in this case)
hes_1112 <- subset(transitions_data_whole, EpiEnd_FY == 1112)

hes_1112 <- select(hes_1112, c("totalLoS_cln","GA_LoS", "cc_LoS", "admimeth_C", "ICD",
                               "agegrp_v3", "cc_start_flg","cohort","epistart_str",
                               "epiend_str","ccstartdate_str"))

## Calculate alt_los as epiend - epistart discard negative values and limit to <= 75
hes_1112$alt_LOS <- as.integer(hes_1112$epiend_str - hes_1112$epistart_str)
hes_1112 <- subset(hes_1112, totalLoS_cln >= 0 | hes_1112$GA_LoS >= 0 | hes_1112$cc_LoS >= 0 | alt_LOS >= 0)
hes_1112 <- subset(hes_1112, alt_LOS >= 0 & alt_LOS <= 75)

## Calculate just GA LoS again, this is entirety of GA stay and assumes one CC stay per episode.

hes_1112$alt_GA <- ifelse(is.na(hes_1112$cc_LoS),hes_1112$alt_LOS, hes_1112$alt_LOS - hes_1112$cc_LoS)
hes_1112 <- subset(hes_1112, alt_GA >= 0 & alt_GA <= 75)

## Calculate CC LoS as CCstart to epiend date for those starting in CC 

hes_1112$alt_CC <- ifelse(is.na(hes_1112$cc_LoS), hes_1112$alt_LOS, hes_1112$epiend_str - hes_1112$ccstartdate_str)
hes_1112 <- subset(hes_1112, alt_CC >= 0 & alt_CC <= 75)

hes_1112$disease <- "ICD"

hes_1112$age <- paste("_AGE", hes_1112$agegrp, sep = "")

hes_1112$ptgrp <- paste0(hes_1112$disease, hes_1112$ICD, hes_1112$age)
hes_1112$a <- ifelse(hes_1112$admimeth == 1, "N", "E")

# Identify patients in Elective vs. Emergency
elective <- subset(hes_1112, admimeth_C == 1)
emergency <- subset(hes_1112, admimeth_C == 2)

# Identify patients who start in GA vs. CC and Elective vs. Emergency
N_ga <- subset(hes_1112, admimeth_C == 1 & cc_start_flg == 0)
E_ga <- subset(hes_1112, admimeth_C == 2 & cc_start_flg == 0)

N_cc <- subset(hes_1112, admimeth_C == 1 & cc_start_flg == 1)
E_cc <- subset(hes_1112, admimeth_C == 2 & cc_start_flg == 1)

# frequency and proportion tables for each category
# All electives
count_elective <- as.data.frame.matrix(table(elective$alt_LOS, elective$ptgrp))
prop_elective <- as.data.frame.matrix(prop.table(table(elective$alt_LOS, elective$ptgrp), 2)*100)

# All emergencies
count_emergency <- as.data.frame.matrix(table(emergency$alt_LOS, emergency$ptgrp))
prop_emergency <- as.data.frame.matrix(prop.table(table(emergency$alt_LOS, emergency$ptgrp), 2)*100)

# Elective patients starting in G&A
count_N_ga <- as.data.frame.matrix(table(N_ga$alt_LOS, N_ga$ptgrp))
prop_N_ga <- as.data.frame.matrix(prop.table(table(N_ga$alt_LOS, N_ga$ptgrp), 2)*100)

# Emergency patients starting in G&A
count_E_ga <- as.data.frame.matrix(table(E_ga$alt_LOS, E_ga$ptgrp))
prop_E_ga <- as.data.frame.matrix(prop.table(table(E_ga$alt_LOS, E_ga$ptgrp), 2)*100)

# Elective patients starting in CC
count_N_cc <- as.data.frame.matrix(table(N_cc$alt_LOS, N_cc$ptgrp))
prop_N_cc <- as.data.frame.matrix(prop.table(table(N_cc$alt_LOS, N_cc$ptgrp), 2)*100)

count_N_cc2 <- as.data.frame.matrix(table(N_cc$cc_LoS, N_cc$ptgrp))
prop_N_cc2 <- as.data.frame.matrix(prop.table(table(N_cc$cc_LoS, N_cc$ptgrp), 2)*100)

# Emergency patients starting in CC
count_E_cc <- as.data.frame.matrix(table(E_cc$alt_LOS, E_cc$ptgrp))
prop_E_cc <- as.data.frame.matrix(prop.table(table(E_cc$alt_LOS, E_cc$ptgrp), 2)*100)

# Write these all out to csvs
# All electives
write.csv(count_elective,
          file = "D:/Dropbox/COVID19/Overflow/JOSH/4_8_data/count_elective_epi.csv",
          row.names = FALSE,
          quote = FALSE)
write.csv(prop_elective,
          file = "D:/Dropbox/COVID19/Overflow/JOSH/4_8_data/prop_elective_epi.csv",
          row.names = FALSE,
          quote = FALSE)
# All emergencies
write.csv(count_emergency,
          file = "D:/Dropbox/COVID19/Overflow/JOSH/4_8_data/count_emergency_epi.csv",
          row.names = FALSE,
          quote = FALSE)
write.csv(prop_emergency,
          file = "D:/Dropbox/COVID19/Overflow/JOSH/4_8_data/prop_emergency_epi.csv",
          row.names = FALSE,
          quote = FALSE)
# Elective patients starting in G&A
write.csv(count_N_ga,
          file = "D:/Dropbox/COVID19/Overflow/JOSH/4_8_data/count_N_ga_epi.csv",
          row.names = FALSE,
          quote = FALSE)
write.csv(prop_N_ga,
          file = "D:/Dropbox/COVID19/Overflow/JOSH/4_8_data/prop_N_ga_epi.csv",
          row.names = FALSE,
          quote = FALSE)
# Emergency patients starting in G&A
write.csv(count_E_ga,
          file = "D:/Dropbox/COVID19/Overflow/JOSH/4_8_data/count_E_ga_epi.csv",
          row.names = FALSE,
          quote = FALSE)
write.csv(prop_E_ga,
          file = "D:/Dropbox/COVID19/Overflow/JOSH/4_8_data/prop_E_ga_epi.csv",
          row.names = FALSE,
          quote = FALSE)
# Elective patients starting in CC
write.csv(count_N_cc,
          file = "D:/Dropbox/COVID19/Overflow/JOSH/4_8_data/count_N_cc_epi.csv",
          row.names = FALSE,
          quote = FALSE)
write.csv(prop_N_cc,
          file = "D:/Dropbox/COVID19/Overflow/JOSH/4_8_data/prop_N_cc_epi.csv",
          row.names = FALSE,
          quote = FALSE)
# Emergency patients starting in CC
write.csv(count_E_cc,
          file = "D:/Dropbox/COVID19/Overflow/JOSH/4_8_data/count_E_cc_epi.csv",
          row.names = FALSE,
          quote = FALSE)
write.csv(prop_E_cc,
          file = "D:/Dropbox/COVID19/Overflow/JOSH/4_8_data/prop_E_cc_epi.csv",
          row.names = FALSE,
          quote = FALSE)



# colnames(count_emergency) <- paste(colnames(count_emergency),"_count",sep = "")
# colnames(prop_emergency) <- paste(colnames(prop_emergency),"_proportion",sep = "")
# 
# colnames(prop_emergency) <- paste(colnames(prop_emergency),"_proportion",sep = "")
# colnames(count_elective) <- paste(colnames(count_elective),"_count",sep = "")
# 
# colnames(count_N_ga) <- paste(colnames(count_N_ga),"_count",sep = "")
# colnames(prop_N_ga) <- paste(colnames(prop_N_ga),"_proportion",sep = "")
# 
# colnames(count_N_cc) <- paste(colnames(count_N_cc),"_count",sep = "")
# colnames(prop_N_cc) <- paste(colnames(prop_N_cc),"_proportion",sep = "")
# 
# colnames(count_E_ga) <- paste(colnames(count_E_ga),"_count",sep = "")
# colnames(prop_E_ga) <- paste(colnames(prop_E_ga),"_proportion",sep = "")
# 
# colnames(count_E_cc) <- paste(colnames(count_E_cc),"_count",sep = "")
# colnames(prop_E_cc) <- paste(colnames(prop_E_cc),"_proportion",sep = "")

count_emergency$LoS <- rownames(count_emergency)
count_emergency <- count_emergency %>% select(LoS, everything())

prop_emergency$LoS <- rownames(prop_emergency)
prop_emergency <- prop_emergency %>% select(LoS, everything())

count_elective$LoS <- rownames(count_elective)
count_elective <- count_elective %>% select(LoS, everything())

prop_elective$LoS <- rownames(prop_elective)
prop_elective <- prop_elective %>% select(LoS, everything())

count_N_ga$LoS <- rownames(count_N_ga)
count_N_ga <- count_N_ga %>% select(LoS, everything())

prop_N_ga$LoS <- rownames(prop_N_ga)
prop_N_ga <- prop_N_ga %>% select(LoS, everything())

count_N_cc$LoS <- rownames(count_N_cc)
count_N_cc <- count_N_cc %>% select(LoS, everything())

prop_N_cc$LoS <- rownames(prop_N_cc)
prop_N_cc <- prop_N_cc %>% select(LoS, everything())

count_E_ga$LoS <- rownames(count_E_ga)
count_E_ga <- count_E_ga %>% select(LoS, everything())

prop_E_ga$LoS <- rownames(prop_E_ga)
prop_E_ga <- prop_E_ga %>% select(LoS, everything())

count_E_cc$LoS <- rownames(count_E_cc)
count_E_cc <- count_E_cc %>% select(LoS, everything())

prop_E_cc$LoS <- rownames(prop_E_cc)
prop_E_cc <- prop_E_cc %>% select(LoS, everything())


# all_emergencies <- dplyr::left_join(count_emergency, prop_emergency, by = c("LoS" = "LoS"))
# all_emergencies <- all_emergencies %>% select(LoS, everything())
# 
# all_elective <- dplyr::left_join(count_elective, prop_elective, by = c("LoS" = "LoS"))
# all_elective <- all_elective %>% select(LoS, everything())
# 
# elective_ga <- dplyr::left_join(count_N_ga, prop_N_ga, by = c("LoS" = "LoS"))
# elective_ga <- elective_ga %>% select(LoS, everything())
# 
# electives_cc <- dplyr::left_join(count_N_cc, prop_N_cc, by = c("LoS" = "LoS"))
# electives_cc <- electives_cc %>% select(LoS, everything())
# 
# emergency_ga <- dplyr::left_join(count_E_ga, prop_E_ga, by = c("LoS" = "LoS"))
# emergency_ga <- emergency_ga %>% select(LoS, everything())
# 
# emergency_cc <- dplyr::left_join(count_E_cc, prop_E_cc, by = c("LoS" = "LoS"))
# emergency_cc <- emergency_cc %>% select(LoS, everything())



df_list <- list("all_emergency_count" = count_emergency, "all_emergency_prop" = prop_emergency,
                "emergency_ga_count" = count_E_ga ,"emergency_ga_prop" = prop_E_ga,
                "emergency_cc_count" = count_E_cc, "emergency_cc_prop" = prop_E_cc,
                "all_elective_count" = count_elective, "all_elective_prop" = prop_elective,
                "elective_ga_count" = count_N_ga, "elective_ga_prop" = prop_N_ga,
                "elective_cc_count" = count_N_cc, "elective_cc_prop" = prop_N_cc)

write.xlsx(df_list, file = "D:/Dropbox/COVID19/Overflow/Distributions/totalLoS_episodes_19_08_2020.xlsx")





