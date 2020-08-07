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

hes_1112 <- select(hes_1112, c("totalLoS_cln","GA_LoS", "cc_LoS", "admimeth_C", "ICD", "agegrp_v3", "cc_start_flg"))

hes_1112 <- subset(hes_1112, totalLoS_cln >= 0 | hes_1112$GA_LoS >= 0 | hes_1112$cc_LoS >= 0)

hes_1112$disease <- "ICD"

hes_1112$age <- paste("_AGE", hes_1112$agegrp, sep = "")

hes_1112$ptgrp <- paste0(hes_1112$disease, hes_1112$ICD, hes_1112$age)
hes_1112$a <- ifelse(hes_1112$admimeth == 1, "N", "E")

# Identify patients in Elective vs. Emergency
elective = subset(hes_1112, admimeth_C == 1)
emergency = subset(hes_1112, admimeth_C == 2)

# Identify patients who start in GA vs. CC and Elective vs. Emergency
N_ga = subset(hes_1112, admimeth_C == 1 & cc_start_flg == 0)
E_ga = subset(hes_1112, admimeth_C == 2 & cc_start_flg == 0)

N_cc = subset(hes_1112, admimeth_C == 1 & cc_start_flg == 1)
E_cc = subset(hes_1112, admimeth_C == 2 & cc_start_flg == 1)

# frequency and proportion tables for each category
# All electives
count_elective <- as.data.frame.matrix(table(elective$totalLoS_cln, elective$ptgrp))
prop_elective <- as.data.frame.matrix(prop.table(table(elective$totalLoS_cln, elective$ptgrp), 2)*100)

# All emergencies
count_emergency <- as.data.frame.matrix(table(emergency$totalLoS_cln, emergency$ptgrp))
prop_emergency <- as.data.frame.matrix(prop.table(table(emergency$totalLoS_cln, emergency$ptgrp), 2)*100)

# Elective patients starting in G&A
count_N_ga <- as.data.frame.matrix(table(N_ga$totalLoS_cln, N_ga$ptgrp))
prop_N_ga <- as.data.frame.matrix(prop.table(table(N_ga$totalLoS_cln, N_ga$ptgrp), 2)*100)

# Emergency patients starting in G&A
count_E_ga <- as.data.frame.matrix(table(E_ga$totalLoS_cln, E_ga$ptgrp))
prop_E_ga <- as.data.frame.matrix(prop.table(table(E_ga$totalLoS_cln, E_ga$ptgrp), 2)*100)

# Elective patients starting in CC
count_N_cc <- as.data.frame.matrix(table(N_cc$totalLoS_cln, N_cc$ptgrp))
prop_N_cc <- as.data.frame.matrix(prop.table(table(N_cc$totalLoS_cln, N_cc$ptgrp), 2)*100)

# Emergency patients starting in CC
count_E_cc <- as.data.frame.matrix(table(E_cc$totalLoS_cln, E_cc$ptgrp))
prop_E_cc <- as.data.frame.matrix(prop.table(table(E_cc$totalLoS_cln, E_cc$ptgrp), 2)*100)

# Write these all out to csvs
# All electives
write.csv(count_elective,
          file = "D:/Overflows/output/count_elective.csv",
          row.names = FALSE,
          quote = FALSE)
write.csv(prop_elective,
          file = "D:/Overflows/output/prop_elective.csv",
          row.names = FALSE,
          quote = FALSE)
# All emergencies
write.csv(count_emergency,
          file = "D:/Overflows/output/count_emergency.csv",
          row.names = FALSE,
          quote = FALSE)
write.csv(prop_emergency,
          file = "D:/Overflows/output/prop_emergency.csv",
          row.names = FALSE,
          quote = FALSE)
# Elective patients starting in G&A
write.csv(count_N_ga,
          file = "D:/Overflows/output/count_N_ga.csv",
          row.names = FALSE,
          quote = FALSE)
write.csv(prop_N_ga,
          file = "D:/Overflows/output/prop_N_ga.csv",
          row.names = FALSE,
          quote = FALSE)
# Emergency patients starting in G&A
write.csv(count_E_ga,
          file = "D:/Overflows/output/count_E_ga.csv",
          row.names = FALSE,
          quote = FALSE)
write.csv(prop_E_ga,
          file = "D:/Overflows/output/prop_E_ga.csv",
          row.names = FALSE,
          quote = FALSE)
# Elective patients starting in CC
write.csv(count_N_cc,
          file = "D:/Overflows/output/count_N_cc.csv",
          row.names = FALSE,
          quote = FALSE)
write.csv(prop_N_cc,
          file = "D:/Overflows/output/prop_N_cc.csv",
          row.names = FALSE,
          quote = FALSE)
# Emergency patients starting in CC
write.csv(count_E_cc,
          file = "D:/Overflows/output/count_E_cc.csv",
          row.names = FALSE,
          quote = FALSE)
write.csv(prop_E_cc,
          file = "D:/Overflows/output/prop_E_cc.csv",
          row.names = FALSE,
          quote = FALSE)

