##########################################################
###### LOS frequency tables ##############################
##########################################################

##############################################################################################################
# Required packages
require(data.table)
require(dplyr)

## Functions ##

aggregate_weekly <- function(count_csv){
  
  ## Function to convert daily counts of LoS into 
  ## 3.5, 10.5, 17.5 ... to 73.5 days of stay 
  ## Assume admitted on Monday:
  ## Mon Tue Wed Thu Fri Sat Sun
  ##  0   1   2   3   4   5   6   days in data
  ##  1   2   3   4   5   6   7   days of stay so split Thursday into 2 (split 3 in data)
  ## count_csv - dataframe with ICD heads and count values, rownames as days of LoS
  

  day_nums <- seq(3.5, 75, 7)
  out_df <- data.frame(matrix(data = 0, ncol = ncol(count_csv) + 1, nrow = length(day_nums)))
  colnames(out_df) <- c("days",colnames(count_csv))
  out_df$days <- day_nums
  out_df_prop <- data.frame(matrix(data = 0, ncol = ncol(count_csv) + 1, nrow = length(day_nums)))
  colnames(out_df_prop) <- c("days",colnames(count_csv))
  out_df_prop$days <- day_nums
  count_csv$los <- as.integer(rownames(count_csv))
  count_csv <- count_csv %>% select(los, everything())
  
  for(k in 2:length(colnames(count_csv))){
    
    base_val <- 0
    for(j in 1:nrow(out_df)){
      
      upper_lim <- floor(out_df$days[j]) - 1
      half_day <- floor(out_df$days[j])
      out_df[j,k] <- sum(count_csv[count_csv$los >= base_val & count_csv$los <= upper_lim, k])
      out_df[j,k] <- out_df[j,k] + floor(count_csv[count_csv$los == half_day, k]/2)
      if(base_val != 0){
        out_df[j,k] <- out_df[j,k] + ceiling(count_csv[count_csv$los == floor(base_val), k]/2)
      }
      
      base_val <- half_day + 0.5
    }
    
    out_df_prop[,k] <- out_df[,k] / sum(out_df[,k])
    
  }
  
  return(list(count = out_df, prop = out_df_prop))
  
}

count_threshold <- function(aggregated_list, cut_off = 10){
  
  ## function to further aggregate groups where values fall under 10 or separate cutoff
  ## Input: Aggregated list: List object from aggregate_weekly function
  ##        cut_off:         Minimum number of values for which to further group (default 10)
  ## Output: List with count and proportion dfs available.
  
  ## Iteratively go through each column checking sums 
  
}

##############################################################################################################
## Load HES data
system.time(transitions_data <- fread("D:/Overflows/data/HES_APC_CC_0913_transitions_all_ICD.csv"))

# Take only 1 year of data (11-12 in this case)
hes_1112 <- subset(transitions_data, EpiEnd_FY == 1112)

hes_1112 <- select(hes_1112, c("totalLoS_cln","GA_LoS", "cc_LoS", "admimeth_C", "ICD",
                               "agegrp_v3", "cc_start_flg","cohort","epistart_str",
                               "epiend_str","ccstartdate_str"))

## Calculate alt_los as epiend - epistart discard negative values and limit to <= 75
hes_1112$alt_LOS <- as.integer(hes_1112$epiend_str - hes_1112$epistart_str)
hes_1112 <- subset(hes_1112, alt_LOS >= 0 & alt_LOS <= 75)

hes_1112$disease <- "ICD"

hes_1112$age <- paste("_AGE", hes_1112$agegrp_v3, sep = "")

hes_1112$ptgrp <- paste0(hes_1112$disease, hes_1112$ICD, hes_1112$age)
hes_1112$a <- ifelse(hes_1112$admimeth_C == 1, "N", "E")

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
aggregated_dfs_elec <- aggregate_weekly(count_elective)

# All emergencies
count_emergency <- as.data.frame.matrix(table(emergency$alt_LOS, emergency$ptgrp))
aggregated_dfs_emerg <- aggregate_weekly(count_emergency)

# Elective patients starting in G&A
count_N_ga <- as.data.frame.matrix(table(N_ga$alt_LOS, N_ga$ptgrp))
aggregated_dfs_elec_GA <- aggregate_weekly(count_N_ga)

# Emergency patients starting in G&A
count_E_ga <- as.data.frame.matrix(table(E_ga$alt_LOS, E_ga$ptgrp))
aggregated_dfs_emerg_GA <- aggregate_weekly(count_E_ga)

# Elective patients starting in CC
count_N_cc <- as.data.frame.matrix(table(N_cc$alt_LOS, N_cc$ptgrp))
aggregated_dfs_elec_CC <- aggregate_weekly(count_N_cc)

# Emergency patients starting in CC
count_E_cc <- as.data.frame.matrix(table(E_cc$alt_LOS, E_cc$ptgrp))
aggregated_dfs_emerg_CC <- aggregate_weekly(count_E_cc)



###############################################################################
## Simplify counts to remove those under 10 ###################################
###############################################################################






# Write these all out to csvs
# All electives



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





