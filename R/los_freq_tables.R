##########################################################
###### LOS frequency tables ##############################
##########################################################

##############################################################################################################
# Required packages
require(data.table)
require(dplyr)
require(openxlsx)

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
  
  ## Iteratively go through each column checking sums above cut_off 
  count_data <- aggregated_list[[1]]
  prop_data <- aggregated_list[[2]]
  
  for(k in 2:ncol(count_data)){
    
    aggregated <- FALSE
    for(j in 1:nrow(count_data)){
      if(!aggregated){
        if(count_data[j,k] <= cut_off){
          count_data[j,k] <- sum(count_data[j:nrow(count_data), k])
          if(j != nrow(count_data)){
            count_data[(j+1):nrow(count_data),k] <- NA
          }
          if(j > 1 & count_data[j,k] <= cut_off){
            count_data[j-1,k] <- count_data[j-1,k] + count_data[j,k]
            count_data[j,k] <- NA
          }else if(j == 1){
            count_data[,k] <- NA
          }
          aggregated <- TRUE
        }
      }
      
    }
    
    prop_data[,k] <- count_data[,k] / sum(count_data[,k], na.rm = TRUE)
    
  }
  
  return(list(count = count_data, prop = prop_data))
  
}

##############################################################################################################
## Load HES data
system.time(transitions_data <- fread("D:/Overflows/data/HES_APC_CC_0913_transitions_all_ICD.csv"))

# Take only 1 year of data (11-12 in this case)
hes_1112 <- subset(transitions_data, EpiEnd_FY == 1112)

hes_1112 <- select(hes_1112, c("totalLoS_cln","GA_LoS", "cc_LoS", "admimeth_C", "ICD",
                               "agegrp_v3", "cc_start_flg","cohort","epistart_str",
                               "epiend_str","ccstartdate_str","totalLoS"))

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

## all electives 
cut_off_elec <- count_threshold(aggregated_dfs_elec, cut_off = 10)

## All emergencies 
cut_off_emerg <- count_threshold(aggregated_dfs_emerg, cut_off = 10)

## Electives starting in G&A 
cut_off_elec_GA <- count_threshold(aggregated_dfs_elec_GA, cut_off = 10)

## Emergencies starting in G&A
cut_off_emerg_GA <- count_threshold(aggregated_dfs_emerg_GA, cut_off = 10)

## Electives starting in CC
cut_off_elec_CC <- count_threshold(aggregated_dfs_elec_CC, cut_off = 10)

## Emergencies starting in CC
cut_off_emerg_CC <- count_threshold(aggregated_dfs_emerg_CC, cut_off = 10)


df_list <- list("all_emergency_count" = cut_off_emerg[[1]], "all_emergency_prop" = cut_off_emerg[[2]],
                "emergency_ga_count" = cut_off_emerg_GA[[1]] ,"emergency_ga_prop" = cut_off_emerg_GA[[2]],
                "emergency_cc_count" = cut_off_emerg_CC[[1]], "emergency_cc_prop" = cut_off_emerg_CC[[2]],
                "all_elective_count" = cut_off_elec[[1]], "all_elective_prop" = cut_off_elec[[2]],
                "elective_ga_count" = cut_off_elec_GA[[1]], "elective_ga_prop" = cut_off_elec_GA[[2]],
                "elective_cc_count" = cut_off_elec_CC[[1]], "elective_cc_prop" = cut_off_elec_CC[[2]])

write.xlsx(df_list, file = "D:/Dropbox/COVID19/Overflow/Distributions/totalLoS_episodes_17_09_2020.xlsx")





