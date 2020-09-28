##########################################################
###### LOS frequency tables ##############################
##########################################################

##############################################################################################################
# Required packages
require(data.table)
require(dplyr)
require(openxlsx)

## Functions ##

make_prop <- function(count_csv){
  
  ## Function to convert daily counts of LoS into proportions 
  ## count_csv - dataframe with ICD heads and count values, rownames as days of LoS
  

  day_nums <- seq(0,4)
  out_df <- data.frame(matrix(data = 0, ncol = ncol(count_csv) + 1, nrow = length(day_nums)))
  colnames(out_df) <- c("days",colnames(count_csv))
  out_df$days <- day_nums
  out_df_prop <- data.frame(matrix(data = 0, ncol = ncol(count_csv) + 1, nrow = length(day_nums)))
  colnames(out_df_prop) <- c("days",colnames(count_csv))
  out_df_prop$days <- day_nums
  count_csv$los <- as.integer(rownames(count_csv))
  count_csv <- count_csv %>% select(los, everything())
  
  for(k in 2:length(colnames(count_csv))){
    
    out_df[,k] <- count_csv[,k]
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

hes_narrow <- select(transitions_data  , c("totalLoS_cln","GA_LoS", "cc_LoS", "admimeth_C", "ICD",
                               "agegrp_v3", "cc_start_flg","cohort","cc"))
hes_narrow$disease <- "ICD"

hes_narrow$age <- paste("_AGE", hes_narrow$agegrp_v3, sep = "")

hes_narrow$ptgrp <- paste0(hes_narrow$disease, hes_narrow$ICD, hes_narrow$age)
hes_narrow$a <- ifelse(hes_narrow$admimeth_C == 1, "N", "E")

# Identify patients who have GA stay and those with CC stay 
N_ga <- subset(hes_narrow, admimeth_C == 1 )
E_ga <- subset(hes_narrow, admimeth_C == 2 )

N_cc <- subset(hes_narrow, admimeth_C == 1 & cc == 1)
E_cc <- subset(hes_narrow, admimeth_C == 2 & cc == 1)

## subset the GA & CC data for first 4 days only
N_ga <- subset(N_ga, GA_LoS >= 0 & GA_LoS <= 4)
E_ga <- subset(E_ga, GA_LoS >= 0 & GA_LoS <= 4)

N_cc <- subset(N_cc, cc_LoS >= 0 & cc_LoS <= 4)
E_cc <- subset(E_cc, cc_LoS >= 0 & cc_LoS <= 4)


# frequency and proportion tables for each category
# All electives

# Elective patients  G&A
count_N_ga <- as.data.frame.matrix(table(N_ga$GA_LoS, N_ga$ptgrp))
aggregated_dfs_elec_GA <- aggregate_weekly(count_N_ga)

# Emergency patients G&A
count_E_ga <- as.data.frame.matrix(table(E_ga$GA_LoS, E_ga$ptgrp))
aggregated_dfs_emerg_GA <- aggregate_weekly(count_E_ga)

# Elective patients CC
count_N_cc <- as.data.frame.matrix(table(N_cc$cc_LoS, N_cc$ptgrp))
aggregated_dfs_elec_CC <- aggregate_weekly(count_N_cc)

# Emergency patients CC
count_E_cc <- as.data.frame.matrix(table(E_cc$cc_LoS, E_cc$ptgrp))
aggregated_dfs_emerg_CC <- aggregate_weekly(count_E_cc)



###############################################################################
## Simplify counts to remove those under 10 ###################################
###############################################################################

## Electives starting in G&A 
cut_off_elec_GA <- count_threshold(aggregated_dfs_elec_GA, cut_off = 10)

## Emergencies starting in G&A
cut_off_emerg_GA <- count_threshold(aggregated_dfs_emerg_GA, cut_off = 10)

## Electives starting in CC
cut_off_elec_CC <- count_threshold(aggregated_dfs_elec_CC, cut_off = 10)

## Emergencies starting in CC
cut_off_emerg_CC <- count_threshold(aggregated_dfs_emerg_CC, cut_off = 10)


df_list <- list("emergency_ga_count" = cut_off_emerg_GA[[1]] ,"emergency_ga_prop" = cut_off_emerg_GA[[2]],
                "emergency_cc_count" = cut_off_emerg_CC[[1]], "emergency_cc_prop" = cut_off_emerg_CC[[2]],
                "elective_ga_count" = cut_off_elec_GA[[1]], "elective_ga_prop" = cut_off_elec_GA[[2]],
                "elective_cc_count" = cut_off_elec_CC[[1]], "elective_cc_prop" = cut_off_elec_CC[[2]])

write.xlsx(df_list, file = "D:/Dropbox/COVID19/Overflow/Distributions/totalLoS_episodes_17_09_2020.xlsx")





