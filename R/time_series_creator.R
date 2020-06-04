###############################################################################
## Narrow down allocation HES data to work in time series level ###############
###############################################################################



###############################################################################
## Load up the data from the cohort allocation work ###########################
###############################################################################


###############################################################################
## Lets create the new df for the emergency admissions during the hes data ####
## period #####################################################################
###############################################################################

row_sumer <- function(row_indy, out_row, hes_df){
  print(row_indy)
  ## This function takes in the output df row
  ## then we subset the hes_df to just this time point, age group and ICD
  ## admissions is nrow of the subset 
  ## prop_frail is the sum(frail)/nrow(subset)
  
  current_icd <- out_row[row_indy, 4]
  current_age <- out_row[row_indy, 5]
  current_week <- out_row[row_indy, 2]
  current_year <- out_row[row_indy, 1]
  
  subset_hes <- hes_df[ hes_df$ICD == current_icd & 
                          hes_df$agegrp_v3 == current_age &
                          hes_df$admidate_YYYY == current_year &
                          hes_df$admidate_week == current_week, ]
  return_row <- out_row[row_indy,]
  
  return_row$Admissions <- nrow(subset_hes)
  
  if(nrow(subset_hes) > 0)
    return_row$prop_Frail <- sum(subset_hes$Frail) / nrow(subset_hes)
  else{
    return_row$prop_Frail <- NA
  }
  
  return(return_row)
  
}


emergency_df_ts <- function(out_indies, Emergency_out_df, emergencies_only){

  ## This function takes in our split emergencies data from the cluster_split function
  ## Then subsets the hes_emergency_Df to just those passed to it from cluster_split
  ## Then creates the df for output and links back to the row_sumer function to 
  ## go through the df appending data
  
  # emergencies_only <- emergencies_only[row_indies,]
  # emergencies_only$week_year <- paste(as.character(emergencies_only$admidate_week),
  #                                     as.character(emergencies_only$admidate_YYYY),
  #                                     sep = "-")
  # icds <- unique(emergencies_only$ICD)
  # age_groupings <- unique(emergencies_only$agegrp_v3)
  # week_year_combos <- unique(emergencies_only$week_year)
  # num_weeks <- length(week_year_combos)
  # num_icd <- length(icds)
  # num_ages <- 3
  # nrows_df <- num_weeks * num_ages * num_icd
  # 
  # emergency_out_Df <- data.frame(matrix(nrow = nrows_df, ncol = 6))
  # colnames(emergency_out_Df) <- c("admidate_YYYY","admidate_week",
  #                                 "Admissions","ICD","agegrp_v3",
  #                                 "prop_Frail")
  # weeks_year_split <- stringr::str_split_fixed(week_year_combos, "-",2)
  # emergency_out_Df$admidate_YYYY <- rep(as.integer(weeks_year_split[,2]),(num_ages * num_icd))
  # emergency_out_Df$admidate_week <- rep(as.integer(weeks_year_split[,1]),(num_ages * num_icd))
  # emergency_out_Df$ICD <- rep(icds, each = num_weeks * num_ages)
  # emergency_out_Df$agegrp_v3 <- rep(rep(age_groupings, each = num_weeks), num_icd)
  # 
  
  
  emergency_out_Df <- lapply(X = out_indies, FUN = row_sumer,
                             out_row = Emergency_out_df,
                             hes_df = emergencies_only)
  emergency_out_Df <- dplyr::bind_rows(emergency_out_Df)
  
  return(emergency_out_Df)
}

running_emergencies_ts_in_parallel <- function(hes_data, num_cores){
  current_env <- environment(running_emergencies_ts_in_parallel)
  
  
  ## Sets up the above functions for running in parallel on windows
  time_start <- Sys.time()
  emergencies_only <- hes_data[hes_data$cohort == 3,]
  
  emergencies_cols <- which(colnames(emergencies_only) %in% c("admidate_week", "admidate_YYYY",
                                                         "Frail","ICD","agegrp_v3"))
  emergencies_only <- emergencies_cols[,emergencies_cols]
  
  
  emergencies_only$week_year <- paste(as.character(emergencies_only$admidate_week),
                                      as.character(emergencies_only$admidate_YYYY),
                                      sep = "-")
  icds <- unique(emergencies_only$ICD)
  age_groupings <- unique(emergencies_only$agegrp_v3)
  week_year_combos <- unique(emergencies_only$week_year)
  num_weeks <- length(week_year_combos)
  num_icd <- length(icds)
  num_ages <- 3
  nrows_df <- num_weeks * num_ages * num_icd
  
  emergency_out_Df <- data.frame(matrix(nrow = nrows_df, ncol = 6))
  colnames(emergency_out_Df) <- c("admidate_YYYY","admidate_week",
                                  "Admissions","ICD","agegrp_v3",
                                  "prop_Frail")
  weeks_year_split <- stringr::str_split_fixed(week_year_combos, "-",2)
  emergency_out_Df$admidate_YYYY <- rep(as.integer(weeks_year_split[,2]),(num_ages * num_icd))
  emergency_out_Df$admidate_week <- rep(as.integer(weeks_year_split[,1]),(num_ages * num_icd))
  emergency_out_Df$ICD <- rep(icds, each = num_weeks * num_ages)
  emergency_out_Df$agegrp_v3 <- rep(rep(age_groupings, each = num_weeks), num_icd)
  emergency_out_Df <- emergency_out_Df
  
  
  emergency_rows <- seq(1, nrows_df)
  print("Setting up the parallel jobs")
  cluster_function <- snow::makeCluster(spec = num_cores)
  function_input <- snow::clusterSplit(cluster_function, emergency_rows)
  snow::clusterExport(cluster_function, "row_sumer")
  snow::clusterExport(cluster_function, "emergency_df_ts")
  print("Loading data onto parallel cores")
  
  snow::clusterExport(cluster_function, "emergencies_only", envir = environment())
  snow::clusterExport(cluster_function, "emergency_out_Df", envir = environment())
  snow::clusterEvalQ(cluster_function, library(dplyr))
  snow::clusterEvalQ(cluster_function, library(stringr))
  print("Running on cores")
  emergency_data_parallel <- snow::clusterApply(cluster_function, function_input,
                                                fun = emergency_df_ts,
                                                emergencies_only = emergencies_only,
                                                Emergency_out_df = emergency_out_Df)
  print("Finished")
  time_end <- Sys.time()
  print((time_end - time_start))
  stopCluster(cluster_function)
  
  emergency_data <- dplyr::bind_rows(emergency_data_parallel)
  emergency_data$date <- as.POSIXct( paste( 1, emergency_data$admidate_week, emergency_data$admidate_YYYY, sep = "-" ), format = "%u-%U-%Y",locale = "UK" ) 
  
  return(emergency_data)
  
}






###############################################################################
## Now for the cohort 1 elective TS creation ##################################
###############################################################################

row_sumer_elective <- function(row_indy, out_row, hes_df){
  
  ## Same as aboe row_sumer just adding in waiting times for Electives
  browser()
  current_icd <- out_row[row_indy, 4]
  current_age <- out_row[row_indy, 5]
  current_week <- out_row[row_indy, 2]
  current_year <- out_row[row_indy, 1]
  
  subset_hes <- hes_df[ hes_df$ICD == current_icd & 
                          hes_df$agegrp_v3 == current_age &
                          hes_df$rttstart_YYYY == current_year &
                          hes_df$rttstart_week == current_week, ]
  return_row <- out_row[row_indy,]
  
  return_row$Admissions <- nrow(subset_hes)
  
  if(nrow(subset_hes) > 0){
    return_row$prop_Frail <- sum(subset_hes$Frail) / nrow(subset_hes)
    return_row$p50_WT_ICDc <- median(subset_hes$WT)
    return_row$mean_WT_ICDc <- mean(subset_hes$WT)
  }else{
    return_row$prop_Frail <- NA
    return_row$p50_WT_ICDc <- NA
    return_row$mean_WT_ICDc <- NA
  }
  
  
  
  
  return(return_row)
  
}


elective_df_ts <- function(out_indies, electives_only_Df, electives_only){
  
  elective_out_Df <- lapply(X = out_indies, FUN = row_sumer_elective,
                             out_row = elective_out_Df,
                             hes_df = electives_only)
  elective_out_Df <- dplyr::bind_rows(elective_out_Df)
  
  return(elective_out_Df)
}

make_wt_variable <- function(elective_df){
  
  one_year_under <- which(as.integer(elective_df$elecdur) <= 365)
  
  elective_df$WT <- NA
  elective_df$WT[one_year_under] <- elective_df$elecdur[one_year_under]
  
  return(elective_df)
  
}


remove_wt_outliers <- function(ts_electives){
  
  keep_same_group <- which(!(ts_electives$ICD %in% c(3,6,7,10,11,12,14,50)))
  keep_same_group <- ts_electives[keep_same_group,]
  
  ## p 95 grouping 
  
  icd_7_group <- ts_electives[ts_electives$ICD == 7,]
  icd_12_group <- ts_electives[ts_electives$ICD == 12,]
  icd_14_group <- ts_electives[ts_electives$ICD == 14,]
  
  p_val_7_med <- quantile(x = icd_7_group$p50_WT_ICDc, probs = 0.95, na.rm = TRUE)
  p_val_7_mean <- quantile(x = icd_7_group$mean_WT_ICDc, probs = 0.95, na.rm = TRUE)
  
  p_val_12_med <- quantile(x = icd_12_group$p50_WT_ICDc, probs = 0.95, na.rm = TRUE)
  p_val_12_mean <- quantile(x = icd_12_group$mean_WT_ICDc, probs = 0.95, na.rm = TRUE)
  
  p_val_14_med <- quantile(x = icd_14_group$p50_WT_ICDc, probs = 0.95, na.rm = TRUE)
  p_val_14_mean <- quantile(x = icd_14_group$mean_WT_ICDc, probs = 0.95, na.rm = TRUE)
  
  icd_7_group$mean_WT_ICDc[icd_7_group$mean_WT_ICDc > p_val_7_mean] <- NA
  icd_7_group$p50_WT_ICDc[icd_7_group$p50_WT_ICDc > p_val_7_med] <- NA
  
  icd_12_group$mean_WT_ICDc[icd_12_group$mean_WT_ICDc > p_val_12_mean] <- NA
  icd_12_group$p50_WT_ICDc[icd_12_group$p50_WT_ICDc > p_val_12_med] <- NA
  
  icd_14_group$mean_WT_ICDc[icd_14_group$mean_WT_ICDc > p_val_14_mean] <- NA
  icd_14_group$p50_WT_ICDc[icd_14_group$p50_WT_ICDc > p_val_14_med] <- NA
  
  ## p96 grouping 
  
  icd_10_group <- ts_electives[ts_electives$ICD == 10,]
  icd_50_group <- ts_electives[ts_electives$ICD == 50,]
  
  p_val_10_med <- quantile(x = icd_10_group$p50_WT_ICDc, probs = 0.96, na.rm = TRUE)
  p_val_10_mean <- quantile(x = icd_10_group$mean_WT_ICDc, probs = 0.96, na.rm = TRUE)
  
  p_val_50_med <- quantile(x = icd_50_group$p50_WT_ICDc, probs = 0.96, na.rm = TRUE)
  p_val_50_mean <- quantile(x = icd_50_group$mean_WT_ICDc, probs = 0.96, na.rm = TRUE)
  
  icd_10_group$mean_WT_ICDc[icd_10_group$mean_WT_ICDc > p_val_10_mean] <- NA
  icd_10_group$p50_WT_ICDc[icd_10_group$p50_WT_ICDc > p_val_10_med] <- NA
  
  icd_50_group$mean_WT_ICDc[icd_50_group$mean_WT_ICDc > p_val_50_mean] <- NA
  icd_50_group$p50_WT_ICDc[icd_50_group$p50_WT_ICDc > p_val_50_med] <- NA
  
  ## p98 grouping 
  
  icd_6_group <- ts_electives[ts_electives$ICD == 6,]
  icd_3_group <- ts_electives[ts_electives$ICD == 3,]
  
  p_val_6_med <- quantile(x = icd_6_group$p50_WT_ICDc, probs = 0.98, na.rm = TRUE)
  p_val_6_mean <- quantile(x = icd_6_group$mean_WT_ICDc, probs = 0.98, na.rm = TRUE)
  
  p_val_3_med <- quantile(x = icd_3_group$p50_WT_ICDc, probs = 0.98, na.rm = TRUE)
  p_val_3_mean <- quantile(x = icd_3_group$mean_WT_ICDc, probs = 0.98, na.rm = TRUE)
  
  icd_6_group$mean_WT_ICDc[icd_6_group$mean_WT_ICDc > p_val_6_mean] <- NA
  icd_6_group$p50_WT_ICDc[icd_6_group$p50_WT_ICDc > p_val_6_med] <- NA
  
  icd_3_group$mean_WT_ICDc[icd_3_group$mean_WT_ICDc > p_val_3_mean] <- NA
  icd_3_group$p50_WT_ICDc[icd_3_group$p50_WT_ICDc > p_val_3_med] <- NA
  
  ## p99 grouping 
  
  icd_11_group <- ts_electives[ts_electives$ICD == 11,]
  
  p_val_11_med <- quantile(x = icd_11_group$p50_WT_ICDc, probs = 0.99, na.rm = TRUE)
  p_val_11_mean <- quantile(x = icd_11_group$mean_WT_ICDc, probs = 0.99, na.rm = TRUE)
  
  icd_11_group$mean_WT_ICDc[icd_11_group$mean_WT_ICDc > p_val_11_mean] <- NA
  icd_11_group$p50_WT_ICDc[icd_11_group$p50_WT_ICDc > p_val_11_med] <- NA
  
  ## Reform the df 
  
  electives_ts_df <- rbind.data.frame(keep_same_group,
                                      icd_7_group, icd_12_group, icd_14_group,
                                      icd_10_group, icd_50_group, 
                                      icd_6_group, icd_3_group,
                                      icd_11_group)
  stopifnot(nrow(electives_ts_df) == nrow(ts_electives))
    
  
  return(electives_ts_df)
  
  
}



running_elective_ts_in_parallel <- function(hes_data, num_cores){

    time_start <- Sys.time()
  electives_only <- hes_data[hes_data$cohort == 1,]
  
  if(!("WT" %in% colnames(hes_data)))
    electives_only <- make_wt_variable(electives_only)
  
  missing_rtt_week <- which(is.na(electives_only$rttstart_week))
  missing_rtt_year <- which(is.na(electives_only$rttstart_YYYY))
  
  electives_only <- electives_only[-c(missing_rtt_week, missing_rtt_year),]
  
  elective_cols <- which(colnames(electives_only) %in% c("rttstart_week", "rttstart_YYYY",
                                                         "WT","Frail","ICD","agegrp_v3"))
  electives_only <- electives_only[,elective_cols]
  
  
  electives_only$week_year <- paste(as.character(electives_only$rttstart_week),
                                    as.character(electives_only$rttstart_YYYY),
                                    sep = "-")
  icds <- unique(electives_only$ICD)
  age_groupings <- unique(electives_only$agegrp_v3)
  week_year_combos <- unique(electives_only$week_year)
  num_weeks <- length(week_year_combos)
  num_icd <- length(icds)
  num_ages <- 3
  nrows_df <- num_weeks * num_ages * num_icd
  
  elective_out_Df <- data.frame(matrix(nrow = nrows_df, ncol = 8))
  colnames(elective_out_Df) <- c("rttstart_YYYY","rttstart_week",
                                  "Admissions","ICD","agegrp_v3",
                                  "prop_Frail","p50_WT_ICDc",
                                  "mean_WT_ICDc")
  weeks_year_split <- str_split_fixed(week_year_combos, "-",2)
  elective_out_Df$rttstart_YYYY <- rep(as.integer(weeks_year_split[,2]),(num_ages * num_icd))
  elective_out_Df$rttstart_week <- rep(as.integer(weeks_year_split[,1]),(num_ages * num_icd))
  elective_out_Df$ICD <- rep(icds, each = num_weeks * num_ages)
  elective_out_Df$agegrp_v3 <- rep(rep(age_groupings, each = num_weeks), num_icd)
  
  print("Creating parallel jobs")
  elective_rows <- seq(1, nrows_df)
  cluster_function <- snow::makeCluster(spec = num_cores)
  function_input <- snow::clusterSplit(cluster_function, elective_rows)
  print("Loading up elective functions")
  snow::clusterExport(cluster_function, "row_sumer_elective")
  snow::clusterExport(cluster_function, "elective_df_ts")
  print("Loading up elective data")
  snow::clusterExport(cluster_function, "electives_only", envir = environment())
  snow::clusterExport(cluster_function, "elective_out_Df", envir = environment())
  
  snow::clusterEvalQ(cluster_function, library(dplyr))
  snow::clusterEvalQ(cluster_function, library(stringr))
  print("Running parallel jobs")
  elective_data_parallel <- snow::clusterApply(cluster_function, function_input,
                                                fun = elective_df_ts,
                                                electives_only = electives_only,
                                               electives_only_Df = elective_out_Df)
  stopCluster(cluster_function)
  
  elective_data <- dplyr::bind_rows(elective_data_parallel)
  
  print("Finished")
  time_end <- Sys.time()
  print((time_end - time_start))
  elective_data$date <- as.POSIXct( paste( 1, elective_data$rttstart_week, elective_data$rttstart_YYYY, sep = "-" ), format = "%u-%U-%Y",locale = "UK" ) 
  
  
  elective_data <- remove_wt_outliers(elective_data)
  return(elective_data)
  
}


time_series_creator <- function(hes_data, num_cores){
  require(stringr)
  require(plyr)
  require(dplyr)
  require(snow)  
  emergency_ts <- running_emergencies_ts_in_parallel(hes_data, num_cores = num_cores)
  electives_ts <- running_elective_ts_in_parallel(hes_data = hes_data,
                                                  num_cores = num_cores)
  
  return(list(emergency_ts, electives_ts))
  
  
}





