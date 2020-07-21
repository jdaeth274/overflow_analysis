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
  
  if(nrow(subset_hes) > 0){
    return_row$prop_Frail <- sum(subset_hes$Frail) / nrow(subset_hes)
    return_row$prop_cc <- sum(subset_hes$cc) / nrow(subset_hes)
    
  }else{
    return_row$prop_Frail <- NA
    return_row$prop_cc <- NA
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
                                                         "Frail","ICD","agegrp_v3","cc"))
  emergencies_only <- emergencies_only[,emergencies_cols]
  
  
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
  
  emergency_out_Df <- data.frame(matrix(nrow = nrows_df, ncol = 7))
  colnames(emergency_out_Df) <- c("admidate_YYYY","admidate_week",
                                  "Admissions","ICD","agegrp_v3",
                                  "prop_Frail","prop_cc")
  weeks_year_split <- stringr::str_split_fixed(week_year_combos, "-",2)
  emergency_out_Df$admidate_YYYY <- rep(as.integer(weeks_year_split[,2]),(num_ages * num_icd))
  emergency_out_Df$admidate_week <- rep(as.integer(weeks_year_split[,1]),(num_ages * num_icd))
  emergency_out_Df$ICD <- rep(icds, each = num_weeks * num_ages)
  emergency_out_Df$agegrp_v3 <- rep(rep(age_groupings, each = num_weeks), num_icd)
  emergency_out_Df <- emergency_out_Df
  
  
  emergency_rows <- seq(1, nrows_df)
  tic("CLuster_Set_up")
  print("Setting up the parallel jobs")
  cluster_function <- snow::makeCluster(spec = num_cores, outfile = "./ts_creat_log.txt")
  function_input <- snow::clusterSplit(cluster_function, emergency_rows)
  toc()
  tic("Function copy")
  snow::clusterExport(cluster_function, "row_sumer")
  snow::clusterExport(cluster_function, "emergency_df_ts")
  toc()
  print("Loading data onto parallel cores")
  tic("Loading data onto parallel cores")
  snow::clusterExport(cluster_function, "emergencies_only", envir = environment())
  snow::clusterExport(cluster_function, "emergency_out_Df", envir = environment())
  snow::clusterEvalQ(cluster_function, library(dplyr))
  snow::clusterEvalQ(cluster_function, library(stringr))
  toc()
  tic("Running on cores")
  print("Running on cores")
  emergency_data_parallel <- snow::clusterApply(cluster_function, function_input,
                                                fun = emergency_df_ts,
                                                emergencies_only = emergencies_only,
                                                Emergency_out_df = emergency_out_Df)
  toc()
  snow::clusterEvalQ(cluster_function, rm(emergencies_only))
  snow::clusterEvalQ(cluster_function, rm(emergencies_only))
  snow::clusterEvalQ(cluster_function, gc())
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
  
  current_icd <- out_row[row_indy, 4]
  current_age <- out_row[row_indy, 5]
  current_week <- out_row[row_indy, 2]
  current_year <- out_row[row_indy, 1]
  
  subset_hes <- hes_df[ hes_df$ICD == current_icd & 
                          hes_df$agegrp_v3 == current_age &
                          hes_df$rttstart_YYYY == current_year &
                          hes_df$rttstart_week == current_week, ]
  return_row <- out_row[row_indy,]
  print(nrow(return_row))
  
  return_row$Admissions <- nrow(subset_hes)
  
  if(nrow(subset_hes) > 0){
    return_row$prop_Frail <- sum(subset_hes$Frail) / nrow(subset_hes)
    return_row$p50_WT_ICDc <- median(subset_hes$WT)
    return_row$mean_WT_ICDc <- mean(subset_hes$WT)
    return_row$prop_cc <- sum(subset_hes$cc) / nrow(subset_hes)
  }else{
    return_row$prop_Frail <- NA
    return_row$p50_WT_ICDc <- NA
    return_row$mean_WT_ICDc <- NA
    return_row$prop_cc <- NA
  }
  
  
  
  
  return(return_row)
  
}


elective_df_ts <- function(out_indies, electives_only_Df, electives_only){
  
  print(out_indies)
  
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


waiting_pool <- function(current_icd, hes_data, end_date){
  ## This function produces our waiting pool, designed to be run with lapply
  current_hes <- hes_data[hes_data$ICD == current_icd,]
  
  out_df <- data.frame(matrix(ncol = 5,nrow = 3))
  colnames(out_df) <- c("age","icd","pool","median_WT","mean_WT")
  out_df$age <- c("<25","25-64", "65+")
  out_df$icd <- current_icd
  
  
  for(age in 1:3){
  
    age_1 <- current_hes[current_hes$agegrp_v3 == age &
                           current_hes$rttstart < end_date &
                           current_hes$admidate_MDY >= end_date,]
    out_df$pool[age] <- nrow(age_1)
    if(nrow(age_1) > 0){
      out_df$median_WT[age] <- median(age_1$WT, na.rm = TRUE)
      out_df$mean_WT[age] <- mean(age_1$WT, na.rm = TRUE)
    }
  
  }

  
  return(out_df)
}


running_elective_ts_in_parallel <- function(hes_data, num_cores, forecast_date, ts_run, waiting_pool_run){

  time_start <- Sys.time()
  electives_only <- hes_data[hes_data$cohort == 1,]
  
  if(!("WT" %in% colnames(hes_data)))
    electives_only <- make_wt_variable(electives_only)
  
  missing_rtt_week <- which(is.na(electives_only$rttstart_week))
  missing_rtt_year <- which(is.na(electives_only$rttstart_YYYY))
  
  if(length(missing_rtt_week) >0)
    electives_only <- electives_only[-missing_rtt_week,]  
  if(length(missing_rtt_year) >0)
    electives_only <- electives_only[-missing_rtt_year,]  
  
  
  elective_cols <- which(colnames(electives_only) %in% c("rttstart_week", "rttstart_YYYY",
                                                         "WT","Frail","ICD","agegrp_v3","cc"))
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
  
  elective_out_Df <- data.frame(matrix(nrow = nrows_df, ncol = 9))
  colnames(elective_out_Df) <- c("rttstart_YYYY","rttstart_week",
                                  "Admissions","ICD","agegrp_v3",
                                  "prop_Frail","p50_WT_ICDc",
                                  "mean_WT_ICDc","prop_cc")
  weeks_year_split <- str_split_fixed(week_year_combos, "-",2)
  elective_out_Df$rttstart_YYYY <- rep(as.integer(weeks_year_split[,2]),(num_ages * num_icd))
  elective_out_Df$rttstart_week <- rep(as.integer(weeks_year_split[,1]),(num_ages * num_icd))
  elective_out_Df$ICD <- rep(icds, each = num_weeks * num_ages)
  elective_out_Df$agegrp_v3 <- rep(rep(age_groupings, each = num_weeks), num_icd)
  
  if(ts_run){
  
  
  tic("Creating parallel jobs")
  print("Creating parallel jobs")
  print(nrows_df)
  elective_rows <- seq(1, nrows_df)
  cluster_function <- snow::makeCluster(spec = num_cores,  outfile = "./cluster_log_file_TS.txt")
  function_input <- snow::clusterSplit(cluster_function, elective_rows)
  toc()
  tic("Loading up elective functions")
  print("Loading up elective functions")
  snow::clusterExport(cluster_function, "row_sumer_elective")
  snow::clusterExport(cluster_function, "elective_df_ts")
  toc()
  tic("Loading up elective data")
  print("Loading up elective data")
  snow::clusterExport(cluster_function, "electives_only", envir = environment())
  snow::clusterExport(cluster_function, "elective_out_Df", envir = environment())
  
  snow::clusterEvalQ(cluster_function, library(dplyr))
  snow::clusterEvalQ(cluster_function, library(stringr))
  toc()
  tic("Running parallel jobs")
  print("Running parallel jobs")
  elective_data_parallel <- snow::clusterApply(cluster_function, function_input,
                                                fun = elective_df_ts,
                                                electives_only = electives_only,
                                               electives_only_Df = elective_out_Df)
  
  snow::clusterEvalQ(cluster_function, rm(electives_only))
  snow::clusterEvalQ(cluster_function, rm(elective_out_Df))
  snow::clusterEvalQ(cluster_function, gc())
    stopCluster(cluster_function)
  elective_data <- dplyr::bind_rows(elective_data_parallel)
  
  print("Finished")
  time_end <- Sys.time()
  print((time_end - time_start))
  elective_data$date <- as.POSIXct( paste( 1, elective_data$rttstart_week, elective_data$rttstart_YYYY, sep = "-" ), format = "%u-%U-%Y",locale = "UK" ) 
  
  
  elective_data <- remove_wt_outliers(elective_data)
  
  }else{
    elective_data <- 0
  }
  ## get waiting pool ##
  
  
  
  if(waiting_pool_run){
  
    if(ts_run){
  
    if(as.Date(forecast_date) %in% elective_data$date){
      end_date <- forecast_date
    }else{
      dates_in_forecast <- as.Date(elective_data$date)
      end_date <- min(dates_in_forecast[dates_in_forecast > forecast_date])
    }
    }else{
      end_date <- as.Date("2012-03-01")
    }
  electives_only <- hes_data[hes_data$cohort == 1,]
  elective_cols <- which(colnames(electives_only) %in% c("rttstart", "admidate_MDY",
                                                         "WT","Frail","ICD","agegrp_v3"))
  electives_only <- electives_only[,elective_cols]
  
    
  waiting_pools <- lapply(icds, FUN = waiting_pool,
                          hes_data = electives_only,
                          end_date)
  waiting_pools <- dplyr::bind_rows(waiting_pools)
  }else{
    waiting_pools <- 0
  }
  
  
  
  
  return(list(elective_data, waiting_pools))
  
}


bundle_props <- function(hes_data){
  
  narrowed_hes <- hes_data[,which(colnames(hes_data) %in% c("admidate_week","admidate_YYYY", "ICD",
                                                            "agegrp_v3","ICD","MainICD10Cat",
                                                            "cohort", "rttstart_week","rttstart_YYYY"))]
  
  ## electives_first 
  elective_data <- narrowed_hes[narrowed_hes$cohort == 1,]
  
  elective_bundle_data <- elective_data[elective_data$ICD == 50,]
  rm(elective_data)
  rm(hes_data)
  elective_bundle_data$non_bundled_icd <- elective_bundle_data$MainICD10Cat
  elective_bundle_data$non_bundled_icd[which(elective_bundle_data$MainICD10Cat %in% c(8,16,17))] <- 50
  elective_bundle_data$one <- 1
  cats_icd <- unique(elective_bundle_data$non_bundled_icd) 
  
  
  
  weekly_ICD_vals_elec <- aggregate(elective_bundle_data, by = list(elective_bundle_data$non_bundled_icd,
                                                               elective_bundle_data$rttstart_week,
                                                               elective_bundle_data$rttstart_YYYY,
                                                               elective_bundle_data$agegrp_v3), FUN = sum)
  tot_vals <- aggregate(elective_bundle_data, by = list(elective_bundle_data$ICD,
                                                               elective_bundle_data$rttstart_week,
                                                               elective_bundle_data$rttstart_YYYY,
                                                               elective_bundle_data$agegrp_v3), FUN = sum)
  tot_vals <- tot_vals[order(tot_vals$Group.4,tot_vals$Group.3,tot_vals$Group.2),]
  weekly_ICD_vals_elec <- weekly_ICD_vals_elec[order(weekly_ICD_vals_elec$Group.1,
                                                     weekly_ICD_vals_elec$Group.4,weekly_ICD_vals_elec$Group.3,
                                                     weekly_ICD_vals_elec$Group.2),]
  
  ## make out df 
  
  proportions_elec_out <- data.frame(matrix(nrow = (length(cats_icd) * nrow(tot_vals)),
                                            ncol = 4))
  colnames(proportions_elec_out) <- c("rttstart_week","rttstart_YYYY","agegrp_v3",
                                      "ICD")
  proportions_elec_out$rttstart_week <- rep(tot_vals$Group.2, length(cats_icd))
  proportions_elec_out$rttstart_YYYY <- rep(tot_vals$Group.3, length(cats_icd))
  proportions_elec_out$agegrp_v3 <- rep(tot_vals$Group.4, length(cats_icd))
  proportions_elec_out$ICD <- rep(cats_icd, each = nrow(tot_vals))
  
  proportions_elec_out <- dplyr::left_join(proportions_elec_out, weekly_ICD_vals_elec,
                                           by = c("rttstart_week" = "Group.2",
                                                  "rttstart_YYYY" = "Group.3",
                                                  "agegrp_v3" = "Group.4",
                                                  "ICD" = "Group.1"))
  
  proportions_elec_out <- proportions_elec_out[,c(1:4,ncol(proportions_elec_out))]
  proportions_elec_out$tot_bundle <- rep(tot_vals$one, length(cats_icd))
  proportions_elec_out$prop_bundle <- proportions_elec_out$one / proportions_elec_out$tot_bundle
  
  rm(elective_bundle_data)
  rm(weekly_ICD_vals_elec)
  rm(tot_vals)
  
  ## Now emergencies 
  emergency_data <- narrowed_hes[narrowed_hes$cohort == 3,]
  
  emergency_bundle_data <- emergency_data[emergency_data$ICD == 51,]
  rm(emergency_data)
  emergency_bundle_data$non_bundled_icd <- emergency_bundle_data$MainICD10Cat
  emergency_bundle_data$non_bundled_icd[which(emergency_bundle_data$MainICD10Cat %in% c(8,16,17))] <- 51
  emergency_bundle_data$one <- 1
  cats_icd <- unique(emergency_bundle_data$non_bundled_icd) 
  
  
  
  weekly_ICD_vals_emerg <- aggregate(emergency_bundle_data, by = list(emergency_bundle_data$non_bundled_icd,
                                                               emergency_bundle_data$admidate_week,
                                                               emergency_bundle_data$admidate_YYYY,
                                                               emergency_bundle_data$agegrp_v3), FUN = sum)
  tot_vals <- aggregate(emergency_bundle_data, by = list(emergency_bundle_data$ICD,
                                                               emergency_bundle_data$admidate_week,
                                                               emergency_bundle_data$admidate_YYYY,
                                                               emergency_bundle_data$agegrp_v3), FUN = sum)
  tot_vals <- tot_vals[order(tot_vals$Group.4,tot_vals$Group.3,tot_vals$Group.2),]
  weekly_ICD_vals_emerg <- weekly_ICD_vals_emerg[order(weekly_ICD_vals_emerg$Group.1,
                                                     weekly_ICD_vals_emerg$Group.4,weekly_ICD_vals_emerg$Group.3,
                                                     weekly_ICD_vals_emerg$Group.2),]
  
  ## make out df 
  
  proportions_emerg_out <- data.frame(matrix(nrow = (length(cats_icd) * nrow(tot_vals)),
                                            ncol = 4))
  colnames(proportions_emerg_out) <- c("admidate_week","admidate_YYYY","agegrp_v3",
                                      "ICD")
  proportions_emerg_out$admidate_week <- rep(tot_vals$Group.2, length(cats_icd))
  proportions_emerg_out$admidate_YYYY <- rep(tot_vals$Group.3, length(cats_icd))
  proportions_emerg_out$agegrp_v3 <- rep(tot_vals$Group.4, length(cats_icd))
  proportions_emerg_out$ICD <- rep(cats_icd, each = nrow(tot_vals))
  
  proportions_emerg_out <- dplyr::left_join(proportions_emerg_out, weekly_ICD_vals_emerg,
                                           by = c("admidate_week" = "Group.2",
                                                  "admidate_YYYY" = "Group.3",
                                                  "agegrp_v3" = "Group.4",
                                                  "ICD" = "Group.1"))
  
  proportions_emerg_out <- proportions_emerg_out[,c(1:4,ncol(proportions_emerg_out))]
  proportions_emerg_out$tot_bundle <- rep(tot_vals$one, length(cats_icd))
  proportions_emerg_out$prop_bundle <- proportions_emerg_out$one / proportions_emerg_out$tot_bundle
  
  
  
  return(list(proportions_elec_out, proportions_emerg_out))
  
}

in_hosp_pool <- function(hes_data, forecast_date){
  narrowed_hes <- hes_data[,which(colnames(hes_data) %in% c("admidate_week","admidate_YYYY","agegrp_v3","ICD",
                                                            "cohort", "admidate_MDY", "disdate_MDY", "cc"))]
  
  narrowed_hes$disdate_MDY <- as.Date(narrowed_hes$disdate_MDY, format = "%d%b%Y")
  ## electives_first 
  elective_data <- narrowed_hes[narrowed_hes$cohort == 1,]
  
  elective_in_hosp_dat <- elective_data[elective_data$admidate_MDY < forecast_date &
                                          elective_data$disdate_MDY >= forecast_date,]
  rm(hes_data)
  
  elective_in_hosp_dat$one <- 1
  num_icds <- unique(elective_data$ICD)
  num_ages <- unique(elective_data$agegrp_v3)
  
  in_hosp_by_ICD <- aggregate(one ~ ICD + agegrp_v3 + cc, elective_in_hosp_dat, sum)
  
  ## make out df 
  
  in_hops_elec_out <- data.frame(matrix(nrow = (length(num_icds) * length(num_ages) * 2),
                                            ncol = 6))
  colnames(in_hops_elec_out) <- c("a","p","s",
                                      "icd","agegrpv3", "cc")
  
  ## get icd in right form 
  in_hops_elec_out$icd <- rep(rep(num_icds, each = length(num_ages)), 2)
  in_hops_elec_out$agegrpv3 <- rep(rep(num_ages, length(num_icds)), 2)
  in_hops_elec_out$cc <- rep(c(0,1), each = (length(num_icds) * length(num_ages)))
  in_hops_elec_out$a <- "N"
  
  for(k in 1:length(num_icds)){
    if(nchar(num_icds[k]) == 1){
      num_icds[k] <- paste("0",num_icds[k], sep = "")
    }
  }
  
  icd_vec <- paste("ICD",num_icds, sep = "")
  age_vec <- paste("_AGE",num_ages, sep = "")
  in_hops_elec_out$p <- paste(rep(rep(icd_vec, each = length(num_ages)), 2),rep(rep(age_vec, length(num_icds)), 2), sep = "")
  in_hops_elec_out$s <- rep(c("G","C"),each = (length(num_ages) * length(num_icds)))
  in_hops_elec_out$y0 <- 0
  
  in_hops_elec_out <- dplyr::left_join(in_hops_elec_out, in_hosp_by_ICD,
                                           by = c("icd" = "ICD",
                                                  "agegrpv3" = "agegrp_v3",
                                                  "cc" = "cc"))
  
  in_hops_elec_out <- in_hops_elec_out[,c(1:3,which(colnames(in_hops_elec_out) == "one"))]
  
  rm(in_hosp_by_ICD)
  rm(icd_vec)
  rm(age_vec)
  
  ## Emergencies now 
  emergency_data <- narrowed_hes[narrowed_hes$cohort == 3,]
  
  emergency_in_hosp_dat <- emergency_data[emergency_data$admidate_MDY < forecast_date &
                                          emergency_data$disdate_MDY >= forecast_date,]
  emergency_in_hosp_dat$one <- 1
  num_icds <- unique(emergency_data$ICD)
  num_ages <- unique(emergency_data$agegrp_v3)
  
  in_hosp_by_ICD <- aggregate(one ~ ICD + agegrp_v3 + cc, emergency_in_hosp_dat, sum)
  
  ## make out df 
  
  in_hops_emerg_out <- data.frame(matrix(nrow = (length(num_icds) * length(num_ages) * 2),
                                            ncol = 6))
  colnames(in_hops_emerg_out) <- c("a","p","s",
                                      "icd","agegrpv3", "cc")
  
  ## get icd in right form 
  in_hops_emerg_out$icd <- rep(rep(num_icds, each = length(num_ages)), 2)
  in_hops_emerg_out$agegrpv3 <- rep(rep(num_ages, length(num_icds)), 2)
  in_hops_emerg_out$cc <- rep(c(0,1), each = (length(num_icds) * length(num_ages)))
  in_hops_emerg_out$a <- "E"
  
  for(k in 1:length(num_icds)){
    if(nchar(num_icds[k]) == 1){
      num_icds[k] <- paste("0",num_icds[k], sep = "")
    }
  }
  
  icd_vec <- paste("ICD",num_icds, sep = "")
  age_vec <- paste("_AGE",num_ages, sep = "")
  in_hops_emerg_out$p <- paste(rep(rep(icd_vec, each = length(num_ages)), 2),rep(rep(age_vec, length(num_icds)), 2), sep = "")
  in_hops_emerg_out$s <- rep(c("G","C"),each = (length(num_ages) * length(num_icds)))
  in_hops_emerg_out$y0 <- 0
  
  in_hops_emerg_out <- dplyr::left_join(in_hops_emerg_out, in_hosp_by_ICD,
                                           by = c("icd" = "ICD",
                                                  "agegrpv3" = "agegrp_v3",
                                                  "cc" = "cc"))
  
  in_hops_emerg_out <- in_hops_emerg_out[,c(1:3,which(colnames(in_hops_emerg_out) == "one"))]
  
  in_hops_tot <- rbind.data.frame(in_hops_elec_out, in_hops_emerg_out)
  
  return(in_hops_tot)
  
}




time_series_creator <- function(hes_data, num_cores, forecast_date, emergency_run = TRUE,
                                elective_res = TRUE, elective_ts = TRUE, waiting_pool = TRUE,
                                forecast_cutoff){
  require(stringr)
  require(plyr)
  require(dplyr)
  require(snow)
  require(tictoc)
  
  if(emergency_run){  
  emergency_ts <- running_emergencies_ts_in_parallel(hes_data, num_cores = num_cores)
  }else{
    emergency_ts <- 0
  }
  if(elective_res){
    electives_res <- running_elective_ts_in_parallel(hes_data = hes_data,
                                                  num_cores = num_cores,
                                                  forecast_date = forecast_date,
                                                  elective_ts, waiting_pool)
    electives_ts <- electives_res[[1]]
    waiting_patient_pool <- electives_res[[2]]
  }else{
    elective_ts <- 0
    waiting_patient_pool <- 0
  }
    
  print("Calculating bundle proportions")

  prop_bundle <- bundle_props(hes_data)
  print("Calculating in hopsital pool")
  in_hosp_pool_tot <- in_hosp_pool(hes_data, forecast_cutoff)
  elec_bundle <- prop_bundle[[1]]
  emerg_bundle <- prop_bundle[[2]]
  
  
  return(list(emergency_ts, electives_ts, waiting_patient_pool,
              elec_bundle, emerg_bundle, in_hosp_pool_tot))
  
  
}





