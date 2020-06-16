###############################################################################
## Split the cohorts into 1,2,3 ###############################################
###############################################################################

# set working directory


# Import HES data

cohort_allocator <- function(hes_data, hes_ids){

  ## WATCH OUT FOR THE DATE FORMAT IN RTTSTART AND ADMIDATE_MDY
  
  #print("Starting cohort identification")
  
    
    
    for(k in 1:nrow(hes_ids)){
     #cat("\r","This percent done", (round((k/nrow(hes_ids)*100), digits = 6)), "%")
    current_id <- as.character(hes_ids[k, 1])
    
    narrowed_df <- hes_data[hes_data$hesid == current_id,]
    index_altered <- NULL
    rtt_start_date_2 <- NULL
    
    if(any(grepl(1,narrowed_df$cohort))){
      different_diags <- plyr::count(narrowed_df$diag_01)
      
      for(j in 1:nrow(different_diags)){
        current_diag <- different_diags[j, 1]
        diag_data <- narrowed_df[narrowed_df$diag_01 == current_diag,]
        if(any(grepl(1,diag_data$cohort)) & any(grepl(3, diag_data$cohort))){
          electives <- which(diag_data$cohort == 1)
          elective_rows <- diag_data[electives,]
          emergencies <- diag_data[-electives,]
          elective_rows <- elective_rows[order(elective_rows$rttstart),]
          
          
          for(elec in 1:nrow(elective_rows)){
            indexers_to_change <- NULL
            current_elective <- elective_rows[elec,]
            start_date <- current_elective$rttstart
            end_date <- current_elective$admidate_MDY
            
            cohort_2ers <- which((emergencies$admidate_MDY >= start_date) &
                                   (emergencies$admidate_MDY <= end_date))
            if(length(cohort_2ers) > 0)
              indexers_to_change <- emergencies$index[cohort_2ers]
            if(length(indexers_to_change) > 0){
              new_indies <- which(!(indexers_to_change %in% index_altered))
              indexers_to_change <- indexers_to_change[new_indies]
              index_altered <- base::append(index_altered, indexers_to_change)
              new_rttstart <- rep(current_elective$rttstart, length(indexers_to_change))
              rtt_start_date_2 <- append(rtt_start_date_2, new_rttstart)
            }
          }
          
          
          
        }
        
        
      }
      
      
    }
    
    if(length(index_altered) > 0){
      
      ## check for index duplicates (where an emergency overlaps two elective rttstarts and ends)
      if(length(index_altered) > 1){
        
        indy_date <- cbind.data.frame(index_altered, rtt_start_date_2)
        colnames(indy_date) <- c("index","rtt")
        
        indy_date <- indy_date[order(indy_date$rtt),]
        
        indy_date <- indy_date[!duplicated(indy_date[,"index"]),]
        
        
        index_altered <- indy_date$index
        rtt_start_date_2 <- indy_date$rtt
        
      }
      
      hes_data[hes_data$index %in% index_altered,"cohort"] <- 2
      
      for(rtt in 1:length(rtt_start_date_2)){
        current_index <- index_altered[rtt]
        current_rtt <- rtt_start_date_2[rtt]
        
        hes_data[hes_data$index == current_index,"rttstart"] <- current_rtt
      
      }
    }
    }
  
  return(hes_data)
  
  
}


cohort_allocator_parallel <- function(unique_rows, hes_dataset, hes_ids){
  #cat("\r","ON this num: ", unique_rows)
  
  input_hes_rows <- hes_ids[unique_rows,]
  input_hes_data <- hes_dataset[hes_dataset$hesid %in% as.character(input_hes_rows[,1]),]
  rm(hes_dataset)
  
  re_shifted_cohort <- cohort_allocator(hes_data = input_hes_data, hes_ids = input_hes_rows)
  
  
  return(re_shifted_cohort)
}


indiv_ICD_para <- function(ICD = "ONE", num_cores, hes_data){
  
  if(ICD != "ONE"){
    
    start_time <- Sys.time()
    current_ICD_data <- hes_data[hes_data$ICD == ICD,]
    current_ICD_data <- as.data.frame(current_ICD_data)
    rm(hes_data)
    
    hes_ids <- plyr::count(current_ICD_data$hesid)
    hesid_rows <- seq(1, nrow(hes_ids))
    
    
    print(paste("Setting up parallel job for ICD:", ICD))
    cluster_function <- snow::makeCluster(spec = num_cores)
    function_input <- snow::clusterSplit(cluster_function, hesid_rows)
    print(paste("Copying over functions for ICD:", ICD))
    snow::clusterExport(cluster_function, "cohort_allocator")
    print(paste("Copying over data for ICD:", ICD))
    copy_start <- Sys.time()
    snow::clusterExport(cluster_function, "hes_ids", envir = environment())
    snow::clusterExport(cluster_function, "current_ICD_data", envir = environment())
    snow::clusterExport(cluster_function, "cohort_allocator_parallel")
    copy_end <- Sys.time()
    print(copy_end - copy_start)
    print(paste("Running cohort allocation jobs for ICD:", ICD))
    jobs_start <- Sys.time()
    hes_data_parallel <- snow::clusterApply(cluster_function, function_input,
                                            fun = cohort_allocator_parallel,
                                            hes_ids = hes_ids,
                                            hes_dataset = current_ICD_data)
    snow::clusterEvalQ(cluster_function, rm(hes_ids, pos = environment()))
    snow::clusterEvalQ(cluster_function, rm(current_ICD_data, pos = environment()))
    snow::clusterEvalQ(cluster_function, gc())
    stopCluster(cluster_function)
    jobs_end <- Sys.time()
    print(jobs_end - jobs_start)
    hes_cohorts_df <- base::as.data.frame(dplyr::bind_rows(hes_data_parallel))
    
    hes_cohorts_df <- hes_cohorts_df[order(hes_cohorts_df$index),]
    current_ICD_data$cohort <- hes_cohorts_df$cohort
    current_ICD_data$rttstart <- hes_cohorts_df$rttstart
    hes_cohorts_df <- current_ICD_data
    
    
    print("Writing out results")
    
    
    print("Finished")
    end_time <- Sys.time()
    
    print((end_time - start_time))
    return(hes_cohorts_df)
  }else{
    start_time <- Sys.time()
    hes_data <- as.data.frame(hes_data)
    hes_ids <- plyr::count(hes_data$hesid)
    hesid_rows <- seq(1, nrow(hes_ids))
    
    
    print(paste("Setting up parallel job for ICD:", ICD))
    cluster_function <- snow::makeCluster(spec = num_cores, outfile = "D:/Overflows/cluster_log_file.txt")
    function_input <- snow::clusterSplit(cluster_function, hesid_rows)
    print(paste("Copying over functions for ICD:", ICD))
    snow::clusterExport(cluster_function, "cohort_allocator")
    print(paste("Copying over data for ICD:", ICD))
    copy_start <- Sys.time()
    snow::clusterExport(cluster_function, "hes_ids", envir = environment())
    snow::clusterExport(cluster_function, "hes_data", envir = environment())
    snow::clusterExport(cluster_function, "cohort_allocator_parallel")
    copy_end <- Sys.time()
    print(copy_end - copy_start)
    print(paste("Running cohort allocation jobs for ICD:", ICD))
    jobs_start <- Sys.time()
    hes_data_parallel <- snow::clusterApply(cluster_function, function_input,
                                            fun = cohort_allocator_parallel,
                                            hes_ids = hes_ids,
                                            hes_dataset = hes_data)
    snow::clusterEvalQ(cluster_function, rm(hes_ids))
    snow::clusterEvalQ(cluster_function, rm(hes_data))
    snow::clusterEvalQ(cluster_function, gc())
    stopCluster(cluster_function)
    
    gc()
    jobs_end <- Sys.time()
    print(jobs_end - jobs_start)
    hes_cohorts_df <- base::as.data.frame(dplyr::bind_rows(hes_data_parallel))
    
    
    hes_cohorts_df <- hes_cohorts_df[order(hes_cohorts_df$index),]
    hes_data$cohort <- hes_cohorts_df$cohort
    hes_data$rttstart <- hes_cohorts_df$rttstart
    
    
    
    print("Writing out results")
    
    
    print("Finished")
    end_time <- Sys.time()
    
    print((end_time - start_time))
    return(hes_data)
  }
  

}


making_WT_variable <- function(hes_data){
  
  hes_data$WaitingTime <- NA
  
  cohort_1 <- hes_data[hes_data$cohort == 1,]
  one_year_under <- which(as.integer(cohort_1$elecdur) <= 365)
  cohort_1$WaitingTime[one_year_under] <- cohort_1$elecdur[one_year_under]
  
  cohort_2 <- hes_data[hes_data$cohort == 2,]
  wait_times <- cohort_2$admidate_MDY - cohort_2$rttstart
  
}



cohort_set_up <- function(num_cores = 0, data_loc = "E:/HES/COVID/HES_APC_CC_0913_TEMP02.csv"){
    
  start_time <- Sys.time()
  require(parallel)
  require(plyr)
  require(snow)
  require(vroom)
  require(pryr)
  print("Loading up data")
  data_start <- Sys.time()
  data <- vroom(data_loc, delim = ",", num_threads = num_cores)
  
  data_end <- Sys.time()
  print(data_end - data_start)

  data$index <- seq(1, nrow(data))
  
  ###############################################################################
  ## Ok so lets first set up our cohort variable on the dataframe:             ##
  ## 1 - Electives, 2 - Electives to Emergencies, 3 - Emergencies              ##
  ###############################################################################
  
  data$cohort <- 3
  
  electives <- which(data$admimeth_C == 1)
  
  data$cohort[electives] <- 1
  
  ###############################################################################
  ## Lets loop through the unique hesids to get to the electives turned        ##
  ## emergencies ################################################################
  ###############################################################################
  
  if(num_cores == 0){
  
    num_cores <- parallel::detectCores()
    
    use_cores <- num_cores - 2
  }
  
  
  
  print("Narrowing data")
  #old_size <- pryr::object_size(data)
  
  if(any(grepl("ICD$", colnames(data), ignore.case = TRUE))){
    
  
  
  parallel_cols <- which(colnames(data) %in% c("hesid","index","diag_01","cohort","rttstart",
                                               "admidate_MDY"))
  icd_col <- base::grep("ICD$", colnames(data), ignore.case = TRUE)
  parallel_cols <- append(parallel_cols, icd_col)
  
  
  parallel_data <- data[,parallel_cols]
  
  
  paralell_admi_nas <- which(is.na(parallel_data$admidate_MDY))
  if(length(paralell_admi_nas) > 0){
    parallel_data <- parallel_data[-paralell_admi_nas,]
  }
  elective_nas <- parallel_data[parallel_data$cohort == 1,]
  elective_na_rows <- which(is.na(elective_nas$rttstart))
  remove_rtt_missing <- elective_nas$index[elective_na_rows]
  if(length(remove_rtt_missing)>0)
    parallel_data <- parallel_data[-remove_rtt_missing,]
  
  
  
  
  parallel_data$rttstart <- as.Date(parallel_data$rttstart, format = "%d%b%Y")
  parallel_data$admidate_MDY <- as.Date(parallel_data$admidate_MDY, format = "%d%b%Y")    
  
  
  #new_size <- pryr::object_size(parallel_data)
  #print(paste("Old size:", old_size, "New size:", new_size))
  
  ## ICD groups ###
  

  
  ICD_groupings <- plyr::count(parallel_data$ICD)
  
  icds <- as.character(ICD_groupings[,1])
  
  
  print("Lapplying through the ICDS")
  lapply_start <- Sys.time()
  
  icds_res <- base::lapply(icds, FUN = indiv_ICD_para,
                           num_cores = num_cores,
                           hes_data = parallel_data)
  lapply_end <- Sys.time()
  print(base::paste("Lapply total time spent", lapply_end - lapply_start))
  
  total_res <- dplyr::bind_rows(icds_res)
  total_res <- total_res[order(total_res$index),]
  missing_indies <- which(!(data$index %in% total_res$index))
  data <- data[-missing_indies,]
  data$cohort_old <- data$cohort
  data$cohort <- total_res$cohort
  data$rttstart <- total_res$rttstart
  data$admidate_MDY <- as.Date(data$admidate_MDY, format = "%d%b%Y")    
  
  }else{
    print("I.m in the one icd section !")
    ## assume you're just running on one icd here 
    parallel_cols <- which(colnames(data) %in% c("hesid","index","diag_01","cohort","rttstart",
                                                 "admidate_MDY"))
    
    parallel_data <- data[,parallel_cols]
   # new_size <- pryr::object_size(parallel_data)
    #print(paste("Old size:", old_size, "New size:", new_size))
    
    paralell_admi_nas <- which(is.na(parallel_data$admidate_MDY))
    if(length(paralell_admi_nas) > 0){
      parallel_data <- parallel_data[-paralell_admi_nas,]
    }
    elective_nas <- parallel_data[parallel_data$cohort == 1,]
    elective_na_rows <- which(is.na(elective_nas$rttstart))
    remove_rtt_missing <- elective_nas$index[elective_na_rows]
    if(length(remove_rtt_missing)>0)
      parallel_data <- parallel_data[-remove_rtt_missing,]
    
    
                              
    
    parallel_data$rttstart <- as.Date(parallel_data$rttstart, format = "%d%b%Y")
    parallel_data$admidate_MDY <- as.Date(parallel_data$admidate_MDY, format = "%d%b%Y")    
    
    total_res <- indiv_ICD_para(num_cores = num_cores,
                              hes_data = parallel_data)
    total_res <- total_res[order(total_res$index),]
    missing_indies <- which(!(data$index %in% total_res$index))
    data <- data[-missing_indies,]
    data$cohort_old <- data$cohort
    
    data$cohort <- total_res$cohort
    data$rttstart <- total_res$rttstart
    data$admidate_MDY <- as.Date(data$admidate_MDY, format = "%d%b%Y")    
    
    
  }
  
  tot_time <- Sys.time()
  print(tot_time - start_time)
  
  return(data)
  
}














