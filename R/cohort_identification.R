###############################################################################
## Split the cohorts into 1,2,3 ###############################################
###############################################################################

# set working directory


# Import HES data

cohort_allocator <- function(hes_data, hes_ids){

  
  
  #print("Starting cohort identification")
  
  
    
    for(k in 1:nrow(hes_ids)){
     # cat("\r","This percent done", (round((k/nrow(hes_ids)*100), digits = 6)), "%")
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
          
          for(elec in 1:nrow(elective_rows)){
            current_elective <- elective_rows[elec,]
            start_date <- current_elective$rttstart
            end_date <- current_elective$admidate_MDY
            
            cohort_2ers <- which((emergencies$admidate_MDY >= start_date) &
                                   (emergencies$admidate_MDY <= end_date))
            indexers_to_change <- emergencies$index[cohort_2ers]
            if(length(indexers_to_change) > 0){
              index_altered <- base::append(index_altered, indexers_to_change)
              new_rttstart <- rep(current_elective$rttstart, length(indexers_to_change))
              rtt_start_date_2 <- append(rtt_start_date_2, new_rttstart)
            }
          }
          
          
          
        }
        
        
      }
      
      
    }
    
    if(length(index_altered) > 0){
      hes_data[hes_data$index %in% index_altered,"cohort"] <- 2
      hes_data[hes_data$index %in% index_altered,"rttstart"] <- rtt_start_date_2
    }
    }
  
  return(hes_data)
  
  
}


cohort_allocator_parallel <- function(unique_rows, hes_dataset, hes_ids){
  cat("\r","ON this num: ", unique_rows)
  input_hes_rows <- hes_ids[unique_rows,]
  input_hes_data <- hes_dataset[hes_dataset$hesid %in% as.character(input_hes_rows[,1]),]
  
  re_shifted_cohort <- cohort_allocator(hes_data = input_hes_data, hes_ids = input_hes_rows)
  
  
  return(re_shifted_cohort)
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
  print("Counting unique hesids")
  count_start <- Sys.time()
  hes_ids <- plyr::count(data$hesid)
  count_end <- Sys.time()
  print(count_end - count_start)
  
  if(num_cores == 0){
  
    num_cores <- parallel::detectCores()
    
    use_cores <- num_cores - 2
  }
  
  hesid_rows <- seq(1, nrow(hes_ids))
  
  print("Narrowing data")
  old_size <- pryr::object_size(data)
  
  parallel_cols <- which(colnames(data) %in% c("hesid","index","diag_01","cohort","rttstart",
                                               "admidate_MDY"))
  
  
  parallel_data <- data[,parallel_cols]
  new_size <- pryr::object_size(parallel_data)
  print(paste("Old size:", old_size, "New size:", new_size))
  
  
  print("Setting up parralel job")
  cluster_function <- snow::makeCluster(spec = num_cores)
  function_input <- snow::clusterSplit(cluster_function, hesid_rows)
  print("COpying over functions")
  snow::clusterExport(cluster_function, "cohort_allocator")
  print("Copying over data")
  copy_start <- Sys.time()
  snow::clusterExport(cluster_function, "hes_ids", envir = environment())
  snow::clusterExport(cluster_function, "parallel_data", envir = environment())
  snow::clusterExport(cluster_function, "cohort_allocator_parallel")
  copy_end <- Sys.time()
  print(copy_end - copy_start)
  print("Running cohort allocation jobs")
  jobs_start <- Sys.time()
  hes_data_parallel <- snow::clusterApply(cluster_function, function_input,
                                          fun = cohort_allocator_parallel,
                                          hes_ids = hes_ids,
                                          hes_dataset = data)
  stopCluster(cluster_function)
  jobs_end <- Sys.time()
  print(jobs_end - jobs_start)
  hes_cohorts_df <- dplyr::bind_rows(hes_data_parallel)
  
  hes_cohorts_df <- hes_cohorts_df[order(hes_cohorts_df$index),]
  data$cohort <- hes_cohorts_df$cohort
  hes_cohorts_df <- data
  
  
  print("Writing out results")
  
  write.csv(hes_cohorts_df,
            file = "E:/temp/hes_cohort_allocations_BIG.csv",
            row.names = FALSE)
  
  print("Finished")
  end_time <- Sys.time()
  
  print((end_time - start_time))
  return(hes_cohorts_df)
  
  
}














