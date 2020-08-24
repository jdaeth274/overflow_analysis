combine_elective_rows <- function(electives_ts_data){
  
  ## Function to combine rows from multiple years of elective time_series data 
  ## Input:
  ## electives_ts_data: Df from multiple years of output from time_series_creator binded together
  
  ## check if date variable present
  if(!("date" %in% colnames(electives_ts_data))){
    ## create quick proxy at the moment of pasted week & year doesn't have to be calendar date
    ## This will though have to be changed before input into the time_series_forecasts
    electives_ts_data$date <- paste(electives_ts_data$rttstart_YYYY, electives_ts_data$rttstart_week,sep = "-")
  }
  
  ## create index to remove duplicate rows with once created sum row 
  electives_ts_data$index <- seq(1, nrow(electives_ts_data))
  
  ## loop over icd and agegrps to get merging rows 
  unique_icds <- unique(electives_ts_data$ICD)
  uniques_ages <- c(1,2,3)
  total_groups <- 3 * length(unique_icds)
  current_group <- 1
  
  out_df <- NULL
  
  ## define summary function:
  admission_mean <- function(admi, mean_var){
    out <- (admi*mean_var)/sum(admi)
    return(out)
  }
  
  
  for(icd in unique_icds){
    for(age in uniques_ages){
      cat("On ICD:",icd, "age group:", age, "total group", current_group, "of", total_groups,"\n")
      current_df <- electives_ts_data[electives_ts_data$ICD == icd &
                                        electives_ts_data$agegrp_v3 == age,]
      ## Just make sure we skip over the 15 - 3 grouping which we remove from our data
      if(nrow(current_df) > 0){
        
        ## lets just take the mean of the medians as our median 
        
        tot_admi_df <- current_df %>% dplyr::group_by(date) %>% dplyr::summarise(admissions_tot = sum(Admissions),
                                                                                 prop_Frail = admission_mean(Admissions, prop_Frail),
                                                                                 prop_cc = admission_mean(Admissions, prop_cc),
                                                                                 mean_WT_ICDc = admission_mean(Admissions, mean_WT_ICDc) ,
                                                                                 p50_WT_ICDc = admission_mean(Admissions, p50_WT_ICDc),
                                                                                 rttstart_YYYY = rttstart_YYYY[1], rttstart_week = rttstart_week[1])
        tot_admi_df$ICD <- icd  
        tot_admi_df$agegrp_v3 <- age
        colnames(tot_admi_df)[2] <- "Admissions"
        
        out_df <- dplyr::bind_rows(out_df, tot_admi_df)
        current_group <- current_group + 1
      }
    }
  }
  
  return(out_df)
  
}
