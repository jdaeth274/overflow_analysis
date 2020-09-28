##########################################################
###### LOS frequency tables ##############################
##########################################################

##############################################################################################################

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
