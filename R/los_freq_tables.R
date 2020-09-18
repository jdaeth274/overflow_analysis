##########################################################
###### LOS frequency tables ##############################
##########################################################

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


