###############################################################################
## regression analyses ########################################################
## elec 2 emerg logit                     #####################################
## emerg transitions multinomial logit    #####################################
## elec transitions multinomial logit     #####################################
###############################################################################
library(foreign)
library(nnet)
# next two are for constructing robust standard errors but do not work with nnet unless we worked them out.
library(sandwich) 
library(lmtest)

elec_to_emergencies <- function(patient_groups, hes_data_orig, start_year, end_year, month_trend = TRUE,
                                time_trend = TRUE){
  out_df <- NULL
  
  for(j in 1:length(patient_groups)){

  print(paste("On patient group:",patient_groups[j],". Number",j, "of",length(patient_groups)))
  tic(paste("Running for patient group:",patient_groups[j]))
  
  split_patient_group <- str_split_fixed(patient_groups[j], "-",2)
  current_icd <- as.integer(split_patient_group[1])
  current_age <- split_patient_group[2]
  
  print("Narrowing hes_data")
  tic("Narrow")
  hes_data <- hes_data_orig[hes_data_orig$elective_icd == current_icd & 
                              hes_data_orig$agegrp_v3 == current_age,]
  toc()
    
    
  year_range <- seq(start_year, end_year)
  year_add <- seq(0,by = 52, length.out = length(year_range))
  
  hes_data$year_add <-  unname(setNames(year_add,year_range)[as.character(hes_data$admidate_YYYY)])
  
  hes_data$reg_week <- hes_data$admidate_week + hes_data$year_add
  
  col_names <- NULL
  
  if(time_trend & month_trend){
    reg_dataset <- hes_data[,which(colnames(hes_data) %in% c(col_names, "Elective2Emergency","WaitingTime","reg_week"))]
  }else if(time_trend & month_trend == FALSE){
    reg_dataset <- hes_data[,which(colnames(hes_data) %in% c("Elective2Emergency","WaitingTime","reg_week"))]
  }else if(month_trend & time_trend == FALSE){
    reg_dataset <- hes_data[,which(colnames(hes_data) %in% c(col_names,"Elective2Emergency","WaitingTime"))]
  }else{
    reg_dataset <- hes_data[,which(colnames(hes_data) %in% c("Elective2Emergency","WaitingTime"))]
  }
  
  print("Checking for switches")
  switches <- plyr::count(reg_dataset$Elective2Emergency)
  
  if(nrow(switches) == 2){
    print("Switches present running glm")
    tic("glm run")
    logit_icd2 <- glm(Elective2Emergency ~ . , family = binomial(link = "logit"), data = reg_dataset)
    
    pred_df <- data.frame(matrix(data = 0, ncol = 1, nrow = 10))
    colnames(pred_df) <- "WaitingTime"
    
    pred_df$WaitingTime <- seq(7,70, by =7)
    
    pred_res <- predict(logit_icd2, newdata = pred_df, "response")
    out_row <- data.frame(matrix(data = 0, ncol = 4, nrow = 1))
    colnames(out_row) <- c("mean_7","age","ICD","patient_group")
    out_row$mean_7 <- mean(pred_res)
    out_row$age <- current_age
    out_row$ICD <- current_icd
    out_row$patient_group <- patient_groups[j]
    toc()
  
  }else{
    
    out_row <- data.frame(matrix(data = 0, ncol = 4, nrow = 1))
    colnames(out_row) <- c("mean_7","age","ICD","patient_group")
    out_row$mean_7 <- 0
    out_row$age <- current_age
    out_row$ICD <- current_icd
    out_row$patient_group <- patient_groups[j]
    
  }
  
  out_df <- dplyr::bind_rows(out_df, out_row)
  
  toc()
  }
  
  return(out_df)
}

make_elective_cohort_variable <- function(cohorts_data){
  ## Making variable to run survival analysis on elective ICD groupings 
  cohorts_data$elective_icd <- cohorts_data$MainICD10Cat
  ## 50 for the elective grouping 
  bundled_group <- c(16,5,8,17,1,15,4)
  rows_bundles <- which(cohorts_data$elective_icd %in% bundled_group)
  cohorts_data$elective_icd[rows_bundles] <- 50
  
  return(cohorts_data)
  
  
}


failure_func_setup <- function(hes_data, start_year = 2009, end_year = 2013, month_trend,
                               time_trend){
  print("Narrowing Failure func data")
  tic("Narrowing dat")
  hes_data <- hes_data[hes_data$cohort != 3,]
  hes_data <- make_elective_cohort_variable(hes_data)
  hes_data <- hes_data[,c("ICD","agegrp_v3","WaitingTime","Elective2Emergency",
                          "admidate_MDY","admidate_MM","admidate_YYYY",
                          "admidate_week","elective_icd")]
  
  toc()
  icd_list <- unique(hes_data$elective_icd)
  icds_to_run <- rep(icd_list, each = 3)
  ages_to_run <- rep(c(1,2,3),length(icd_list))
  patient_groups_to_run <- paste(icds_to_run,ages_to_run, sep = "-")
  tot_patient_groups <- c(patient_groups_to_run, patient_groups_to_run_cc)
  
  failure_func <- elec_to_emergencies(patient_groups_to_run, hes_data_orig = hes_data,
                                      start_year = start_year, end_year = end_year,
                                      month_trend = month_trend, time_trend = time_trend)
  
  
  
}


week_num_func <- function(current_week,start_week){
  
  weeknum <- difftime(strptime(current_week, format = "%Y-%m-%d"),
                      strptime(start_week, format = "%Y-%m-%d"), units = "weeks")
  weeknum <- floor(weeknum)
  
  return(weeknum)
  
}




elective_regression <- function(patient_group,hes_data_orig, start_date, forecast_length, forecast_start,
                                results_pdf = "D:/Dropbox/COVID19/Overflow/regressions_ga_output.pdf"){
  ## Input ICD hes data for one ICD one age group one patient group (cc/ga)
  ## to be run in parrallel
  tic("Whole process")
  whole_df <- NULL
  forecast_seq <- seq(as.Date(forecast_start), by = "week", length.out = 52)
  require(stringr)
  require(lubridate)  
  require(reshape2)
  require(tictoc)
  
  pdf(file = results_pdf, paper = "A4r", width = 10, height = 7)  
  
  for(j in 1:length(patient_group)){
    print(paste("On patient group:",patient_group[j],". Number",j, "of",length(patient_group)))
    tic(paste("Running for patient group:",patient_group[j]))
    
    split_patient_group <- str_split_fixed(patient_group[j], "-",3)
    current_icd <- as.integer(split_patient_group[1])
    current_ward <- split_patient_group[2]
    current_age <- as.integer(split_patient_group[3])
    
    hes_data <- hes_data_orig[hes_data_orig$ICD == current_icd & 
                                hes_data_orig$agegrp_v3 == current_age,]
    
    if(current_ward == "cc"){
      
      tic("Narrowing to CC only")
      hes_data <- hes_data[hes_data$cc == 1,]
      toc()
      ## remove na trans
      tic("Removing NAs")
      na_rows <- which(is.na(hes_data$cc_transitions))
      if(length(na_rows) > 0)
        hes_data <- hes_data[-na_rows,]
      
      na_ga_los <- which(is.na(hes_data$cc_LoS))
      if(length(na_ga_los) > 0)
        hes_data <- hes_data[-na_ga_los,]
      toc()
      
      ## set up week and year fixed effects
      
      tic("Getting the reg week data")
      
      hes_data$reg_week <- sapply(hes_data$admidate_MDY,week_num_func,start_week = start_date)
      toc()
      col_names <- NULL
      
  
      
      for(k in 2:length(month.name)){
        col_name <- paste(month.name[k], "fixed_effect",sep = "_")
        hes_data[,col_name] <- 0
        hes_data[hes_data$admidate_MM == k,col_name] <- 1
        col_names <- append(col_names, col_name)
        
      }
      
      
      ## Set up the outcome variable 
      tic("Set up outcome variable")
      hes_data$outcome <- NA
      
      hes_data[hes_data$cc_LoS >= 7 ,"outcome"]<-"CC"
      hes_data[hes_data$cc_LoS < 7 & hes_data$cc_transitions == 3,"outcome"] <- "Dead"
      hes_data[hes_data$cc_LoS < 7 & hes_data$cc_transitions == 2,"outcome"] <- "GA"
      hes_data[hes_data$cc_LoS < 7 & hes_data$cc_transitions == 1,"outcome"] <- "Discharged"
      
      
      
      hes_data$stay <- hes_data$cc_LoS
      hes_data[hes_data$stay > 7,"stay"]<-7
      
      
      
      hes_data$outcome <- factor(hes_data$outcome, levels = c("CC","Dead","GA","Discharged"))
      hes_data$outcome <- relevel(hes_data$outcome, ref = "GA")
      
      
      toc()
      ## set up the reg data 
      tic("Reg set up")
      
      reg_data_no_stay <- hes_data[,which(colnames(hes_data) %in% c(col_names, "outcome","WaitingTime","reg_week"))]
  
      actual_prop_dat <- reg_data_no_stay
      actual_prop_dat$one <- 1
      actual_prop_dat_agg <- aggregate(one ~ reg_week, data = actual_prop_dat, FUN = sum)
      actual_prop_dat_outcome <- aggregate(one ~ reg_week + outcome, data = actual_prop_dat, FUN = sum)
      
      toc()
      
      
      outcomes_in_dat <- plyr::count(reg_data_no_stay$outcome)
      ## Checking if enough outcomes for model else will return one for probs
      
      if(nrow(outcomes_in_dat) < 2){
        out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
        colnames(out_df) <- levels(hes_data$outcome)
        out_df[,outcomes_in_dat[1,1]] <- 1
        mean_wt_pred <- as.data.frame(out_df)
        mean_wt_pred$patient_group <- patient_group[j]
        mean_wt_pred$ICD <- current_icd
        mean_wt_pred$age <- current_age
        mean_wt_pred$WT <- "mean"
        median_wt_pred <- mean_wt_pred
        
        
      }else{
        
        
        
        # ml depending on WT only
        tic("Regreesion model fitting")
        ml_stay <- multinom(outcome ~ ., data = reg_data_no_stay)
        toc()
        
        # we can average the probability over WT
        
        
        
        # Or calculate it at the mean WT
        tic("Probability creation")
        means <- data.frame(matrix(data = 0,ncol = ncol(reg_data_no_stay) - 1, nrow = forecast_length))
        forecast_week_start <- week_num_func(forecast_start,"2009-01-01")
        forecast_week_end <- week_num_func(forecast_seq[forecast_length],"2009-01-01")
        forecast_week_nums <- seq(forecast_week_start, forecast_week_end)
        
        ## make the actual data df to compare to predictions 
        
        actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
        colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
        actual_dat$reg_week <- forecast_week_nums
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_agg, by = c("reg_week"="reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "tot"
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "GA_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "CC_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Dead_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Discharged_tot"
        
        actual_dat$GA <- actual_dat$GA_tot / actual_dat$tot 
        actual_dat$CC <- actual_dat$CC_tot / actual_dat$tot 
        actual_dat$Dead <- actual_dat$Dead_tot / actual_dat$tot 
        actual_dat$Discharged <- actual_dat$Discharged_tot / actual_dat$tot 
        actual_dat$WT <- "actual"
        colnames(means) <- colnames(reg_data_no_stay)[-which(colnames(reg_data_no_stay) == "outcome")]
        
        means$reg_week <- forecast_week_nums
        month_nums <- month(forecast_seq)
        for(k in 1:length(month_nums)){
          current_month <- month_nums[k]
          if(current_month != 1)
            means[k,current_month + 1] <- 1
          
        }
        
        
        means[,1] <- mean(reg_data_no_stay$WaitingTime)
        
        
        
        mean_wt_pred <- predict(ml_stay, newdata = means, "probs")
        if (!is.matrix(mean_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- mean_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          mean_wt_pred <- out_df
        }else if(ncol(mean_wt_pred) != 4){
          
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(mean_wt_pred)){
            out_df[,column] <- mean_wt_pred[,column]
            
          }
          
          mean_wt_pred <- out_df
        }
        
        mean_wt_pred <- as.data.frame(mean_wt_pred)
        mean_wt_pred$patient_group <- patient_group[j]
        mean_wt_pred$ICD <- current_icd
        mean_wt_pred$age <- current_age
        mean_wt_pred$WT <- "mean"
        # Or calculate it at the median
        medians <- means
        medians$WaitingTime <- median(reg_data_no_stay$WaitingTime)
        
        
        median_wt_pred <- predict(ml_stay, newdata = medians, "probs")
        if (!is.matrix(median_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- median_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          median_wt_pred <- out_df
        }else if(ncol(median_wt_pred) != 4){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(median_wt_pred)){
            out_df[,column] <- median_wt_pred[,column]
            
          }
          
          median_wt_pred <- out_df
        }
        
        median_wt_pred <- as.data.frame(median_wt_pred)
        
        median_wt_pred$patient_group <- patient_group[j]
        median_wt_pred$ICD <- current_icd
        median_wt_pred$age <- current_age
        median_wt_pred$WT <- "median"
        
        toc()
        
      }
    } else if(current_ward == "ga"){
      
      tic("Narrowing to GA")
      hes_data <- hes_data[hes_data$cc == 0,]
      toc()
      ## remove na trans
      tic("Removing NAs")
      na_rows <- which(is.na(hes_data$ga_transitions))
      if(length(na_rows) > 0)
        hes_data <- hes_data[-na_rows,]
      
      na_ga_los <- which(is.na(hes_data$GA_LoS))
      if(length(na_ga_los) > 0)
        hes_data <- hes_data[-na_ga_los,]
      toc()
      ## set up week and year fixed effects
      
      tic("Getting the reg weeks")
      hes_data$reg_week <- sapply(hes_data$admidate_MDY,week_num_func,start_week = start_date)
      toc()
      col_names <- NULL
      
  
      tic("Setting the month fixed effects")
      for(k in 2:length(month.name)){
        col_name <- paste(month.name[k], "fixed_effect",sep = "_")
        hes_data[,col_name] <- 0
        hes_data[hes_data$admidate_MM == k,col_name] <- 1
        col_names <- append(col_names, col_name)
        
      }
      toc()
      hes_data <- as.data.frame(hes_data)
      
      ## Set up the outcome variable 
      tic("Setting up the outcome variable")
      hes_data$outcome <- "NA"
      
      #hes_data$outcome <- ifelse(!is.na(hes_data$outcome[hes_data$GA_LoS >= 7]) & hes_data$GA_LoS >= 7, "GA","other")
      hes_data[!is.na(hes_data$GA_LoS) & hes_data$GA_LoS >= 7,]$outcome <- "GA"
      hes_data[hes_data$GA_LoS < 7 & hes_data$ga_transitions == 3,"outcome"] <- "Dead"
      hes_data[hes_data$GA_LoS < 7 & hes_data$ga_transitions == 2,"outcome"] <- "CC"
      hes_data[hes_data$GA_LoS < 7 & hes_data$ga_transitions == 1,"outcome"] <- "Discharged"
      
      
      
      hes_data$stay <- hes_data$GA_LoS
      hes_data[hes_data$stay > 7,"stay"]<-7
      
      
      
      hes_data$outcome <- factor(hes_data$outcome, levels = c("CC","GA","Dead","Discharged"))
      hes_data$outcome <- relevel(hes_data$outcome, ref = "CC")
      
      toc()
      
      ## set up the reg data 
      
      tic("Narrowing reg df")
      reg_data_no_stay <- hes_data[,which(colnames(hes_data) %in% c(col_names, "outcome","WaitingTime","reg_week"))]
      toc()
      
      ## set up the actual df 
      
      tic("Setting up the actual df")
      actual_prop_dat <- reg_data_no_stay
      actual_prop_dat$one <- 1
      actual_prop_dat_agg <- aggregate(one ~ reg_week, data = actual_prop_dat, FUN = sum)
      actual_prop_dat_outcome <- aggregate(one ~ reg_week + outcome, data = actual_prop_dat, FUN = sum)
      toc()
      
      
      # ml depending on WT only
      tic("Fitting Regression")
      ml_stay <- multinom(outcome ~ ., data = reg_data_no_stay)
      toc()
      
      # we can average the probability over WT
      
      tic("Probability forecasting")
      
      # Or calculate it at the mean WT
      means <- data.frame(matrix(data = 0,ncol = ncol(reg_data_no_stay) - 1, nrow = forecast_length))
      forecast_week_start <- week_num_func(forecast_start,"2009-01-01")
      forecast_week_end <- week_num_func(forecast_seq[forecast_length],"2009-01-01")
      forecast_week_nums <- seq(forecast_week_start, forecast_week_end)
      
      
      ## make the actual data df to compare to predictions 
      
      actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
      colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
      actual_dat$reg_week <- forecast_week_nums
      actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_agg, by = c("reg_week"="reg_week"))
      colnames(actual_dat)[ncol(actual_dat)] <- "tot"
      actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA",c("reg_week","one")],
                                     by = c("reg_week" = "reg_week"))
      colnames(actual_dat)[ncol(actual_dat)] <- "GA_tot"
      
      actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC",c("reg_week","one")],
                                     by = c("reg_week" = "reg_week"))
      colnames(actual_dat)[ncol(actual_dat)] <- "CC_tot"
      
      actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead",c("reg_week","one")],
                                     by = c("reg_week" = "reg_week"))
      colnames(actual_dat)[ncol(actual_dat)] <- "Dead_tot"
      
      actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged",c("reg_week","one")],
                                     by = c("reg_week" = "reg_week"))
      colnames(actual_dat)[ncol(actual_dat)] <- "Discharged_tot"
      
      actual_dat$GA <- actual_dat$GA_tot / actual_dat$tot 
      actual_dat$CC <- actual_dat$CC_tot / actual_dat$tot 
      actual_dat$Dead <- actual_dat$Dead_tot / actual_dat$tot 
      actual_dat$Discharged <- actual_dat$Discharged_tot / actual_dat$tot 
      actual_dat$WT <- "actual"
      
      ## back to the mean df
      
      colnames(means) <- colnames(reg_data_no_stay)[-which(colnames(reg_data_no_stay) == "outcome")]
      
      means$reg_week <- forecast_week_nums
      month_nums <- month(forecast_seq)
      for(k in 1:length(month_nums)){
        current_month <- month_nums[k]
        if(current_month != 1)
          means[k,current_month + 1] <- 1
        
      }
      
      
      means[,1] <- mean(reg_data_no_stay$WaitingTime)
      
  
      
      mean_wt_pred <- predict(ml_stay, newdata = means, "probs")
      if (!is.matrix(mean_wt_pred)){
        out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
        colnames(out_df) <- levels(hes_data$outcome)
        single_val <- plyr::count(predict(ml_stay, means))
        out_df[,as.character(single_val[1,1])] <- mean_wt_pred
        if(length(ml_stay$lev) > 1){
          other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
          out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
        }
        
        mean_wt_pred <- out_df
      }else if(ncol(mean_wt_pred) != 4){
        
        out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
        colnames(out_df) <- levels(hes_data$outcome)
        for(column in colnames(mean_wt_pred)){
          out_df[,column] <- mean_wt_pred[,column]
          
        }
        
        mean_wt_pred <- out_df
      }
      
      
      mean_wt_pred <- as.data.frame(mean_wt_pred)
      mean_wt_pred$patient_group <- patient_group[j]
      mean_wt_pred$ICD <- current_icd
      mean_wt_pred$age <- current_age
      mean_wt_pred$WT <- "mean"
      # Or calculate it at the median
      medians <- means
      medians$WaitingTime <- median(reg_data_no_stay$WaitingTime)
      
      
      median_wt_pred <- predict(ml_stay, newdata = medians, "probs")
      if (!is.matrix(median_wt_pred)){
        out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
        colnames(out_df) <- levels(hes_data$outcome)
        single_val <- plyr::count(predict(ml_stay, means))
        out_df[,as.character(single_val[1,1])] <- median_wt_pred
        if(length(ml_stay$lev) > 1){
          other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
          out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
        }
        
        median_wt_pred <- out_df
      }else if(ncol(median_wt_pred) != 4){
        out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
        colnames(out_df) <- levels(hes_data$outcome)
        for(column in colnames(median_wt_pred)){
          out_df[,column] <- median_wt_pred[,column]
          
        }
        
        median_wt_pred <- out_df
      }
      
      
      
      median_wt_pred <- as.data.frame(median_wt_pred)
      
      median_wt_pred$patient_group <- patient_group[j]
      median_wt_pred$ICD <- current_icd
      median_wt_pred$age <- current_age
      median_wt_pred$WT <- "median"
      toc()
      
    }
    
    
    
    mean_wt_pred$reg_week <- forecast_week_nums
    median_wt_pred$reg_week <- forecast_week_nums
    
    tot_df <- dplyr::bind_rows(mean_wt_pred, median_wt_pred, actual_dat)
    
    graphing_df <- melt(tot_df, id.vars = colnames(tot_df)[5:14])
    graphing_df$line_group <- paste(graphing_df$variable, graphing_df$WT, sep = "-")
    
    graph_plot <- ggplot(data = graphing_df, aes(x = reg_week, y = value, group = line_group)) +
      geom_line(aes(color = variable, linetype = WT)) + theme_bw() +
      ggtitle(paste(patient_group[j],"month fixed and trend times"))
    
    print(graph_plot)
    toc()
    
    whole_df <- dplyr::bind_rows(whole_df, graphing_df)
  
  }
  dev.off()  
  
  toc()
  return(whole_df)

}


emergency_regression <- function(patient_group,hes_data_orig, start_date, forecast_length, forecast_start,
                                results_pdf = "D:/Dropbox/COVID19/Overflow/regressions_ga_output.pdf",
                                time_trend = TRUE, month_trend = TRUE){
  ## Input ICD hes data for one ICD one age group one patient group (cc/ga)
  ## to be run in parrallel
  tic("Whole process")
  whole_df <- NULL
  forecast_seq <- seq(as.Date(forecast_start), by = "week", length.out = 52)
  require(stringr)
  require(lubridate)  
  require(reshape2)
  require(tictoc)
  
  pdf(file = results_pdf, paper = "A4r", width = 10, height = 7)  
  
  for(j in 1:length(patient_group)){
    print(paste("On patient group:",patient_group[j],". Number",j, "of",length(patient_group)))
    tic(paste("Running for patient group:",patient_group[j]))
    
    split_patient_group <- str_split_fixed(patient_group[j], "-",3)
    current_icd <- as.integer(split_patient_group[1])
    current_ward <- split_patient_group[2]
    current_age <- as.integer(split_patient_group[3])
    
    hes_data <- hes_data_orig[hes_data_orig$ICD == current_icd & 
                                hes_data_orig$agegrp_v3 == current_age,]
    
    if(current_ward == "cc"){
      
      tic("Narrowing to CC only")
      hes_data <- hes_data[hes_data$cc == 1,]
      toc()
      ## remove na trans
      tic("Removing NAs")
      na_rows <- which(is.na(hes_data$cc_transitions))
      if(length(na_rows) > 0)
        hes_data <- hes_data[-na_rows,]
      
      na_ga_los <- which(is.na(hes_data$cc_LoS))
      if(length(na_ga_los) > 0)
        hes_data <- hes_data[-na_ga_los,]
      toc()
      
      ## set up transitions 
      tic("Set up outcome variable")
      hes_data$outcome <- NA
      
      hes_data[hes_data$cc_LoS >= 7 ,"outcome"]<-"CC"
      hes_data[hes_data$cc_LoS < 7 & hes_data$cc_transitions == 3,"outcome"] <- "Dead"
      hes_data[hes_data$cc_LoS < 7 & hes_data$cc_transitions == 2,"outcome"] <- "GA"
      hes_data[hes_data$cc_LoS < 7 & hes_data$cc_transitions == 1,"outcome"] <- "Discharged"
      
      
      
      hes_data$stay <- hes_data$cc_LoS
      hes_data[hes_data$stay > 7,"stay"]<-7
      
      
      
      hes_data$outcome <- factor(hes_data$outcome, levels = c("CC","Dead","GA","Discharged"))
      hes_data$outcome <- relevel(hes_data$outcome, ref = "GA")
      
      
      toc()
      
    }else{
      tic("Narrowing to GA only")
      hes_data <- hes_data[hes_data$cc == 0,]
      toc()
      ## remove na trans
      tic("Removing NAs")
      na_rows <- which(is.na(hes_data$ga_transitions))
      if(length(na_rows) > 0)
        hes_data <- hes_data[-na_rows,]
      
      na_ga_los <- which(is.na(hes_data$GA_LoS))
      if(length(na_ga_los) > 0)
        hes_data <- hes_data[-na_ga_los,]
      toc()
      
      ## set up transitions 
      tic("Set up outcome variable")
      hes_data$outcome <- NA
      
      hes_data[hes_data$GA_LoS >= 7 ,"outcome"]<-"GA"
      hes_data[hes_data$GA_LoS < 7 & hes_data$ga_transitions == 3,"outcome"] <- "Dead"
      hes_data[hes_data$GA_LoS < 7 & hes_data$ga_transitions == 2,"outcome"] <- "CC"
      hes_data[hes_data$GA_LoS < 7 & hes_data$ga_transitions == 1,"outcome"] <- "Discharged"
      
      
      
      hes_data$stay <- hes_data$GA_LoS
      hes_data[hes_data$stay > 7,"stay"]<-7
      
      
      
      hes_data$outcome <- factor(hes_data$outcome, levels = c("CC","Dead","GA","Discharged"))
      hes_data$outcome <- relevel(hes_data$outcome, ref = "GA")
      
      
      toc()
    }
      ## set up week and year fixed effects
      
      forecast_week_start <- week_num_func(forecast_start,"2009-01-01")
      forecast_week_end <- week_num_func(forecast_seq[forecast_length],"2009-01-01")
      forecast_week_nums <- seq(forecast_week_start, forecast_week_end)
      
      
      
      
        tic("Getting the reg week data")
        hes_data$reg_week <- sapply(hes_data$admidate_MDY,week_num_func,start_week = start_date)
        toc()
        
      
      if(month_trend == TRUE){
        col_names <- NULL
        
        
        
        for(k in 2:length(month.name)){
          col_name <- paste(month.name[k], "fixed_effect",sep = "_")
          hes_data[,col_name] <- 0
          hes_data[hes_data$admidate_MM == k,col_name] <- 1
          col_names <- append(col_names, col_name)
          
        }
        
      }
      ## Set up the outcome variable 
      ## set up the reg data 
      tic("Reg set up")
      
      if(time_trend == TRUE & month_trend == TRUE){
        reg_data_no_stay <- hes_data[,which(colnames(hes_data) %in% c(col_names, "outcome","reg_week"))]
      }else if(time_trend == TRUE & month_trend == FALSE){
        reg_data_no_stay <- hes_data[,which(colnames(hes_data) %in% c("outcome","reg_week"))]
      }else if(time_trend == FALSE & month_trend == TRUE){
        reg_data_no_stay <- hes_data[,which(colnames(hes_data) %in% c(col_names, "outcome"))]
      }else{
        reg_data_no_stay <- hes_data[,which(colnames(hes_data) %in% c("outcome"))]
      }
      
      
      
      
      
      actual_prop_dat <- hes_data[,which(colnames(hes_data) %in% c("outcome","reg_week"))]
      actual_prop_dat$one <- 1
      actual_prop_dat_agg <- aggregate(one ~ reg_week, data = actual_prop_dat, FUN = sum)
      actual_prop_dat_outcome <- aggregate(one ~ reg_week + outcome, data = actual_prop_dat, FUN = sum)
      
      
      toc()
      
      
      
      outcomes_in_dat <- plyr::count(reg_data_no_stay$outcome)
      ## Checking if enough outcomes for model else will return one for probs
      
      if(nrow(outcomes_in_dat) < 2){
        out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
        colnames(out_df) <- levels(hes_data$outcome)
        out_df[,outcomes_in_dat[1,1]] <- 1
        mean_wt_pred <- as.data.frame(out_df)
        mean_wt_pred$patient_group <- patient_group[j]
        mean_wt_pred$ICD <- current_icd
        mean_wt_pred$age <- current_age
        mean_wt_pred$WT <- "mean"
        median_wt_pred <- mean_wt_pred
        median_wt_pred$WT <- "median"
        actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
        colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
        actual_dat$reg_week <- forecast_week_nums
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_agg, by = c("reg_week"="reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "tot"
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "GA_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "CC_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Dead_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Discharged_tot"
        
        actual_dat$GA <- actual_dat$GA_tot / actual_dat$tot 
        actual_dat$CC <- actual_dat$CC_tot / actual_dat$tot 
        actual_dat$Dead <- actual_dat$Dead_tot / actual_dat$tot 
        actual_dat$Discharged <- actual_dat$Discharged_tot / actual_dat$tot 
        actual_dat$WT <- "actual"
        
        
      }else{
        
        
        
        # ml depending on WT only
        
        if(time_trend | month_trend){
        
          tic("Regreesion model fitting")
          ml_stay <- multinom(outcome ~ ., data = reg_data_no_stay)
          toc()
        }else{
          mean_wt_pred <- data.frame(matrix(data = 0, nrow = forecast_length, ncol = 4))
          colnames(mean_wt_pred) <- levels(hes_data$outcome)
          denominator <- sum(actual_prop_dat_agg$one)
          GA_num <- sum(actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA", "one"])
          cc_num <- sum(actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC", "one"])
          dead_num <- sum(actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead", "one"])
          dis_num <- sum(actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged", "one"])
          
          mean_wt_pred$GA <- GA_num / denominator
          mean_wt_pred$CC <- cc_num / denominator 
          mean_wt_pred$Dead <- dead_num / denominator
          mean_wt_pred$Discharged <- dis_num / denominator
          mean_wt_pred$patient_group <- patient_group[j]
          mean_wt_pred$ICD <- current_icd
          mean_wt_pred$age <- current_age
          mean_wt_pred$WT <- "mean"
          
          median_wt_pred <- mean_wt_pred
          median_wt_pred$WT <- "median"
          
          actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
          colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
          actual_dat$reg_week <- forecast_week_nums
          actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_agg, by = c("reg_week"="reg_week"))
          colnames(actual_dat)[ncol(actual_dat)] <- "tot"
          actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA",c("reg_week","one")],
                                         by = c("reg_week" = "reg_week"))
          colnames(actual_dat)[ncol(actual_dat)] <- "GA_tot"
          
          actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC",c("reg_week","one")],
                                         by = c("reg_week" = "reg_week"))
          colnames(actual_dat)[ncol(actual_dat)] <- "CC_tot"
          
          actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead",c("reg_week","one")],
                                         by = c("reg_week" = "reg_week"))
          colnames(actual_dat)[ncol(actual_dat)] <- "Dead_tot"
          
          actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged",c("reg_week","one")],
                                         by = c("reg_week" = "reg_week"))
          colnames(actual_dat)[ncol(actual_dat)] <- "Discharged_tot"
          
          actual_dat$GA <- actual_dat$GA_tot / actual_dat$tot 
          actual_dat$CC <- actual_dat$CC_tot / actual_dat$tot 
          actual_dat$Dead <- actual_dat$Dead_tot / actual_dat$tot 
          actual_dat$Discharged <- actual_dat$Discharged_tot / actual_dat$tot 
          actual_dat$WT <- "actual"
          
          
          
        }
        if(time_trend & month_trend){
        # we can average the probability over WT
        
        
        
        # Or calculate it at the mean WT
        tic("Probability creation")
        means <- data.frame(matrix(data = 0,ncol = ncol(reg_data_no_stay) - 1, nrow = forecast_length))
        
        ## make the actual data df to compare to predictions 
        
        actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
        colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
        actual_dat$reg_week <- forecast_week_nums
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_agg, by = c("reg_week"="reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "tot"
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "GA_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "CC_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Dead_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Discharged_tot"
        
        actual_dat$GA <- actual_dat$GA_tot / actual_dat$tot 
        actual_dat$CC <- actual_dat$CC_tot / actual_dat$tot 
        actual_dat$Dead <- actual_dat$Dead_tot / actual_dat$tot 
        actual_dat$Discharged <- actual_dat$Discharged_tot / actual_dat$tot 
        actual_dat$WT <- "actual"
        colnames(means) <- colnames(reg_data_no_stay)[-which(colnames(reg_data_no_stay) == "outcome")]
        
        means$reg_week <- forecast_week_nums
        month_nums <- month(forecast_seq)
        for(k in 1:length(month_nums)){
          current_month <- month_nums[k]
          if(current_month != 1)
            means[k,current_month + 1] <- 1
          
        }
        
        
        mean_wt_pred <- predict(ml_stay, newdata = means, "probs")
        if (!is.matrix(mean_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- mean_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          mean_wt_pred <- out_df
        }else if(ncol(mean_wt_pred) != 4){
          
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(mean_wt_pred)){
            out_df[,column] <- mean_wt_pred[,column]
            
          }
          
          mean_wt_pred <- out_df
        }
        
        mean_wt_pred <- as.data.frame(mean_wt_pred)
        mean_wt_pred$patient_group <- patient_group[j]
        mean_wt_pred$ICD <- current_icd
        mean_wt_pred$age <- current_age
        mean_wt_pred$WT <- "mean"
        # Or calculate it at the median
        medians <- means
        
        median_wt_pred <- predict(ml_stay, newdata = medians, "probs")
        if (!is.matrix(median_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- median_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          median_wt_pred <- out_df
        }else if(ncol(median_wt_pred) != 4){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(median_wt_pred)){
            out_df[,column] <- median_wt_pred[,column]
            
          }
          
          median_wt_pred <- out_df
        }
        
        median_wt_pred <- as.data.frame(median_wt_pred)
        
        median_wt_pred$patient_group <- patient_group[j]
        median_wt_pred$ICD <- current_icd
        median_wt_pred$age <- current_age
        median_wt_pred$WT <- "median"
        
        toc()
        }else if(time_trend & month_trend == FALSE){
          tic("Probability creation")
          means <- data.frame(matrix(data = 0,ncol = ncol(reg_data_no_stay) - 1, nrow = forecast_length))
          
          ## make the actual data df to compare to predictions 
          
          actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
          colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
          actual_dat$reg_week <- forecast_week_nums
          actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_agg, by = c("reg_week"="reg_week"))
          colnames(actual_dat)[ncol(actual_dat)] <- "tot"
          actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA",c("reg_week","one")],
                                         by = c("reg_week" = "reg_week"))
          colnames(actual_dat)[ncol(actual_dat)] <- "GA_tot"
          
          actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC",c("reg_week","one")],
                                         by = c("reg_week" = "reg_week"))
          colnames(actual_dat)[ncol(actual_dat)] <- "CC_tot"
          
          actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead",c("reg_week","one")],
                                         by = c("reg_week" = "reg_week"))
          colnames(actual_dat)[ncol(actual_dat)] <- "Dead_tot"
          
          actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged",c("reg_week","one")],
                                         by = c("reg_week" = "reg_week"))
          colnames(actual_dat)[ncol(actual_dat)] <- "Discharged_tot"
          
          actual_dat$GA <- actual_dat$GA_tot / actual_dat$tot 
          actual_dat$CC <- actual_dat$CC_tot / actual_dat$tot 
          actual_dat$Dead <- actual_dat$Dead_tot / actual_dat$tot 
          actual_dat$Discharged <- actual_dat$Discharged_tot / actual_dat$tot 
          actual_dat$WT <- "actual"
          colnames(means) <- colnames(reg_data_no_stay)[-which(colnames(reg_data_no_stay) == "outcome")]
          
          means$reg_week <- forecast_week_nums
          
          mean_wt_pred <- predict(ml_stay, newdata = means, "probs")
          if (!is.matrix(mean_wt_pred)){
            out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
            colnames(out_df) <- levels(hes_data$outcome)
            single_val <- plyr::count(predict(ml_stay, means))
            out_df[,as.character(single_val[1,1])] <- mean_wt_pred
            if(length(ml_stay$lev) > 1){
              other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
              out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
            }
            
            mean_wt_pred <- out_df
          }else if(ncol(mean_wt_pred) != 4){
            
            out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
            colnames(out_df) <- levels(hes_data$outcome)
            for(column in colnames(mean_wt_pred)){
              out_df[,column] <- mean_wt_pred[,column]
              
            }
            
            mean_wt_pred <- out_df
          }
          
          mean_wt_pred <- as.data.frame(mean_wt_pred)
          mean_wt_pred$patient_group <- patient_group[j]
          mean_wt_pred$ICD <- current_icd
          mean_wt_pred$age <- current_age
          mean_wt_pred$WT <- "mean"
          # Or calculate it at the median
          medians <- means
          
          median_wt_pred <- predict(ml_stay, newdata = medians, "probs")
          if (!is.matrix(median_wt_pred)){
            out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
            colnames(out_df) <- levels(hes_data$outcome)
            single_val <- plyr::count(predict(ml_stay, means))
            out_df[,as.character(single_val[1,1])] <- median_wt_pred
            if(length(ml_stay$lev) > 1){
              other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
              out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
            }
            
            median_wt_pred <- out_df
          }else if(ncol(median_wt_pred) != 4){
            out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
            colnames(out_df) <- levels(hes_data$outcome)
            for(column in colnames(median_wt_pred)){
              out_df[,column] <- median_wt_pred[,column]
              
            }
            
            median_wt_pred <- out_df
          }
          
          median_wt_pred <- as.data.frame(median_wt_pred)
          
          median_wt_pred$patient_group <- patient_group[j]
          median_wt_pred$ICD <- current_icd
          median_wt_pred$age <- current_age
          median_wt_pred$WT <- "median"
          
          toc()
          
          
          
          
        }else if(time_trend == FALSE & month_trend){
          tic("Probability creation")
          means <- data.frame(matrix(data = 0,ncol = ncol(reg_data_no_stay) - 1, nrow = forecast_length))
          
          ## make the actual data df to compare to predictions 
          
          actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
          colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
          actual_dat$reg_week <- forecast_week_nums
          actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_agg, by = c("reg_week"="reg_week"))
          colnames(actual_dat)[ncol(actual_dat)] <- "tot"
          actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA",c("reg_week","one")],
                                         by = c("reg_week" = "reg_week"))
          colnames(actual_dat)[ncol(actual_dat)] <- "GA_tot"
          
          actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC",c("reg_week","one")],
                                         by = c("reg_week" = "reg_week"))
          colnames(actual_dat)[ncol(actual_dat)] <- "CC_tot"
          
          actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead",c("reg_week","one")],
                                         by = c("reg_week" = "reg_week"))
          colnames(actual_dat)[ncol(actual_dat)] <- "Dead_tot"
          
          actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged",c("reg_week","one")],
                                         by = c("reg_week" = "reg_week"))
          colnames(actual_dat)[ncol(actual_dat)] <- "Discharged_tot"
          
          actual_dat$GA <- actual_dat$GA_tot / actual_dat$tot 
          actual_dat$CC <- actual_dat$CC_tot / actual_dat$tot 
          actual_dat$Dead <- actual_dat$Dead_tot / actual_dat$tot 
          actual_dat$Discharged <- actual_dat$Discharged_tot / actual_dat$tot 
          actual_dat$WT <- "actual"
          colnames(means) <- colnames(reg_data_no_stay)[-which(colnames(reg_data_no_stay) == "outcome")]
          
          
          month_nums <- month(forecast_seq)
          for(k in 1:length(month_nums)){
            current_month <- month_nums[k]
            if(current_month != 1)
              means[k,current_month + 1] <- 1
            
          }
          
          
          mean_wt_pred <- predict(ml_stay, newdata = means, "probs")
          if (!is.matrix(mean_wt_pred)){
            out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
            colnames(out_df) <- levels(hes_data$outcome)
            single_val <- plyr::count(predict(ml_stay, means))
            out_df[,as.character(single_val[1,1])] <- mean_wt_pred
            if(length(ml_stay$lev) > 1){
              other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
              out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
            }
            
            mean_wt_pred <- out_df
          }else if(ncol(mean_wt_pred) != 4){
            
            out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
            colnames(out_df) <- levels(hes_data$outcome)
            for(column in colnames(mean_wt_pred)){
              out_df[,column] <- mean_wt_pred[,column]
              
            }
            
            mean_wt_pred <- out_df
          }
          
          mean_wt_pred <- as.data.frame(mean_wt_pred)
          mean_wt_pred$patient_group <- patient_group[j]
          mean_wt_pred$ICD <- current_icd
          mean_wt_pred$age <- current_age
          mean_wt_pred$WT <- "mean"
          # Or calculate it at the median
          medians <- means
          
          median_wt_pred <- predict(ml_stay, newdata = medians, "probs")
          if (!is.matrix(median_wt_pred)){
            out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
            colnames(out_df) <- levels(hes_data$outcome)
            single_val <- plyr::count(predict(ml_stay, means))
            out_df[,as.character(single_val[1,1])] <- median_wt_pred
            if(length(ml_stay$lev) > 1){
              other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
              out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
            }
            
            median_wt_pred <- out_df
          }else if(ncol(median_wt_pred) != 4){
            out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
            colnames(out_df) <- levels(hes_data$outcome)
            for(column in colnames(median_wt_pred)){
              out_df[,column] <- median_wt_pred[,column]
              
            }
            
            median_wt_pred <- out_df
          }
          
          median_wt_pred <- as.data.frame(median_wt_pred)
          
          median_wt_pred$patient_group <- patient_group[j]
          median_wt_pred$ICD <- current_icd
          median_wt_pred$age <- current_age
          median_wt_pred$WT <- "median"
          
          toc()
          
          
        }
        
      }
    
    mean_wt_pred$reg_week <- forecast_week_nums
    median_wt_pred$reg_week <- forecast_week_nums
    
    tot_df <- dplyr::bind_rows(mean_wt_pred, median_wt_pred, actual_dat)
    
    graphing_df <- melt(tot_df, id.vars = colnames(tot_df)[5:14])
    graphing_df$line_group <- paste(graphing_df$variable, graphing_df$WT, sep = "-")
    
    graph_plot <- ggplot(data = graphing_df, aes(x = reg_week, y = value, group = line_group)) +
      geom_line(aes(color = variable, linetype = WT)) + theme_bw() +
      ggtitle(paste(patient_group[j],"month fixed and trend times"))
    
    print(graph_plot)
    toc()
    
    whole_df <- dplyr::bind_rows(whole_df, graphing_df)
    
  }
  dev.off()  
  
  toc()
  return(whole_df)
  
}




elective_regression <- function(patient_group,hes_data_orig, start_date, forecast_length, forecast_start,
                                 results_pdf = "D:/Dropbox/COVID19/Overflow/regressions_ga_output.pdf",
                                 time_trend = TRUE, month_trend = TRUE){
  ## Input ICD hes data for one ICD one age group one patient group (cc/ga)
  ## to be run in parrallel
  tic("Whole process")
  whole_df <- NULL
  forecast_seq <- seq(as.Date(forecast_start), by = "week", length.out = 52)
  require(stringr)
  require(lubridate)  
  require(reshape2)
  require(tictoc)
  
  pdf(file = results_pdf, paper = "A4r", width = 10, height = 7)  

  for(j in 1:length(patient_group)){
    print(paste("On patient group:",patient_group[j],". Number",j, "of",length(patient_group)))
    tic(paste("Running for patient group:",patient_group[j]))
    
    split_patient_group <- str_split_fixed(patient_group[j], "-",3)
    current_icd <- as.integer(split_patient_group[1])
    current_ward <- split_patient_group[2]
    current_age <- as.integer(split_patient_group[3])
    
    hes_data <- hes_data_orig[hes_data_orig$ICD == current_icd & 
                                hes_data_orig$agegrp_v3 == current_age,]
    
    if(current_ward == "cc"){
      
      tic("Narrowing to CC only")
      hes_data <- hes_data[hes_data$cc == 1,]
      toc()
      ## remove na trans
      tic("Removing NAs")
      na_rows <- which(is.na(hes_data$cc_transitions))
      if(length(na_rows) > 0)
        hes_data <- hes_data[-na_rows,]
      
      na_cc_los <- which(is.na(hes_data$cc_LoS))
      if(length(na_cc_los) > 0)
        hes_data <- hes_data[-na_cc_los,]
      
      na_cc_WT <- which(is.na(hes_data$WaitingTime))
      if(length(na_cc_WT) > 0)
        hes_data <- hes_data[-na_cc_WT,]
      
      
      toc()
      
      ## set up transitions 
      tic("Set up outcome variable")
      hes_data$outcome <- NA
      
      hes_data[hes_data$cc_LoS >= 7 ,"outcome"]<-"CC"
      hes_data[hes_data$cc_LoS < 7 & hes_data$cc_transitions == 3,"outcome"] <- "Dead"
      hes_data[hes_data$cc_LoS < 7 & hes_data$cc_transitions == 2,"outcome"] <- "GA"
      hes_data[hes_data$cc_LoS < 7 & hes_data$cc_transitions == 1,"outcome"] <- "Discharged"
      
      
      
      hes_data$stay <- hes_data$cc_LoS
      hes_data[hes_data$stay > 7,"stay"]<-7
      
      
      
      hes_data$outcome <- factor(hes_data$outcome, levels = c("CC","Dead","GA","Discharged"))
      hes_data$outcome <- relevel(hes_data$outcome, ref = "GA")
      
      
      toc()
      
    }else{
      ## remove na trans
      tic("Removing NAs")
      na_rows <- which(is.na(hes_data$ga_transitions))
      if(length(na_rows) > 0)
        hes_data <- hes_data[-na_rows,]
      
      na_ga_los <- which(is.na(hes_data$GA_LoS))
      if(length(na_ga_los) > 0)
        hes_data <- hes_data[-na_ga_los,]
      
      na_ga_WT <- which(is.na(hes_data$WaitingTime))
      if(length(na_ga_WT) > 0)
        hes_data <- hes_data[-na_ga_WT,]
      
      
      toc()
      
      ## set up transitions 
      tic("Set up outcome variable")
      hes_data$outcome <- NA
      
      hes_data[hes_data$GA_LoS >= 7 ,"outcome"]<-"GA"
      hes_data[hes_data$GA_LoS < 7 & hes_data$ga_transitions == 3,"outcome"] <- "Dead"
      hes_data[hes_data$GA_LoS < 7 & hes_data$ga_transitions == 2,"outcome"] <- "CC"
      hes_data[hes_data$GA_LoS < 7 & hes_data$ga_transitions == 1,"outcome"] <- "Discharged"
      
      
      
      hes_data$stay <- hes_data$GA_LoS
      hes_data[hes_data$stay > 7,"stay"]<-7
      
      
      
      hes_data$outcome <- factor(hes_data$outcome, levels = c("CC","Dead","GA","Discharged"))
      hes_data$outcome <- relevel(hes_data$outcome, ref = "GA")
      
      
      toc()
    }
    ## set up week and year fixed effects
    
    forecast_week_start <- week_num_func(forecast_start,"2009-01-01")
    forecast_week_end <- week_num_func(forecast_seq[forecast_length],"2009-01-01")
    forecast_week_nums <- seq(forecast_week_start, forecast_week_end)
    
    
    
    
    tic("Getting the reg week data")
    hes_data$reg_week <- sapply(hes_data$admidate_MDY,week_num_func,start_week = start_date)
    toc()
    
    
    if(month_trend == TRUE){
      col_names <- NULL
      
      
      
      for(k in 2:length(month.name)){
        col_name <- paste(month.name[k], "fixed_effect",sep = "_")
        hes_data[,col_name] <- 0
        hes_data[hes_data$admidate_MM == k,col_name] <- 1
        col_names <- append(col_names, col_name)
        
      }
      
    }
    ## Set up the outcome variable 
    ## set up the reg data 
    tic("Reg set up")
    
    if(time_trend == TRUE & month_trend == TRUE){
      reg_data_no_stay <- hes_data[,which(colnames(hes_data) %in% c(col_names, "outcome","reg_week","WaitingTime"))]
    }else if(time_trend == TRUE & month_trend == FALSE){
      reg_data_no_stay <- hes_data[,which(colnames(hes_data) %in% c("outcome","reg_week","WaitingTime"))]
    }else if(time_trend == FALSE & month_trend == TRUE){
      reg_data_no_stay <- hes_data[,which(colnames(hes_data) %in% c(col_names, "outcome","WaitingTime"))]
    }else{
      reg_data_no_stay <- hes_data[,which(colnames(hes_data) %in% c("outcome","WaitingTime"))]
    }
    
    
    
    
    
    actual_prop_dat <- hes_data[,which(colnames(hes_data) %in% c("outcome","reg_week"))]
    actual_prop_dat$one <- 1
    actual_prop_dat_agg <- aggregate(one ~ reg_week, data = actual_prop_dat, FUN = sum)
    actual_prop_dat_outcome <- aggregate(one ~ reg_week + outcome, data = actual_prop_dat, FUN = sum)
    
    
    toc()
    
    
    
    outcomes_in_dat <- plyr::count(reg_data_no_stay$outcome)
    ## Checking if enough outcomes for model else will return one for probs
    
    if(nrow(outcomes_in_dat) < 2){
      out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
      colnames(out_df) <- levels(hes_data$outcome)
      out_df[,outcomes_in_dat[1,1]] <- 1
      mean_wt_pred <- as.data.frame(out_df)
      mean_wt_pred$patient_group <- patient_group[j]
      mean_wt_pred$ICD <- current_icd
      mean_wt_pred$age <- current_age
      mean_wt_pred$WT <- "mean"
      median_wt_pred <- mean_wt_pred
      median_wt_pred$WT <- "median"
      actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
      colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
      actual_dat$reg_week <- forecast_week_nums
      actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_agg, by = c("reg_week"="reg_week"))
      colnames(actual_dat)[ncol(actual_dat)] <- "tot"
      actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA",c("reg_week","one")],
                                     by = c("reg_week" = "reg_week"))
      colnames(actual_dat)[ncol(actual_dat)] <- "GA_tot"
      
      actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC",c("reg_week","one")],
                                     by = c("reg_week" = "reg_week"))
      colnames(actual_dat)[ncol(actual_dat)] <- "CC_tot"
      
      actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead",c("reg_week","one")],
                                     by = c("reg_week" = "reg_week"))
      colnames(actual_dat)[ncol(actual_dat)] <- "Dead_tot"
      
      actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged",c("reg_week","one")],
                                     by = c("reg_week" = "reg_week"))
      colnames(actual_dat)[ncol(actual_dat)] <- "Discharged_tot"
      
      actual_dat$GA <- actual_dat$GA_tot / actual_dat$tot 
      actual_dat$CC <- actual_dat$CC_tot / actual_dat$tot 
      actual_dat$Dead <- actual_dat$Dead_tot / actual_dat$tot 
      actual_dat$Discharged <- actual_dat$Discharged_tot / actual_dat$tot 
      actual_dat$WT <- "actual"
      
      
    }else{
      
      
      
      # ml depending on WT only
      
      tic("Regreesion model fitting")
      ml_stay <- multinom(outcome ~ ., data = reg_data_no_stay)
      

      wt_0_dat <- data.frame(matrix(0, ncol = ncol(reg_data_no_stay) - 1, nrow = 1))
      colnames(wt_0_dat) <- colnames(reg_data_no_stay)[-which(colnames(reg_data_no_stay) == "outcome")]
      wt_1_dat <- wt_0_dat
      wt_1_dat$WaitingTime <- 1
      
      wt_0_pred <- predict(ml_stay, newdata = wt_0_dat, "probs")
      wt_1_pred <- predict(ml_stay, newdata = wt_1_dat, "probs")
      
      wt_diff <- wt_1_pred - wt_0_pred
      
      if(length(wt_diff) != 4){
        
        missing_factors <- which(!(levels(reg_data_no_stay$outcome) %in% names(wt_diff)))
        old_names <- names(wt_diff)
        wt_diff <- append(wt_diff,rep(0, length(missing_factors)))
        names(wt_diff) <- c(old_names, levels(reg_data_no_stay$outcome)[missing_factors])
        
      }
      current_coef_df <- data.frame(matrix(0, ncol = 5, nrow = 4))
      colnames(current_coef_df) <- c("transition","coef","patient_group","ICD","age")
      current_coef_df$transition <- names(wt_diff)
      current_coef_df$coef <- wt_diff
      current_coef_df$patient_group <- patient_group[j]
      current_coef_df$ICD <- current_icd
      current_coef_df$age <- current_age
      
      toc()
      
      if(time_trend & month_trend){
        # we can average the probability over WT
        
        
        
        # Or calculate it at the mean WT
        tic("Probability creation")
        means <- data.frame(matrix(data = 0,ncol = ncol(reg_data_no_stay) - 1, nrow = forecast_length))
        
        ## make the actual data df to compare to predictions 
        
        actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
        colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
        actual_dat$reg_week <- forecast_week_nums
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_agg, by = c("reg_week"="reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "tot"
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "GA_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "CC_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Dead_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Discharged_tot"
        
        actual_dat$GA <- actual_dat$GA_tot / actual_dat$tot 
        actual_dat$CC <- actual_dat$CC_tot / actual_dat$tot 
        actual_dat$Dead <- actual_dat$Dead_tot / actual_dat$tot 
        actual_dat$Discharged <- actual_dat$Discharged_tot / actual_dat$tot 
        actual_dat$WT <- "actual"
        colnames(means) <- colnames(reg_data_no_stay)[-which(colnames(reg_data_no_stay) == "outcome")]
        
        means$WaitingTime <- mean(reg_data_no_stay$WaitingTime, na.rm = TRUE)
        means$reg_week <- forecast_week_nums
        month_nums <- month(forecast_seq)
        for(k in 1:length(month_nums)){
          current_month <- month_nums[k]
          if(current_month != 1)
            means[k,current_month + 1] <- 1
          
        }
        
        
        mean_wt_pred <- predict(ml_stay, newdata = means, "probs")
        if (!is.matrix(mean_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- mean_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          mean_wt_pred <- out_df
        }else if(ncol(mean_wt_pred) != 4){
          
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(mean_wt_pred)){
            out_df[,column] <- mean_wt_pred[,column]
            
          }
          
          mean_wt_pred <- out_df
        }
        
        mean_wt_pred <- as.data.frame(mean_wt_pred)
        mean_wt_pred$patient_group <- patient_group[j]
        mean_wt_pred$ICD <- current_icd
        mean_wt_pred$age <- current_age
        mean_wt_pred$WT <- "mean"
        # Or calculate it at the median
        medians <- means
        medians$WaitingTime <- median(reg_data_no_stay$WaitingTime)
        
        median_wt_pred <- predict(ml_stay, newdata = medians, "probs")
        if (!is.matrix(median_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- median_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          median_wt_pred <- out_df
        }else if(ncol(median_wt_pred) != 4){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(median_wt_pred)){
            out_df[,column] <- median_wt_pred[,column]
            
          }
          
          median_wt_pred <- out_df
        }
        
        median_wt_pred <- as.data.frame(median_wt_pred)
        
        median_wt_pred$patient_group <- patient_group[j]
        median_wt_pred$ICD <- current_icd
        median_wt_pred$age <- current_age
        median_wt_pred$WT <- "median"
        
        toc()
      }else if(time_trend & month_trend == FALSE){
        tic("Probability creation")
        means <- data.frame(matrix(data = 0,ncol = ncol(reg_data_no_stay) - 1, nrow = forecast_length))
        
        ## make the actual data df to compare to predictions 
        
        actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
        colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
        actual_dat$reg_week <- forecast_week_nums
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_agg, by = c("reg_week"="reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "tot"
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "GA_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "CC_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Dead_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Discharged_tot"
        
        actual_dat$GA <- actual_dat$GA_tot / actual_dat$tot 
        actual_dat$CC <- actual_dat$CC_tot / actual_dat$tot 
        actual_dat$Dead <- actual_dat$Dead_tot / actual_dat$tot 
        actual_dat$Discharged <- actual_dat$Discharged_tot / actual_dat$tot 
        actual_dat$WT <- "actual"
        colnames(means) <- colnames(reg_data_no_stay)[-which(colnames(reg_data_no_stay) == "outcome")]
        
        means$reg_week <- forecast_week_nums
        
        mean_wt_pred <- predict(ml_stay, newdata = means, "probs")
        if (!is.matrix(mean_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- mean_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          mean_wt_pred <- out_df
        }else if(ncol(mean_wt_pred) != 4){
          
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(mean_wt_pred)){
            out_df[,column] <- mean_wt_pred[,column]
            
          }
          
          mean_wt_pred <- out_df
        }
        
        mean_wt_pred <- as.data.frame(mean_wt_pred)
        mean_wt_pred$patient_group <- patient_group[j]
        mean_wt_pred$ICD <- current_icd
        mean_wt_pred$age <- current_age
        mean_wt_pred$WT <- "mean"
        # Or calculate it at the median
        medians <- means
        
        median_wt_pred <- predict(ml_stay, newdata = medians, "probs")
        if (!is.matrix(median_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- median_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          median_wt_pred <- out_df
        }else if(ncol(median_wt_pred) != 4){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(median_wt_pred)){
            out_df[,column] <- median_wt_pred[,column]
            
          }
          
          median_wt_pred <- out_df
        }
        
        median_wt_pred <- as.data.frame(median_wt_pred)
        
        median_wt_pred$patient_group <- patient_group[j]
        median_wt_pred$ICD <- current_icd
        median_wt_pred$age <- current_age
        median_wt_pred$WT <- "median"
        
        toc()
        
        
        
        
      }else if(time_trend == FALSE & month_trend){
        tic("Probability creation")
        means <- data.frame(matrix(data = 0,ncol = ncol(reg_data_no_stay) - 1, nrow = forecast_length))
        
        ## make the actual data df to compare to predictions 
        
        actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
        colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
        actual_dat$reg_week <- forecast_week_nums
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_agg, by = c("reg_week"="reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "tot"
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "GA_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "CC_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Dead_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Discharged_tot"
        
        actual_dat$GA <- actual_dat$GA_tot / actual_dat$tot 
        actual_dat$CC <- actual_dat$CC_tot / actual_dat$tot 
        actual_dat$Dead <- actual_dat$Dead_tot / actual_dat$tot 
        actual_dat$Discharged <- actual_dat$Discharged_tot / actual_dat$tot 
        actual_dat$WT <- "actual"
        colnames(means) <- colnames(reg_data_no_stay)[-which(colnames(reg_data_no_stay) == "outcome")]
        
        
        month_nums <- month(forecast_seq)
        for(k in 1:length(month_nums)){
          current_month <- month_nums[k]
          if(current_month != 1)
            means[k,current_month + 1] <- 1
          
        }
        
        
        mean_wt_pred <- predict(ml_stay, newdata = means, "probs")
        if (!is.matrix(mean_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- mean_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          mean_wt_pred <- out_df
        }else if(ncol(mean_wt_pred) != 4){
          
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(mean_wt_pred)){
            out_df[,column] <- mean_wt_pred[,column]
            
          }
          
          mean_wt_pred <- out_df
        }
        
        mean_wt_pred <- as.data.frame(mean_wt_pred)
        mean_wt_pred$patient_group <- patient_group[j]
        mean_wt_pred$ICD <- current_icd
        mean_wt_pred$age <- current_age
        mean_wt_pred$WT <- "mean"
        # Or calculate it at the median
        medians <- means
        
        median_wt_pred <- predict(ml_stay, newdata = medians, "probs")
        if (!is.matrix(median_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- median_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          median_wt_pred <- out_df
        }else if(ncol(median_wt_pred) != 4){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(median_wt_pred)){
            out_df[,column] <- median_wt_pred[,column]
            
          }
          
          median_wt_pred <- out_df
        }
        
        median_wt_pred <- as.data.frame(median_wt_pred)
        
        median_wt_pred$patient_group <- patient_group[j]
        median_wt_pred$ICD <- current_icd
        median_wt_pred$age <- current_age
        median_wt_pred$WT <- "median"
        
        toc()
        
        
      }else if( month_trend == FALSE & time_trend == FALSE){
        tic("Probability creation")
        means <- data.frame(matrix(data = 0,ncol = ncol(reg_data_no_stay) - 1, nrow = forecast_length))
        
        ## make the actual data df to compare to predictions 
        
        actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
        colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
        actual_dat$reg_week <- forecast_week_nums
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_agg, by = c("reg_week"="reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "tot"
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "GA_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "CC_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Dead_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Discharged_tot"
        
        actual_dat$GA <- actual_dat$GA_tot / actual_dat$tot 
        actual_dat$CC <- actual_dat$CC_tot / actual_dat$tot 
        actual_dat$Dead <- actual_dat$Dead_tot / actual_dat$tot 
        actual_dat$Discharged <- actual_dat$Discharged_tot / actual_dat$tot 
        actual_dat$WT <- "actual"
        colnames(means) <- colnames(reg_data_no_stay)[-which(colnames(reg_data_no_stay) == "outcome")]
        
        mean_wt_pred <- predict(ml_stay, newdata = means, "probs")
        if (!is.matrix(mean_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- mean_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          mean_wt_pred <- out_df
        }else if(ncol(mean_wt_pred) != 4){
          
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(mean_wt_pred)){
            out_df[,column] <- mean_wt_pred[,column]
            
          }
          
          mean_wt_pred <- out_df
        }
        
        mean_wt_pred <- as.data.frame(mean_wt_pred)
        mean_wt_pred$patient_group <- patient_group[j]
        mean_wt_pred$ICD <- current_icd
        mean_wt_pred$age <- current_age
        mean_wt_pred$WT <- "mean"
        # Or calculate it at the median
        medians <- means
        
        median_wt_pred <- predict(ml_stay, newdata = medians, "probs")
        if (!is.matrix(median_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- median_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          median_wt_pred <- out_df
        }else if(ncol(median_wt_pred) != 4){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(median_wt_pred)){
            out_df[,column] <- median_wt_pred[,column]
            
          }
          
          median_wt_pred <- out_df
        }
        
        median_wt_pred <- as.data.frame(median_wt_pred)
        
        median_wt_pred$patient_group <- patient_group[j]
        median_wt_pred$ICD <- current_icd
        median_wt_pred$age <- current_age
        median_wt_pred$WT <- "median"
        
        toc()
        
        
        
        
        
      }
      
    }
    
    mean_wt_pred$reg_week <- forecast_week_nums
    median_wt_pred$reg_week <- forecast_week_nums
    
    tot_df <- dplyr::bind_rows(mean_wt_pred, median_wt_pred, actual_dat)
    
    graphing_df <- melt(tot_df, id.vars = colnames(tot_df)[5:14])
    graphing_df$line_group <- paste(graphing_df$variable, graphing_df$WT, sep = "-")
    
    graph_plot <- ggplot(data = graphing_df, aes(x = reg_week, y = value, group = line_group)) +
      geom_line(aes(color = variable, linetype = WT)) + theme_bw() +
      ggtitle(paste(patient_group[j],"month fixed and trend times"))
    
    print(graph_plot)
    toc()
    
    whole_df <- dplyr::bind_rows(whole_df, graphing_df)
    
    
  }
  dev.off()  
  
  toc()
  return(whole_df)
  
}


###############################################################################
## Write the functions for clustering #########################################
###############################################################################

elective_regression_cluster <- function(patient_group,hes_data_orig, start_date, forecast_length, forecast_start,
                                                                time_trend = TRUE, month_trend = TRUE, wt_variable = "linear"){
  ## Input ICD hes data for one ICD one age group one patient group (cc/ga)
  ## to be run in parrallel
  tic("Whole process")
  whole_df <- NULL
  coef_df <- NULL
  forecast_seq <- seq(as.Date(forecast_start), by = "week", length.out = 52)
  require(stringr)
  require(lubridate)  
  require(reshape2)
  require(tictoc)
  

  
  for(j in 1:length(patient_group)){
    print(paste("On patient group:",patient_group[j],". Number",j, "of",length(patient_group)))
    tic(paste("Running for patient group:",patient_group[j]))
    
    split_patient_group <- str_split_fixed(patient_group[j], "-",3)
    current_icd <- as.integer(split_patient_group[1])
    current_ward <- split_patient_group[2]
    current_age <- as.integer(split_patient_group[3])
    
    hes_data <- hes_data_orig[hes_data_orig$ICD == current_icd & 
                                hes_data_orig$agegrp_v3 == current_age,]
    
    if(current_ward == "cc"){
      
      tic("Narrowing to CC only")
      hes_data <- hes_data[hes_data$cc == 1,]
      toc()
      ## remove na trans
      tic("Removing NAs")
      na_rows <- which(is.na(hes_data$cc_transitions))
      if(length(na_rows) > 0)
        hes_data <- hes_data[-na_rows,]
      
      na_cc_los <- which(is.na(hes_data$cc_LoS))
      if(length(na_cc_los) > 0)
        hes_data <- hes_data[-na_cc_los,]
      
      na_cc_WT <- which(is.na(hes_data$WaitingTime))
      if(length(na_cc_WT) > 0)
        hes_data <- hes_data[-na_cc_WT,]
      
      
      toc()
      
      ## set up transitions 
      tic("Set up outcome variable")
      hes_data$outcome <- NA
      
      hes_data[hes_data$cc_LoS >= 7 ,"outcome"]<-"CC"
      hes_data[hes_data$cc_LoS < 7 & hes_data$cc_transitions == 3,"outcome"] <- "Dead"
      hes_data[hes_data$cc_LoS < 7 & hes_data$cc_transitions == 2,"outcome"] <- "GA"
      hes_data[hes_data$cc_LoS < 7 & hes_data$cc_transitions == 1,"outcome"] <- "Discharged"
      
      
      
      hes_data$stay <- hes_data$cc_LoS
      hes_data[hes_data$stay > 7,"stay"]<-7
      
      
      
      hes_data$outcome <- factor(hes_data$outcome, levels = c("CC","Dead","GA","Discharged"))
      hes_data$outcome <- relevel(hes_data$outcome, ref = "GA")
      
      
      toc()
      
    }else{
      ## remove na trans
      tic("Removing NAs")
      na_rows <- which(is.na(hes_data$ga_transitions))
      if(length(na_rows) > 0)
        hes_data <- hes_data[-na_rows,]
      
      na_ga_los <- which(is.na(hes_data$GA_LoS))
      if(length(na_ga_los) > 0)
        hes_data <- hes_data[-na_ga_los,]
      
      na_ga_WT <- which(is.na(hes_data$WaitingTime))
      if(length(na_ga_WT) > 0)
        hes_data <- hes_data[-na_ga_WT,]
      
      
      toc()
      
      ## set up transitions 
      tic("Set up outcome variable")
      hes_data$outcome <- NA
      
      hes_data[hes_data$GA_LoS >= 7 ,"outcome"]<-"GA"
      hes_data[hes_data$GA_LoS < 7 & hes_data$ga_transitions == 3,"outcome"] <- "Dead"
      hes_data[hes_data$GA_LoS < 7 & hes_data$ga_transitions == 2,"outcome"] <- "CC"
      hes_data[hes_data$GA_LoS < 7 & hes_data$ga_transitions == 1,"outcome"] <- "Discharged"
      
      
      
      hes_data$stay <- hes_data$GA_LoS
      hes_data[hes_data$stay > 7,"stay"]<-7
      
      
      
      hes_data$outcome <- factor(hes_data$outcome, levels = c("CC","Dead","GA","Discharged"))
      hes_data$outcome <- relevel(hes_data$outcome, ref = "GA")
      
      
      toc()
    }
    ## set up week and year fixed effects
    
    forecast_week_start <- week_num_func(forecast_start,"2009-01-01")
    forecast_week_end <- week_num_func(forecast_seq[forecast_length],"2009-01-01")
    forecast_week_nums <- seq(forecast_week_start, forecast_week_end)
    
    
    
    
    tic("Getting the reg week data")
    hes_data$reg_week <- sapply(hes_data$admidate_MDY,week_num_func,start_week = start_date)
    toc()
    
    
    if(month_trend == TRUE){
      col_names <- NULL
      
      
      
      for(k in 2:length(month.name)){
        col_name <- paste(month.name[k], "fixed_effect",sep = "_")
        hes_data[,col_name] <- 0
        hes_data[hes_data$admidate_MM == k,col_name] <- 1
        col_names <- append(col_names, col_name)
        
      }
      
    }
    ## Set up the outcome variable 
    ## set up the reg data 
    tic("Reg set up")
    
    if(time_trend == TRUE & month_trend == TRUE){
      reg_data_no_stay <- hes_data[,which(colnames(hes_data) %in% c(col_names, "outcome","reg_week","WaitingTime"))]
    }else if(time_trend == TRUE & month_trend == FALSE){
      reg_data_no_stay <- hes_data[,which(colnames(hes_data) %in% c("outcome","reg_week","WaitingTime"))]
    }else if(time_trend == FALSE & month_trend == TRUE){
      reg_data_no_stay <- hes_data[,which(colnames(hes_data) %in% c(col_names, "outcome","WaitingTime"))]
    }else{
      reg_data_no_stay <- hes_data[,which(colnames(hes_data) %in% c("outcome","WaitingTime"))]
    }
    
    if(wt_variable == "squared"){
      reg_data_no_stay$WaitingTime <- reg_data_no_stay$WaitingTime^2 
    }else if(wt_variable == "log"){
      reg_data_no_stay$WaitingTime <- log(reg_data_no_stay$WaitingTime)
    }
      
    
    
    
    actual_prop_dat <- hes_data[,which(colnames(hes_data) %in% c("outcome","reg_week"))]
    actual_prop_dat$one <- 1
    actual_prop_dat_agg <- aggregate(one ~ reg_week, data = actual_prop_dat, FUN = sum)
    actual_prop_dat_outcome <- aggregate(one ~ reg_week + outcome, data = actual_prop_dat, FUN = sum)
    
    
    toc()
    
    
    
    outcomes_in_dat <- plyr::count(reg_data_no_stay$outcome)
    ## Checking if enough outcomes for model else will return one for probs
    
    if(nrow(outcomes_in_dat) < 2){
      out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
      colnames(out_df) <- levels(hes_data$outcome)
      out_df[,outcomes_in_dat[1,1]] <- 1
      mean_wt_pred <- as.data.frame(out_df)
      mean_wt_pred$patient_group <- patient_group[j]
      mean_wt_pred$ICD <- current_icd
      mean_wt_pred$age <- current_age
      mean_wt_pred$WT <- "mean"
      median_wt_pred <- mean_wt_pred
      median_wt_pred$WT <- "median"
      seven_wt_pred <- mean_wt_pred
      seven_wt_pred$WT <- "seven"
      
      actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
      colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
      actual_dat$reg_week <- forecast_week_nums
      actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_agg, by = c("reg_week"="reg_week"))
      colnames(actual_dat)[ncol(actual_dat)] <- "tot"
      actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA",c("reg_week","one")],
                                     by = c("reg_week" = "reg_week"))
      colnames(actual_dat)[ncol(actual_dat)] <- "GA_tot"
      
      actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC",c("reg_week","one")],
                                     by = c("reg_week" = "reg_week"))
      colnames(actual_dat)[ncol(actual_dat)] <- "CC_tot"
      
      actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead",c("reg_week","one")],
                                     by = c("reg_week" = "reg_week"))
      colnames(actual_dat)[ncol(actual_dat)] <- "Dead_tot"
      
      actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged",c("reg_week","one")],
                                     by = c("reg_week" = "reg_week"))
      colnames(actual_dat)[ncol(actual_dat)] <- "Discharged_tot"
      
      actual_dat$GA <- actual_dat$GA_tot / actual_dat$tot 
      actual_dat$CC <- actual_dat$CC_tot / actual_dat$tot 
      actual_dat$Dead <- actual_dat$Dead_tot / actual_dat$tot 
      actual_dat$Discharged <- actual_dat$Discharged_tot / actual_dat$tot 
      actual_dat$WT <- "actual"
      actual_dat$patient_group <- patient_group[j]
      
      current_coef_df <- data.frame(matrix(0, ncol = 5, nrow = 4))
      colnames(current_coef_df) <- c("transition","coef","patient_group","ICD","age")
      current_coef_df$transition <- c("CC","GA","Dead","Discharged")
      current_coef_df$patient_group <- patient_group[j]
      current_coef_df$ICD <- current_icd
      current_coef_df$age <- current_age
      
      
      
    }else{
      # ml depending on WT only
      
      tic("Regreesion model fitting")
      ml_stay <- multinom(outcome ~ ., data = reg_data_no_stay)
      wt_0_dat <- data.frame(matrix(0, ncol = ncol(reg_data_no_stay) - 1, nrow = 1))
      colnames(wt_0_dat) <- colnames(reg_data_no_stay)[-which(colnames(reg_data_no_stay) == "outcome")]
      wt_1_dat <- wt_0_dat
      wt_1_dat$WaitingTime <- 1
      
      wt_0_pred <- predict(ml_stay, newdata = wt_0_dat, "probs")
      wt_1_pred <- predict(ml_stay, newdata = wt_1_dat, "probs")
      
      wt_diff <- wt_1_pred - wt_0_pred
      
      if(length(wt_diff) != 4){
        if(length(wt_diff) == 1){
          current_lev <- ml_stay$lev[1]
          other_lev <- ml_stay$lev[2]
          wt_diff <- c(wt_diff, ((1- wt_1_pred) - (1- wt_0_pred)))
          names(wt_diff) <- c(current_lev, other_lev)
        }
          missing_factors <- which(!(levels(reg_data_no_stay$outcome) %in% names(wt_diff)))
          old_names <- names(wt_diff)
          wt_diff <- append(wt_diff,rep(0, length(missing_factors)))
          names(wt_diff) <- c(old_names, levels(reg_data_no_stay$outcome)[missing_factors])
          
      }
      current_coef_df <- data.frame(matrix(0, ncol = 5, nrow = 4))
      colnames(current_coef_df) <- c("transition","coef","patient_group","ICD","age")
      current_coef_df$transition <- names(wt_diff)
      current_coef_df$coef <- wt_diff
      current_coef_df$patient_group <- patient_group[j]
      current_coef_df$ICD <- current_icd
      current_coef_df$age <- current_age
      
      
      toc()
      
      if(time_trend & month_trend){
        # we can average the probability over WT
        
        
        
        # Or calculate it at the mean WT
        tic("Probability creation")
        means <- data.frame(matrix(data = 0,ncol = ncol(reg_data_no_stay) - 1, nrow = forecast_length))
        
        ## make the actual data df to compare to predictions 
        
        actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
        colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
        actual_dat$reg_week <- forecast_week_nums
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_agg, by = c("reg_week"="reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "tot"
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "GA_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "CC_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Dead_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Discharged_tot"
        
        actual_dat$GA <- actual_dat$GA_tot / actual_dat$tot 
        actual_dat$CC <- actual_dat$CC_tot / actual_dat$tot 
        actual_dat$Dead <- actual_dat$Dead_tot / actual_dat$tot 
        actual_dat$Discharged <- actual_dat$Discharged_tot / actual_dat$tot 
        actual_dat$WT <- "actual"
        actual_dat$patient_group <- patient_group[j]
        colnames(means) <- colnames(reg_data_no_stay)[-which(colnames(reg_data_no_stay) == "outcome")]
        
        means$WaitingTime <- mean(reg_data_no_stay$WaitingTime, na.rm = TRUE)
        means$reg_week <- forecast_week_nums
        month_nums <- month(forecast_seq)
        for(k in 1:length(month_nums)){
          current_month <- month_nums[k]
          if(current_month != 1)
            means[k,current_month + 1] <- 1
          
        }
        
        
        mean_wt_pred <- predict(ml_stay, newdata = means, "probs")
        if (!is.matrix(mean_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- mean_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          mean_wt_pred <- out_df
        }else if(ncol(mean_wt_pred) != 4){
          
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(mean_wt_pred)){
            out_df[,column] <- mean_wt_pred[,column]
            
          }
          
          mean_wt_pred <- out_df
        }
        
        mean_wt_pred <- as.data.frame(mean_wt_pred)
        mean_wt_pred$patient_group <- patient_group[j]
        mean_wt_pred$ICD <- current_icd
        mean_wt_pred$age <- current_age
        mean_wt_pred$WT <- "mean"
        # Or calculate it at the median
        medians <- means
        medians$WaitingTime <- median(reg_data_no_stay$WaitingTime)
        
        median_wt_pred <- predict(ml_stay, newdata = medians, "probs")
        if (!is.matrix(median_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, medians))
          out_df[,as.character(single_val[1,1])] <- median_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          median_wt_pred <- out_df
        }else if(ncol(median_wt_pred) != 4){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(median_wt_pred)){
            out_df[,column] <- median_wt_pred[,column]
            
          }
          
          median_wt_pred <- out_df
        }
        
        median_wt_pred <- as.data.frame(median_wt_pred)
        
        median_wt_pred$patient_group <- patient_group[j]
        median_wt_pred$ICD <- current_icd
        median_wt_pred$age <- current_age
        median_wt_pred$WT <- "median"
        
        ## 7 day average
        
        seven_days <- means
        seven_days$WaitingTime <- 7
        
        seven_wt_pred <- predict(ml_stay, newdata = seven_days, "probs")
        if (!is.matrix(seven_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, seven_days))
          out_df[,as.character(single_val[1,1])] <- seven_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          seven_wt_pred <- out_df
        }else if(ncol(seven_wt_pred) != 4){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(seven_wt_pred)){
            out_df[,column] <- seven_wt_pred[,column]
            
          }
          
          seven_wt_pred <- out_df
        }
        
        seven_wt_pred <- as.data.frame(seven_wt_pred)
        
        seven_wt_pred$patient_group <- patient_group[j]
        seven_wt_pred$ICD <- current_icd
        seven_wt_pred$age <- current_age
        seven_wt_pred$WT <- "seven"
        
        
        toc()
      }else if(time_trend & month_trend == FALSE){
        tic("Probability creation")
        means <- data.frame(matrix(data = 0,ncol = ncol(reg_data_no_stay) - 1, nrow = forecast_length))
        
        ## make the actual data df to compare to predictions 
        
        actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
        colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
        actual_dat$reg_week <- forecast_week_nums
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_agg, by = c("reg_week"="reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "tot"
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "GA_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "CC_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Dead_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Discharged_tot"
        
        actual_dat$GA <- actual_dat$GA_tot / actual_dat$tot 
        actual_dat$CC <- actual_dat$CC_tot / actual_dat$tot 
        actual_dat$Dead <- actual_dat$Dead_tot / actual_dat$tot 
        actual_dat$Discharged <- actual_dat$Discharged_tot / actual_dat$tot 
        actual_dat$WT <- "actual"
        actual_dat$patient_group <- patient_group[j]
        colnames(means) <- colnames(reg_data_no_stay)[-which(colnames(reg_data_no_stay) == "outcome")]
        
        means$reg_week <- forecast_week_nums
        
        mean_wt_pred <- predict(ml_stay, newdata = means, "probs")
        if (!is.matrix(mean_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- mean_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          mean_wt_pred <- out_df
        }else if(ncol(mean_wt_pred) != 4){
          
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(mean_wt_pred)){
            out_df[,column] <- mean_wt_pred[,column]
            
          }
          
          mean_wt_pred <- out_df
        }
        
        mean_wt_pred <- as.data.frame(mean_wt_pred)
        mean_wt_pred$patient_group <- patient_group[j]
        mean_wt_pred$ICD <- current_icd
        mean_wt_pred$age <- current_age
        mean_wt_pred$WT <- "mean"
        # Or calculate it at the median
        medians <- means
        
        median_wt_pred <- predict(ml_stay, newdata = medians, "probs")
        if (!is.matrix(median_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, medians))
          out_df[,as.character(single_val[1,1])] <- median_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          median_wt_pred <- out_df
        }else if(ncol(median_wt_pred) != 4){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(median_wt_pred)){
            out_df[,column] <- median_wt_pred[,column]
            
          }
          
          median_wt_pred <- out_df
        }
        
        median_wt_pred <- as.data.frame(median_wt_pred)
        
        median_wt_pred$patient_group <- patient_group[j]
        median_wt_pred$ICD <- current_icd
        median_wt_pred$age <- current_age
        median_wt_pred$WT <- "median"
        
        ## 7 day average
        
        seven_days <- means
        seven_days$WaitingTime <- 7
        
        seven_wt_pred <- predict(ml_stay, newdata = seven_days, "probs")
        if (!is.matrix(seven_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, seven_days))
          out_df[,as.character(single_val[1,1])] <- seven_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          seven_wt_pred <- out_df
        }else if(ncol(seven_wt_pred) != 4){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(seven_wt_pred)){
            out_df[,column] <- seven_wt_pred[,column]
            
          }
          
          seven_wt_pred <- out_df
        }
        
        seven_wt_pred <- as.data.frame(seven_wt_pred)
        
        seven_wt_pred$patient_group <- patient_group[j]
        seven_wt_pred$ICD <- current_icd
        seven_wt_pred$age <- current_age
        seven_wt_pred$WT <- "seven"
        
        toc()
        
        
        
        
      }else if(time_trend == FALSE & month_trend){
        tic("Probability creation")
        means <- data.frame(matrix(data = 0,ncol = ncol(reg_data_no_stay) - 1, nrow = forecast_length))
        
        ## make the actual data df to compare to predictions 
        
        actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
        colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
        actual_dat$reg_week <- forecast_week_nums
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_agg, by = c("reg_week"="reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "tot"
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "GA_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "CC_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Dead_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Discharged_tot"
        
        actual_dat$GA <- actual_dat$GA_tot / actual_dat$tot 
        actual_dat$CC <- actual_dat$CC_tot / actual_dat$tot 
        actual_dat$Dead <- actual_dat$Dead_tot / actual_dat$tot 
        actual_dat$Discharged <- actual_dat$Discharged_tot / actual_dat$tot 
        actual_dat$WT <- "actual"
        actual_dat$patient_group <- patient_group[j]
        colnames(means) <- colnames(reg_data_no_stay)[-which(colnames(reg_data_no_stay) == "outcome")]
        
        
        month_nums <- month(forecast_seq)
        for(k in 1:length(month_nums)){
          current_month <- month_nums[k]
          if(current_month != 1)
            means[k,current_month + 1] <- 1
          
        }
        
        
        mean_wt_pred <- predict(ml_stay, newdata = means, "probs")
        if (!is.matrix(mean_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- mean_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          mean_wt_pred <- out_df
        }else if(ncol(mean_wt_pred) != 4){
          
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(mean_wt_pred)){
            out_df[,column] <- mean_wt_pred[,column]
            
          }
          
          mean_wt_pred <- out_df
        }
        
        mean_wt_pred <- as.data.frame(mean_wt_pred)
        mean_wt_pred$patient_group <- patient_group[j]
        mean_wt_pred$ICD <- current_icd
        mean_wt_pred$age <- current_age
        mean_wt_pred$WT <- "mean"
        # Or calculate it at the median
        medians <- means
        
        median_wt_pred <- predict(ml_stay, newdata = medians, "probs")
        if (!is.matrix(median_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- median_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          median_wt_pred <- out_df
        }else if(ncol(median_wt_pred) != 4){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(median_wt_pred)){
            out_df[,column] <- median_wt_pred[,column]
            
          }
          
          median_wt_pred <- out_df
        }
        
        median_wt_pred <- as.data.frame(median_wt_pred)
        
        median_wt_pred$patient_group <- patient_group[j]
        median_wt_pred$ICD <- current_icd
        median_wt_pred$age <- current_age
        median_wt_pred$WT <- "median"
        
        ## 7 day average
        
        seven_days <- means
        seven_days$WaitingTime <- 7
        
        seven_wt_pred <- predict(ml_stay, newdata = seven_days, "probs")
        if (!is.matrix(seven_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, seven_days))
          out_df[,as.character(single_val[1,1])] <- seven_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          seven_wt_pred <- out_df
        }else if(ncol(seven_wt_pred) != 4){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(seven_wt_pred)){
            out_df[,column] <- seven_wt_pred[,column]
            
          }
          
          seven_wt_pred <- out_df
        }
        
        seven_wt_pred <- as.data.frame(seven_wt_pred)
        
        seven_wt_pred$patient_group <- patient_group[j]
        seven_wt_pred$ICD <- current_icd
        seven_wt_pred$age <- current_age
        seven_wt_pred$WT <- "seven"
        
        
        toc()
        
        
      }else if( month_trend == FALSE & time_trend == FALSE){
        tic("Probability creation")
        means <- data.frame(matrix(data = 0,ncol = ncol(reg_data_no_stay) - 1, nrow = forecast_length))
        
        ## make the actual data df to compare to predictions 
        
        actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
        colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
        actual_dat$reg_week <- forecast_week_nums
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_agg, by = c("reg_week"="reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "tot"
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "GA_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "CC_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Dead_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Discharged_tot"
        
        actual_dat$GA <- actual_dat$GA_tot / actual_dat$tot 
        actual_dat$CC <- actual_dat$CC_tot / actual_dat$tot 
        actual_dat$Dead <- actual_dat$Dead_tot / actual_dat$tot 
        actual_dat$Discharged <- actual_dat$Discharged_tot / actual_dat$tot 
        actual_dat$WT <- "actual"
        actual_dat$patient_group <- patient_group[j]
        colnames(means) <- colnames(reg_data_no_stay)[-which(colnames(reg_data_no_stay) == "outcome")]
        
        mean_wt_pred <- predict(ml_stay, newdata = means, "probs")
        if (!is.matrix(mean_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- mean_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          mean_wt_pred <- out_df
        }else if(ncol(mean_wt_pred) != 4){
          
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(mean_wt_pred)){
            out_df[,column] <- mean_wt_pred[,column]
            
          }
          
          mean_wt_pred <- out_df
        }
        
        mean_wt_pred <- as.data.frame(mean_wt_pred)
        mean_wt_pred$patient_group <- patient_group[j]
        mean_wt_pred$ICD <- current_icd
        mean_wt_pred$age <- current_age
        mean_wt_pred$WT <- "mean"
        # Or calculate it at the median
        medians <- means
        
        median_wt_pred <- predict(ml_stay, newdata = medians, "probs")
        if (!is.matrix(median_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- median_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          median_wt_pred <- out_df
        }else if(ncol(median_wt_pred) != 4){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(median_wt_pred)){
            out_df[,column] <- median_wt_pred[,column]
            
          }
          
          median_wt_pred <- out_df
        }
        
        median_wt_pred <- as.data.frame(median_wt_pred)
        
        median_wt_pred$patient_group <- patient_group[j]
        median_wt_pred$ICD <- current_icd
        median_wt_pred$age <- current_age
        median_wt_pred$WT <- "median"
        
        ## 7 day average
        
        seven_days <- means
        seven_days$WaitingTime <- 7
        
        seven_wt_pred <- predict(ml_stay, newdata = seven_days, "probs")
        if (!is.matrix(seven_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, seven_days))
          out_df[,as.character(single_val[1,1])] <- seven_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          seven_wt_pred <- out_df
        }else if(ncol(seven_wt_pred) != 4){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(seven_wt_pred)){
            out_df[,column] <- seven_wt_pred[,column]
            
          }
          
          seven_wt_pred <- out_df
        }
        
        seven_wt_pred <- as.data.frame(seven_wt_pred)
        
        seven_wt_pred$patient_group <- patient_group[j]
        seven_wt_pred$ICD <- current_icd
        seven_wt_pred$age <- current_age
        seven_wt_pred$WT <- "seven"
        
        toc()
        
        
        
        
        
      }
      
    }
    
    mean_wt_pred$reg_week <- forecast_week_nums
    median_wt_pred$reg_week <- forecast_week_nums
    seven_wt_pred$reg_week <- forecast_week_nums
    
    tot_df <- dplyr::bind_rows(mean_wt_pred, median_wt_pred, seven_wt_pred, actual_dat)
    
    graphing_df <- melt(tot_df, id.vars = colnames(tot_df)[5:14])
    graphing_df$line_group <- paste(graphing_df$variable, graphing_df$WT, sep = "-")
    
    
    whole_df <- dplyr::bind_rows(whole_df, graphing_df)
    coef_df <- dplyr::bind_rows(coef_df, current_coef_df)
    
  }
  
  
  toc()
  return(list(whole_df, coef_df))
  
}

emergency_regression_cluster <- function(patient_group,hes_data_orig, start_date, forecast_length, forecast_start,
                                         time_trend = TRUE, month_trend = TRUE){
  ## Input ICD hes data for one ICD one age group one patient group (cc/ga)
  ## to be run in parrallel
  tic("Whole process")
  whole_df <- NULL
  forecast_seq <- seq(as.Date(forecast_start), by = "week", length.out = 52)
  require(stringr)
  require(lubridate)  
  require(reshape2)
  require(tictoc)
  
  
  
  for(j in 1:length(patient_group)){
    print(paste("On patient group:",patient_group[j],". Number",j, "of",length(patient_group)))
    tic(paste("Running for patient group:",patient_group[j]))
    
    split_patient_group <- str_split_fixed(patient_group[j], "-",3)
    current_icd <- as.integer(split_patient_group[1])
    current_ward <- split_patient_group[2]
    current_age <- as.integer(split_patient_group[3])
    
    hes_data <- hes_data_orig[hes_data_orig$ICD == current_icd & 
                                hes_data_orig$agegrp_v3 == current_age,]
    if(current_ward == "cc"){
      
      tic("Narrowing to CC only")
      hes_data <- hes_data[hes_data$cc == 1,]
      toc()
      ## remove na trans
      tic("Removing NAs")
      na_rows <- which(is.na(hes_data$cc_transitions))
      if(length(na_rows) > 0)
        hes_data <- hes_data[-na_rows,]
      
      na_ga_los <- which(is.na(hes_data$cc_LoS))
      if(length(na_ga_los) > 0)
        hes_data <- hes_data[-na_ga_los,]
      toc()
      
      
      hes_data_orig
      
      
      ## set up transitions 
      tic("Set up outcome variable")
      hes_data$outcome <- NA
      
      hes_data[hes_data$cc_LoS >= 7 ,"outcome"]<-"CC"
      hes_data[hes_data$cc_LoS < 7 & hes_data$cc_transitions == 3,"outcome"] <- "Dead"
      hes_data[hes_data$cc_LoS < 7 & hes_data$cc_transitions == 2,"outcome"] <- "GA"
      hes_data[hes_data$cc_LoS < 7 & hes_data$cc_transitions == 1,"outcome"] <- "Discharged"
      
      
      
      hes_data$stay <- hes_data$cc_LoS
      hes_data[hes_data$stay > 7,"stay"]<-7
      
      
      
      hes_data$outcome <- factor(hes_data$outcome, levels = c("CC","Dead","GA","Discharged"))
      hes_data$outcome <- relevel(hes_data$outcome, ref = "GA")
      
      
      toc()
      
    }else{
      tic("Narrowing to GA only")
      hes_data <- hes_data[hes_data$cc == 0,]
      toc()
      ## remove na trans
      tic("Removing NAs")
      na_rows <- which(is.na(hes_data$ga_transitions))
      if(length(na_rows) > 0)
        hes_data <- hes_data[-na_rows,]
      
      na_ga_los <- which(is.na(hes_data$GA_LoS))
      if(length(na_ga_los) > 0)
        hes_data <- hes_data[-na_ga_los,]
      toc()
      
      ## set up transitions 
      tic("Set up outcome variable")
      hes_data$outcome <- NA
      
      hes_data[hes_data$GA_LoS >= 7 ,"outcome"]<-"GA"
      hes_data[hes_data$GA_LoS < 7 & hes_data$ga_transitions == 3,"outcome"] <- "Dead"
      hes_data[hes_data$GA_LoS < 7 & hes_data$ga_transitions == 2,"outcome"] <- "CC"
      hes_data[hes_data$GA_LoS < 7 & hes_data$ga_transitions == 1,"outcome"] <- "Discharged"
      
      
      
      hes_data$stay <- hes_data$GA_LoS
      hes_data[hes_data$stay > 7,"stay"]<-7
      
      
      
      hes_data$outcome <- factor(hes_data$outcome, levels = c("CC","Dead","GA","Discharged"))
      hes_data$outcome <- relevel(hes_data$outcome, ref = "GA")
      
      
      toc()
    }
    ## set up week and year fixed effects
    
    forecast_week_start <- week_num_func(forecast_start,"2009-01-01")
    forecast_week_end <- week_num_func(forecast_seq[forecast_length],"2009-01-01")
    forecast_week_nums <- seq(forecast_week_start, forecast_week_end)
    
    
    
    
    tic("Getting the reg week data")
    hes_data$reg_week <- sapply(hes_data$admidate_MDY,week_num_func,start_week = start_date)
    toc()
    
    
    if(month_trend == TRUE){
      col_names <- NULL
      
      
      
      for(k in 2:length(month.name)){
        col_name <- paste(month.name[k], "fixed_effect",sep = "_")
        hes_data[,col_name] <- 0
        hes_data[hes_data$admidate_MM == k,col_name] <- 1
        col_names <- append(col_names, col_name)
        
      }
      
    }
    ## Set up the outcome variable 
    ## set up the reg data 
    tic("Reg set up")
    
    if(time_trend == TRUE & month_trend == TRUE){
      reg_data_no_stay <- hes_data[,which(colnames(hes_data) %in% c(col_names, "outcome","reg_week"))]
    }else if(time_trend == TRUE & month_trend == FALSE){
      reg_data_no_stay <- hes_data[,which(colnames(hes_data) %in% c("outcome","reg_week"))]
    }else if(time_trend == FALSE & month_trend == TRUE){
      reg_data_no_stay <- hes_data[,which(colnames(hes_data) %in% c(col_names, "outcome"))]
    }else{
      reg_data_no_stay <- hes_data[,which(colnames(hes_data) %in% c("outcome"))]
    }
    
    
    
    
    
    actual_prop_dat <- hes_data[,which(colnames(hes_data) %in% c("outcome","reg_week"))]
    actual_prop_dat$one <- 1
    actual_prop_dat_agg <- aggregate(one ~ reg_week, data = actual_prop_dat, FUN = sum)
    actual_prop_dat_outcome <- aggregate(one ~ reg_week + outcome, data = actual_prop_dat, FUN = sum)
    
    
    toc()
    
    
    
    outcomes_in_dat <- plyr::count(reg_data_no_stay$outcome)
    ## Checking if enough outcomes for model else will return one for probs
    
    if(nrow(outcomes_in_dat) < 2){
      out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
      colnames(out_df) <- levels(hes_data$outcome)
      out_df[,outcomes_in_dat[1,1]] <- 1
      mean_wt_pred <- as.data.frame(out_df)
      mean_wt_pred$patient_group <- patient_group[j]
      mean_wt_pred$ICD <- current_icd
      mean_wt_pred$age <- current_age
      mean_wt_pred$WT <- "mean"
      median_wt_pred <- mean_wt_pred
      median_wt_pred$WT <- "median"
      actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
      colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
      actual_dat$reg_week <- forecast_week_nums
      actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_agg, by = c("reg_week"="reg_week"))
      colnames(actual_dat)[ncol(actual_dat)] <- "tot"
      actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA",c("reg_week","one")],
                                     by = c("reg_week" = "reg_week"))
      colnames(actual_dat)[ncol(actual_dat)] <- "GA_tot"
      
      actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC",c("reg_week","one")],
                                     by = c("reg_week" = "reg_week"))
      colnames(actual_dat)[ncol(actual_dat)] <- "CC_tot"
      
      actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead",c("reg_week","one")],
                                     by = c("reg_week" = "reg_week"))
      colnames(actual_dat)[ncol(actual_dat)] <- "Dead_tot"
      
      actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged",c("reg_week","one")],
                                     by = c("reg_week" = "reg_week"))
      colnames(actual_dat)[ncol(actual_dat)] <- "Discharged_tot"
      
      actual_dat$GA <- actual_dat$GA_tot / actual_dat$tot 
      actual_dat$CC <- actual_dat$CC_tot / actual_dat$tot 
      actual_dat$Dead <- actual_dat$Dead_tot / actual_dat$tot 
      actual_dat$Discharged <- actual_dat$Discharged_tot / actual_dat$tot 
      actual_dat$patient_group <- patient_group[j]
      actual_dat$WT <- "actual"
      
      
    }else{
      
      
      
      # ml depending on WT only
      
      if(time_trend | month_trend){
        
        tic("Regreesion model fitting")
        ml_stay <- multinom(outcome ~ ., data = reg_data_no_stay)
        toc()
      }else{
        mean_wt_pred <- data.frame(matrix(data = 0, nrow = forecast_length, ncol = 4))
        colnames(mean_wt_pred) <- levels(hes_data$outcome)
        denominator <- sum(actual_prop_dat_agg$one)
        GA_num <- sum(actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA", "one"])
        cc_num <- sum(actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC", "one"])
        dead_num <- sum(actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead", "one"])
        dis_num <- sum(actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged", "one"])
        
        mean_wt_pred$GA <- GA_num / denominator
        mean_wt_pred$CC <- cc_num / denominator 
        mean_wt_pred$Dead <- dead_num / denominator
        mean_wt_pred$Discharged <- dis_num / denominator
        mean_wt_pred$patient_group <- patient_group[j]
        mean_wt_pred$ICD <- current_icd
        mean_wt_pred$age <- current_age
        mean_wt_pred$WT <- "mean"
        
        median_wt_pred <- mean_wt_pred
        median_wt_pred$WT <- "median"
        
        actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
        colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
        actual_dat$reg_week <- forecast_week_nums
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_agg, by = c("reg_week"="reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "tot"
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "GA_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "CC_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Dead_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Discharged_tot"
        
        actual_dat$GA <- actual_dat$GA_tot / actual_dat$tot 
        actual_dat$CC <- actual_dat$CC_tot / actual_dat$tot 
        actual_dat$Dead <- actual_dat$Dead_tot / actual_dat$tot 
        actual_dat$Discharged <- actual_dat$Discharged_tot / actual_dat$tot 
        actual_dat$WT <- "actual"
        
        
        
      }
      if(time_trend & month_trend){
        # we can average the probability over WT
        
        
        
        # Or calculate it at the mean WT
        tic("Probability creation")
        means <- data.frame(matrix(data = 0,ncol = ncol(reg_data_no_stay) - 1, nrow = forecast_length))
        
        ## make the actual data df to compare to predictions 
        
        actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
        colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
        actual_dat$reg_week <- forecast_week_nums
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_agg, by = c("reg_week"="reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "tot"
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "GA_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "CC_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Dead_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Discharged_tot"
        
        actual_dat$GA <- actual_dat$GA_tot / actual_dat$tot 
        actual_dat$CC <- actual_dat$CC_tot / actual_dat$tot 
        actual_dat$Dead <- actual_dat$Dead_tot / actual_dat$tot 
        actual_dat$Discharged <- actual_dat$Discharged_tot / actual_dat$tot 
        actual_dat$WT <- "actual"
        actual_dat$patient_group <- patient_group[j]
        colnames(means) <- colnames(reg_data_no_stay)[-which(colnames(reg_data_no_stay) == "outcome")]
        
        means$reg_week <- forecast_week_nums
        month_nums <- month(forecast_seq)
        for(k in 1:length(month_nums)){
          current_month <- month_nums[k]
          if(current_month != 1)
            means[k,current_month + 1] <- 1
          
        }
        
        
        mean_wt_pred <- predict(ml_stay, newdata = means, "probs")
        if (!is.matrix(mean_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- mean_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          mean_wt_pred <- out_df
        }else if(ncol(mean_wt_pred) != 4){
          
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(mean_wt_pred)){
            out_df[,column] <- mean_wt_pred[,column]
            
          }
          
          mean_wt_pred <- out_df
        }
        
        mean_wt_pred <- as.data.frame(mean_wt_pred)
        mean_wt_pred$patient_group <- patient_group[j]
        mean_wt_pred$ICD <- current_icd
        mean_wt_pred$age <- current_age
        mean_wt_pred$WT <- "mean"
        # Or calculate it at the median
        medians <- means
        
        median_wt_pred <- predict(ml_stay, newdata = medians, "probs")
        if (!is.matrix(median_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- median_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          median_wt_pred <- out_df
        }else if(ncol(median_wt_pred) != 4){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(median_wt_pred)){
            out_df[,column] <- median_wt_pred[,column]
            
          }
          
          median_wt_pred <- out_df
        }
        
        median_wt_pred <- as.data.frame(median_wt_pred)
        
        median_wt_pred$patient_group <- patient_group[j]
        median_wt_pred$ICD <- current_icd
        median_wt_pred$age <- current_age
        median_wt_pred$WT <- "median"
        
        toc()
      }else if(time_trend & month_trend == FALSE){
        tic("Probability creation")
        means <- data.frame(matrix(data = 0,ncol = ncol(reg_data_no_stay) - 1, nrow = forecast_length))
        
        ## make the actual data df to compare to predictions 
        
        actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
        colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
        actual_dat$reg_week <- forecast_week_nums
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_agg, by = c("reg_week"="reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "tot"
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "GA_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "CC_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Dead_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Discharged_tot"
        
        actual_dat$GA <- actual_dat$GA_tot / actual_dat$tot 
        actual_dat$CC <- actual_dat$CC_tot / actual_dat$tot 
        actual_dat$Dead <- actual_dat$Dead_tot / actual_dat$tot 
        actual_dat$Discharged <- actual_dat$Discharged_tot / actual_dat$tot 
        actual_dat$WT <- "actual"
        actual_dat$patient_group <- patient_group[j]
        colnames(means) <- colnames(reg_data_no_stay)[-which(colnames(reg_data_no_stay) == "outcome")]
        
        means$reg_week <- forecast_week_nums
        
        mean_wt_pred <- predict(ml_stay, newdata = means, "probs")
        if (!is.matrix(mean_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- mean_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          mean_wt_pred <- out_df
        }else if(ncol(mean_wt_pred) != 4){
          
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(mean_wt_pred)){
            out_df[,column] <- mean_wt_pred[,column]
            
          }
          
          mean_wt_pred <- out_df
        }
        
        mean_wt_pred <- as.data.frame(mean_wt_pred)
        mean_wt_pred$patient_group <- patient_group[j]
        mean_wt_pred$ICD <- current_icd
        mean_wt_pred$age <- current_age
        mean_wt_pred$WT <- "mean"
        # Or calculate it at the median
        medians <- means
        
        median_wt_pred <- predict(ml_stay, newdata = medians, "probs")
        if (!is.matrix(median_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- median_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          median_wt_pred <- out_df
        }else if(ncol(median_wt_pred) != 4){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(median_wt_pred)){
            out_df[,column] <- median_wt_pred[,column]
            
          }
          
          median_wt_pred <- out_df
        }
        
        median_wt_pred <- as.data.frame(median_wt_pred)
        
        median_wt_pred$patient_group <- patient_group[j]
        median_wt_pred$ICD <- current_icd
        median_wt_pred$age <- current_age
        median_wt_pred$WT <- "median"
        
        toc()
        
        
        
        
      }else if(time_trend == FALSE & month_trend){
        tic("Probability creation")
        means <- data.frame(matrix(data = 0,ncol = ncol(reg_data_no_stay) - 1, nrow = forecast_length))
        
        ## make the actual data df to compare to predictions 
        
        actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
        colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
        actual_dat$reg_week <- forecast_week_nums
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_agg, by = c("reg_week"="reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "tot"
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "GA",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "GA_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "CC",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "CC_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Dead",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Dead_tot"
        
        actual_dat <- dplyr::left_join(actual_dat, actual_prop_dat_outcome[actual_prop_dat_outcome$outcome == "Discharged",c("reg_week","one")],
                                       by = c("reg_week" = "reg_week"))
        colnames(actual_dat)[ncol(actual_dat)] <- "Discharged_tot"
        
        actual_dat$GA <- actual_dat$GA_tot / actual_dat$tot 
        actual_dat$CC <- actual_dat$CC_tot / actual_dat$tot 
        actual_dat$Dead <- actual_dat$Dead_tot / actual_dat$tot 
        actual_dat$Discharged <- actual_dat$Discharged_tot / actual_dat$tot 
        actual_dat$WT <- "actual"
        actual_dat$patient_group <- patient_group[j]
        colnames(means) <- colnames(reg_data_no_stay)[-which(colnames(reg_data_no_stay) == "outcome")]
        
        
        month_nums <- month(forecast_seq)
        for(k in 1:length(month_nums)){
          current_month <- month_nums[k]
          if(current_month != 1)
            means[k,current_month + 1] <- 1
          
        }
        
        
        mean_wt_pred <- predict(ml_stay, newdata = means, "probs")
        if (!is.matrix(mean_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- mean_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          mean_wt_pred <- out_df
        }else if(ncol(mean_wt_pred) != 4){
          
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(mean_wt_pred)){
            out_df[,column] <- mean_wt_pred[,column]
            
          }
          
          mean_wt_pred <- out_df
        }
        
        mean_wt_pred <- as.data.frame(mean_wt_pred)
        mean_wt_pred$patient_group <- patient_group[j]
        mean_wt_pred$ICD <- current_icd
        mean_wt_pred$age <- current_age
        mean_wt_pred$WT <- "mean"
        # Or calculate it at the median
        medians <- means
        
        median_wt_pred <- predict(ml_stay, newdata = medians, "probs")
        if (!is.matrix(median_wt_pred)){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          single_val <- plyr::count(predict(ml_stay, means))
          out_df[,as.character(single_val[1,1])] <- median_wt_pred
          if(length(ml_stay$lev) > 1){
            other_level <- ml_stay$lev[which(ml_stay$lev != single_val[1,1])]
            out_df[,other_level] <- 1 - out_df[,as.character(single_val[1,1])]
          }
          
          median_wt_pred <- out_df
        }else if(ncol(median_wt_pred) != 4){
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in colnames(median_wt_pred)){
            out_df[,column] <- median_wt_pred[,column]
            
          }
          
          median_wt_pred <- out_df
        }
        
        median_wt_pred <- as.data.frame(median_wt_pred)
        
        median_wt_pred$patient_group <- patient_group[j]
        median_wt_pred$ICD <- current_icd
        median_wt_pred$age <- current_age
        median_wt_pred$WT <- "median"
        
        toc()
        
        
      }
      
    }
    
    mean_wt_pred$reg_week <- forecast_week_nums
    median_wt_pred$reg_week <- forecast_week_nums
    seven
    
    tot_df <- dplyr::bind_rows(mean_wt_pred, median_wt_pred, actual_dat)
    
    graphing_df <- melt(tot_df, id.vars = colnames(tot_df)[5:14])
    graphing_df$line_group <- paste(graphing_df$variable, graphing_df$WT, sep = "-")
    
    toc()
    
    whole_df <- dplyr::bind_rows(whole_df, graphing_df)
    
  }
  
  
  toc()
  return(whole_df)
  
}


regression_cluster_set_up <- function(patient_group, hes_data, forecast_length, forecast_start = "2012-03-05",
                                      start_date = "2009-01-01", time_trend = TRUE, month_trend = TRUE, wt_variable = "linear"){

  ## Takes in data and the patient group as either "elective" or "emergency"
  print("Listing the icds now") 
  tic("Total run through:")
  
  whole_graph_df <- NULL
  coefdf_tot <- NULL
  
  if(patient_group == "elective"){
    tic("Narrowing down the data")
    
    hes_data <- hes_data[hes_data$cohort == 1,]
    hes_data <- hes_data[,c("ICD","agegrp_v3","cc","WaitingTime","cc_LoS",
                            "GA_LoS","ga_transitions","admidate_MDY","admidate_MM","admidate_YYYY",
                            "admidate_week","cc_transitions")]
    
    toc()
    icd_list <- unique(hes_data$ICD)
    icds_to_run <- rep(icd_list, each = 3)
    ages_to_run <- rep(c(1,2,3),length(icd_list))
    patient_groups_to_run <- paste(icds_to_run,"ga",ages_to_run, sep = "-")
    patient_groups_to_run_cc <- paste(icds_to_run,"cc",ages_to_run, sep = "-")
    tot_patient_groups <- c(patient_groups_to_run, patient_groups_to_run_cc)
    
    tic("Setting up elective cluster")
    regression_cluster <- snow::makeCluster(spec = 14, outfile = "./crr_cluster_log.txt")
    
    regression_input <- snow::clusterSplit(regression_cluster, tot_patient_groups)
    toc()
    snow::clusterExport(regression_cluster, "elective_regression_cluster")
    snow::clusterExport(regression_cluster, "week_num_func")
    print(paste("Copying over data for ICD"))
    copy_start <- Sys.time()
    snow::clusterExport(regression_cluster, "hes_data", envir = environment())
    snow::clusterExport(regression_cluster, "forecast_length", envir = environment())
    snow::clusterExport(regression_cluster, "forecast_start", envir = environment())
    snow::clusterExport(regression_cluster, "start_date", envir = environment())
    snow::clusterExport(regression_cluster, "time_trend", envir = environment())
    snow::clusterExport(regression_cluster, "month_trend", envir = environment())
    snow::clusterExport(regression_cluster, "wt_variable", envir = environment())
    copy_end <- Sys.time()
    print(copy_end - copy_start)
    snow::clusterEvalQ(regression_cluster, require(dplyr))
    snow::clusterEvalQ(regression_cluster, require(lubridate))
    snow::clusterEvalQ(regression_cluster, require(stringr))
    snow::clusterEvalQ(regression_cluster, require(tictoc))
    snow::clusterEvalQ(regression_cluster, require(nnet))
    snow::clusterEvalQ(regression_cluster, require(foreign))
    snow::clusterEvalQ(regression_cluster, require(reshape2))
    
    print(paste("Running Elective jobs for ICD"))
    jobs_start <- Sys.time()
    regression_jobs_parallel <- snow::clusterApply(regression_cluster, regression_input,
                                            fun = elective_regression_cluster,
                                            hes_data_orig = hes_data,
                                            forecast_length = forecast_length, forecast_start = forecast_start,
                                            start_date = start_date, time_trend = time_trend,
                                            month_trend = month_trend,
                                            wt_variable = wt_variable)
    jobs_end <- Sys.time()
    print(jobs_end - jobs_start)
    snow::clusterEvalQ(regression_cluster, rm(hes_data))
    snow::clusterEvalQ(regression_cluster, rm(hes_data_orig))
    snow::clusterEvalQ(regression_cluster, gc())
    stopCluster(regression_cluster)
    
    for(k in 1:length(regression_jobs_parallel)){
      whole_df <- regression_jobs_parallel[[k]][[1]]
      coef_df <- regression_jobs_parallel[[k]][[2]]
      
      whole_graph_df <- dplyr::bind_rows(whole_graph_df, whole_df)
      coefdf_tot <- dplyr::bind_rows(coefdf_tot, coef_df)
      
    }
    
    
  }else if(patient_group == "emergency"){
    tic("Narrowing down the data")
    hes_data <- hes_data[hes_data$cohort == 3,]
    
    hes_data <- hes_data[,c("ICD","agegrp_v3","cc","cc_LoS",
                            "GA_LoS","ga_transitions","admidate_MDY","admidate_MM","admidate_YYYY",
                            "admidate_week","cc_transitions")]
    icd_list <- unique(hes_data$ICD)
    icds_to_run <- rep(icd_list, each = 3)
    ages_to_run <- rep(c(1,2,3),length(icd_list))
    patient_groups_to_run <- paste(icds_to_run,"ga",ages_to_run, sep = "-")
    patient_groups_to_run_cc <- paste(icds_to_run,"cc",ages_to_run, sep = "-")
    tot_patient_groups <- c(patient_groups_to_run, patient_groups_to_run_cc)
    icd_15_age_3_rm <- which(tot_patient_groups %in% c("15-ga-3","15-cc-3"))
    tot_patient_groups <- tot_patient_groups[-icd_15_age_3_rm]
    print(length(tot_patient_groups))

    toc()
    tic("Setting up Emergency cluster")
    regression_cluster <- snow::makeCluster(spec = 15, outfile = "./crr_cluster_log.txt")
    
    regression_input <- snow::clusterSplit(regression_cluster, tot_patient_groups)
    toc()
    snow::clusterExport(regression_cluster, "emergency_regression_cluster")
    snow::clusterExport(regression_cluster, "week_num_func")
    print(paste("Copying over data for ICD"))
    copy_start <- Sys.time()
    snow::clusterExport(regression_cluster, "hes_data", envir = environment())
    snow::clusterExport(regression_cluster, "forecast_length", envir = environment())
    snow::clusterExport(regression_cluster, "forecast_start", envir = environment())
    snow::clusterExport(regression_cluster, "start_date", envir = environment())
    snow::clusterExport(regression_cluster, "time_trend", envir = environment())
    snow::clusterExport(regression_cluster, "month_trend", envir = environment())
    copy_end <- Sys.time()
    print(copy_end - copy_start)
    snow::clusterEvalQ(regression_cluster, require(dplyr))
    snow::clusterEvalQ(regression_cluster, require(lubridate))
    snow::clusterEvalQ(regression_cluster, require(stringr))
    snow::clusterEvalQ(regression_cluster, require(tictoc))
    snow::clusterEvalQ(regression_cluster, require(nnet))
    snow::clusterEvalQ(regression_cluster, require(foreign))
    snow::clusterEvalQ(regression_cluster, require(reshape2))
    
    print(paste("Running Emergency jobs for ICD"))
    jobs_start <- Sys.time()
    regression_jobs_parallel <- snow::clusterApply(regression_cluster, regression_input,
                                                   fun = emergency_regression_cluster,
                                                   hes_data_orig = hes_data,
                                                   forecast_length = forecast_length, forecast_start = forecast_start,
                                                   start_date = start_date, time_trend = time_trend,
                                                   month_trend = month_trend)
    jobs_end <- Sys.time()
    print(jobs_end - jobs_start)
    snow::clusterEvalQ(regression_cluster, rm(hes_data))
    snow::clusterEvalQ(regression_cluster, rm(hes_data_orig))
    snow::clusterEvalQ(regression_cluster, gc())
    stopCluster(regression_cluster)
    
    whole_graph_df <- dplyr::bind_rows(regression_jobs_parallel)
    
  }
  
  toc()  
  return(list(whole_graph_df, coefdf_tot))  
  
  
  
}





multi_graph_plotter <- function(reg_res, out_file, trend_types = "month and time", patient_admi = "elective"){
  
  pdf(file = out_file, paper = "A4r", width = 10, height = 7)  
  
  patient_groupings <- unique(reg_res$patient_group)
  
  for(k in 1:length(patient_groupings)){
    current_grouping <- patient_groupings[k]
    current_graph_df <- reg_res[reg_res$patient_group == current_grouping,]
    
    graph_plot <- ggplot(data = current_graph_df, aes(x = reg_week, y = value, group = line_group)) +
      geom_line(aes(color = variable, linetype = WT)) + theme_bw() +
      ggtitle(paste(current_grouping,trend_types, patient_admi))
    
    print(graph_plot)
    
      
    
    
    
  }
  
  dev.off()
  
}





















