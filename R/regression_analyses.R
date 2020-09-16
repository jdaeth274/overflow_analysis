###############################################################################
## regression analyses ########################################################
## elec 2 emerg logit                     #####################################
## emerg transitions multinomial logit    #####################################
## elec transitions multinomial logit     #####################################
###############################################################################

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

covid_regression <- function(covid_data, week_num = c(1,2,3,4,5), half_week = TRUE){
  ## Assume COVID data has the following variables:
  ## GA_LoS - Length of stay in days for G&A patients
  ## cc_LoS - Length of stay in days for CC patients
  ## cc     - Flag if record for patient in CC, 1 = Yes, 0 = No
  ## cc_transitions - Transitions for patients in cc 1 = Discharged straight from cc, 2 = to G&A, 3 = Death, NA for GA patients
  ## ga_transitions - Transitions for patients in GA 1 = Discharged, 2 = to CC, 3 = death, NA for cc patients
  ## agegrp_v3 - what agegrp patient in, 1 = 0-24, 2 = 25-64, 3 = 65+
  
  ## Set up index for keeping track of half days
  
  covid_data$index <- seq(1,nrow(covid_data))
  
  index_df <- NULL
  
  ## Create output table 
  
  tot_out_df <- NULL
  
  for(j in 1:length(week_num)){
    
    out_df <- data.frame(matrix(ncol = 7, nrow = 24))
    colnames(out_df) <- c("a", "p",	"s",	"sbar",	"pi_y",	"coeff","	variance")
    out_df$a <- "N"
    out_df$p <- paste("COVID_",rep(c("AGE1","AGE2","AGE3"),each = 8), sep = "")
    out_df$s <- rep(rep(c("G","C"), each = 4), 3)
    out_df$sbar <- rep(c("H","C","D","G","H","G","D","C"), 3)
    
    
    
    ## create count variable
    
    covid_data$one <- 1
    
    
    ## Create the separate cc and ga datasets and outcome variable by day 7
    
    covid_cc <- covid_data[covid_data$cc == 1,]
    
    covid_ga <- covid_data[-which(is.na(covid_data$ga_transitions)),]
    
    covid_cc$outcome <- "NA"
    skip_cc <- FALSE
    
    if(half_week){
      ## To get the half day vals first look at 3 days 
      
      days_of_stay_variable <- seq(3,length.out = length(week_num), by = 7)
      days_of_stay_variable <- c(-1,days_of_stay_variable)
      lower_end_stay <- days_of_stay_variable[current_week] + 1
      upper_end_stay <- days_of_stay_variable[current_week + 1]
      
      ## set upper day as half day, so 3.5, 10.5, 17.5 ...
      half_day <- upper_end_stay
      
      ## subset only half days 
      half_dayers <- covid_cc[covid_cc$cc_LoS == half_day,]
      if(nrow(half_dayers) >0){
        ## randomly sample half the values (round down if odd number of rows)
        adders <- sample(1:nrow(half_dayers),floor(nrow(half_dayers)/2), replace = F)
        half_row_to_add <- half_dayers[adders,]
        ## These rows to lose will be added to the index df to keep track of
        half_row_to_lose <- half_dayers[-adders,]
        
        patient_group_indy <- rep("cc", nrow(half_row_to_lose))
        indies <- half_row_to_lose$index
        new_indies <- cbind.data.frame(patient_group_indy, indies)
        colnames(new_indies) <- c("patient_group","index")
        
        
      }else{
        
        half_row_to_add <- NULL
        patient_group_indy <- "cc"
        indies <- 0
        new_indies <- cbind.data.frame(patient_group_indy, indies)
        colnames(new_indies) <- c("patient_group","index")
        
      }
      
      
      ## go through the whole days first, set any above whole day cutoff as GA initially
      
      covid_cc[covid_cc$cc_LoS >= upper_end_stay ,"outcome"]<-"CC"
      covid_cc[covid_cc$cc_LoS >= lower_end_stay & covid_cc$cc_LoS < upper_end_stay & covid_cc$cc_transitions == 3,"outcome"] <- "Dead"
      covid_cc[covid_cc$cc_LoS >= lower_end_stay & covid_cc$cc_LoS < upper_end_stay & covid_cc$cc_transitions == 2,"outcome"] <- "GA"
      covid_cc[covid_cc$cc_LoS >= lower_end_stay & covid_cc$cc_LoS < upper_end_stay & covid_cc$cc_transitions == 1,"outcome"] <- "Discharged"
      
      if(!is.null(half_row_to_add)){
        ## now we reallocate the half day values so the 3-3.5 values in
        covid_cc[covid_cc$index %in% half_row_to_add$index & covid_cc$cc_transitions == 3,"outcome"] <- "Dead"
        covid_cc[covid_cc$index %in% half_row_to_add$index & covid_cc$cc_transitions == 2,"outcome"] <- "GA"
        covid_cc[covid_cc$index %in% half_row_to_add$index & covid_cc$cc_transitions == 1,"outcome"] <- "Discharged"
      }
      
      if(!is.null(index_df)){
        ## now we check if previous runs left a half out and then add that in to this, so the 3.5-4 runs 
        patient_df <- index_df[index_df$patient_group == "cc",]
        if(nrow(patient_df) >0){
          if(patient_df$index[1] != 0){
            
            covid_cc[covid_cc$index %in% patient_df$index & covid_cc$cc_transitions == 3,"outcome"] <- "Dead"
            covid_cc[covid_cc$index %in% patient_df$index & covid_cc$cc_transitions == 2,"outcome"] <- "GA"
            covid_cc[covid_cc$index %in% patient_df$index & covid_cc$cc_transitions == 1,"outcome"] <- "Discharged"
            
            index_df <- index_df[-which(index_df$patient_group == "cc"),]
            index_df <- dplyr::bind_rows(index_df, new_indies)
          }else{
            index_df <- index_df[-which(index_df$patient_group == "cc"),]
            index_df <- dplyr::bind_rows(index_df, new_indies)
          }
          
          
          
          
        }else{
          
          index_df <- dplyr::bind_rows(index_df, new_indies)
        }
        
        
      }else{
        
        index_df <- dplyr::bind_rows(index_df, new_indies)
      }
    }else{
      seven_seq <- seq(0, (7*week_num[length(week_num)]),7)
      lower_day_bound <- seven_seq[j]
      upper_day_bound <- seven_seq[j + 1]
      
      
      
      covid_cc[covid_cc$cc_LoS >= upper_day_bound ,"outcome"]<-"CC"
      covid_cc[covid_cc$cc_LoS >= lower_day_bound & covid_cc$cc_LoS < upper_day_bound & covid_cc$cc_transitions == 3,"outcome"] <- "Dead"
      covid_cc[covid_cc$cc_LoS >= lower_day_bound & covid_cc$cc_LoS < upper_day_bound & covid_cc$cc_transitions == 2,"outcome"] <- "GA"
      covid_cc[covid_cc$cc_LoS >= lower_day_bound & covid_cc$cc_LoS < upper_day_bound & covid_cc$cc_transitions == 1,"outcome"] <- "Discharged"
      
    }
    
    
    ## Check there's any data to get proportions from 
    
    
    if(length(which(covid_cc$outcome == "NA")) > 0)
      covid_cc <- covid_cc[-which(covid_cc$outcome == "NA"),]
    
    
    ## If there is no data then we'll set all the probs to 0
    
    if(nrow(covid_cc) == 0){
      
      skip_cc <- TRUE
      out_df[c(5:8,13:16,21:24),"coeff"] <- 0
      out_df$week <- j
      
    }
    
    
    skip_ga <- FALSE
    
    covid_ga$outcome <- "NA"
    if(half_week){
      ## To get the half day vals first look at 3 days 
      
      days_of_stay_variable <- seq(3,length.out = length(week_num), by = 7)
      days_of_stay_variable <- c(-1,days_of_stay_variable)
      lower_end_stay <- days_of_stay_variable[current_week] + 1
      upper_end_stay <- days_of_stay_variable[current_week + 1]
      
      ## set upper day as half day, so 3.5, 10.5, 17.5 ...
      half_day <- upper_end_stay
      
      ## subset only half days 
      half_dayers <- covid_ga[covid_ga$cc_LoS == half_day,]
      if(nrow(half_dayers) >0){
        ## randomly sample half the values (round down if odd number of rows)
        adders <- sample(1:nrow(half_dayers),floor(nrow(half_dayers)/2), replace = F)
        half_row_to_add <- half_dayers[adders,]
        ## These rows to lose will be added to the index df to keep track of
        half_row_to_lose <- half_dayers[-adders,]
        
        patient_group_indy <- rep("GA", nrow(half_row_to_lose))
        indies <- half_row_to_lose$index
        new_indies <- cbind.data.frame(patient_group_indy, indies)
        colnames(new_indies) <- c("patient_group","index")
        
        
      }else{
        
        half_row_to_add <- NULL
        patient_group_indy <- "GA"
        indies <- 0
        new_indies <- cbind.data.frame(patient_group_indy, indies)
        colnames(new_indies) <- c("patient_group","index")
        
      }
      
      
      ## go through the whole days first, set any above whole day cutoff as GA initially
      
      covid_ga[covid_ga$cc_LoS >= upper_end_stay ,"outcome"]<-"GA"
      covid_ga[covid_ga$cc_LoS >= lower_end_stay & covid_ga$cc_LoS < upper_end_stay & covid_ga$cc_transitions == 3,"outcome"] <- "Dead"
      covid_ga[covid_ga$cc_LoS >= lower_end_stay & covid_ga$cc_LoS < upper_end_stay & covid_ga$cc_transitions == 2,"outcome"] <- "CC"
      covid_ga[covid_ga$cc_LoS >= lower_end_stay & covid_ga$cc_LoS < upper_end_stay & covid_ga$cc_transitions == 1,"outcome"] <- "Discharged"
      
      if(!is.null(half_row_to_add)){
        ## now we reallocate the half day values so the 3-3.5 values in
        covid_ga[covid_ga$index %in% half_row_to_add$index & covid_ga$cc_transitions == 3,"outcome"] <- "Dead"
        covid_ga[covid_ga$index %in% half_row_to_add$index & covid_ga$cc_transitions == 2,"outcome"] <- "CC"
        covid_ga[covid_ga$index %in% half_row_to_add$index & covid_ga$cc_transitions == 1,"outcome"] <- "Discharged"
      }
      
      if(!is.null(index_df)){
        ## now we check if previous runs left a half out and then add that in to this, so the 3.5-4 runs 
        patient_df <- index_df[index_df$patient_group == "GA",]
        if(nrow(patient_df) >0){
          if(patient_df$index[1] != 0){
            
            covid_ga[covid_ga$index %in% patient_df$index & covid_ga$cc_transitions == 3,"outcome"] <- "Dead"
            covid_ga[covid_ga$index %in% patient_df$index & covid_ga$cc_transitions == 2,"outcome"] <- "CC"
            covid_ga[covid_ga$index %in% patient_df$index & covid_ga$cc_transitions == 1,"outcome"] <- "Discharged"
            
            index_df <- index_df[-which(index_df$patient_group == "GA"),]
            index_df <- dplyr::bind_rows(index_df, new_indies)
          }else{
            index_df <- index_df[-which(index_df$patient_group == "GA"),]
            index_df <- dplyr::bind_rows(index_df, new_indies)
          }
          
          
          
          
        }else{
          
          index_df <- dplyr::bind_rows(index_df, new_indies)
        }
        
        
      }else{
        
        index_df <- dplyr::bind_rows(index_df, new_indies)
      }
    }else{
      seven_seq <- seq(0, (7*week_num[length(week_num)]),7)
      lower_day_bound <- seven_seq[j]
      upper_day_bound <- seven_seq[j + 1]
      
      
      
      covid_ga[covid_ga$cc_LoS >= upper_day_bound ,"outcome"]<-"GA"
      covid_ga[covid_ga$cc_LoS >= lower_day_bound & covid_ga$cc_LoS < upper_day_bound & covid_ga$cc_transitions == 3,"outcome"] <- "Dead"
      covid_ga[covid_ga$cc_LoS >= lower_day_bound & covid_ga$cc_LoS < upper_day_bound & covid_ga$cc_transitions == 2,"outcome"] <- "CC"
      covid_ga[covid_ga$cc_LoS >= lower_day_bound & covid_ga$cc_LoS < upper_day_bound & covid_ga$cc_transitions == 1,"outcome"] <- "Discharged"
      
    }
    if(length(which(covid_ga$outcome == "NA")) > 0)
      covid_ga <- covid_ga[-which(covid_ga$outcome == "NA"),]
    
    
    ## If there is no data then we'll set all the probs to 0
    
    if(nrow(covid_ga) == 0){
      
      skip_ga <- TRUE
      out_df[c(1:4,9:12,17:20),"coeff"] <- 0
      out_df$week <- j
      
      skip_ga <- TRUE
      
    }
    
    if(skip_ga == FALSE){
    
      covid_ga$death <- ifelse(covid_ga$outcome == "Dead" , 1,0)
      covid_ga$discharges <- ifelse(covid_ga$outcome == "Discharged", 1, 0)
      covid_ga$ward_switch <- ifelse(covid_ga$outcome == "CC", 1, 0)
      covid_ga$remain <- ifelse(covid_ga$outcome == "GA", 1, 0)
      
      
      ## GA transitions 
      
      tot_ages <- aggregate(one ~ agegrp_v3, covid_ga, sum)$one
      
      if(length(tot_ages) != 3){
        missing_age <- c(1,2,3)[which(!(c(1,2,3) %in% covid_ga$agegrp_v3))]
        present_age <- c(1,2,3)[which((c(1,2,3) %in% covid_ga$agegrp_v3))]
        ## get index of missing age to replace with 0 
        
        tot_ages_agg <- tot_ages
        tot_ages <- rep(0,3)
        tot_ages[present_age] <- tot_ages_agg
        
        ga_deaths <- rep(0,3)
        ga_deaths_agg <- aggregate(death ~ agegrp_v3, covid_ga, sum)$death
        ga_deaths[present_age] <- ga_deaths_agg
        ga_discharges <- rep(0,3)
        ga_discharges_agg <- aggregate(discharges ~ agegrp_v3, covid_ga, sum)$discharges
        ga_discharges[present_age] <- ga_discharges_agg
        ga_ga <- rep(0,3)
        ga_ga_agg <- aggregate(remain ~ agegrp_v3, covid_ga, sum)$remain
        ga_ga[present_age] <- ga_ga_agg
        ga_cc <- rep(0,3)
        ga_cc_agg <- aggregate(ward_switch ~ agegrp_v3, covid_ga, sum)$ward_switch
        ga_cc[present_age] <- ga_cc_agg
        
        
        ga_dat <- NULL
        for(k in 1:3){
          if(k %in% missing_age){
            age_dat <- rep(0, 4)
          }else{
            age_dat <- c(ga_discharges[k] / tot_ages[k], ga_cc[k] / tot_ages[k], 
                         ga_deaths[k] / tot_ages[k], ga_ga[k] / tot_ages[k])
          }
          
          ga_dat <- append(ga_dat, age_dat)
        }
        
        
      }else{
        
        
        ga_deaths <- aggregate(death ~ agegrp_v3, covid_ga, sum)$death
        ga_discharges <- aggregate(discharges ~ agegrp_v3, covid_ga, sum)$discharges
        ga_ga <- aggregate(remain ~ agegrp_v3, covid_ga, sum)$remain
        ga_cc <- aggregate(ward_switch ~ agegrp_v3, covid_ga, sum)$ward_switch
        
        age_1_dat <- c(ga_discharges[1] / tot_ages[1], ga_cc[1] / tot_ages[1], 
                       ga_deaths[1] / tot_ages[1], ga_ga[1] / tot_ages[1])
        age_2_dat <- c(ga_discharges[2] / tot_ages[2], ga_cc[2] / tot_ages[2], 
                       ga_deaths[2] / tot_ages[2], ga_ga[2] / tot_ages[2])
        age_3_dat <- c(ga_discharges[3] / tot_ages[3], ga_cc[3] / tot_ages[3], 
                       ga_deaths[3] / tot_ages[3], ga_ga[3] / tot_ages[3])
        ga_dat <- c(age_1_dat, age_2_dat, age_3_dat)
      }
      
      out_df[c(1:4,9:12,17:20),"coeff"] <- ga_dat
    }
    ## CC transitions
    
    if(skip_cc == FALSE){
      
      covid_cc$death <- ifelse(covid_cc$outcome == "Dead" , 1,0)
      covid_cc$discharges <- ifelse(covid_cc$outcome == "Discharged", 1, 0)
      covid_cc$ward_switch <- ifelse(covid_cc$outcome == "GA", 1, 0)
      covid_cc$remain <- ifelse(covid_cc$outcome == "CC", 1, 0)
      
      
      
      tot_ages <- aggregate(one ~ agegrp_v3, covid_cc, sum)$one
      
      ## Check if any missing ages 
      
      if(length(tot_ages) != 3){
        missing_age <- c(1,2,3)[which(!(c(1,2,3) %in% covid_cc$agegrp_v3))]
        present_age <- c(1,2,3)[which((c(1,2,3) %in% covid_cc$agegrp_v3))]
        ## get index of missing age to replace with 0 
        
        tot_ages_agg <- tot_ages
        tot_ages <- rep(0,3)
        tot_ages[present_age] <- tot_ages_agg
        
        cc_deaths <- rep(0,3)
        cc_deaths_agg <- aggregate(death ~ agegrp_v3, covid_cc, sum)$death
        cc_deaths[present_age] <- cc_deaths_agg
        cc_discharges <- rep(0,3)
        cc_discharges_agg <- aggregate(discharges ~ agegrp_v3, covid_cc, sum)$discharges
        cc_discharges[present_age] <- cc_discharges_agg
        cc_cc <- rep(0,3)
        cc_cc_agg <- aggregate(remain ~ agegrp_v3, covid_cc, sum)$remain
        cc_cc[present_age] <- cc_cc_agg
        cc_ga <- rep(0,3)
        cc_ga_agg <- aggregate(ward_switch ~ agegrp_v3, covid_cc, sum)$ward_switch
        cc_ga[present_age] <- cc_ga_agg
        
        
        cc_dat <- NULL
        for(k in 1:3){
          if(k %in% missing_age){
            age_dat <- rep(0, 4)
          }else{
            age_dat <- c(cc_discharges[k] / tot_ages[k], cc_ga[k] / tot_ages[k], 
                         cc_deaths[k] / tot_ages[k], cc_ga[k] / tot_ages[k])
          }
          
          cc_dat <- append(cc_dat, age_dat)
        }
        
                            
      }else{
      
      
        cc_deaths <- aggregate(death ~ agegrp_v3, covid_cc, sum)$death
        cc_discharges <- aggregate(discharges ~ agegrp_v3, covid_cc, sum)$discharges
        cc_cc <- aggregate(remain ~ agegrp_v3, covid_cc, sum)$remain
        cc_ga <- aggregate(ward_switch ~ agegrp_v3, covid_cc, sum)$ward_switch
        
        age_1_dat <- c(cc_discharges[1] / tot_ages[1], cc_ga[1] / tot_ages[1], 
                       cc_deaths[1] / tot_ages[1], cc_ga[1] / tot_ages[1])
        age_2_dat <- c(cc_discharges[2] / tot_ages[2], cc_ga[2] / tot_ages[2], 
                       cc_deaths[2] / tot_ages[2], cc_ga[2] / tot_ages[2])
        age_3_dat <- c(cc_discharges[3] / tot_ages[3], cc_ga[3] / tot_ages[3], 
                       cc_deaths[3] / tot_ages[3], cc_ga[3] / tot_ages[3])
        cc_dat <- c(age_1_dat, age_2_dat, age_3_dat)
      }
      
      
      
      out_df[c(5:8,13:16,21:24),"coeff"] <- cc_dat
    }
    
    tot_out_df <- dplyr::bind_rows(tot_out_df, out_df)
  }
  
  
  return(tot_out_df)
  
  
}







###############################################################################
## Write the functions for clustering #########################################
###############################################################################

elective_regression_cluster <- function(patient_group,hes_data_orig, start_date, forecast_length, forecast_start,
                                                                time_trend = TRUE, month_trend = TRUE, wt_variable = "linear",
                                        week_num, half_week){
  ## Input ICD hes data for one ICD one age group one patient group (cc/ga)
  ## to be run in parrallel
  tic("Whole process")
  tot_whole_df <- NULL
  tot_coef_df <- NULL
  forecast_seq <- seq(as.Date(forecast_start), by = "week", length.out = 52)
  index_df <- NULL
  hes_data_orig$index <- seq(1,nrow(hes_data_orig),1)
  require(stringr)
  require(lubridate)  
  require(reshape2)
  require(tictoc)
    
  for(current_week in week_num){
    whole_df <- NULL
    coef_df <- NULL
    if(half_week){
      append_week <- (current_week + (current_week - 1)) /2 
    }else{
      append_week <- current_week
    }
    
    
    for(j in 1:length(patient_group)){

      print(paste("On patient group:",patient_group[j],". Number",j, "of",length(patient_group)))
      tic(paste("Running for patient group:",patient_group[j]))
      
          
      split_patient_group <- str_split_fixed(patient_group[j], "-",3)
      current_icd <- as.integer(split_patient_group[1])
      current_ward <- split_patient_group[2]
      current_age <- as.integer(split_patient_group[3])
      
      hes_data <- hes_data_orig[hes_data_orig$ICD == current_icd & 
                                  hes_data_orig$agegrp_v3 == current_age,]
      
      forecast_week_start <- week_num_func(forecast_start, start_date)
      forecast_week_end <- week_num_func(forecast_seq[length(forecast_seq)], start_date)
      forecast_week_nums <- seq(forecast_week_start, forecast_week_end)
      
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
        hes_data$outcome <- "NA"
        
        if(half_week){
          ## To get the half day vals first look at 3 days 
                
          
          days_of_stay_variable <- seq(3,length.out = length(week_num), by = 7)
          days_of_stay_variable <- c(-1,days_of_stay_variable)
          lower_end_stay <- days_of_stay_variable[current_week] + 1
          upper_end_stay <- days_of_stay_variable[current_week + 1]
          half_day <- upper_end_stay
          
          half_dayers <- hes_data[hes_data$cc_LoS == half_day,]
          if(nrow(half_dayers) >0){
            adders <- sample(1:nrow(half_dayers),floor(nrow(half_dayers)/2), replace = F)
            half_row_to_add <- half_dayers[adders,]
            half_row_to_lose <- half_dayers[-adders,]
            
            patient_group_indy <- rep(patient_group[j], nrow(half_row_to_lose))
            indies <- half_row_to_lose$index
            new_indies <- cbind.data.frame(patient_group_indy, indies)
            colnames(new_indies) <- c("patient_group","index")
            
            
          }else{
            
            half_row_to_add <- NULL
            patient_group_indy <- patient_group[j]
            indies <- 0
            new_indies <- cbind.data.frame(patient_group_indy, indies)
            colnames(new_indies) <- c("patient_group","index")
            
          }
          
          
          ## go through the whole days first, set any above whole day cutoff as GA initially
          
          hes_data[hes_data$cc_LoS >= upper_end_stay ,"outcome"]<-"CC"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 3,"outcome"] <- "Dead"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 2,"outcome"] <- "GA"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 1,"outcome"] <- "Discharged"
          
          if(!is.null(half_row_to_add)){
            ## now we reallocate the half day values so the 3-3.5 values in
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$cc_transitions == 3,"outcome"] <- "Dead"
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$cc_transitions == 2,"outcome"] <- "GA"
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$cc_transitions == 1,"outcome"] <- "Discharged"
          }
          
          if(!is.null(index_df)){
            ## now we check if previous runs left a half out and then add that in to this so the 3.5-4 runs 
            patient_df <- index_df[index_df$patient_group == patient_group[j],]
            if(nrow(patient_df) >0){
              if(patient_df$index[1] != 0){
                
        
                hes_data[hes_data$index %in% patient_df$index & hes_data$cc_transitions == 3,"outcome"] <- "Dead"
                hes_data[hes_data$index %in% patient_df$index & hes_data$cc_transitions == 2,"outcome"] <- "GA"
                hes_data[hes_data$index %in% patient_df$index & hes_data$cc_transitions == 1,"outcome"] <- "Discharged"
                
                index_df <- index_df[-which(index_df$patient_group == patient_group[j]),]
                index_df <- dplyr::bind_rows(index_df, new_indies)
              }else{
                index_df <- index_df[-which(index_df$patient_group == patient_group[j]),]
                index_df <- dplyr::bind_rows(index_df, new_indies)
              }
              
              
              
              
            }else{
              
              index_df <- dplyr::bind_rows(index_df, new_indies)
            }
            
            
          }else{
            
            index_df <- dplyr::bind_rows(index_df, new_indies)
          }
        }else{
          days_of_stay_variable <- seq(0, 42, 7)
          lower_end_stay <- days_of_stay_variable[current_week]
          upper_end_stay <- days_of_stay_variable[current_week + 1]
          
          
          hes_data[hes_data$cc_LoS >= upper_end_stay ,"outcome"]<-"CC"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 3,"outcome"] <- "Dead"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 2,"outcome"] <- "GA"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 1,"outcome"] <- "Discharged"
        }
        if(length(which(is.na(hes_data$outcome))) > 0)
          hes_data <- hes_data[-which(is.na(hes_data$outcome)),]
        
        
        if(nrow(hes_data) == 0){
          actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
          colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
          actual_dat$reg_week <- forecast_week_nums
          actual_dat[,2:5] <- 0
          actual_dat$WT <- "actual"
          actual_dat$patient_group <- patient_group[j]
          actual_dat$ICD <- current_icd
          actual_dat$age <- current_age
          
          current_coef_df <- data.frame(matrix(0, ncol = 5, nrow = 4))
          colnames(current_coef_df) <- c("transition","coef","patient_group","ICD","age")
          current_coef_df$transition <- c("CC","GA","Dead","Discharged")
          current_coef_df$patient_group <- patient_group[j]
          current_coef_df$ICD <- current_icd
          current_coef_df$age <- current_age
          
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- c("GA","Dead","CC","Discharged")
          mean_wt_pred <- as.data.frame(out_df)
          mean_wt_pred$patient_group <- patient_group[j]
          mean_wt_pred$ICD <- current_icd
          mean_wt_pred$age <- current_age
          mean_wt_pred$WT <- "mean"
          median_wt_pred <- mean_wt_pred
          median_wt_pred$WT <- "median"
          seven_wt_pred <- mean_wt_pred
          seven_wt_pred$WT <- "seven"
          mean_wt_pred$reg_week <- forecast_week_nums
          median_wt_pred$reg_week <- forecast_week_nums
          seven_wt_pred$reg_week <- forecast_week_nums
          
          tot_df <- dplyr::bind_rows(mean_wt_pred, median_wt_pred, seven_wt_pred, actual_dat)
          
          graphing_df <- melt(tot_df, id.vars = colnames(tot_df)[5:9])
          graphing_df$line_group <- paste(graphing_df$variable, graphing_df$WT, sep = "-")
          
          
          whole_df <- dplyr::bind_rows(whole_df, graphing_df)
          coef_df <- dplyr::bind_rows(coef_df, current_coef_df)
          
          next()
          
        }
        
        
        
        
        
        hes_data$outcome <- factor(hes_data$outcome, levels = c("CC","Dead","GA","Discharged"))
        hes_data$outcome <- relevel(hes_data$outcome, ref = "GA")
        
        
        toc()
        
      }else{
        ## remove na trans
        hes_data <- hes_data[hes_data$cc_start_flg == 0,]
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
        hes_data$outcome <- "NA"
        
        if(half_week){
          ## To get the half day vals first look at 3 days 
          
          
          days_of_stay_variable <- seq(3,length.out = length(week_num), by = 7)
          days_of_stay_variable <- c(-1,days_of_stay_variable)
          lower_end_stay <- days_of_stay_variable[current_week] + 1
          upper_end_stay <- days_of_stay_variable[current_week + 1]
          half_day <- upper_end_stay
          
          half_dayers <- hes_data[hes_data$GA_LoS == half_day,]
          if(nrow(half_dayers) >0){
            adders <- sample(1:nrow(half_dayers),floor(nrow(half_dayers)/2), replace = F)
            half_row_to_add <- half_dayers[adders,]
            half_row_to_lose <- half_dayers[-adders,]
            
            patient_group_indy <- rep(patient_group[j], nrow(half_row_to_lose))
            indies <- half_row_to_lose$index
            new_indies <- cbind.data.frame(patient_group_indy, indies)
            colnames(new_indies) <- c("patient_group","index")
            
            
          }else{
            
            half_row_to_add <- NULL
            patient_group_indy <- patient_group[j]
            indies <- 0
            new_indies <- cbind.data.frame(patient_group_indy, indies)
            colnames(new_indies) <- c("patient_group","index")
            
          }
          
          
          ## go through the whole days first, set any above whole day cutoff as GA initially
          
          hes_data[hes_data$GA_LoS >= upper_end_stay ,"outcome"]<-"GA"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 3,"outcome"] <- "Dead"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 2,"outcome"] <- "CC"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 1,"outcome"] <- "Discharged"
          
          if(!is.null(half_row_to_add)){
            ## now we reallocate the half day values so the 3-3.5 values in
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$ga_transitions == 3,"outcome"] <- "Dead"
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$ga_transitions == 2,"outcome"] <- "CC"
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$ga_transitions == 1,"outcome"] <- "Discharged"
          }
          
          if(!is.null(index_df)){
            ## now we check if previous runs left a half out and then add that in to this so the 3.5-4 runs 
            patient_df <- index_df[index_df$patient_group == patient_group[j],]
            if(nrow(patient_df) >0){
              if(patient_df$index[1] != 0){
                hes_data[hes_data$index %in% patient_df$index & hes_data$ga_transitions == 3,"outcome"] <- "Dead"
                hes_data[hes_data$index %in% patient_df$index & hes_data$ga_transitions == 2,"outcome"] <- "CC"
                hes_data[hes_data$index %in% patient_df$index & hes_data$ga_transitions == 1,"outcome"] <- "Discharged"
                
                index_df <- index_df[-which(index_df$patient_group == patient_group[j]),]
                index_df <- dplyr::bind_rows(index_df, new_indies)
              }else{
                index_df <- index_df[-which(index_df$patient_group == patient_group[j]),]
                index_df <- dplyr::bind_rows(index_df, new_indies)
              }
              
              
              
              
            }else{
              
              index_df <- dplyr::bind_rows(index_df, new_indies)
            }
            
            
          }else{
            
            index_df <- dplyr::bind_rows(index_df, new_indies)
          }
          
          
          
          
        }else{
          
          
          days_of_stay_variable <- seq(0, 42, 7)
          lower_end_stay <- days_of_stay_variable[current_week]
          upper_end_stay <- days_of_stay_variable[current_week + 1]
          
          
          hes_data[hes_data$GA_LoS >= upper_end_stay ,"outcome"]<-"GA"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 3,"outcome"] <- "Dead"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 2,"outcome"] <- "CC"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 1,"outcome"] <- "Discharged"
        }
        if(length(which(is.na(hes_data$outcome))) > 0)
          hes_data <- hes_data[-which(is.na(hes_data$outcome)),]
        
        if(nrow(hes_data) == 0){
          actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
          colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
          actual_dat$reg_week <- forecast_week_nums
          actual_dat[,2:5] <- 0
          actual_dat$WT <- "actual"
          actual_dat$patient_group <- patient_group[j]
          actual_dat$ICD <- current_icd
          actual_dat$age <- current_age
          
          
          current_coef_df <- data.frame(matrix(0, ncol = 5, nrow = 4))
          colnames(current_coef_df) <- c("transition","coef","patient_group","ICD","age")
          current_coef_df$transition <- c("CC","GA","Dead","Discharged")
          current_coef_df$patient_group <- patient_group[j]
          current_coef_df$ICD <- current_icd
          current_coef_df$age <- current_age
          
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- c("GA","Dead","CC","Discharged")
          mean_wt_pred <- as.data.frame(out_df)
          mean_wt_pred$patient_group <- patient_group[j]
          mean_wt_pred$ICD <- current_icd
          mean_wt_pred$age <- current_age
          mean_wt_pred$WT <- "mean"
          median_wt_pred <- mean_wt_pred
          median_wt_pred$WT <- "median"
          seven_wt_pred <- mean_wt_pred
          seven_wt_pred$WT <- "seven"
          mean_wt_pred$reg_week <- forecast_week_nums
          median_wt_pred$reg_week <- forecast_week_nums
          seven_wt_pred$reg_week <- forecast_week_nums
          
          tot_df <- dplyr::bind_rows(mean_wt_pred, median_wt_pred, seven_wt_pred, actual_dat)
          
          graphing_df <- melt(tot_df, id.vars = colnames(tot_df)[5:9])
          graphing_df$line_group <- paste(graphing_df$variable, graphing_df$WT, sep = "-")
          
          
          whole_df <- dplyr::bind_rows(whole_df, graphing_df)
          coef_df <- dplyr::bind_rows(coef_df, current_coef_df)
          
          next()
            
        }
        
        
        
        
        
        
        hes_data$outcome <- factor(hes_data$outcome, levels = c("CC","Dead","GA","Discharged"))
        hes_data$outcome <- relevel(hes_data$outcome, ref = "GA")
        
        
        toc()
      }
      ## set up week and year fixed effects
      
      forecast_week_start <- week_num_func(forecast_start,start_date)
      forecast_week_end <- week_num_func(forecast_seq[forecast_length],start_date)
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
      
      transformation <- function(x){x}
      if(wt_variable == "squared"){
        reg_data_no_stay$WaitingTime <- reg_data_no_stay$WaitingTime^2 
        transformation <- function(x){x^2}
      }else if(wt_variable == "log"){
        reg_data_no_stay$WaitingTime <- log(reg_data_no_stay$WaitingTime)
        transformation <- function(x){log(x)}
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
        actual_dat$ICD <- current_icd
        actual_dat$age <- current_age
        
        
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
          actual_dat$ICD <- current_icd
          actual_dat$age <- current_age
          
          
          means$WaitingTime <- transformation(mean(hes_data$WaitingTime, na.rm = TRUE))
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
          medians$WaitingTime <- transformation(median(hes_data$WaitingTime))
          
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
          seven_days$WaitingTime <- transformation(7)
          
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
          actual_dat$ICD <- current_icd
          actual_dat$age <- current_age
          
          colnames(means) <- colnames(reg_data_no_stay)[-which(colnames(reg_data_no_stay) == "outcome")]
          
          means$reg_week <- forecast_week_nums
          means$WaitingTime <- transformation(mean(hes_data$WaitingTime, na.rm = TRUE))
          
          
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
          medians$WaitingTime <- transformation(median(hes_data$WaitingTime, na.rm = TRUE))
          
          
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
          seven_days$WaitingTime <- transformation(7)
          
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
          
          colnames(means) <- colnames(reg_data_no_stay)[-which(colnames(reg_data_no_stay) == "outcome")]
          
          means$WaitingTime <- transformation(mean(hes_data$WaitingTime, na.rm = TRUE))
          
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
          medians$WaitingTime <- transformation(median(hes_data$WaitingTime, na.rm = TRUE))
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
          seven_days$WaitingTime <- transformation(7)
          
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
          actual_dat$ICD <- current_icd
          actual_dat$age <- current_age
          actual_dat$patient_group <- patient_group[j]
          colnames(means) <- colnames(reg_data_no_stay)[-which(colnames(reg_data_no_stay) == "outcome")]
          
          means$WaitingTime <- transformation(mean(hes_data$WaitingTime, na.rm = TRUE))
          
          
          
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
          medians$WaitingTime <- transformation(median(hes_data$WaitingTime, na.rm = TRUE))
          
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
          seven_days$WaitingTime <- transformation(7)
          
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
    whole_df$week <- append_week
    coef_df$week <- append_week
    
    tot_whole_df <- dplyr::bind_rows(tot_whole_df, whole_df)
    tot_coef_df <- dplyr::bind_rows(tot_coef_df, coef_df)
  }
  
  
  
  return(list(tot_whole_df, tot_coef_df))
  
}

emergency_regression_cluster <- function(patient_group,hes_data_orig, start_date, forecast_length, forecast_start,
                                         time_trend = TRUE, month_trend = TRUE, week_num, half_week){
  ## Input ICD hes data for one ICD one age group one patient group (cc/ga)
  ## to be run in parrallel
  tic("Whole process")
  tot_whole_df <- NULL  
  forecast_seq <- seq(as.Date(forecast_start), by = "week", length.out = 52)
  hes_data_orig$index <- seq(1,nrow(hes_data_orig),1)
  require(stringr)
  require(lubridate)  
  require(reshape2)
  require(tictoc)
  
  index_df <- NULL  
  
  for(current_week in week_num){
    whole_df <- NULL
    if(half_week){
      append_week <- (current_week + (current_week - 1)) /2 
    }else{
      append_week <- current_week
    }
    
    
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
        hes_data$outcome <- "NA"
        if(half_week){
          ## To get the half day vals first look at 3 days 
          
          
          days_of_stay_variable <- seq(3,length.out = length(week_num), by = 7)
          days_of_stay_variable <- c(-1,days_of_stay_variable)
          lower_end_stay <- days_of_stay_variable[current_week] + 1
          upper_end_stay <- days_of_stay_variable[current_week + 1]
          half_day <- upper_end_stay
          
          half_dayers <- hes_data[hes_data$cc_LoS == half_day,]
          if(nrow(half_dayers) >0){
            adders <- sample(1:nrow(half_dayers),floor(nrow(half_dayers)/2), replace = F)
            half_row_to_add <- half_dayers[adders,]
            half_row_to_lose <- half_dayers[-adders,]
            
            patient_group_indy <- rep(patient_group[j], nrow(half_row_to_lose))
            indies <- half_row_to_lose$index
            new_indies <- cbind.data.frame(patient_group_indy, indies)
            colnames(new_indies) <- c("patient_group","index")
            
            
          }else{
            
            half_row_to_add <- NULL
            patient_group_indy <- patient_group[j]
            indies <- 0
            new_indies <- cbind.data.frame(patient_group_indy, indies)
            colnames(new_indies) <- c("patient_group","index")
            
          }
          
          
          ## go through the whole days first, set any above whole day cutoff as GA initially
          
          hes_data[hes_data$cc_LoS >= upper_end_stay ,"outcome"]<-"CC"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 3,"outcome"] <- "Dead"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 2,"outcome"] <- "GA"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 1,"outcome"] <- "Discharged"
          
          if(!is.null(half_row_to_add)){
            ## now we reallocate the half day values so the 3-3.5 values in
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$cc_transitions == 3,"outcome"] <- "Dead"
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$cc_transitions == 2,"outcome"] <- "GA"
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$cc_transitions == 1,"outcome"] <- "Discharged"
          }
          
          if(!is.null(index_df)){
            ## now we check if previous runs left a half out and then add that in to this so the 3.5-4 runs 
            patient_df <- index_df[index_df$patient_group == patient_group[j],]
            if(nrow(patient_df) >0){
              if(patient_df$index[1] != 0){
                
                
                hes_data[hes_data$index %in% patient_df$index & hes_data$cc_transitions == 3,"outcome"] <- "Dead"
                hes_data[hes_data$index %in% patient_df$index & hes_data$cc_transitions == 2,"outcome"] <- "GA"
                hes_data[hes_data$index %in% patient_df$index & hes_data$cc_transitions == 1,"outcome"] <- "Discharged"
                
                index_df <- index_df[-which(index_df$patient_group == patient_group[j]),]
                index_df <- dplyr::bind_rows(index_df, new_indies)
              }else{
                index_df <- index_df[-which(index_df$patient_group == patient_group[j]),]
                index_df <- dplyr::bind_rows(index_df, new_indies)
              }
              
              
              
              
            }else{
              
              index_df <- dplyr::bind_rows(index_df, new_indies)
            }
            
            
          }else{
            
            index_df <- dplyr::bind_rows(index_df, new_indies)
          }
        }else{
          days_of_stay_variable <- seq(0, 42, 7)
          lower_end_stay <- days_of_stay_variable[current_week]
          upper_end_stay <- days_of_stay_variable[current_week + 1]
          
          
          hes_data[hes_data$cc_LoS >= upper_end_stay ,"outcome"]<-"CC"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 3,"outcome"] <- "Dead"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 2,"outcome"] <- "GA"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 1,"outcome"] <- "Discharged"
        }
        if(length(which(is.na(hes_data$outcome))) > 0)
          hes_data <- hes_data[-which(is.na(hes_data$outcome)),]
        
        if(nrow(hes_data) == 0){
          actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
          colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
          actual_dat$reg_week <- forecast_week_nums
          actual_dat[,2:5] <- 0
          actual_dat$WT <- "actual"
          actual_dat$patient_group <- patient_group[j]
          
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- c("GA","Dead","CC","Discharged")
          mean_wt_pred <- as.data.frame(out_df)
          mean_wt_pred$patient_group <- patient_group[j]
          mean_wt_pred$ICD <- current_icd
          mean_wt_pred$age <- current_age
          mean_wt_pred$WT <- "mean"
          median_wt_pred <- mean_wt_pred
          median_wt_pred$WT <- "median"
          seven_wt_pred <- mean_wt_pred
          seven_wt_pred$WT <- "seven"
          mean_wt_pred$reg_week <- forecast_week_nums
          median_wt_pred$reg_week <- forecast_week_nums
          seven_wt_pred$reg_week <- forecast_week_nums
          
          tot_df <- dplyr::bind_rows(mean_wt_pred, median_wt_pred, seven_wt_pred, actual_dat)
          
          graphing_df <- melt(tot_df, id.vars = colnames(tot_df)[5:9])
          graphing_df$line_group <- paste(graphing_df$variable, graphing_df$WT, sep = "-")
          
          
          whole_df <- dplyr::bind_rows(whole_df, graphing_df)
          
          
          next()
          
        }
        
        
        
        
        
        hes_data$outcome <- factor(hes_data$outcome, levels = c("CC","Dead","GA","Discharged"))
        hes_data$outcome <- relevel(hes_data$outcome, ref = "GA")
        
        
        toc()
        
      }else{
        tic("Narrowing to GA only")
        hes_data <- hes_data[hes_data$cc_start_flg == 0,]
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
        hes_data$outcome <- "NA"
        
        if(half_week){
          ## To get the half day vals first look at 3 days 
          
          
          days_of_stay_variable <- seq(3,length.out = length(week_num), by = 7)
          days_of_stay_variable <- c(-1,days_of_stay_variable)
          lower_end_stay <- days_of_stay_variable[current_week] + 1
          upper_end_stay <- days_of_stay_variable[current_week + 1]
          half_day <- upper_end_stay
          
          half_dayers <- hes_data[hes_data$GA_LoS == half_day,]
          if(nrow(half_dayers) >0){
            adders <- sample(1:nrow(half_dayers),floor(nrow(half_dayers)/2), replace = F)
            half_row_to_add <- half_dayers[adders,]
            half_row_to_lose <- half_dayers[-adders,]
            
            patient_group_indy <- rep(patient_group[j], nrow(half_row_to_lose))
            indies <- half_row_to_lose$index
            new_indies <- cbind.data.frame(patient_group_indy, indies)
            colnames(new_indies) <- c("patient_group","index")
            
            
          }else{
            
            half_row_to_add <- NULL
            patient_group_indy <- patient_group[j]
            indies <- 0
            new_indies <- cbind.data.frame(patient_group_indy, indies)
            colnames(new_indies) <- c("patient_group","index")
            
          }
          
          
          ## go through the whole days first, set any above whole day cutoff as GA initially
          
          hes_data[hes_data$GA_LoS >= upper_end_stay ,"outcome"]<-"GA"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 3,"outcome"] <- "Dead"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 2,"outcome"] <- "CC"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 1,"outcome"] <- "Discharged"
          
          if(!is.null(half_row_to_add)){
            ## now we reallocate the half day values so the 3-3.5 values in
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$ga_transitions == 3,"outcome"] <- "Dead"
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$ga_transitions == 2,"outcome"] <- "CC"
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$ga_transitions == 1,"outcome"] <- "Discharged"
          }
          
          if(!is.null(index_df)){
            ## now we check if previous runs left a half out and then add that in to this so the 3.5-4 runs 
            patient_df <- index_df[index_df$patient_group == patient_group[j],]
            if(nrow(patient_df) >0){
              if(patient_df$index[1] != 0){
                hes_data[hes_data$index %in% patient_df$index & hes_data$ga_transitions == 3,"outcome"] <- "Dead"
                hes_data[hes_data$index %in% patient_df$index & hes_data$ga_transitions == 2,"outcome"] <- "CC"
                hes_data[hes_data$index %in% patient_df$index & hes_data$ga_transitions == 1,"outcome"] <- "Discharged"
                
                index_df <- index_df[-which(index_df$patient_group == patient_group[j]),]
                index_df <- dplyr::bind_rows(index_df, new_indies)
              }else{
                index_df <- index_df[-which(index_df$patient_group == patient_group[j]),]
                index_df <- dplyr::bind_rows(index_df, new_indies)
              }
              
              
              
              
            }else{
              
              index_df <- dplyr::bind_rows(index_df, new_indies)
            }
            
            
          }else{
            
            index_df <- dplyr::bind_rows(index_df, new_indies)
          }
          
          
          
          
        }else{
        
        
          days_of_stay_variable <- seq(0, 42, 7)
          lower_end_stay <- days_of_stay_variable[current_week]
          upper_end_stay <- days_of_stay_variable[current_week + 1]
          
          
          hes_data[hes_data$GA_LoS >= upper_end_stay ,"outcome"]<-"GA"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 3,"outcome"] <- "Dead"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 2,"outcome"] <- "CC"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 1,"outcome"] <- "Discharged"
        }
          
        if(length(which(is.na(hes_data$outcome))) > 0)
          hes_data <- hes_data[-which(is.na(hes_data$outcome)),]
        
        if(nrow(hes_data) == 0){
          actual_dat <- data.frame(matrix(data = NA, nrow = forecast_length, ncol = 5))
          colnames(actual_dat) <- c("reg_week","GA","Dead","CC","Discharged")
          actual_dat$reg_week <- forecast_week_nums
          actual_dat[,2:5] <- 0
          actual_dat$WT <- "actual"
          actual_dat$patient_group <- patient_group[j]
          
          
          out_df <- matrix(data = 0,ncol = 4, nrow = forecast_length)
          colnames(out_df) <- c("GA","Dead","CC","Discharged")
          mean_wt_pred <- as.data.frame(out_df)
          mean_wt_pred$patient_group <- patient_group[j]
          mean_wt_pred$ICD <- current_icd
          mean_wt_pred$age <- current_age
          mean_wt_pred$WT <- "mean"
          median_wt_pred <- mean_wt_pred
          median_wt_pred$WT <- "median"
          seven_wt_pred <- mean_wt_pred
          seven_wt_pred$WT <- "seven"
          mean_wt_pred$reg_week <- forecast_week_nums
          median_wt_pred$reg_week <- forecast_week_nums
          seven_wt_pred$reg_week <- forecast_week_nums
          
          tot_df <- dplyr::bind_rows(mean_wt_pred, median_wt_pred, seven_wt_pred, actual_dat)
          
          graphing_df <- melt(tot_df, id.vars = colnames(tot_df)[5:9])
          graphing_df$line_group <- paste(graphing_df$variable, graphing_df$WT, sep = "-")
          
          
          whole_df <- dplyr::bind_rows(whole_df, graphing_df)
          
          
          next()
          
        }
        
          
        hes_data$outcome <- factor(hes_data$outcome, levels = c("CC","Dead","GA","Discharged"))
        hes_data$outcome <- relevel(hes_data$outcome, ref = "GA")
        
        
        toc()
      }
      ## set up week and year fixed effects
      
      forecast_week_start <- week_num_func(forecast_start,start_date)
      forecast_week_end <- week_num_func(forecast_seq[forecast_length],start_date)
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
      
      
      
      if(month_trend | time_trend)
        outcomes_in_dat <- plyr::count(reg_data_no_stay$outcome)
      else
        outcomes_in_dat <- plyr::count(reg_data_no_stay)
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
          actual_dat$patient_group <- patient_group[j]
          
          
          
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
      
      tot_df <- dplyr::bind_rows(mean_wt_pred, median_wt_pred, actual_dat)
      
      graphing_df <- melt(tot_df, id.vars = colnames(tot_df)[5:14])
      graphing_df$line_group <- paste(graphing_df$variable, graphing_df$WT, sep = "-")
      
      toc()
      
      whole_df <- dplyr::bind_rows(whole_df, graphing_df)
      
      
        
    }
    
    whole_df$week <- append_week
    tot_whole_df <- dplyr::bind_rows(tot_whole_df, whole_df)    
    
  }
  toc()
  return(tot_whole_df)
  
}

elective_regression2 <- function(patient_group,hes_data_orig, start_date, forecast_length, forecast_start,
                                 time_trend = TRUE, month_trend = TRUE, wt_variable = "linear",
                                 week_num, half_week){
  ## Function to creater mulitnomial logit on waiting time for transitions within hospital for patients
  ## Works as a function in SNOW clustering.
  ## Input: Patient group, strings in the form of ICD-ward-AGE (02-ga-1)
  ##        hes_data_orig - transitions data frame narrowed form cluster set up function
  ##        week_num - vector of weeks to calculate at 
  ##        half_week - TRUE calculate values at 3.5 10.5 etc 
  

  ## Set up the outdfs 
  tic("Whole process")
  tot_whole_df <- NULL
  tot_coef_df <- NULL
  index_df <- NULL
  hes_data_orig$index <- seq(1,nrow(hes_data_orig),1)
  require(stringr)
  require(lubridate)  
  require(reshape2)
  require(tictoc)
  

  ## First loop through the week_numbers input 
  
  for(current_week in week_num){
    whole_df <- data.frame(matrix(ncol = 9, nrow = length(patient_group)))
    colnames(whole_df) <- c("GA","CC","Dead","Discharged","patient_group","ICD","age","WT", "group_size")
    coef_df <- NULL
    if(half_week){
      append_week <- (current_week + (current_week - 1)) /2 
    }else{
      append_week <- current_week
    }
    
    ## For each week now loop through the patient groups 
    for(j in 1:length(patient_group)){
      
      print(paste("On patient group:",patient_group[j],". Number",j, "of",length(patient_group)))
      tic(paste("Running for patient group:",patient_group[j]))
      
      ## Get the patient group charecteristics 
      split_patient_group <- str_split_fixed(patient_group[j], "-",3)
      current_icd <- as.integer(split_patient_group[1])
      current_ward <- split_patient_group[2]
      current_age <- as.integer(split_patient_group[3])
      ## Narrow down the HES data
      hes_data <- hes_data_orig[hes_data_orig$ICD == current_icd & 
                                  hes_data_orig$agegrp_v3 == current_age,]
      

      if(current_ward == "cc"){
        
        tic("Narrowing to CC only")
        hes_data <- hes_data[hes_data$cc == 1,]
        toc()
        ## remove na trans
        tic("Removing NAs")
        na_rows <- which(is.na(hes_data$cc_transitions))
        if(length(na_rows) > 0){
          hes_data <- hes_data[-na_rows,]
        }
        na_cc_los <- which(is.na(hes_data$cc_LoS))
        if(length(na_cc_los) > 0){
          hes_data <- hes_data[-na_cc_los,]
        }
        na_cc_WT <- which(is.na(hes_data$WaitingTime))
        if(length(na_cc_WT) > 0){
          hes_data <- hes_data[-na_cc_WT,]
        }
        
        toc()
        
        ## set up transitions 
        tic("Set up outcome variable")
        hes_data$outcome <- NA_character_
        
        if(half_week){
          ## To get the half day vals first look at 3 days 
          ## First set up the sequence of days to draw from 
          days_of_stay_variable <- seq(3,length.out = length(week_num), by = 7)
          days_of_stay_variable <- c(-1,days_of_stay_variable)
          ## Now we set an upper and lower initial bound ( -1 is to ensure continuity of the sequence)
          ## upper is at day 3 10 etc, lower is at 0, 4, 11, upper is the half day so 3.5 
          lower_end_stay <- days_of_stay_variable[current_week] + 1
          upper_end_stay <- days_of_stay_variable[current_week + 1]
          half_day <- upper_end_stay
          
          half_dayers <- hes_data[hes_data$cc_LoS == half_day,]
          if(nrow(half_dayers) >0){
            ## If there is any LoS at half day then we look at randomly sampling 
            ## half these to take as if they transitioned by midpoint of the day. 
            adders <- sample(1:nrow(half_dayers),floor(nrow(half_dayers)/2), replace = F)
            half_row_to_add <- half_dayers[adders,]
            half_row_to_lose <- half_dayers[-adders,]
            
            ## Keep track of those lost so we can add them back in for the next week's 
            ## values 
            patient_group_indy <- rep(patient_group[j], nrow(half_row_to_lose))
            indies <- half_row_to_lose$index
            new_indies <- cbind.data.frame(patient_group_indy, indies)
            colnames(new_indies) <- c("patient_group","index")
            
            
          }else{
            
            half_row_to_add <- NULL
            patient_group_indy <- patient_group[j]
            indies <- 0
            new_indies <- cbind.data.frame(patient_group_indy, indies)
            colnames(new_indies) <- c("patient_group","index")
            
          }
          
          
          ## go through the whole days first, set any above whole day cutoff as GA initially
          
          hes_data[hes_data$cc_LoS >= upper_end_stay ,"outcome"]<-"CC"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 3,"outcome"] <- "Dead"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 2,"outcome"] <- "GA"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 1,"outcome"] <- "Discharged"
          
          if(!is.null(half_row_to_add)){
            ## now we reallocate the half day values so the 3-3.5 values in
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$cc_transitions == 3,"outcome"] <- "Dead"
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$cc_transitions == 2,"outcome"] <- "GA"
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$cc_transitions == 1,"outcome"] <- "Discharged"
          }
          
          if(!is.null(index_df)){
            ## now we check if previous runs left a half out and then add that in to this so the 3.5-4 runs 
            patient_df <- index_df[index_df$patient_group == patient_group[j],]
            if(nrow(patient_df) >0){
              if(patient_df$index[1] != 0){
                
                
                hes_data[hes_data$index %in% patient_df$index & hes_data$cc_transitions == 3,"outcome"] <- "Dead"
                hes_data[hes_data$index %in% patient_df$index & hes_data$cc_transitions == 2,"outcome"] <- "GA"
                hes_data[hes_data$index %in% patient_df$index & hes_data$cc_transitions == 1,"outcome"] <- "Discharged"
                
                index_df <- index_df[-which(index_df$patient_group == patient_group[j]),]
                index_df <- dplyr::bind_rows(index_df, new_indies)
              }else{
                index_df <- index_df[-which(index_df$patient_group == patient_group[j]),]
                index_df <- dplyr::bind_rows(index_df, new_indies)
              }
              
              
              
              
            }else{
              
              index_df <- dplyr::bind_rows(index_df, new_indies)
            }
            
            
          }else{
            
            index_df <- dplyr::bind_rows(index_df, new_indies)
          }
        }else{
          ## Simpler for whole week values, just set the upper and lower bounds on the multiples of 7
          
          days_of_stay_variable <- seq(0, 42, 7)
          lower_end_stay <- days_of_stay_variable[current_week]
          upper_end_stay <- days_of_stay_variable[current_week + 1]
          
          
          hes_data[hes_data$cc_LoS >= upper_end_stay ,"outcome"]<-"CC"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 3,"outcome"] <- "Dead"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 2,"outcome"] <- "GA"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 1,"outcome"] <- "Discharged"
        }
        if(length(which(is.na(hes_data$outcome))) > 0){
          hes_data <- hes_data[-which(is.na(hes_data$outcome)),]
        }
        
        if(nrow(hes_data) == 0){
          ## If there is no data for this week of stay we just set the probabilities to 0
          current_coef_df <- data.frame(matrix(0, ncol = 5, nrow = 4))
          colnames(current_coef_df) <- c("transition","coef","patient_group","ICD","age")
          current_coef_df$transition <- c("CC","GA","Dead","Discharged")
          current_coef_df$patient_group <- patient_group[j]
          current_coef_df$ICD <- current_icd
          current_coef_df$age <- current_age
          
          out_df <- as.data.frame(matrix(data = 0,ncol = 4, nrow = 1))
          colnames(out_df) <- c("GA","Dead","CC","Discharged")
          seven_wt_pred <- as.data.frame(out_df)
          seven_wt_pred$patient_group <- patient_group[j]
          seven_wt_pred$ICD <- current_icd
          seven_wt_pred$age <- current_age
          seven_wt_pred$WT <- "seven"
          seven_wt_pred$group_size <- 0
          
          seven_wt_pred <- seven_wt_pred %>% select(GA,CC,Dead,Discharged,patient_group,ICD,age,WT, group_size)
          whole_df[j,] <- seven_wt_pred 
          
          
          coef_df <- dplyr::bind_rows(coef_df, current_coef_df)
          
          next()
          
        }
        
        
        
        
        
        hes_data$outcome <- factor(hes_data$outcome, levels = c("CC","Dead","GA","Discharged"))
        hes_data$outcome <- relevel(hes_data$outcome, ref = "GA")
        
        
        toc()
        
      }else{
        ## remove na trans
        hes_data <- hes_data[hes_data$cc_start_flg == 0,]
        tic("Removing NAs")
        na_rows <- which(is.na(hes_data$ga_transitions))
        if(length(na_rows) > 0){
          hes_data <- hes_data[-na_rows,]
        }
        na_ga_los <- which(is.na(hes_data$GA_LoS))
        if(length(na_ga_los) > 0)
          hes_data <- hes_data[-na_ga_los,]
        
        na_ga_WT <- which(is.na(hes_data$WaitingTime))
        if(length(na_ga_WT) > 0)
          hes_data <- hes_data[-na_ga_WT,]
        
        
        toc()
        
        ## set up transitions 
        tic("Set up outcome variable")
        hes_data$outcome <- NA_character_
        
        if(half_week){
          ## To get the half day vals first look at 3 days
          ## Same methodology as above for CC patients 
          
          days_of_stay_variable <- seq(3,length.out = length(week_num), by = 7)
          days_of_stay_variable <- c(-1,days_of_stay_variable)
          lower_end_stay <- days_of_stay_variable[current_week] + 1
          upper_end_stay <- days_of_stay_variable[current_week + 1]
          half_day <- upper_end_stay
          
          half_dayers <- hes_data[hes_data$GA_LoS == half_day,]
          if(nrow(half_dayers) >0){
            adders <- sample(1:nrow(half_dayers),floor(nrow(half_dayers)/2), replace = F)
            half_row_to_add <- half_dayers[adders,]
            half_row_to_lose <- half_dayers[-adders,]
            
            patient_group_indy <- rep(patient_group[j], nrow(half_row_to_lose))
            indies <- half_row_to_lose$index
            new_indies <- cbind.data.frame(patient_group_indy, indies)
            colnames(new_indies) <- c("patient_group","index")
            
            
          }else{
            
            half_row_to_add <- NULL
            patient_group_indy <- patient_group[j]
            indies <- 0
            new_indies <- cbind.data.frame(patient_group_indy, indies)
            colnames(new_indies) <- c("patient_group","index")
            
          }
          
          
          ## go through the whole days first, set any above whole day cutoff as GA initially
          
          hes_data[hes_data$GA_LoS >= upper_end_stay ,"outcome"]<-"GA"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 3,"outcome"] <- "Dead"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 2,"outcome"] <- "CC"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 1,"outcome"] <- "Discharged"
          
          if(!is.null(half_row_to_add)){
            ## now we reallocate the half day values so the 3-3.5 values in
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$ga_transitions == 3,"outcome"] <- "Dead"
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$ga_transitions == 2,"outcome"] <- "CC"
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$ga_transitions == 1,"outcome"] <- "Discharged"
          }
          
          if(!is.null(index_df)){
            ## now we check if previous runs left a half out and then add that in to this so the 3.5-4 runs 
            patient_df <- index_df[index_df$patient_group == patient_group[j],]
            if(nrow(patient_df) >0){
              if(patient_df$index[1] != 0){
                hes_data[hes_data$index %in% patient_df$index & hes_data$ga_transitions == 3,"outcome"] <- "Dead"
                hes_data[hes_data$index %in% patient_df$index & hes_data$ga_transitions == 2,"outcome"] <- "CC"
                hes_data[hes_data$index %in% patient_df$index & hes_data$ga_transitions == 1,"outcome"] <- "Discharged"
                
                index_df <- index_df[-which(index_df$patient_group == patient_group[j]),]
                index_df <- dplyr::bind_rows(index_df, new_indies)
              }else{
                index_df <- index_df[-which(index_df$patient_group == patient_group[j]),]
                index_df <- dplyr::bind_rows(index_df, new_indies)
              }
              
            }else{
              
              index_df <- dplyr::bind_rows(index_df, new_indies)
            }
            
            
          }else{
            
            index_df <- dplyr::bind_rows(index_df, new_indies)
          }
          
          
          
          
        }else{
          
          
          days_of_stay_variable <- seq(0, 42, 7)
          lower_end_stay <- days_of_stay_variable[current_week]
          upper_end_stay <- days_of_stay_variable[current_week + 1]
          
          
          hes_data[hes_data$GA_LoS >= upper_end_stay ,"outcome"]<-"GA"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 3,"outcome"] <- "Dead"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 2,"outcome"] <- "CC"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 1,"outcome"] <- "Discharged"
        }
        if(length(which(is.na(hes_data$outcome))) > 0)
          hes_data <- hes_data[-which(is.na(hes_data$outcome)),]
        
        if(nrow(hes_data) == 0){
          current_coef_df <- data.frame(matrix(0, ncol = 5, nrow = 4))
          colnames(current_coef_df) <- c("transition","coef","patient_group","ICD","age")
          current_coef_df$transition <- c("CC","GA","Dead","Discharged")
          current_coef_df$patient_group <- patient_group[j]
          current_coef_df$ICD <- current_icd
          current_coef_df$age <- current_age
          
          out_df <- as.data.frame(matrix(data = 0,ncol = 4, nrow = 1))
          colnames(out_df) <- c("GA","Dead","CC","Discharged")
          seven_wt_pred <- as.data.frame(out_df)
          seven_wt_pred$patient_group <- patient_group[j]
          seven_wt_pred$ICD <- current_icd
          seven_wt_pred$age <- current_age
          seven_wt_pred$WT <- "seven"
          
          seven_wt_pred <- seven_wt_pred %>% select(GA,CC,Dead,Discharged,patient_group,ICD,age,WT)
          whole_df[j,] <- seven_wt_pred 
          
          coef_df <- dplyr::bind_rows(coef_df, current_coef_df)
          
          next()
          
        }
        
        hes_data$outcome <- factor(hes_data$outcome, levels = c("CC","Dead","GA","Discharged"))
        hes_data$outcome <- relevel(hes_data$outcome, ref = "GA")
        
      }
        
      toc()
      reg_data_no_stay <- hes_data[,which(colnames(hes_data) %in% c("outcome","WaitingTime"))]

      outcomes_in_dat <- plyr::count(reg_data_no_stay$outcome)
      
      if(nrow(outcomes_in_dat) < 2){
        
        ## Need at least two separate outcomes to fit the model 
        ## If only 1 set the transitions as one for that outcome 
        
        out_df <- as.data.frame(matrix(data = 0,ncol = 4, nrow = 1))
        colnames(out_df) <- c("GA","CC","Dead","Discharged")
        out_df[,as.character(outcomes_in_dat[1,1])] <- 1
        mean_wt_pred <- as.data.frame(out_df)
        mean_wt_pred$patient_group <- patient_group[j]
        mean_wt_pred$ICD <- current_icd
        mean_wt_pred$age <- current_age
        mean_wt_pred$WT <- "seven"
        mean_wt_pred$group_size <- nrow(reg_data_no_stay)
        
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
        
        ## Create the coefs by looking at a one day increase in Waiting time and how that affects
        ## outcomes 
        
        wt_0_dat <- as.data.frame(matrix(0, ncol = ncol(reg_data_no_stay) - 1, nrow = 1))
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
        
        
        tic("Probability creation")
        ## Create one row df of Waiting time at 7 days for the transitions probabilities 
        means <- as.data.frame(matrix(data = 0,ncol = ncol(reg_data_no_stay) - 1, nrow = 1))
          
        ## make the actual data df to compare to predictions 
        means$WaitingTime <- 7
        ## Now predict 
          
        mean_wt_pred <- predict(ml_stay, newdata = means, "probs")
        ## Now go through outcomes to ensure that prediction returns values for all transitions
        ## If two outcomes probably won't be a matirx 
        print("Done the predicting")
        if(length(mean_wt_pred) == 1){
          print("In length 1")
          
          out_df <- matrix(data = 0,ncol = 4, nrow = 1)
          colnames(out_df) <- c("GA","CC","Dead","Discharged")
          single_val <- mean_wt_pred
          max_count <- which.is.max(outcomes_in_dat[,2])
          if(single_val > 0.5){
            out_df[,outcomes_in_dat[max_count, 1]] <- mean_wt_pred
            out_df[,outcomes_in_dat[-max_count,1]] <- 1 - mean_wt_pred
          }else{
            out_df[,outcomes_in_dat[max_count, 1]] <- 1 - mean_wt_pred
            out_df[,outcomes_in_dat[-max_count,1]] <- mean_wt_pred
          }
            
          mean_wt_pred <- out_df
        }else if(length(mean_wt_pred) > 1 ){
          print("In the length > 1 section")
          out_df <- matrix(data = 0,ncol = 4, nrow = 1)
          colnames(out_df) <- levels(hes_data$outcome)
          for(column in names(mean_wt_pred)){
            out_df[,column] <- mean_wt_pred[column]
          
            }
            
            mean_wt_pred <- out_df
          }
          
          mean_wt_pred <- as.data.frame(mean_wt_pred)
          mean_wt_pred$patient_group <- patient_group[j]
          mean_wt_pred$ICD <- current_icd
          mean_wt_pred$age <- current_age
          mean_wt_pred$WT <- "seven"
          mean_wt_pred$group_size <- nrow(reg_data_no_stay)
                    
          toc()
          
          
      }
      
      
      
      mean_wt_pred <- mean_wt_pred %>% select(GA,CC,Dead,Discharged,patient_group,ICD,age,WT, group_size)
      whole_df[j,] <- mean_wt_pred 
      coef_df <- dplyr::bind_rows(coef_df, current_coef_df)
      
      }
      
      
      toc()
      whole_df$week <- append_week
      coef_df$week <- append_week
      
      tot_whole_df <- dplyr::bind_rows(tot_whole_df, whole_df)
      tot_coef_df <- dplyr::bind_rows(tot_coef_df, coef_df)
    }
  

    
    return(list(tot_whole_df, tot_coef_df))
    
}
      
emergency_regression2 <- function(patient_group,hes_data_orig, start_date, forecast_length, forecast_start,
                                 time_trend = TRUE, month_trend = TRUE, wt_variable = "linear",
                                 week_num, half_week){
  ## Function to get mean proportions for transitions within hospital for patients as Electives
  ## Works as a function in SNOW clustering.
  ## Input: Patient group, strings in the form of ICD-ward-AGE (02-ga-1)
  ##        hes_data_orig - transitions data frame narrowed form cluster set up function
  ##        week_num - vector of weeks to calculate at 
  ##        half_week - TRUE calculate values at 3.5 10.5 etc 
  

  ## Set up the outdfs 
  tic("Whole process")
  tot_whole_df <- NULL
  tot_coef_df <- NULL
  index_df <- NULL
  hes_data_orig$index <- seq(1,nrow(hes_data_orig),1)
  require(stringr)
  require(lubridate)  
  require(reshape2)
  require(tictoc)
  
  ## First loop through the week_numbers input 
  for(current_week in week_num){
    whole_df <- data.frame(matrix(ncol = 9, nrow = length(patient_group)))
    colnames(whole_df) <- c("GA","CC","Dead","Discharged","patient_group","ICD","age","WT", "group_size")
    if(half_week){
      append_week <- (current_week + (current_week - 1)) /2 
    }else{
      append_week <- current_week
    }
    
    ## For each week now loop through the patient groups 
    for(j in 1:length(patient_group)){
      
      print(paste("On patient group:",patient_group[j],". Number",j, "of",length(patient_group)))
      tic(paste("Running for patient group:",patient_group[j]))
      
      ## Get the patient group charecteristics 
      split_patient_group <- str_split_fixed(patient_group[j], "-",3)
      current_icd <- as.integer(split_patient_group[1])
      current_ward <- split_patient_group[2]
      current_age <- as.integer(split_patient_group[3])
      ## Narrow down the HES data
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
        
        
        toc()
        ## set up transitions 
        tic("Set up outcome variable")
        hes_data$outcome <- NA_character_
        
        if(half_week){
          ## To get the half day vals first look at 3 days 
          ## First set up the sequence of days to draw from 
          days_of_stay_variable <- seq(3,length.out = length(week_num), by = 7)
          days_of_stay_variable <- c(-1,days_of_stay_variable)
          ## Now we set an upper and lower initial bound ( -1 is to ensure continuity of the sequence)
          ## upper is at day 3 10 etc, lower is at 0, 4, 11, upper is the half day so 3.5 
          lower_end_stay <- days_of_stay_variable[current_week] + 1
          upper_end_stay <- days_of_stay_variable[current_week + 1]
          half_day <- upper_end_stay
          
          half_dayers <- hes_data[hes_data$cc_LoS == half_day,]
          if(nrow(half_dayers) >0){
            ## If there is any LoS at half day then we look at randomly sampling 
            ## half these to take as if they transitioned by midpoint of the day. 
            adders <- sample(1:nrow(half_dayers),floor(nrow(half_dayers)/2), replace = F)
            half_row_to_add <- half_dayers[adders,]
            half_row_to_lose <- half_dayers[-adders,]
            
            ## Keep track of those lost so we can add them back in for the next week's 
            ## values 
            patient_group_indy <- rep(patient_group[j], nrow(half_row_to_lose))
            indies <- half_row_to_lose$index
            new_indies <- cbind.data.frame(patient_group_indy, indies)
            colnames(new_indies) <- c("patient_group","index")
            
            
          }else{
            
            half_row_to_add <- NULL
            patient_group_indy <- patient_group[j]
            indies <- 0
            new_indies <- cbind.data.frame(patient_group_indy, indies)
            colnames(new_indies) <- c("patient_group","index")
            
          }
          
          
          ## go through the whole days first, set any above whole day cutoff as GA initially
          
          hes_data[hes_data$cc_LoS >= upper_end_stay ,"outcome"]<-"CC"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 3,"outcome"] <- "Dead"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 2,"outcome"] <- "GA"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 1,"outcome"] <- "Discharged"
          
          if(!is.null(half_row_to_add)){
            ## now we reallocate the half day values so the 3-3.5 values in
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$cc_transitions == 3,"outcome"] <- "Dead"
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$cc_transitions == 2,"outcome"] <- "GA"
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$cc_transitions == 1,"outcome"] <- "Discharged"
          }
          
          if(!is.null(index_df)){
            ## now we check if previous runs left a half out and then add that in to this so the 3.5-4 runs 
            patient_df <- index_df[index_df$patient_group == patient_group[j],]
            if(nrow(patient_df) >0){
              if(patient_df$index[1] != 0){
                
                
                hes_data[hes_data$index %in% patient_df$index & hes_data$cc_transitions == 3,"outcome"] <- "Dead"
                hes_data[hes_data$index %in% patient_df$index & hes_data$cc_transitions == 2,"outcome"] <- "GA"
                hes_data[hes_data$index %in% patient_df$index & hes_data$cc_transitions == 1,"outcome"] <- "Discharged"
                
                index_df <- index_df[-which(index_df$patient_group == patient_group[j]),]
                index_df <- dplyr::bind_rows(index_df, new_indies)
              }else{
                index_df <- index_df[-which(index_df$patient_group == patient_group[j]),]
                index_df <- dplyr::bind_rows(index_df, new_indies)
              }
              
              
              
              
            }else{
              
              index_df <- dplyr::bind_rows(index_df, new_indies)
            }
            
            
          }else{
            
            index_df <- dplyr::bind_rows(index_df, new_indies)
          }
        }else{
          ## Simpler for whole week values, just set the upper and lower bounds on the multiples of 7
          
          days_of_stay_variable <- seq(0, 42, 7)
          lower_end_stay <- days_of_stay_variable[current_week]
          upper_end_stay <- days_of_stay_variable[current_week + 1]
          
          
          hes_data[hes_data$cc_LoS >= upper_end_stay ,"outcome"]<-"CC"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 3,"outcome"] <- "Dead"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 2,"outcome"] <- "GA"
          hes_data[hes_data$cc_LoS >= lower_end_stay & hes_data$cc_LoS < upper_end_stay & hes_data$cc_transitions == 1,"outcome"] <- "Discharged"
        }
        if(length(which(is.na(hes_data$outcome))) > 0){
          hes_data <- hes_data[-which(is.na(hes_data$outcome)),]
        }
        
        if(nrow(hes_data) == 0){
          ## If there is no data for this week of stay we just set the probabilities to 0
          out_df <- matrix(data = 0,ncol = 4, nrow = 1)
          colnames(out_df) <- c("GA","Dead","CC","Discharged")
          seven_wt_pred <- as.data.frame(out_df)
          seven_wt_pred$patient_group <- patient_group[j]
          seven_wt_pred$ICD <- current_icd
          seven_wt_pred$age <- current_age
          seven_wt_pred$WT <- "seven"
          seven_wt_pred$group_size <- 0
          
          seven_wt_pred <- seven_wt_pred %>% select(GA,CC,Dead,Discharged,patient_group,ICD,age,WT, group_size)
          whole_df[j,] <- seven_wt_pred 
          
          
          next()
          
        }
        
        toc()
        
      }else{
        ## remove na trans
        hes_data <- hes_data[hes_data$cc_start_flg == 0,]
        tic("Removing NAs")
        na_rows <- which(is.na(hes_data$ga_transitions))
        if(length(na_rows) > 0){
          hes_data <- hes_data[-na_rows,]
        }
        
        na_ga_los <- which(is.na(hes_data$GA_LoS))
        if(length(na_ga_los) > 0){
          hes_data <- hes_data[-na_ga_los,]
        }

        
        toc()
        
        ## set up transitions 
        tic("Set up outcome variable")
        hes_data$outcome <- NA_character_
        
        if(half_week){
          ## To get the half day vals first look at 3 days
          ## Same methodology as above for CC patients 
          
          days_of_stay_variable <- seq(3,length.out = length(week_num), by = 7)
          days_of_stay_variable <- c(-1,days_of_stay_variable)
          lower_end_stay <- days_of_stay_variable[current_week] + 1
          upper_end_stay <- days_of_stay_variable[current_week + 1]
          half_day <- upper_end_stay
          
          half_dayers <- hes_data[hes_data$GA_LoS == half_day,]
          if(nrow(half_dayers) >0){
            adders <- sample(1:nrow(half_dayers),floor(nrow(half_dayers)/2), replace = F)
            half_row_to_add <- half_dayers[adders,]
            half_row_to_lose <- half_dayers[-adders,]
            
            patient_group_indy <- rep(patient_group[j], nrow(half_row_to_lose))
            indies <- half_row_to_lose$index
            new_indies <- cbind.data.frame(patient_group_indy, indies)
            colnames(new_indies) <- c("patient_group","index")
            
            
          }else{
            
            half_row_to_add <- NULL
            patient_group_indy <- patient_group[j]
            indies <- 0
            new_indies <- cbind.data.frame(patient_group_indy, indies)
            colnames(new_indies) <- c("patient_group","index")
            
          }
          
          
          ## go through the whole days first, set any above whole day cutoff as GA initially
          
          hes_data[hes_data$GA_LoS >= upper_end_stay ,"outcome"]<-"GA"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 3,"outcome"] <- "Dead"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 2,"outcome"] <- "CC"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 1,"outcome"] <- "Discharged"
          
          if(!is.null(half_row_to_add)){
            ## now we reallocate the half day values so the 3-3.5 values in
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$ga_transitions == 3,"outcome"] <- "Dead"
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$ga_transitions == 2,"outcome"] <- "CC"
            hes_data[hes_data$index %in% half_row_to_add$index & hes_data$ga_transitions == 1,"outcome"] <- "Discharged"
          }
          
          if(!is.null(index_df)){
            ## now we check if previous runs left a half out and then add that in to this so the 3.5-4 runs 
            patient_df <- index_df[index_df$patient_group == patient_group[j],]
            if(nrow(patient_df) >0){
              if(patient_df$index[1] != 0){
                hes_data[hes_data$index %in% patient_df$index & hes_data$ga_transitions == 3,"outcome"] <- "Dead"
                hes_data[hes_data$index %in% patient_df$index & hes_data$ga_transitions == 2,"outcome"] <- "CC"
                hes_data[hes_data$index %in% patient_df$index & hes_data$ga_transitions == 1,"outcome"] <- "Discharged"
                
                index_df <- index_df[-which(index_df$patient_group == patient_group[j]),]
                index_df <- dplyr::bind_rows(index_df, new_indies)
              }else{
                index_df <- index_df[-which(index_df$patient_group == patient_group[j]),]
                index_df <- dplyr::bind_rows(index_df, new_indies)
              }
              
            }else{
              
              index_df <- dplyr::bind_rows(index_df, new_indies)
            }
            
            
          }else{
            
            index_df <- dplyr::bind_rows(index_df, new_indies)
          }
          
          
          
          
        }else{
          
          
          days_of_stay_variable <- seq(0, 42, 7)
          lower_end_stay <- days_of_stay_variable[current_week]
          upper_end_stay <- days_of_stay_variable[current_week + 1]
          
          
          hes_data[hes_data$GA_LoS >= upper_end_stay ,"outcome"]<-"GA"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 3,"outcome"] <- "Dead"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 2,"outcome"] <- "CC"
          hes_data[hes_data$GA_LoS >= lower_end_stay & hes_data$GA_LoS < upper_end_stay & hes_data$ga_transitions == 1,"outcome"] <- "Discharged"
        }
        if(length(which(is.na(hes_data$outcome))) > 0){
          hes_data <- hes_data[-which(is.na(hes_data$outcome)),]
        }
        
        if(nrow(hes_data) == 0){
          out_df <- matrix(data = 0,ncol = 4, nrow = 1)
          colnames(out_df) <- c("GA","Dead","CC","Discharged")
          seven_wt_pred <- as.data.frame(out_df)
          seven_wt_pred$patient_group <- patient_group[j]
          seven_wt_pred$ICD <- current_icd
          seven_wt_pred$age <- current_age
          seven_wt_pred$WT <- "seven"
          seven_wt_pred$group_size <- 0
          
          seven_wt_pred <- seven_wt_pred %>% select(GA,CC,Dead,Discharged,patient_group,ICD,age,WT, group_size)
          whole_df[j,] <- seven_wt_pred 
          
          
          
          next()
          
        }
        
        
      }
      
      toc()
      reg_data_no_stay <- hes_data[,which(colnames(hes_data) %in% c("outcome"))]
      
      tic("Probability creation")
      
      mean_wt_pred <- as.data.frame(matrix(data = 0, nrow = 1, ncol = 4))
      colnames(mean_wt_pred) <- c("GA","CC","Dead","Discharged")
      denominator <- nrow(reg_data_no_stay)
      GA_num <- length(reg_data_no_stay[reg_data_no_stay == "GA"])
      cc_num <- length(reg_data_no_stay[reg_data_no_stay == "CC"])
      dead_num <- length(reg_data_no_stay[reg_data_no_stay == "Dead"])
      dis_num <- length(reg_data_no_stay[reg_data_no_stay == "Discharged"])
        
      mean_wt_pred$GA <- GA_num / denominator
      mean_wt_pred$CC <- cc_num / denominator 
      mean_wt_pred$Dead <- dead_num / denominator
      mean_wt_pred$Discharged <- dis_num / denominator
      mean_wt_pred$patient_group <- patient_group[j]
      mean_wt_pred$ICD <- current_icd
      mean_wt_pred$age <- current_age
      mean_wt_pred$WT <- "seven"
      mean_wt_pred$group_size <- denominator
      mean_wt_pred <- mean_wt_pred %>% select(GA,CC,Dead,Discharged,patient_group,ICD,age,WT, group_size)
        
      toc()
      
      
      whole_df[j,] <- mean_wt_pred 
      
    }
    
    
    toc()
    whole_df$week <- append_week
    
    tot_whole_df <- dplyr::bind_rows(tot_whole_df, whole_df)
    
  }
  
  
  
  return(tot_whole_df)
  
}




regression_cluster_set_up <- function(patient_group, hes_data, forecast_length = 52, forecast_start = "2012-03-05",
                                      start_date = "2009-01-01", time_trend = TRUE, month_trend = TRUE, wt_variable = "linear",
                                      week_num = 1, failure_function_run = TRUE, half_week = TRUE, num_cores = 8){

  ## Takes in data and the patient group as either "elective" or "emergency"
  print("Listing the icds now") 
  tic("Total run through:")
  
  whole_graph_df <- NULL
  coefdf_tot <- NULL
  failure_df <- NULL
  actual_out_df <- NULL
  
  if(failure_function_run){
    ## use the try catch function to run through all the icds for the failure function to
    ## run on 
    print("Narrowing down to cohort 1 & 2")
    tic("Narrow down HES")
    cohort_12 <- hes_data[hes_data$cohort != 3,]
    toc()
    cohort_12 <- make_elective_cohort_variable(cohort_12)
    elective_groupings <- unique(cohort_12$elective_icd)
    print(elective_groupings)
    
    for(k in 1:length(elective_groupings)){
    
      base_dir <- "./"
      current_icd <- elective_groupings[k]
      cohort12_icd <- cohort_12[cohort_12$elective_icd == current_icd,]
      
      failure_csv <- paste(base_dir, "survival/surival_cohort_12_ICD",current_icd,".csv", sep = "")
      failure_pdf <- paste(base_dir, "survival/surival_cohort_12_ICD",current_icd,".pdf", sep = "")
      
      ## try catch for survival analysis 
      ## failure function 
      
      
      cohort12_failure_func <- failure_func_try(cohort23 =  cohort12_icd,
                                                csv_survival_name = failure_csv,
                                                pdf_to_export = failure_pdf, current_ICD = current_icd)   
      if(length(cohort12_failure_func) != 0){
        failure_df <- rbind.data.frame(failure_df, cohort12_failure_func[[2]])
      }
    }
    
    
  }
  
  
  
    
  if(patient_group == "elective"){
    tic("Narrowing down the data")
    
    hes_data <- hes_data[hes_data$cohort == 1,]
    hes_data <- hes_data[,c("ICD","agegrp_v3","cc","WaitingTime","cc_LoS",
                            "GA_LoS","ga_transitions","admidate_MDY","admidate_YYYY",
                            "admidate_week","cc_transitions","cc_start_flg")]
    
    toc()
    icd_list <- unique(hes_data$ICD)
    icds_to_run <- rep(icd_list, each = 3)
    ages_to_run <- rep(c(1,2,3),length(icd_list))
    patient_groups_to_run <- paste(icds_to_run,"ga",ages_to_run, sep = "-")
    patient_groups_to_run_cc <- paste(icds_to_run,"cc",ages_to_run, sep = "-")
    tot_patient_groups <- c(patient_groups_to_run, patient_groups_to_run_cc)
    tic("Setting up elective cluster")
    regression_cluster <- snow::makeCluster(spec = num_cores, outfile = "./crr_cluster_log.txt", type = "SOCK")
    regression_input <- snow::clusterSplit(regression_cluster, tot_patient_groups)
    toc()
    snow::clusterExport(regression_cluster, "elective_regression2")
    print(paste("Copying over data for ICD"))
    copy_start <- Sys.time()
    snow::clusterExport(regression_cluster, "hes_data", envir = environment())
    snow::clusterExport(regression_cluster, "week_num", envir = environment())
    snow::clusterExport(regression_cluster, "half_week", envir = environment())
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
                                            fun = elective_regression2,
                                            hes_data_orig = hes_data,
                                            week_num = week_num,
                                            half_week = half_week)
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
    
    if(half_week){
      week_num <- seq(0.5,length.out = length(week_num), by = 1)
    }
    
    for(j in week_num){
      current_whole_graph_df <- whole_graph_df[whole_graph_df$week == j,]
      current_out <- reg_data_summariser2(current_whole_graph_df, admi_type = "N",week_num = j)
      
      actual_out_df <- dplyr::bind_rows(actual_out_df, current_out)
      
    }
    
    
  }else if(patient_group == "emergency"){
    tic("Narrowing down the data")
    hes_data <- hes_data[hes_data$cohort != 1,]
    
    hes_data <- hes_data[,c("ICD","agegrp_v3","cc","cc_LoS",
                            "GA_LoS","ga_transitions","admidate_MDY","admidate_YYYY",
                            "admidate_week","cc_transitions", "cc_start_flg")]
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
    regression_cluster <- snow::makeCluster(spec = num_cores/2, outfile = "./crr_cluster_log.txt", type = "SOCK")
    
    regression_input <- snow::clusterSplit(regression_cluster, tot_patient_groups)
    toc()
    snow::clusterExport(regression_cluster, "emergency_regression_cluster")
    print(paste("Copying over data for ICD"))
    copy_start <- Sys.time()
    snow::clusterExport(regression_cluster, "hes_data", envir = environment())
    snow::clusterExport(regression_cluster, "wt_variable", envir = environment())
    snow::clusterExport(regression_cluster, "half_week", envir = environment())
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
                                                   fun = emergency_regression2,
                                                   hes_data_orig = hes_data,
                                                   week_num = week_num,
                                                   half_week = half_week)
    jobs_end <- Sys.time()
    print(jobs_end - jobs_start)
    snow::clusterEvalQ(regression_cluster, rm(hes_data))
    snow::clusterEvalQ(regression_cluster, rm(hes_data_orig))
    snow::clusterEvalQ(regression_cluster, gc())
    stopCluster(regression_cluster)
    
    whole_graph_df <- dplyr::bind_rows(regression_jobs_parallel)
    
    actual_out_df <- NULL
    if(half_week){
      week_num <- seq(0.5,length.out = length(week_num), by = 1)
    }
    
    
    for(j in week_num){
      current_whole_graph_df <- whole_graph_df[whole_graph_df$week == j,]
      current_out <- reg_data_summariser2(current_whole_graph_df, admi_type = "E",week_num = j)
      
      actual_out_df <- dplyr::bind_rows(actual_out_df, current_out)
      
    }
    
  }
  
  toc()  
  return(list(actual_out_df, coefdf_tot, failure_df))   
  
  
  
}





multi_graph_plotter <- function(reg_res, out_file, trend_types = "month and time", patient_admi = "elective"){
  
  pdf(file = out_file, paper = "A4r", width = 10, height = 7)  
  
  patient_groupings <- unique(reg_res$patient_group)
  if(length(which(is.na(patient_groupings))) > 0)
    patient_groupings <- patient_groupings[-which(is.na(patient_groupings))]
  
  for(k in 1:length(patient_groupings)){
    current_grouping <- patient_groupings[k]
    current_graph_df <- reg_res[reg_res$patient_group == current_grouping,]
    if(length(which(is.na(current_graph_df$patient_group))) > 0)
      current_graph_df <- current_graph_df[-which(is.na(current_graph_df$patient_group)),]
    
    
    graph_plot <- ggplot(data = current_graph_df, aes(x = reg_week, y = value, group = line_group)) +
      geom_line(aes(color = variable, linetype = WT)) + theme_bw() +
      ggtitle(paste(current_grouping,trend_types, patient_admi))
    
    print(graph_plot)
    
      
    
    
    
  }
  
  dev.off()
  
}


memory_less_observations <- function(hes_data){
  start_time <- Sys.time()
  ## Need to output the following:
  ## - average LoS for G&A and CC for all patient groups
  ## - Overall admissions in the data for CC and G&A all patient groups
  ## - Overall deaths in the data for CC and G&A for all patient groups 
  ## - Overall discharges in the data for CC and G&A for all patient groups 
  
  ## Narrow HES 
  print("Narrowing HES data")
  tic()
  hes_data <- hes_data[,c("GA_LoS","cc_LoS","ga_transitions","cc_transitions","cc",
                          "ICD","cohort","agegrp_v3")]
  
  hes_data$admi <- ifelse(hes_data$cohort == 1,1,3)
  hes_data$ga_2 <- hes_data$ga_transitions
  hes_data$cc_2 <- hes_data$cc_transitions
  hes_data$ga_2[which(is.na(hes_data$ga_2))] <- 0
  hes_data$cc_2[which(is.na(hes_data$cc_2))] <- 0
  
  hes_data$death <- ifelse(hes_data$ga_2 == 3 | hes_data$cc_2 == 3, 1,0)
  hes_data$discharges <- ifelse(hes_data$ga_2 == 1 | hes_data$cc_2 == 1, 1, 0)
  hes_data$one <- 1
  toc()
  ## Los Averages first 
  print("Setting up separate groupings")
  tic()
  emerg_rows <- hes_data[hes_data$admi == 3,]
  elec_rows <- hes_data[hes_data$admi == 1,]
  
  emerg_cc <- emerg_rows[emerg_rows$cc == 1,]
  elec_cc <- elec_rows[elec_rows$cc == 1,]
  emerg_cc <- emerg_cc[-which(is.na(emerg_cc$cc_transitions)),]
  elec_cc <- elec_cc[-which(is.na(elec_cc$cc_transitions)),]
  
  emerg_ga <- emerg_rows[-which(is.na(emerg_rows$ga_transitions)),]
  elec_ga <- elec_rows[-which(is.na(elec_rows$ga_transitions)),]
  
  emerg_icds <- unique(emerg_rows$ICD)
  elec_icds <- unique(elec_rows$ICD)
  
  emerg_icds_seq <- rep(emerg_icds, each = 3*2)
  elec_icds_seq <- rep(elec_icds, each = 3*2)
  ages <- rep(c(1,2,3), 2 *(length(emerg_icds) + length(elec_icds)))
  ward <- rep(rep(c("GA","CC"), each = 3), (length(emerg_icds) + length(elec_icds)))
  admis <- c(rep("emergency",length(emerg_icds) * 3 * 2), rep("elective",length(elec_icds) * 2 * 3))
  tot_icds <- c(emerg_icds_seq, elec_icds_seq)
  toc()
  print("LoS calculations")
  tic()
  LoS_df <- data.frame(matrix(data = 0, ncol = 4, nrow = length(tot_icds)))
  colnames(LoS_df) <- c("ICD","agegrp_v3","admi","ward")
  LoS_df[,1] <- tot_icds
  LoS_df[,2] <- ages
  LoS_df[,3] <- admis
  LoS_df[,4] <- ward
  
  emerg_ga_Los <- aggregate(GA_LoS ~ ICD + agegrp_v3, emerg_ga, FUN = mean, na.rm = TRUE)
  colnames(emerg_ga_Los)[3] <- "LoS"
  emerg_ga_Los$admi <- "emergency"
  emerg_ga_Los$ward <- "GA"
  
  elec_ga_Los <- aggregate(GA_LoS ~ ICD + agegrp_v3, elec_ga, FUN = mean, na.rm = TRUE)
  colnames(elec_ga_Los)[3] <- "LoS"
  elec_ga_Los$admi <-"elective"
  elec_ga_Los$ward <- "GA"
  
  emerg_cc_Los <- aggregate(cc_LoS ~ ICD + agegrp_v3, emerg_cc, FUN = mean, na.rm = TRUE)
  colnames(emerg_cc_Los)[3] <- "LoS"
  emerg_cc_Los$admi <- "emergency"
  emerg_cc_Los$ward <- "CC"
  
  elec_cc_Los <- aggregate(cc_LoS ~ ICD + agegrp_v3, elec_cc, FUN = mean, na.rm = TRUE)
  colnames(elec_cc_Los)[3] <- "LoS"
  elec_cc_Los$admi <- "elective"
  elec_cc_Los$ward <- "CC"
  
  tot_df <- dplyr::bind_rows(emerg_ga_Los, emerg_cc_Los, elec_ga_Los, elec_cc_Los)
  LoS_df <- dplyr::left_join(LoS_df, tot_df)
  toc()
  ## admissions ##
  print("Admission calculations")
  tic()
  admi_df <- data.frame(matrix(data = 0, ncol = 4, nrow = length(tot_icds)))
  colnames(admi_df) <- c("ICD","agegrp_v3","admi","ward")
  admi_df[,1] <- tot_icds
  admi_df[,2] <- ages
  admi_df[,3] <- admis
  admi_df[,4] <- ward
  
  emerg_ga_admi <- aggregate(one ~ ICD + agegrp_v3, emerg_ga, FUN = sum)
  colnames(emerg_ga_admi)[3] <- "Admissions"
  emerg_ga_admi$admi <- "emergency"
  emerg_ga_admi$ward <- "GA"
  
  elec_ga_admi <- aggregate(one ~ ICD + agegrp_v3, elec_ga, FUN = sum)
  colnames(elec_ga_admi)[3] <- "Admissions"
  elec_ga_admi$admi <-"elective"
  elec_ga_admi$ward <- "GA"
  
  emerg_cc_admi <- aggregate(one ~ ICD + agegrp_v3, emerg_cc,FUN = sum)
  colnames(emerg_cc_admi)[3] <- "Admissions"
  emerg_cc_admi$admi <- "emergency"
  emerg_cc_admi$ward <- "CC"
  
  elec_cc_admi <- aggregate(one ~ ICD + agegrp_v3, elec_cc, FUN = sum)
  colnames(elec_cc_admi)[3] <- "Admissions"
  elec_cc_admi$admi <- "elective"
  elec_cc_admi$ward <- "CC"
  
  tot_df <- dplyr::bind_rows(emerg_ga_admi, emerg_cc_admi, elec_ga_admi, elec_cc_admi)
  admi_df <- dplyr::left_join(admi_df, tot_df)
  toc()
  ## deaths ##
  
  print("Deaths calculations")
  tic()
  deaths_df <- data.frame(matrix(data = 0, ncol = 4, nrow = length(tot_icds)))
  colnames(deaths_df) <- c("ICD","agegrp_v3","admi","ward")
  deaths_df[,1] <- tot_icds
  deaths_df[,2] <- ages
  deaths_df[,3] <- admis
  deaths_df[,4] <- ward
  
  emerg_ga_def <- aggregate(death ~ ICD + agegrp_v3, emerg_ga, FUN = sum, na.rm = TRUE)
  colnames(emerg_ga_def)[3] <- "death"
  emerg_ga_def$admi <- "emergency"
  emerg_ga_def$ward <- "GA"
  
  elec_ga_def <- aggregate(death ~ ICD + agegrp_v3, elec_ga, FUN = sum, na.rm = TRUE)
  colnames(elec_ga_def)[3] <- "death"
  elec_ga_def$admi <-"elective"
  elec_ga_def$ward <- "GA"
  
  emerg_cc_def <- aggregate(death ~ ICD + agegrp_v3, emerg_cc,FUN = sum, na.rm = TRUE)
  colnames(emerg_cc_def)[3] <- "death"
  emerg_cc_def$admi <- "emergency"
  emerg_cc_def$ward <- "CC"
  
  elec_cc_def <- aggregate(death ~ ICD + agegrp_v3, elec_cc, FUN = sum, na.rm = TRUE)
  colnames(elec_cc_admi)[3] <- "death"
  elec_cc_def$admi <- "elective"
  elec_cc_def$ward <- "CC"
  
  tot_df <- dplyr::bind_rows(emerg_ga_def, emerg_cc_def, elec_ga_def, elec_cc_def)
  deaths_df <- dplyr::left_join(deaths_df, tot_df)
  toc()
  ## Discharges
  print("Discharges now")
  tic()
  dis_df <- data.frame(matrix(data = 0, ncol = 4, nrow = length(tot_icds)))
  colnames(dis_df) <- c("ICD","agegrp_v3","admi","ward")
  dis_df[,1] <- tot_icds
  dis_df[,2] <- ages
  dis_df[,3] <- admis
  dis_df[,4] <- ward
  
  emerg_ga_dis <- aggregate(discharges ~ ICD + agegrp_v3, emerg_ga, FUN = sum, na.rm = TRUE)
  colnames(emerg_ga_dis)[3] <- "discharges"
  emerg_ga_dis$admi <- "emergency"
  emerg_ga_dis$ward <- "GA"
  
  elec_ga_dis <- aggregate(discharges ~ ICD + agegrp_v3, elec_ga, FUN = sum, na.rm = TRUE)
  colnames(elec_ga_dis)[3] <- "discharges"
  elec_ga_dis$admi <-"elective"
  elec_ga_dis$ward <- "GA"
  
  emerg_cc_dis <- aggregate(discharges ~ ICD + agegrp_v3, emerg_cc,FUN = sum, na.rm = TRUE)
  colnames(emerg_cc_dis)[3] <- "discharges"
  emerg_cc_dis$admi <- "emergency"
  emerg_cc_dis$ward <- "CC"
  
  elec_cc_dis <- aggregate(discharges ~ ICD + agegrp_v3, elec_cc, FUN = sum, na.rm = TRUE)
  colnames(elec_cc_dis)[3] <- "discharges"
  elec_cc_dis$admi <- "elective"
  elec_cc_dis$ward <- "CC"
  
  tot_df <- dplyr::bind_rows(emerg_ga_dis, emerg_cc_dis, elec_ga_dis, elec_cc_dis)
  dis_df <- dplyr::left_join(dis_df, tot_df)
  toc()
  
  end_time <- Sys.time()
  
  print("Total time spent:")
  print(end_time - start_time)
  return(list(LoS_df, admi_df, deaths_df, dis_df))
  
}



failure_func_try <- function(cohort23, csv_survival_name, pdf_to_export, current_ICD){
  
  out <- tryCatch({
    message("Trying out the failure function")
    
    failure_function(cohort23 = cohort23, csv_survial_name = csv_survival_name,
                     pdf_to_export = pdf_to_export, 
                     current_ICD = current_ICD)   
  },
  error = function(cond){
    message(paste("Cohort 1 to 2 failure function not working for ICD:", current_ICD))
    message(cond)
    return(NULL)
  },
  warning = function(cond){
    message(paste("Following warning message for ICD:", current_ICD))
    message(cond)
    return(out)
    
  },
  finally = {
    message(paste("Processed cohort 1 to 2 for ICD:", current_ICD))
  }
  
  
  )
  
  
  return(out)
}

failure_function <- function(cohort23, csv_survial_name,
                             pdf_to_export = "electives_to_emergencies_survival_plot",
                             current_ICD){
  ## Need to double check whether Waiting time is in there!
  start_time <- Sys.time()
  
  print("Failure function for electives to emergencies")
  if(!("WaitingTime" %in% colnames(cohort23))){
    cohort23 <- make_wt_variable(cohort23)
    WT_colname <- which(colnames(cohort23) == "WT")
    colnames(cohort23)[WT_colname] <- "WaitingTime"
  }
  
  # 1.1 Fit survival function
  
  survival_df <- NULL
  ages <- unique(cohort23$agegrp_v3)
  
  for(k in 1:length(ages)){
    current_age_dat <- cohort23[cohort23$agegrp_v3 == ages[k],]
    surv_object <- Surv(time = current_age_dat$WaitingTime, event = current_age_dat$Elective2Emergency)
    e2e <- survfit(surv_object ~ 1, data = current_age_dat, type = "kaplan-meier")
    
    age_rows <- cbind.data.frame(e2e$time, e2e$surv, rep(ages[k], length(e2e$surv)), e2e$cumhaz,
                                 e2e$std.err)
    colnames(age_rows) <- c("time","survival","age_group", "cumulative_hazard", "std_err_hazard")
    survival_df <- rbind.data.frame(survival_df, age_rows)
    colnames(survival_df) <- c("time","survival","age_group", "cumulative_hazard", "std_err_hazard")
    
  }
  
  survival_df$ICD <- current_ICD
  # 1.2 Plot survival function
  
  graph_df <- survival_df
  graph_df$age_group <- as.character(graph_df$age_group)
  
  survival_plot_2 <- ggplot(data = graph_df, aes(x = time, y = survival, group = age_group, colour = age_group)) +
    geom_line() + xlab("Waiting time (days)") + ylab("Survival probability") + ylim(c(0,1)) + labs(colour = "Age group") +
    ggtitle("Survival function")
  
  
  # 1.3 Plot failure function
  cumhaz_plot_2 <- ggplot(data = graph_df, aes(x = time, y = cumulative_hazard, group = age_group, colour = age_group)) +
    geom_line() + xlab("Waiting time (days)") + ylab("Failure probability") + ylim(c(0,1)) + labs(colour = "Age group") +
    ggtitle("Failure function (cumulative hazard)")
  
  
  # 1.4 Put failure function and standard errors into a df --> export to csv
  #### This needs to be in the correct dir, from the setwd from above ####
  
  pdf(file = pdf_to_export, paper = "a4r")
  
  print(survival_plot_2)
  print(cumhaz_plot_2)
  
  dev.off()
  
  write.csv(survival_df,
            file = csv_survial_name,
            row.names = FALSE)
  
  ## get the 7 day multiples ##
  
  out_df <- data.frame(matrix(ncol = 7, nrow = 3))
  colnames(out_df) <- c("age","ICD","day_7","mean_7","median_7", "num_weeks", "num_switches")
  out_df$age <- ages
  out_df$ICD <- current_ICD

  for(k in 1:length(ages)){
    
    ## get day 7 vals first 
    day_7s <- survival_df[survival_df$time == 7 &
                            survival_df$age_group == out_df$age[k],]
    if(nrow(day_7s) == 1){
      out_df$day_7[k] <- day_7s$cumulative_hazard
    }else{
      day_7s <- survival_df[survival_df$time < 7 &
                              survival_df$age_group == out_df$age[k],]
      if(nrow(day_7s) == 0){
        out_df$day_7[k] <- 0
      }else{
        out_df$day_7[k] <- max(day_7s$cumulative_hazard)
      }
      
    }
    ## now multiples 
    seven_multiples <- c(70, max(survival_df$time))
    seven_mults <- seq(7, max(survival_df$time), by = 7)
    
    actual_vec <- out_df$day_7[k]
    multi_df <- out_df$day_7[k]
    age_surve <- survival_df[survival_df$age_group == out_df$age[k],]
    
    out_df$num_weeks[k] <- length(seven_mults)
    current_age_dat <- cohort23[cohort23$agegrp_v3 == ages[k] & cohort23$WaitingTime <= min(seven_multiples),]
    out_df$num_switches[k] <- nrow(current_age_dat[current_age_dat$cohort == 2,])
    
    for(j in 2:length(seven_mults)){
      current_mult_vals <- age_surve[age_surve$time > seven_mults[j-1] &
                                       age_surve$time <= seven_mults[j],]
      if(nrow(current_mult_vals) == 0){
        multi_df <- append(multi_df, 0)
        actual_vec <- append(actual_vec, actual_vec[length(actual_vec)])
      }else{
        diff_val <- max(current_mult_vals$cumulative_hazard) - actual_vec[length(actual_vec)]
        multi_df <- append(multi_df, diff_val)
        actual_vec <- append(actual_vec, max(current_mult_vals$cumulative_hazard))
      }
      
    }
    out_df$mean_7[k] <- mean(actual_vec)
    out_df$median_7[k] <- median(actual_vec)

  }
  
  end_time <- Sys.time()
  print(end_time - start_time)
  return(list(survival_df, out_df))
  
}


reg_data_summariser <- function(reg_data,admi_type, week_num){

  
  icds_in_play <- unique(reg_data[reg_data$WT == "seven","ICD"])
  if(length(which(is.na(icds_in_play))) > 0)
    icds_in_play <- icds_in_play[-which(is.na(icds_in_play))]
    
  row_nums <- length(icds_in_play) * 8 * 3
  
  icd_names <- icds_in_play
  for(k in 1:length(icd_names)){
    if(nchar(icd_names[k]) == 1)
      icd_names[k] <- paste("0",icd_names[k], sep = "")
  }
  
  
  
  out_df <- data.frame(matrix(ncol = 7, nrow = row_nums))
  colnames(out_df) <- c("a", "p",	"s",	"sbar",	"pi_y",	"coeff","variance")
  out_df$a <- admi_type
  out_df$p <- paste("ICD",rep(icd_names, each = 24),"_AGE",rep(rep(c(1,2,3),each = 8), length(icd_names)), sep = "")
  out_df$s <- rep(rep(c("G","C"), each = 4), 3*length(icd_names))
  out_df$sbar <- rep(c("H","C","D","G","H","G","D","C"), 3*length(icd_names))
  
  if(length(which(is.na(reg_data$WT))) > 0)
    reg_data <- reg_data[-which(is.na(reg_data$WT))]
  
  if(admi_type == "N"){
    icd_age_groups <- paste(rep(icds_in_play,each = 3),rep(c(1,2,3), length(icds_in_play)), sep = "-")
    current_row_set <- 1
    
    
    for(k in 1:length(icd_age_groups)){
      current_icd <- as.integer(str_split_fixed(icd_age_groups[k],"-",2)[1])
      current_age <- as.integer(str_split_fixed(icd_age_groups[k],"-",2)[2])
      ga_grouping <- paste(current_icd,"ga",current_age,sep = "-")
      cc_grouping <- paste(current_icd,"cc",current_age,sep = "-")
      
      current_df <- reg_data[reg_data$ICD == current_icd & reg_data$age == current_age, ]
      current_df <- current_df[current_df$WT == "seven",]
      if(length(which(is.na(current_df$ICD))) > 0)
        current_df <- current_df[-which(is.na(current_df$ICD)),]
      
      ga_df <- current_df[current_df$patient_group == ga_grouping,]
      cc_df <- current_df[current_df$patient_group == cc_grouping,]
      
      ga_ga <- ga_df[ga_df$variable == "GA","value"][1]
      ga_cc <- ga_df[ga_df$variable == "CC","value"][1]
      ga_dead <- ga_df[ga_df$variable == "Dead","value"][1]
      ga_dis <- ga_df[ga_df$variable == "Discharged","value"][1]
      ## cc
      cc_ga <- cc_df[cc_df$variable == "GA","value"][1]
      cc_cc <- cc_df[cc_df$variable == "CC","value"][1]
      cc_dead <- cc_df[cc_df$variable == "Dead","value"][1]
      cc_dis <- cc_df[cc_df$variable == "Discharged","value"][1]
      
      row_vals <- seq((current_row_set * 8) - 7,(current_row_set * 8))
      
      out_df$pi_y[row_vals] <- c(ga_dis, ga_cc, ga_dead, ga_ga,cc_dis, cc_ga, cc_dead, cc_cc) 
      
      print(out_df$p[row_vals[1]])
      print(paste(current_icd, current_age,sep = "--"))
      current_row_set <- current_row_set + 1
      
    }
    
    
    
    
    
  }else if(admi_type == "E"){
    icd_age_groups <- paste(rep(icds_in_play,each = 3),rep(c(1,2,3), length(icds_in_play)), sep = "-")
    current_row_set <- 1
    
    
    for(k in 1:length(icd_age_groups)){
      current_icd <- as.integer(str_split_fixed(icd_age_groups[k],"-",2)[1])
      current_age <- as.integer(str_split_fixed(icd_age_groups[k],"-",2)[2])
      ga_grouping <- paste(current_icd,"ga",current_age,sep = "-")
      cc_grouping <- paste(current_icd,"cc",current_age,sep = "-")
      
      current_df <- reg_data[reg_data$ICD == current_icd & reg_data$age == current_age, ]
      if(length(which(is.na(current_df$ICD))) > 0)
        current_df <- current_df[-which(is.na(current_df$ICD)),]
      
      current_df <- current_df[current_df$WT == "mean",]
      ga_df <- current_df[current_df$patient_group == ga_grouping,]
      cc_df <- current_df[current_df$patient_group == cc_grouping,]
      
      ga_ga <- ga_df[ga_df$variable == "GA","value"][1]
      ga_cc <- ga_df[ga_df$variable == "CC","value"][1]
      ga_dead <- ga_df[ga_df$variable == "Dead","value"][1]
      ga_dis <- ga_df[ga_df$variable == "Discharged","value"][1]
      ## cc
      cc_ga <- cc_df[cc_df$variable == "GA","value"][1]
      cc_cc <- cc_df[cc_df$variable == "CC","value"][1]
      cc_dead <- cc_df[cc_df$variable == "Dead","value"][1]
      cc_dis <- cc_df[cc_df$variable == "Discharged","value"][1]
      
      row_vals <- seq((current_row_set * 8) - 7,(current_row_set * 8))
      
      out_df$pi_y[row_vals] <- c(ga_dis, ga_cc, ga_dead, ga_ga,cc_dis, cc_ga, cc_dead, cc_cc) 
      
      print(out_df$p[row_vals[1]])
      print(paste(current_icd, current_age,sep = "--"))
      current_row_set <- current_row_set + 1
      
    }
    
    
  }
    
  out_df$week <- week_num
  return(out_df)
}

reg_data_summariser2 <- function(reg_data,admi_type, week_num){
  
  
  icds_in_play <- unique(reg_data[reg_data$WT == "seven","ICD"])
  if(length(which(is.na(icds_in_play))) > 0){
    icds_in_play <- icds_in_play[-which(is.na(icds_in_play))]
  }
  row_nums <- length(icds_in_play) * 8 * 3
  
  icd_names <- icds_in_play
  for(k in 1:length(icd_names)){
    if(nchar(icd_names[k]) == 1)
      icd_names[k] <- paste("0",icd_names[k], sep = "")
  }
  
  
  
  out_df <- data.frame(matrix(ncol = 8, nrow = row_nums))
  colnames(out_df) <- c("a", "p",	"s",	"sbar",	"pi_y",	"coeff","variance", "group_size")
  out_df$group_size <- 0
  out_df$a <- admi_type
  out_df$p <- paste("ICD",rep(icd_names, each = 24),"_AGE",rep(rep(c(1,2,3),each = 8), length(icd_names)), sep = "")
  out_df$s <- rep(rep(c("G","C"), each = 4), 3*length(icd_names))
  out_df$sbar <- rep(c("H","C","D","G","H","G","D","C"), 3*length(icd_names))
  
  if(length(which(is.na(reg_data$WT))) > 0){
    reg_data <- reg_data[-which(is.na(reg_data$WT))]
  }
  if(admi_type == "N"){
    icd_age_groups <- paste(rep(icds_in_play,each = 3),rep(c(1,2,3), length(icds_in_play)), sep = "-")
    current_row_set <- 1
    
    
    for(k in 1:length(icd_age_groups)){
      current_icd <- as.integer(str_split_fixed(icd_age_groups[k],"-",2)[1])
      current_age <- as.integer(str_split_fixed(icd_age_groups[k],"-",2)[2])
      ga_grouping <- paste(current_icd,"ga",current_age,sep = "-")
      cc_grouping <- paste(current_icd,"cc",current_age,sep = "-")
      
      current_df <- reg_data[reg_data$ICD == current_icd & reg_data$age == current_age, ]
      current_df <- current_df[current_df$WT == "seven",]
      if(length(which(is.na(current_df$ICD))) > 0){
        current_df <- current_df[-which(is.na(current_df$ICD)),]
      }
      
      ga_df <- current_df[current_df$patient_group == ga_grouping,]
      cc_df <- current_df[current_df$patient_group == cc_grouping,]
      
      ga_ga <- ga_df[1,"GA"]
      ga_cc <- ga_df[1,"CC"][1]
      ga_dead <- ga_df[1,"Dead"][1]
      ga_dis <- ga_df[1,"Discharged"][1]
      ## cc
      cc_ga <- cc_df[1,"GA"][1]
      cc_cc <- cc_df[1,"CC"][1]
      cc_dead <- cc_df[1,"Dead"][1]
      cc_dis <- cc_df[1,"Discharged"][1]
      
      row_vals <- seq((current_row_set * 8) - 7,(current_row_set * 8))
      out_df$group_size[row_vals[1:4]] <- ga_df[1, "group_size"]
      out_df$group_size[row_vals[5:8]] <- cc_df[1, "group_size"]
      out_df$pi_y[row_vals] <- c(ga_dis, ga_cc, ga_dead, ga_ga,cc_dis, cc_ga, cc_dead, cc_cc) 
      
      print(out_df$p[row_vals[1]])
      print(paste(current_icd, current_age,sep = "--"))
      current_row_set <- current_row_set + 1
      
    }
    
    
    
    
    
  }else if(admi_type == "E"){
    icd_age_groups <- paste(rep(icds_in_play,each = 3),rep(c(1,2,3), length(icds_in_play)), sep = "-")
    current_row_set <- 1
    
    
    for(k in 1:length(icd_age_groups)){
      current_icd <- as.integer(str_split_fixed(icd_age_groups[k],"-",2)[1])
      current_age <- as.integer(str_split_fixed(icd_age_groups[k],"-",2)[2])
      ga_grouping <- paste(current_icd,"ga",current_age,sep = "-")
      cc_grouping <- paste(current_icd,"cc",current_age,sep = "-")
      
      current_df <- reg_data[reg_data$ICD == current_icd & reg_data$age == current_age, ]
      if(length(which(is.na(current_df$ICD))) > 0){
        current_df <- current_df[-which(is.na(current_df$ICD)),]
      }
      
      current_df <- current_df[current_df$WT == "seven",]
      ga_df <- current_df[current_df$patient_group == ga_grouping,]
      cc_df <- current_df[current_df$patient_group == cc_grouping,]
      
      ga_ga <- ga_df[1,"GA"]
      ga_cc <- ga_df[1,"CC"][1]
      ga_dead <- ga_df[1,"Dead"][1]
      ga_dis <- ga_df[1,"Discharged"][1]
      ## cc
      cc_ga <- cc_df[1,"GA"][1]
      cc_cc <- cc_df[1,"CC"][1]
      cc_dead <- cc_df[1,"Dead"][1]
      cc_dis <- cc_df[1,"Discharged"][1]
      row_vals <- seq((current_row_set * 8) - 7,(current_row_set * 8))
      out_df$group_size[row_vals[1:4]] <- ga_df[1, "group_size"]
      out_df$group_size[row_vals[5:8]] <- cc_df[1, "group_size"]
      
      out_df$pi_y[row_vals] <- c(ga_dis, ga_cc, ga_dead, ga_ga,cc_dis, cc_ga, cc_dead, cc_cc) 
      
      print(out_df$p[row_vals[1]])
      print(paste(current_icd, current_age,sep = "--"))
      current_row_set <- current_row_set + 1
      
    }
    
    
  }
  
  out_df$week <- week_num
  return(out_df)
}












