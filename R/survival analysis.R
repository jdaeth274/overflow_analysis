# Title: Survival Analysis for OverFlow Deaths Project
# Last updated: 4 June 2020
###########################################################################################################################
# install necessary packages

###########################################################################################################################
# set working directory
###########################################################################################################################
###########################################################################################################################
# 1. Failure function for Electives to Emergencies (Cohort = 2 & 3, event = Elective2Emergency == 1, time = WaitingTime)

make_wt_variable <- function(elective_df){
  
  one_year_under <- which(as.integer(elective_df$elecdur) <= 365)
  
  elective_df$WT <- NA
  elective_df$WT[one_year_under] <- elective_df$elecdur[one_year_under]
  
  return(elective_df)
  
}




failure_function <- function(cohort23,csv_survial_name,
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
  
  
  
  survival_plot_2 <- ggplot(data = survival_df, aes(x = time, y = survival, group = age_group, colour = age_group)) +
          geom_line() + xlab("Waiting time (days)") + ylab("Survival probability") + ylim(c(0,1)) + labs(colour = "Age group") +
          ggtitle("Survival function")
  
  
  # 1.3 Plot failure function
  cumhaz_plot_2 <- ggplot(data = survival_df, aes(x = time, y = cumulative_hazard, group = age_group, colour = age_group)) +
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
  
  out_df <- data.frame(matrix(ncol = 5, nrow = 3))
  colnames(out_df) <- c("age","ICD","day_7","mean_7","median_7")
  out_df$age <- ages
  out_df$ICD <- current_ICD
  
  for(k in 1:length(ages)){
    
    ## get day 7 vals first 
    day_7s <- survival_df[survival_df$time == 7 &
                            survival_df$age == out_df$age[k],]
    if(nrow(day_7s) == 1){
      out_df$day_7[k] <- day_7s$cumulative_hazard
    }else{
      day_7s <- survival_df[survival_df$time < 7 &
                                  survival_df$age == out_df$age[k],]
      if(nrow(day_7s) == 0){
        out_df$day_7[k] <- 0
      }else{
        out_df$day_7[k] <- max(day_7s$cumulative_hazard)
      }
      
    }
    ## now multiples 
    seven_mults <- seq(7, max(survival_df$time), by = 7)
    
    actual_vec <- out_df$day_7[k]
    multi_df <- out_df$day_7[k]
    age_surve <- survival_df[survival_df$age == out_df$age[k],]
    
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
###########################################################################################################################
###########################################################################################################################
# 2. Competing Risks - Cohort 1 (True Emergencies)
# use "cmprsk" package since no WT covariate
# NEED TO CODE LOOPING: by ICD, agegrp, transition number, vector of GA_LoS & ga_transitions vs. vector of cc_LoS & cc_transitions)

# 2.1 GA transitions (Cohort = 1, event = ga_transitions == 1, 2, or 3, time = GA_LoS)

cuminc_list_to_df <- function(CI_dataset, group_names = c("Recovered, 1", "Recovered, 2", "Recovered, 3",
                                                          "CC, 1", "CC, 2", "CC, 3","Died, 1",
                                                          "Died, 2","Died, 3"),
                              ICD, unique_only = FALSE, missing_age = 0){
  ## This function takes the output of the competing risks model for emergency admissions and 
  ## then plots our cumulative incidence curves and then returns (in LOS) the transitions probs
  ## at day 7 and multiples of 7
  
  tests_pos <- grep("Tests", names(CI_dataset), ignore.case = TRUE)
  
  cuminc_df <- NULL
  
  ## Just check we've got all the age groups we need 
  
  
  
  for(k in 1:(tests_pos - 1)){
    current_title <- group_names[k]
    current_age <- stringr::str_split_fixed(current_title, ",", 2)[2]
    current_age <- gsub("\\s","", current_age)
    current_transition <- stringr::str_split_fixed(current_title, ",", 2)[1]
    current_list_out <- CI_dataset[[k]]
    
    new_rowz <- cbind.data.frame(current_list_out$time, current_list_out$est,
                                 rep(group_names[k], length(current_list_out$time)),
                                 rep(current_age, length(current_list_out$time)),
                                 rep(current_transition, length(current_list_out$time)),
                                 current_list_out$var, rep(ICD, length(current_list_out$time)))
    colnames(new_rowz) <- c("time","CIF",'age_tran', "age","transition","variance", "ICD_num")
    
    
    
    cuminc_df <- rbind.data.frame(cuminc_df, new_rowz)
    
  }
  

  
  max_time <- max(cuminc_df$time)
  multiplesofseven <- seq(7,max_time, by = 7)
  seventh_day_timies <- timepoints(CI_dataset,multiplesofseven) 
  
  
  
  
  
  
  ## Now we get the summarised CIF for 7 day Length of stay numbers for each age and 
  ## transitions group.
  
  ## First we transform the cuminc data to get the cumulative prob at each time point, 
  ## so no multiples.
  
  los_data <- NULL
  trans_age_groups <- plyr::count(cuminc_df$age_tran)
  
  for(k in 1:nrow(trans_age_groups)){
    current_group_data <- cuminc_df[cuminc_df$age_tran == trans_age_groups[k, 1],]
    ordered_current_age_data <- current_group_data[base::order(current_group_data$CIF, decreasing = TRUE),]
    current_data_max <- ordered_current_age_data[!base::duplicated(ordered_current_age_data$time),]
    current_data_max <- current_data_max[base::order(current_data_max$time ),]
    los_data <- rbind.data.frame(los_data, current_data_max)
    
  }
  
  ## Sometimes the data will have no transitions, we'll just set these to 0
  
  if(nrow(trans_age_groups) != length(group_names)){
    missing_age_trans <- group_names[which(!(group_names %in% trans_age_groups[,1]))]
    missing_times <- rep(0, length(missing_age_trans))
    missing_CIF <- rep(0, length(missing_age_trans))
    missing_age <- str_split_fixed(missing_age_trans, ",\\s",2)[,2]
    missing_transition <- str_split_fixed(missing_age_trans, ",\\s",2)[,1]
    missing_variance <- rep(0, length(missing_age_trans)) 
    missing_ICD <- rep(ICD, length(missing_age_trans))
    
    missing_rows <- cbind.data.frame(missing_times, missing_CIF, missing_age_trans,
                                     missing_age, missing_transition, missing_variance,
                                     missing_ICD)
    colnames(missing_rows) <- c("time","CIF",'age_tran', "age","transition","variance", "ICD_num")
    los_data <- rbind.data.frame(los_data, missing_rows)
  }
  
  
  
  age_groupings <- plyr::count(los_data$age)
  transitions <- plyr::count(los_data$transition)
  los_df <- data.frame(matrix(ncol = 5, nrow = 3*nrow(age_groupings)))
  
  
  for(k in 1:nrow(age_groupings)){
    age_rows <- c((k*3) -2, (k*3) -1, (k*3))
    current_age_data <- los_data[los_data$age == age_groupings[k,1],]
    
    if(max(current_age_data$time) > 7)
      multi_7 <- seq(7,max(current_age_data$time), by = 7)
    
    
    colnames(los_df) <- c("Transition","Day_7","mean_7_multiple","median_7_multiple", "age")
    los_df$age[age_rows] <- rep(as.character(age_groupings[k,1]), 3)
    
    los_df$Transition[age_rows] <- as.character(transitions[,1])
    
    current_7_data <- current_age_data[current_age_data$time == 7,]
    
    
    if(nrow(current_7_data) != 3){
      missing_trans <- as.character(transitions[which(!(transitions[,1] %in% current_7_data$transition)) ,1])
      for( j in 1:length(missing_trans)){
        current_missing_dat_7 <- current_age_data[current_age_data$transition == missing_trans[j] &
                                                    current_age_data$time <= 7, ]
        if(nrow(current_missing_dat_7) < 1){
          ## IF none at 7 or less days, take first and assume linear change in prob
          ## such that at CIF 7 is CIF x / (time_x/7)
          new_rows <- current_age_data[current_age_data$transition == missing_trans[j],]
          if(nrow(new_rows) < 1){
            new_row <- cuminc_df[1,]
            new_row$time <- 7
            new_row$CIF <- NaN
            new_row$age_tran <- "x"
            new_row$transition <- missing_trans[j]
            new_row$variance <- NaN
            
          }else{
            day_diff <- new_row$time[1] / 7
            new_row <- new_rows[1,]
            new_row$CIF <- new_row$CIF / day_diff
          }
        }else{
          new_rows <- current_missing_dat_7[order(current_missing_dat_7$time, decreasing = TRUE),]
          new_row <- new_rows[1,]
          
          
        }
        current_7_data <- rbind.data.frame(current_7_data, new_row)
      }
    }
    
    
    los_df$Day_7[age_rows[1]] <- current_7_data[current_7_data$transition == transitions[1,1],2]
    los_df$Day_7[age_rows[2]] <- current_7_data[current_7_data$transition == transitions[2,1],2]
    los_df$Day_7[age_rows[3]] <- current_7_data[current_7_data$transition == transitions[3,1],2]
    
    trans_nights <- c(k, k+nrow(age_groupings), k+(2*nrow(age_groupings)))
    average_df <- as.data.frame(seventh_day_timies$est)
    na_1 <- which(is.na(average_df[trans_nights,1]))
    if(length(na_1) > 0)
      average_df[trans_nights[na_1],1] <- los_df$Day_7[age_rows[na_1]] 
    
    
    
    
    
    
    tran_1 <- as.numeric(average_df[trans_nights[1],])
    tran_2 <- as.numeric(average_df[trans_nights[2],])
    tran_3 <- as.numeric(average_df[trans_nights[3],])
    
    
    
    
    los_df$mean_7_multiple[age_rows[1]] <- base::mean(tran_1,na.rm = TRUE)
    los_df$mean_7_multiple[age_rows[2]] <- base::mean(tran_2,na.rm = TRUE)
    los_df$mean_7_multiple[age_rows[3]] <- base::mean(tran_3,na.rm = TRUE)
    los_df$median_7_multiple[age_rows[1]] <- median(tran_1, na.rm = TRUE)
    los_df$median_7_multiple[age_rows[2]] <- median(tran_2, na.rm = TRUE)
    los_df$median_7_multiple[age_rows[3]] <- median(tran_3, na.rm = TRUE)
    
    
  }
  
  ## cuminc_df with duplicates
  ## los_df 7th day values
  ## los_data is cuminc_df without steps for plotting. 
  
  los_df$ICD <- ICD
  
  return(list(cuminc_df, los_df, los_data))
}




cohort_3_comp_risks <- function(cohort1, ICD_group, results_pdf){
  library(cmprsk)


  # 2.1.1 Give frequency of different GA transitions
  
  
  # 2.1.2 Set survival datatset (grouped by age groups)
  print("On GA transitions")
  ga_start <- Sys.time()
  ages <- unique(cohort1$agegrp_v3)
  CI.byagegrp <- cuminc(ftime = cohort1$GA_LoS, fstatus = cohort1$ga_transitions, group = cohort1$agegrp_v3)
  
  # THIS GIVES US THE CIF AT THE SPECIFIC TIME POINT (I DON'T KNOW HOW TO GET A SPECIFIC RISK/AGE-GROUP, IT JUST RETURNS ALL OF THEM)
  if(length(unique(cohort1$agegrp_v3)) == 3){
    cuminc_df_list <- cuminc_list_to_df(CI.byagegrp, ICD = ICD_group, unique_only = TRUE) ## JD to export to CSV
  }else{
    ages_present <- sort(unique(cohort1$agegrp_v3))
    transitions <- rep(c("Recovered,","CC,","Died,"), each = length(ages_present))
    ages_to_use <- rep(ages_present, 3)
    transitions_to_use <- paste(transitions, ages_to_use, sep = " ")
    cuminc_df_list <- cuminc_list_to_df(CI.byagegrp, ICD = ICD_group, unique_only = TRUE,
                                        group_names = transitions_to_use) ## JD to export to CSV
    
  }
  los_df <- cuminc_df_list[[2]]
  cuminc_df <- cuminc_df_list[[1]]
  
  # 2.1.4 Plot CIFs
  cuminc_plot <- ggplot(data = cuminc_df, aes(x = time, y = CIF, group = age_tran)) +
    geom_line(aes(linetype = transition, colour = age)) + xlab("Length of stay (days)") +
    ylab("Cumulative transition probability") +ggtitle("Cumualtive transition prob, across 3 transitions")
  
  
  ## Just recovered now ##
  
  recovered_df <- cuminc_df[grep("recovered", cuminc_df$transition, ignore.case = TRUE),]
  
  recovered_plot <- ggplot(data = recovered_df, aes(x = time, y = CIF, group = age_tran)) +
    geom_line(aes(colour = age)) + xlab("Length of stay (days)") +
    ylab("Cumulative Recovery probability") + ggtitle("Recovery probabilities")
  
  
  ## Just CC transition ##
  
  cc_df <- cuminc_df[grep("CC", cuminc_df$transition, ignore.case = TRUE),]
  
  cc_plot <- ggplot(data = cc_df, aes(x = time, y = CIF, group = age_tran)) +
    geom_line(aes(colour = age)) + xlab("Length of stay (days)") +
    ylab("Cumulative move to CC probability") + ggtitle("G&A > CC probabilities")
  
  
  ## Just Died ##
  
  died_df <- cuminc_df[grep("died", cuminc_df$transition, ignore.case = TRUE),]
  
  died_plot <- ggplot(data = died_df, aes(x = time, y = CIF, group = age_tran)) +
    geom_line(aes(colour = age)) + xlab("Length of stay (days)") +
    ylab("Cumulative Death probability") + ggtitle("Death probabilities")
  
  
  # 2.1.5 Test for significant differences in CIF curves between the age-groups at each competing risk (recover, cc, dying)
  CI.byagegrp$Tests
  
  end_time <- Sys.time()
  print(end_time - ga_start)
  print("Now on CC transitions")
  
  cc_start <- Sys.time()
  pdf(file = results_pdf, paper = "a4r")
  
  cuminc_plot
  recovered_plot
  cc_plot
  died_plot
  
  
  
  ###########################################################################################################################
  # 2.2 DO THE SAME FOR CC TRANSITIONS FOR ALL AGE-GROUPS
  
  
  CI.byagegrp_CC <- cuminc(ftime = cohort1$cc_LoS, fstatus = cohort1$cc_transitions, group = cohort1$agegrp_v3)
  
  if(length(CI.byagegrp_CC) != ((3*length(ages)) + 1)){
    
    cc_cuminc_df <- 0
    cc_los <- los_df
    cc_los[which(los_df$Transition == "CC"),1] <- "GA"
    cc_los[,c(2,3,4)] <- NA
   
   }else{
      if(length(unique(cohort1$agegrp_v3)) == 3){
        cc_df_list <- cuminc_list_to_df(CI_dataset = CI.byagegrp_CC, group_names = c("Recovered, <25", "Recovered, 25-64", "Recovered, 65+",
                                                                                 "GA, <25", "GA, 25-64", "GA, 65+","Died, <25",
                                                                                 "Died, 25-64","Died, 65+"), ICD = ICD_group)
    
      }else{
        if(ICD_group == 15)
          
        ages_present <- sort(unique(cohort1$agegrp_v3))
        transitions <- rep(c("Recovered,","GA,","Died,"), each = length(ages_present))
        ages_to_use <- rep(ages_present, 3)
        transitions_to_use <- paste(transitions, ages_to_use, sep = " ")
        cc_df_list <- cuminc_list_to_df(CI.byagegrp_CC, ICD = ICD_group, unique_only = TRUE,
                                            group_names = transitions_to_use) ## JD to export to CSV
      
    }    
  
    
    cc_cuminc_df <- cc_df_list[[1]]
    cc_los <- cc_df_list[[2]]
    
    cuminc_plot_cc <- ggplot(data = cc_cuminc_df, aes(x = time, y = CIF, group = age_tran)) +
      geom_line(aes(linetype = transition, colour = age)) + xlab("Length of stay (days)") +
      ylab("Cumulative transition probability") +ggtitle("Cumualtive transition prob, across 3 transitions")
    
    
    ## Just recovered now ##
    
    recovered_df_cc <- cc_cuminc_df[grep("recovered", cc_cuminc_df$transition, ignore.case = TRUE),]
    
    recovered_plot_cc <- ggplot(data = recovered_df_cc, aes(x = time, y = CIF, group = age_tran)) +
      geom_line(aes(colour = age)) + xlab("Length of stay (days)") +
      ylab("Cumulative Recovery probability") + ggtitle("Recovery probabilities")
    
    
    ## Just CC transition ##
    
    ga_df_cc <- cc_cuminc_df[grep("GA", cc_cuminc_df$transition, ignore.case = TRUE),]
    
    ga_plot_cc <- ggplot(data = ga_df_cc, aes(x = time, y = CIF, group = age_tran)) +
      geom_line(aes(colour = age)) + xlab("Length of stay (days)") +
      ylab("Cumulative move to CC probability") + ggtitle("CC > G&A probabilities")
    
    
    ## Just Died ##
    
    died_df_cc <- cc_cuminc_df[grep("died", cc_cuminc_df$transition, ignore.case = TRUE),]
    
    died_plot_cc <- ggplot(data = died_df_cc, aes(x = time, y = CIF, group = age_tran)) +
      geom_line(aes(colour = age)) + xlab("Length of stay (days)") +
      ylab("Cumulative Death probability") + ggtitle("Death probabilities")
    
    
    cuminc_plot_cc
    recovered_plot_cc
    ga_plot_cc
    died_plot_cc
  }
  dev.off()
  
  
  end_time <- Sys.time()
  print(end_time - cc_start)
  
  return(list(cuminc_df, los_df, cc_cuminc_df,cc_los))
  
}
# 2.1.3 Extract cumulative incidence functions (CIFs) and variance by transition type and age-group --> put in a df and export to csv


###########################################################################################################################
###########################################################################################################################
# 3. Competing Risks - Cohort 3 (Electives staying Electives)
# use "crr" package since we have WaitingTime as a covariate)
# 3.1 GA transitions (Cohort = 3, event = ga_transitions == 1, 2, or 3, time = GA_LoS)


crr_cluster_run <- function(fail_type, cohort_data){
  current_fail <- as.integer(str_split_fixed(fail_type, "-",2)[1])
  current_type <- as.character(str_split_fixed(fail_type,"-",2)[2])
  
  
  
  if(current_type == "ga"){
    if(nrow(cohort_data[cohort_data$ga_transitions == current_fail,]) <= 1){
      return(list(0,0))
    }else{
        agegrp_3_2_time <- cohort_data$GA_LoS
        agegrp_3_2_trans <- cohort_data$ga_transitions
        agegrp_3_2_WT <- cohort_data$WaitingTime
        
        na_tims <- which(is.na(agegrp_3_2_time))
        if(length(na_tims) > 0){
          agegrp_3_2_time <- agegrp_3_2_time[-na_tims]
          agegrp_3_2_trans <- agegrp_3_2_trans[-na_tims]
          agegrp_3_2_WT <- agegrp_3_2_WT[-na_tims]
        }
        old_na <- which(is.na(agegrp_3_2_trans))
        if(length(old_na) > 0){
          agegrp_3_2_time <- agegrp_3_2_time[-old_na]
          agegrp_3_2_trans <- agegrp_3_2_trans[-old_na]
          agegrp_3_2_WT <- agegrp_3_2_WT[-old_na]
        }
        if(current_fail != 1){
          old_2ers <- which(agegrp_3_2_trans == current_fail)
          old_1ers <- which(agegrp_3_2_trans == 1)
          
          agegrp_3_2_trans[old_2ers] <- 1
          agegrp_3_2_trans[old_1ers] <- 2
          
        }
        
        
        CI.agegrp1_t1 <- fastCrr(Crisk(agegrp_3_2_time, agegrp_3_2_trans, failcode = 1) ~ agegrp_3_2_WT,
                                          variance = TRUE, returnDataFrame = TRUE)
        
        # cif1_pred <- predict(CI.agegrp1_t1, newdata = 0, tL = 1)
        # cif1_pred <- cbind(cif1_pred$ftime, cif1_pred$CIF)
        cif1_pred <- fastcmprsk:::predict.fcrr(CI.agegrp1_t1, newdata = 0, tL = 1)
        cif1_pred <- cbind(cif1_pred$ftime, cif1_pred$CIF)
    }
  }else if(current_type == "cc"){
    if(nrow(cohort_data[cohort_data$cc_transitions == current_fail,]) <= 1){
      return(list(0,0))
    }else{
      CI.agegrp1_t1 <- crr(ftime = cohort_data$cc_LoS, fstatus = cohort_data$cc_transitions,
                           cov1 = cohort_data$WaitingTime, failcode = current_fail)
      cif1_pred <- predict.crr(CI.agegrp1_t1, cov1 = 0)
    }
  }
  
  return(list(CI.agegrp1_t1, cif1_pred))
    
}

crr_cluster_run_18 <- function(fail_type, cohort_data){
  ## fail_type is a string containing the current fail code to look at
  ## the the ward type (ga or cc)
  ## and the current age 

  current_fail <- as.integer(str_split_fixed(fail_type, "-",3)[1])
  current_type <- as.character(str_split_fixed(fail_type,"-",3)[2])
  current_age <- as.character(str_split_fixed(fail_type,"-",3)[3])
  
  cohort_data <- cohort_data[cohort_data$agegrp_v3 == as.integer(current_age),]
  
  ## check transitions exists 
  
  
  if(current_type == "ga"){
    agegrp_3_2_time <- cohort_data$GA_LoS
    agegrp_3_2_trans <- cohort_data$ga_transitions
    agegrp_3_2_WT <- cohort_data$WaitingTime
    
    if(!(current_fail %in% agegrp_3_2_trans)){
      return(list(0,0))
      
    }else{
    
    na_tims <- which(is.na(agegrp_3_2_time))
    if(length(na_tims) > 0){
      agegrp_3_2_time <- agegrp_3_2_time[-na_tims]
      agegrp_3_2_trans <- agegrp_3_2_trans[-na_tims]
      agegrp_3_2_WT <- agegrp_3_2_WT[-na_tims]
    }
    old_na <- which(is.na(agegrp_3_2_trans))
    if(length(old_na) > 0){
      agegrp_3_2_time <- agegrp_3_2_time[-old_na]
      agegrp_3_2_trans <- agegrp_3_2_trans[-old_na]
      agegrp_3_2_WT <- agegrp_3_2_WT[-old_na]
    }
    if(current_fail != 1){
      old_2ers <- which(agegrp_3_2_trans == current_fail)
      old_1ers <- which(agegrp_3_2_trans == 1)
      
      agegrp_3_2_trans[old_2ers] <- 1
      agegrp_3_2_trans[old_1ers] <- 2
      
      if(!(1 %in% unique(agegrp_3_2_trans)))
        stop("1 not in agegrp 3_2 trans for some reason")
      
    }
    
    CI.agegrp1_t1 <- NULL
    print(unique(agegrp_3_2_trans))
    try(CI.agegrp1_t1 <- fastCrr(Crisk(agegrp_3_2_time, agegrp_3_2_trans, failcode = 1) ~ agegrp_3_2_WT,
                                 variance = TRUE, returnDataFrame = TRUE))
    if(length(CI.agegrp1_t1) != 0){
      
      ## Sometimes the bootstrapping in the predict.fcrr will mean with low event numbers it can't 
      ## estimate them due to using Crisk again, we'll try turning this off if we get an error
      cif1_pred <- NULL
      no_var <- FALSE
      try(cif1_pred <- fastcmprsk:::predict.fcrr(CI.agegrp1_t1, newdata = 0, tL = 1, type = "interval"), silent = TRUE)
      if(length(cif1_pred) == 0){
        try(cif1_pred <- fastcmprsk:::predict.fcrr(CI.agegrp1_t1, newdata = 0, tL = 1, type = "interval",
                                                   getBootstrapVariance = FALSE), silent = TRUE)
        no_var <- TRUE
        if(length(cif1_pred) == 0){
          cif1_pred <- 0
          CI.agegrp1_t1 <- 0
        }
      } 
      if(length(cif1_pred) > 0){
        if(no_var)
          cif1_pred <- cbind(cif1_pred$ftime, cif1_pred$CIF)
        else
          cif1_pred <- cbind(cif1_pred$ftime, cif1_pred$CIF, cif1_pred$lower, cif1_pred$upper )
      }
    }else{
      cif1_pred <- 0
      CI.agegrp1_t1 <- 0
    }
    
    
    # cif1_pred <- predict(CI.agegrp1_t1, newdata = 0, tL = 1)
    # cif1_pred <- cbind(cif1_pred$ftime, cif1_pred$CIF)
    
    }
    
  }else if(current_type == "cc"){
    
    if(!(current_fail %in% cohort_data$cc_transitions)){
      return(list(0,0))
    }else{
    
      CI.agegrp1_t1 <- NULL
      try(CI.agegrp1_t1 <- crr(ftime = cohort_data$cc_LoS, fstatus = cohort_data$cc_transitions,
                               cov1 = cohort_data$WaitingTime, failcode = current_fail))
    
      if(length(CI.agegrp1_t1) != 0){
        cif1_pred <- predict.crr(CI.agegrp1_t1, cov1 = 0)
      }else{
        cif1_pred <- 0
        CI.agegrp1_t1 <- 0
      }
      
    
    }
  }
  
  return(list(CI.agegrp1_t1, cif1_pred))
  
}


cluster_sum_up <- function(cluster_res, out_dfs, current_k){
  
  
    
  
  transitions_state <- seq(((current_k*6) - 5), current_k*6)
  ## Above loops through the transitions states for an age given by current k
  ## tells us the position in the transitions so we can be cc or ga
  counter <- 0
  ## tells us the row of the outdf to use basically the transitions state (not age!)
  out_df_trans <- 0
  ## tells us which of the out dfs to use
  ward <- 1
  ## gives us all the rows for this age group 
  age_rows <- c((current_k*3) -2, (current_k*3) -1, (current_k*3))
  
  
  for(trans in transitions_state){
    counter <- counter + 1
    out_df_trans <- out_df_trans + 1
    if(counter == 4){
      out_df_trans <- 1
      ward <- 2
    }
    
    if(length(cluster_res[[trans]][[1]]) > 1){
    
      CI.agegrp1_t1 <- cluster_res[[trans]][[1]]
      
      tic("Predicting trans 1")
      cif1_pred <- cluster_res[[trans]][[2]]
      toc()
      
      if(counter < 4){
        cif1 <- cluster_res[[trans]][[2]]
      }else{
        cif1 <- cluster_res[[trans]][[2]]
      }
      
      
      
      
      if(7 %in% cif1[,1]){
        seven_t1 <- cif1[cif1[,1] == 7, 2]
      }else if(nrow(cif1[cif1[,1] <7, ,drop = FALSE]) > 0){
        seven_t1 <- cif1[cif1[,1] < 7, 2]
        seven_t1 <- max(seven_t1)
      }else{
        seven_t1 <- 0
      }
      
      max_time <- max(cif1[,1])
      if(max_time <= 7){
        multi_df <- seven_t1
      }else{
      
        seven_mults <- seq(7, max_time, by = 7)
        start_sevens <- seven_t1
        
        diff_vec <- start_sevens
        actual_vec <- start_sevens
          
          for(j in 2:length(seven_mults)){
            current_mult_vals <- cif1[cif1[,1] > seven_mults[j-1] &
                                        cif1[,1] <= seven_mults[j],,drop = FALSE]
            if(nrow(current_mult_vals) == 0){
              diff_vec <- append(diff_vec, 0)
              actual_vec <- append(actual_vec, actual_vec[length(actual_vec)])
            }else{
              diff_val <- max(current_mult_vals[,2]) - actual_vec[length(actual_vec)]
              diff_vec <- append(diff_vec, diff_val)
              actual_vec <- append(actual_vec, max(current_mult_vals[,2]))
            }
            
          }
        multi_df <- actual_vec
      }
      out_dfs[[ward]]$coeff[age_rows[out_df_trans]] <- CI.agegrp1_t1$coef
      out_dfs[[ward]]$variance[age_rows[out_df_trans]] <- CI.agegrp1_t1$var[1,1]
      out_dfs[[ward]]$day_7[age_rows[out_df_trans]] <- seven_t1
      out_dfs[[ward]]$mean_7[age_rows[out_df_trans]] <- mean(multi_df, na.rm = TRUE)
      out_dfs[[ward]]$median_7[age_rows[out_df_trans]] <- median(multi_df, na.rm = TRUE)
    }else{
      out_dfs[[ward]]$coeff[age_rows[out_df_trans]] <- 0
      out_dfs[[ward]]$variance[age_rows[out_df_trans]] <- 0
      out_dfs[[ward]]$day_7[age_rows[out_df_trans]] <- 0
      out_dfs[[ward]]$mean_7[age_rows[out_df_trans]] <- 0
      out_dfs[[ward]]$median_7[age_rows[out_df_trans]] <- 0
    }
  }
  
  return(out_dfs)
  
}

cluster_crr_try <- function(ftime , fstatus , fail_code, wt_variable){
  print("I'm here at the start GA")
  crr_res <- tryCatch({
    fastCrr(Crisk(ftime, fstatus, failcode = fail_code) ~ wt_variable,
            variance = TRUE, returnDataFrame = TRUE)
    
  },
  error = function(cond){
    message(paste("Cohort 1 CRR not working for Failcode:", fail_code))
    return(NULL)
  },
  warning = function(cond){
    
  },
  finally = {
    message(paste("Processed cohort 1 CRR for Failcode:", fail_code))
  }
  
  
  )
  
  print("I'm here now")
  
  return(crr_res)
}

cluster_crr_try_cc <- function(ftime , fstatus , fail_code, wt_variable){
  print("I'm here at the start CC")
  crr_res <- tryCatch({
    message("Trying out the CRR function")
    crr(ftime = ftime, fstatus = fstatus,
        cov1 = wt_variable, failcode = fail_code)
  },
  error = function(cond){
    message(paste("Cohort 1 CRR not working for Failcode:", fail_code))
    
    return(NULL)
  },
  finally = {
    message(paste("Processed cohort 1 CRR for Failcode:", fail_code))
  }
  
  
  )
  
  
  return(crr_res)
}



cohort_1_competing_risk <- function(cohort1, current_ICD,core_18 = FALSE){
  
  if(!("WaitingTime" %in% colnames(cohort1))){
    cohort23 <- make_wt_variable(cohort1)
    WT_colname <- which(colnames(cohort1) == "WT")
    colnames(cohort1)[WT_colname] <- "WaitingTime"
  }
  
  print("Now on Competing risk cohort 1")
  
  # 3.1.1 Subset age-groups within cohort 3
  out_df <- data.frame(matrix(ncol = 8, nrow = 9))
  colnames(out_df) <- c("age","transitions", "ICD","coeff","variance","day_7","mean_7","median_7")
  out_df$age <- rep(c(1,2,3), each = 3)
  out_df$transition <- rep(c(1,2,3), by = 3)
  out_df$ICD <- current_ICD
  
  cc_out_df <- data.frame(matrix(ncol = 8, nrow = 9))
  colnames(cc_out_df) <- c("age","transitions", "ICD","coeff","variance","day_7","mean_7","median_7")
  cc_out_df$age <- rep(c(1,2,3), each = 3)
  cc_out_df$transition <- rep(c(1,2,3), by = 3)
  cc_out_df$ICD <- current_ICD
  
  
  # 3.1.2 Give frequency of different GA transitions for each age-group
  if(core_18 == FALSE){
      transitions_type <- rep(c("1","2","3"), 2)
      area_codes <- rep(c("ga","cc"), each = 3)
      typies <- paste(transitions_type, area_codes, sep = "-")
  }else{
    transitions_type <- rep(c("1","2","3"), 6)
    area_codes <- rep(rep(c("ga","cc"), each = 3), 3)
    age_codes <- rep(c(1,2,3), each = 6)
    age_typies <- paste(transitions_type,area_codes,age_codes, sep = "-")
  }
  
  if(core_18){
    print("Running all ICD at once")
    print("Setting up cluster run")
    tic("Cluster set up:")
    
    crr_cluster <- snow::makeCluster(spec = 18, outfile = "./crr_cluster_log.txt")
    toc()
    crr_input <- snow::clusterSplit(crr_cluster, age_typies)
    
    snow::clusterExport(crr_cluster, "crr_cluster_run_18")
    snow::clusterExport(crr_cluster, "cluster_crr_try")
    snow::clusterExport(crr_cluster, "cluster_crr_try_cc")
    print(paste("Copying over data for ICD"))
    copy_start <- Sys.time()
    snow::clusterExport(crr_cluster, "cohort1", envir = environment())
    copy_end <- Sys.time()
    print(copy_end - copy_start)
    snow::clusterEvalQ(crr_cluster, require(cmprsk))
    snow::clusterEvalQ(crr_cluster, require(stringr))
    snow::clusterEvalQ(crr_cluster, require(fastcmprsk, lib.loc = "C:/R-3.6.2/library"))
    print(paste("Running CRR jobs for ICD"))
    jobs_start <- Sys.time()
    crr_data_parallel <- snow::clusterApply(crr_cluster, crr_input,
                                            fun = crr_cluster_run_18,
                                            cohort_data = cohort1)
    snow::clusterEvalQ(crr_cluster, rm(cohort1))
    snow::clusterEvalQ(crr_cluster, gc())
    stopCluster(crr_cluster)
    jobs_end <- Sys.time()
    print(jobs_end - jobs_start)
    
    joint_ccr_out <- list(out_df, cc_out_df)
    ghost_list <- list(out_df, cc_out_df)
    ## This loop takes the crr results and sums up our results for output
    
    for(k in 1:3){

      if(k == 1){
      joint_ccr_out <- cluster_sum_up(cluster_res = crr_data_parallel, 
                                      current_k = k, 
                                      out_dfs = ghost_list)
      }else{
        joint_ccr_out <- cluster_sum_up(cluster_res = crr_data_parallel, 
                                        current_k = k, 
                                        out_dfs = joint_ccr_out)
      }
      
    }
    
    out_df <- joint_ccr_out[[1]]
    cc_out_df <- joint_ccr_out[[2]]
    
  }else{
  
  for(k in 1:3){
    print(paste("On GA age group:", k))
    agegrp1 <- subset(cohort1, subset = agegrp_v3 == k)
    
    
    
    crr_cluster <- snow::makeCluster(spec = 6, outfile = "./crr_cluster_log.txt")
    crr_input <- snow::clusterSplit(crr_cluster, typies)
    
    snow::clusterExport(crr_cluster, "crr_cluster_run")
    print(paste("Copying over data for CRR age:", k))
    copy_start <- Sys.time()
    snow::clusterExport(crr_cluster, "agegrp1", envir = environment())
    copy_end <- Sys.time()
    print(copy_end - copy_start)
    snow::clusterEvalQ(crr_cluster, require(cmprsk))
    snow::clusterEvalQ(crr_cluster, require(stringr))
    snow::clusterEvalQ(crr_cluster, require(fastcmprsk, lib.loc = "C:/R-3.6.2/library"))
    print(paste("Running CRR jobs for age:", k))
    jobs_start <- Sys.time()
    crr_data_parallel <- snow::clusterApply(crr_cluster, crr_input,
                                            fun = crr_cluster_run,
                                            cohort_data = agegrp1)
    snow::clusterEvalQ(crr_cluster, rm(agegrp1))
    snow::clusterEvalQ(crr_cluster, gc())
    snow::stopCluster(crr_cluster)
    jobs_end <- Sys.time()
    print(jobs_end - jobs_start)
    
    
    CI.agegrp1_t1 <- crr_data_parallel[[1]][[1]]
    CI.agegrp1_t2 <- crr_data_parallel[[2]][[1]]
    CI.agegrp1_t3 <- crr_data_parallel[[3]][[1]]
    
    CI.agegrp1_t1_cc <- crr_data_parallel[[4]][[1]]
    CI.agegrp1_t2_cc <- crr_data_parallel[[5]][[1]]
    CI.agegrp1_t3_cc <- crr_data_parallel[[6]][[1]]
    
    
    
    age_rows <- c((k*3) -2, (k*3) -1, (k*3))
    tic("Predicting trans 1")
    cif1_pred <- crr_data_parallel[[1]][[2]]
    toc()
    tic("Predicting trans 2")
    cif2_pred <- crr_data_parallel[[2]][[2]]
    toc()
    tic("Predicting trans 3")
    cif3_pred <- crr_data_parallel[[3]][[2]]
    toc()
    
    cif1 <- cbind(cif1_pred$ftime, cif1_pred$CIF)
    cif2 <- cbind(cif2_pred$ftime, cif2_pred$CIF)
    cif3 <- cbind(cif3_pred$ftime, cif3_pred$CIF)
    
    
    if(7 %in% cif1[,1]){
      seven_t1 <- cif1[cif1[,1] == 7, 2]
    }else if(nrow(cif1[cif1[,1] <7, ,drop = FALSE]) > 0){
      seven_t1 <- cif1[cif1[,1] < 7, 2]
      seven_t1 <- max(seven_t1)
    }else{
      seven_t1 <- 0
    }
    if(7 %in% cif2[,1]){
      seven_t2 <- cif2[cif2[,1] == 7, 2]
    }else if(nrow(cif2[cif2[,1] <7, ,drop = FALSE]) > 0){
      seven_t2 <- cif2[cif2[,1] < 7, 2]
      seven_t2 <- max(seven_t2)
    }else{
      seven_t2 <- 0
    }
    if(7 %in% cif3[,1]){
      seven_t3 <- cif3[cif3[,1] == 7, 2]
    }else if(nrow(cif3[cif3[,1] <7, ,drop = FALSE]) > 0){
      seven_t3 <- cif3[cif3[,1] < 7, 2]
      seven_t3 <- max(seven_t3)
    }else{
      seven_t3 <- 0
    }
    
    max_time <- max(c(cif1[,1], cif2[,1], cif3[,1]))
    seven_mults <- seq(7, max_time, by = 7)
    start_sevens <- c(seven_t1, sevent_t2,seven_t3)
    cifs <- list(cif1, cif2, cif3)
    multi_outs <- data.frame(matrix(ncol = 3, nrow = length(seven_mults)))
    
    for(trans in 1:3){
      
      multi_df <- start_sevens[trans]
      
      for(j in 2:length(seven_mults)){
        current_mult_vals <- cifs[[trans]][cifs[[trans]][,1] > seven_mults[k-1] &
                                             cifs[[trans]][,1] <= seven_mults[k],]
        if(nrow(current_mult_vals) == 0){
          multi_df <- append(multi_df, 0)
        }else{
          diff_val <- max(current_mult_vals$cumulative_hazard) - multi_df[length(multi_df)]
          multi_df <- append(multi_df, diff_val)
        }
        
      }
      
      multi_outs[,k] <- multi_df
      
    }
    
    
    out_df$coeff[age_rows[1]] <- crr_data_parallel[[1]]$coef
    out_df$variance[age_rows[1]] <- crr_data_parallel[[1]]$var[1,1]
    out_df$day_7[age_rows[1]] <- seven_t1
    out_df$mean_7[age_rows[1]] <- mean(multi_outs[,1])
    out_df$median_7[age_rows[1]] <- median(multi_outs[,1])
    
    
    
    out_df$coeff[age_rows[2]] <- crr_data_parallel[[2]]$coef
    out_df$variance[age_rows[2]] <- crr_data_parallel[[2]]$var[1,1]
    out_df$day_7[age_rows[2]] <- seven_t2
    out_df$mean_7[age_rows[2]] <- mean(multi_outs[,2])
    out_df$median_7[age_rows[2]] <- median(multi_outs[,2])
    
    
    out_df$coeff[age_rows[3]] <- CI.agegrp1_t3$coef
    out_df$variance[age_rows[3]] <- CI.agegrp1_t3$var[1,1]
    out_df$day_7[age_rows[3]] <- seven_t3
    out_df$mean_7[age_rows[3]] <- mean(multi_outs[,3])
    out_df$median_7[age_rows[3]] <- median(multi_outs[,3])
    
    
    age_rows <- c((k*3) -2, (k*3) -1, (k*3))
    
    cif1_cc <- crr_data_parallel[[4]][[2]]
    cif2_cc <- crr_data_parallel[[5]][[2]]
    cif3_cc <- crr_data_parallel[[6]][[2]]
    
    
    if(7 %in% cif1_cc[,1]){
      seven_t1_cc <- cif1_cc[cif1_cc[,1] == 7, 2]
    }else if(nrow(cif1_cc[cif1_cc[,1] <7,,drop = FALSE ]) > 0){
      seven_t1_cc <- cif1_cc[cif1_cc[,1] < 7, 2]
      seven_t1_cc <- max(seven_t1_cc)
    }else{
      seven_t1_cc <- 0
    }
    if(7 %in% cif2_cc[,1]){
      seven_t2_cc <- cif2_cc[cif2_cc[,1] == 7, 2]
    }else if(nrow(cif2_cc[cif2_cc[,1] <7, ,drop = FALSE]) > 0){
      seven_t2_cc <- cif2_cc[cif2_cc[,1] < 7, 2]
      seven_t2_cc <- max(seven_t2_cc)
    }else{
      seven_t2_cc <- 0
    }
    if(nrow(cif3_cc)>0){
      
      if(7 %in% cif3_cc[,1]){
        seven_t3_cc <- cif3_cc[cif3_cc[,1] == 7, 2]
      }else if(nrow(cif3_cc[cif3_cc[,1] <7,,drop = FALSE ]) > 0){
        seven_t3_cc <- cif3_cc[cif3_cc[,1] < 7, 2]
        seven_t3_cc <- max(seven_t3_cc)
      }else{
        seven_t3_cc <- 0
      }
    }else{
      seven_t3_cc <- 0
    }
    max_time <- max(c(cif1_cc[,1], cif2_cc[,1], cif3_cc[,1]))
    
    seven_mults <- seq(7, max_time, by = 7)
    start_sevens <- c(seven_t1_cc, sevent_t2_cc,seven_t3_cc)
    cifs <- list(cif1_cc, cif2_cc, cif3_cc)
    multi_outs <- data.frame(matrix(ncol = 3, nrow = length(seven_mults)))
    
    for(trans in 1:3){
      
      multi_df <- start_sevens[trans]
      
      for(j in 2:length(seven_mults)){
        current_mult_vals <- cifs[[trans]][cifs[[trans]][,1] > seven_mults[k-1] &
                                             cifs[[trans]][,1] <= seven_mults[k],]
        if(nrow(current_mult_vals) == 0){
          multi_df <- append(multi_df, 0)
        }else{
          diff_val <- max(current_mult_vals$cumulative_hazard) - multi_df[length(multi_df)]
          multi_df <- append(multi_df, diff_val)
        }
        
      }
      
      multi_outs[,k] <- multi_df
      
    }
    
    
    
    
    cc_out_df$coeff[age_rows[1]] <- CI.agegrp1_t1_cc$coef
    cc_out_df$variance[age_rows[1]] <- CI.agegrp1_t1_cc$var[1,1]
    cc_out_df$day_7[age_rows[1]] <- seven_t1_cc
    cc_out_df$mean_7[age_rows[1]] <- mean(multi_outs[,1])
    cc_out_df$median_7[age_rows[1]] <- median(multi_outs[,1])
    
    cc_out_df$coeff[age_rows[2]] <- CI.agegrp1_t2_cc$coef
    cc_out_df$variance[age_rows[2]] <- CI.agegrp1_t2_cc$var[1,1]
    cc_out_df$day_7[age_rows[2]] <- seven_t2_cc
    cc_out_df$mean_7[age_rows[2]] <- mean(multi_outs[,2])
    cc_out_df$median_7[age_rows[2]] <- median(multi_outs[,2])
    
    
    cc_out_df$coeff[age_rows[3]] <- CI.agegrp1_t3_cc$coef
    cc_out_df$variance[age_rows[3]] <- CI.agegrp1_t3_cc$var[1,1]
    cc_out_df$day_7[age_rows[3]] <- seven_t3_cc
    cc_out_df$mean_7[age_rows[3]] <- mean(multi_outs[,3])
    cc_out_df$median_7[age_rows[3]] <- median(multi_outs[,3])
    
    
    
  }
  }
  
  # run this one as an example because it takes less time than the above
  
  ## what we need from this output: coef, var
  
  
  return(list(out_df, cc_out_df))
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

crr_func_try <- function(cohort1 , current_ICD ,
                         core_18){
  out <- tryCatch({
    message("Trying out the CRR function")
    cohort_1_competing_risk(cohort1 = cohort1 , current_ICD = current_ICD,
                            core_18 = core_18)  
    
  },
  error = function(cond){
    message(paste("Cohort 1 CRR not working for ICD:", current_ICD))
    message(cond)
    return(NULL)
  },
  warning = function(cond){
    message(paste("Following warning message for ICD:", current_ICD))
    message(cond)
    return(out)
    
  },
  finally = {
    message(paste("Processed cohort 1 CRR for ICD:", current_ICD))
  }
  
  
  )
  
  
  return(out)
}


cohort_3_func_try <- function(cohort3, ICD_group, results_pdf){
  out <- tryCatch({
    message("Trying out the Cohort 3 competing risks")
    
    cohort_3_comp_risks(cohort1 = cohort3,
                        ICD_group = ICD_group,
                        results_pdf = results_pdf)   
  },
  error = function(cond){
    message(paste("Cohort 3 competing risk not working for ICD:", ICD_group))
    message(cond)
    return(NULL)
    
  # },
  # warning = function(cond){
  #   message(paste("Following warning message for ICD:", ICD_group))
  #   message(cond)
  #   return(out)
  #   
  },
  finally = {
    message(paste("Processed Cohort 3 competing risks for ICD:", ICD_group))
  }
  
  
  )
  
  
  return(out)
}


survival_analysis_set_up <- function(cohorts_data, single_ICD = TRUE, base_dir, single_icd, core_18 = FALSE,
                                     elective_run = TRUE, emergency_run = TRUE, crr_try = TRUE){
  require(tidyverse)
  require(survival)
  require(plyr)
  require(timereg)
  require(ggplot2)
  require(ggpubr) ## need rtools installed for multi-plots
  require(dplyr)
  require(tictoc)
  require(cmprsk)
  require(fastcmprsk, lib.loc = "C:/R-3.6.2/library")
  
  # Import survival analysis data
  
  
  # Subset cohorts
  
  if(single_ICD){
    
  
    cols_needed <- c("WaitingTime","Elective2Emergency","GA_LoS","ga_transitions",
                     "agegrp_v3","cc_LoS","cc_transitions")
    cols_to_keep <- which(colnames(cohorts_data) %in% cols_needed)
    
    cohort1 <- cohorts_data[cohorts_data$cohort == 1, cols_to_keep]
    cohort2 <- cohorts_data[cohorts_data$cohort == 2, cols_to_keep]
    cohort3 <- cohorts_data[cohorts_data$cohort == 3, cols_to_keep]
    cohort12 <- cohorts_data[cohorts_data$cohort != 3, cols_to_keep]
    
    
    failure_pdf <- paste(base_dir, "/survival_plot.pdf", sep = "")
    failure_csv <- paste(base_dir, "/survival_for_groups.csv", sep = "")
    transitions_pdf <- paste(base_dir, "/one_ICD_results.pdf", sep = "")
    
    if(elective_run){
    cohort12_failure_func <- failure_function(cohort23 = cohort12, csv_survial_name = failure_csv,
                                              pdf_to_export = failure_pdf, current_ICD = single_icd)   
    cohort1_transitions <- cohort_1_competing_risk(cohort1, current_ICD = single_icd, core_18 = core_18)
    }else{
      cohort12_failure_func <- 0
      cohort1_transitions <- list(0,0)
      }
    if(emergency_run){
      cohort_3_transitions <- cohort_3_comp_risks(cohort3, ICD_group = single_icd, results_pdf = transitions_pdf)
    }else{
      cohort_3_transitions <- list(0,0,0,0)
    }
    
    failure_res <- cohort12_failure_func
    cohort_3_ga <- cohort_3_transitions[[2]]
    cohort_3_cc <- cohort_3_transitions[[4]]
    cohort_1_res_ga <- cohort1_transitions[[1]]
    cohort_1_res_cc <- cohort1_transitions[[2]]
    
  
  }else{
    cols_needed <- c("WaitingTime","Elective2Emergency","GA_LoS","ga_transitions",
                     "agegrp_v3","cc_LoS","cc_transitions","MainICD10Cat", "ICD")
    cols_to_keep <- which(colnames(cohorts_data) %in% cols_needed)
    
    
    
    cohort1 <- cohorts_data[which(cohorts_data$cohort == 1),]
    cohort2 <- cohorts_data[which(cohorts_data$cohort == 2),]
    cohort3 <- cohorts_data[which(cohorts_data$cohort == 3),]
    cohort12 <- cohorts_data[-which(cohorts_data$cohort == 3),]
    
    cohort12 <- make_elective_cohort_variable(cohort12)
    elective_groupings <- c(2,3,6,7,9,10,11,12,13,14,18,19,21,50)
    emergency_groupings <- c(1,2,4,5,6,9,10,11,12,13,14,15,18,19,51)
    
    failure_res <- NULL
    cohort_1_res_ga <- NULL
    cohort_1_res_cc <- NULL
    cohort_3_ga <- NULL
    cohort_3_cc <- NULL
    if(elective_run == TRUE){
    
      for(k in 1:length(elective_groupings)){
        
        current_icd <- elective_groupings[k]
        cohort12_icd <- cohort12[cohort12$elective_icd == current_icd,]
        cohort1_icd <- cohort1[cohort1$ICD == current_icd,]
        
        failure_csv <- paste(base_dir, "./surival_cohort_12_ICD",current_icd,".csv", sep = "")
        failure_pdf <- paste(base_dir, "./surival_cohort_12_ICD",current_icd,".pdf", sep = "")
        
        ## try catch for survival analysis 
        ## failure function 
        
        
        cohort12_failure_func <- failure_func_try(cohort23 =  cohort12_icd,
                                                  csv_survival_name = failure_csv,
                                                  pdf_to_export = failure_pdf, current_ICD = current_icd)   
        
        if(length(cohort12_failure_func) != 0){
          failure_res <- rbind.data.frame(failure_res, cohort12_failure_func[[2]])
        }
        
        if(crr_try){
        
           cohort1_transitions <- crr_func_try(cohort1 = cohort1_icd , current_ICD = current_icd,
                                                         core_18 = core_18)
          if(length(cohort12_failure_func) != 0){
            cohort_1_res_ga <- rbind.data.frame(cohort_1_res_ga, cohort1_transitions[[1]])
            cohort_1_res_cc <- rbind.data.frame(cohort_1_res_cc, cohort1_transitions[[2]])
          }
    }else{
      cohort_1_res_ga <- 0
      cohort_1_res_cc <- 0
    }
    }  
    }
    if(emergency_run == TRUE){
    for(k in 1:length(emergency_groupings)){
      current_icd <- emergency_groupings[k]
      current_cohort_3 <- cohort3[cohort3$ICD == current_icd,]
      
      transitions_pdf <- paste(base_dir, "/cohort_3_ICD",current_icd,"_results.pdf", sep = "")
      
      
      current_3_transitions <- cohort_3_func_try(cohort3 = current_cohort_3,ICD_group = current_icd,
                                                 results_pdf = transitions_pdf)
      
      if(length(current_3_transitions) != 0){
        cohort_3_ga <- rbind.data.frame(cohort_3_ga, current_3_transitions[[2]])
        cohort_3_cc <- rbind.data.frame(cohort_3_cc, current_3_transitions[[4]])
      }
      
      
      
      
    }
    
    }
    
  }
  
  ## Lets write these dfs out 
  if(substr(base_dir, nchar(base_dir),nchar(base_dir)) == "/")
    base_dir <- substr(base_dir,1,nchar(base_dir) - 1 )
  
  failure_csv <- paste(base_dir, "/summary_cohort_1_to_3.csv", sep = "")
  elective_ga_csv <- paste(base_dir, "/cohort_1_ga_transitions.csv", sep = "")
  elective_cc_csv <- paste(base_dir, "/cohort_1_cc_transitions.csv", sep = "")
  emergency_ga_csv <- paste(base_dir, "/cohort_3_ga_transitions.csv", sep = "")
  emergency_cc_csv <- paste(base_dir, "/cohort_3_cc_transitions.csv", sep = "")
  
  if(elective_run){
    write.csv(failure_res, file = failure_csv, row.names = FALSE, quote = FALSE)
    if(crr_try){
      write.csv(cohort_1_res_ga, file = elective_ga_csv, row.names = FALSE, quote = FALSE)
      write.csv(cohort_1_res_cc, file = elective_cc_csv, row.names = FALSE, quote = FALSE)
    }
      
  }
  if(emergency_run){
    print(emergency_ga_csv)
    print(emergency_cc_csv)
    print(cohort_3_cc)
    write.csv(cohort_3_ga, file = emergency_ga_csv, row.names = FALSE, quote = FALSE)
    write.csv(cohort_3_cc, file = emergency_cc_csv, row.names = FALSE, quote = FALSE)
  }
  
  
  return(list(failure_res, cohort_3_ga, cohort_3_cc, cohort_1_res_ga, cohort_1_res_cc))
}



