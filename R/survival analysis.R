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
  browser()
  one_year_under <- which(as.integer(elective_df$elecdur) <= 365)
  
  elective_df$WT <- NA
  elective_df$WT[one_year_under] <- elective_df$elecdur[one_year_under]
  
  return(elective_df)
  
}




failure_function <- function(cohort23, pdf_to_export = "electives_to_emergencies_survival_plot",
                             csv_survial_name){
  ## Need to double check whether Waiting time is in there!
  start_time <- Sys.time()
  print("Failure function for electives to emergencies")
  if(!("WaitingTime" %in% colnames(cohort23))){
    cohort23 <- make_wt_variable(cohort23)
    WT_colname <- which(colnames(cohort23) == "WT")
    colnames(cohort23)[WT_colname] <- "WaitingTime"
  }
  
  # 1.1 Fit survival function
  surv_object <- Surv(time = cohort23$WaitingTime, event = cohort23$Elective2Emergency)
  e2e <- survfit(surv_object ~ agegrp_v3, data = cohort23, type = "kaplan-meier")
  
  # 1.2 Plot survival function
  survival_df <- cbind.data.frame(e2e$time, e2e$surv, rep(c("<24", "25-64", "65+"), each = length(e2e$surv)/3),
                                  e2e$cumhaz, e2e$std.err)
  colnames(survival_df) <- c("time","survival","age_group", "cumulative_hazard", "std_err_hazard")
  survival_df$ICD <- 2
  
  survival_plot_2 <- ggplot(data = survival_df, aes(x = time, y = survival, group = age_group, colour = age_group)) +
          geom_line() + xlab("Waiting time (days)") + ylab("Survival probability") + ylim(c(0,1)) + labs(colour = "Age group") +
          ggtitle("Survival function")
  
  
  # 1.3 Plot failure function
  cumhaz_plot_2 <- ggplot(data = survival_df, aes(x = time, y = cumulative_hazard, group = age_group, colour = age_group)) +
          geom_line() + xlab("Waiting time (days)") + ylab("Failure probability") + ylim(c(0,1)) + labs(colour = "Age group") +
          ggtitle("Failure function (cumulative hazard)")
  
  
  # 1.4 Put failure function and standard errors into a df --> export to csv
  #### This needs to be in the correct dir, from the setwd from above ####
  csv_survival_name <- paste("survival_electives_to_emergencies_ICD",as.character(2),".csv",sep = "")
  
  pdf(file = pdf_to_export, paper = "a4r")
  
  print(survival_plot_2)
  print(cumhaz_plot_2)
  
  dev.off()
  
  write.csv(survival_df,
            file = csv_survival_name,
            row.names = FALSE)
  
  end_time <- Sys.time()
  print(end_time - start_time)
  
  return(survival_df)

}
###########################################################################################################################
###########################################################################################################################
# 2. Competing Risks - Cohort 1 (True Emergencies)
# use "cmprsk" package since no WT covariate
# NEED TO CODE LOOPING: by ICD, agegrp, transition number, vector of GA_LoS & ga_transitions vs. vector of cc_LoS & cc_transitions)

# 2.1 GA transitions (Cohort = 1, event = ga_transitions == 1, 2, or 3, time = GA_LoS)

cuminc_list_to_df <- function(CI_dataset, group_names = c("Recovered, <25", "Recovered, 25-64", "Recovered, 65+",
                                                          "CC, <25", "CC, 25-64", "CC, 65+","Died, <25",
                                                          "Died, 25-64","Died, 65+"),
                              ICD, unique_only = FALSE){
  ## This function takes the output of the competing risks model for emergency admissions and 
  ## then plots our cumulative incidence curves and then returns (in LOS) the transitions probs
  ## at day 7 and multiples of 7

  tests_pos <- grep("Tests", names(CI_dataset), ignore.case = TRUE)
  
  
  
  cuminc_df <- NULL
  
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
  
  
  
  
  age_groupings <- plyr::count(los_data$age)
  transitions <- plyr::count(los_data$transition)
  
  los_df <- data.frame(matrix(ncol = 5, nrow = 3*nrow(age_groupings)))
  
  
  for(k in 1:nrow(age_groupings)){
    age_rows <- c((k*3) -2, (k*3) -1, (k*3))
    current_age_data <- los_data[los_data$age == age_groupings[k,1],]
    
    multi_7 <- seq(7,max(current_age_data$time), by = 7)
    
    
    colnames(los_df) <- c("Transition","Day_7","mean_7_multiple","median_7_multiple", "age")
    los_df$age[age_rows] <- rep(as.character(age_groupings[k,1]), 3)
    los_df$Transition[age_rows] <- as.character(transitions[,1])
    
    current_7_data <- current_age_data[current_age_data$time == 7,]
    
    
    if(nrow(current_7_data) != 3){
      missing_trans <- as.character(transitions[-which(transitions[,1] %in% current_7_data$transition) ,1])
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
    
    
    average_df <- as.data.frame(seventh_day_timies$est)
    trans_nights <- c(k, k+3, k+6)
    
    average_df$zero <- 0
    average_df <- average_df[,c(ncol(average_df),1:(ncol(average_df) - 1))]
    diff_df <- data.frame(matrix(ncol = (ncol(average_df) - 1), nrow = nrow(average_df)))
    for(p in 2:ncol(average_df)){
      
      diff_df[,p-1] <- average_df[,p] - average_df[,(p-1)]
    }
    
    
    tran_1 <- as.numeric(diff_df[trans_nights[1],])
    tran_2 <- as.numeric(diff_df[trans_nights[2],])
    tran_3 <- as.numeric(diff_df[trans_nights[3],])
    
    
    
    
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
  library(msSurv)
    
  # 2.1.1 Give frequency of different GA transitions
  
  
  # 2.1.2 Set survival datatset (grouped by age groups)
  print("On GA transitions")
  ga_start <- Sys.time()
  CI.byagegrp <- cuminc(ftime = cohort1$GA_LoS, fstatus = cohort1$ga_transitions, group = cohort1$agegrp_v3)
  
  # THIS GIVES US THE CIF AT THE SPECIFIC TIME POINT (I DON'T KNOW HOW TO GET A SPECIFIC RISK/AGE-GROUP, IT JUST RETURNS ALL OF THEM)
  
  cuminc_df_list <- cuminc_list_to_df(CI.byagegrp, ICD = ICD_group, unique_only = TRUE) ## JD to export to CSV
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
  
  cc_df_list <- cuminc_list_to_df(CI_dataset = CI.byagegrp_CC, group_names = c("Recovered, <25", "Recovered, 25-64", "Recovered, 65+",
                                                                               "GA, <25", "GA, 25-64", "GA, 65+","Died, <25",
                                                                               "Died, 25-64","Died, 65+"), ICD = 2)
  
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


cohort_1_competing_risk <- function(cohort1, current_ICD){
  
  if(!("WaitingTime" %in% colnames(cohort1))){
    cohort23 <- make_wt_variable(cohort1)
    WT_colname <- which(colnames(cohort1) == "WT")
    colnames(cohort1)[WT_colname] <- "WaitingTime"
  }
  
  print("Now on Competing risk cohort 1")
  
  # 3.1.1 Subset age-groups within cohort 3
  out_df <- data.frame(matrix(ncol = 8, nrow = 9))
  colnames(out_df) <- c("age","transitions", "ICD","coeff","variance","day_7","mean_7","median_7")
  out_df$age <- rep(c("<25","25-64","65+"), each = 3)
  out_df$transition <- rep(c(1,2,3), by = 3)
  out_df$ICD <- current_ICD
  
  cc_out_df <- data.frame(matrix(ncol = 8, nrow = 9))
  colnames(cc_out_df) <- c("age","transitions", "ICD","coeff","variance","day_7","mean_7","median_7")
  cc_out_df$age <- rep(c("<25","25-64","65+"), each = 3)
  cc_out_df$transition <- rep(c(1,2,3), by = 3)
  cc_out_df$ICD <- current_ICD
  
  
  # 3.1.2 Give frequency of different GA transitions for each age-group
  
  
  
  
  
  for(k in 1:3){
    print(paste("On GA age group:", k))
    agegrp1 <- subset(cohort1, subset = agegrp_v3 == k)
    tic("crr runtime")
    print(paste("start running crrs agegrp:", k))
    CI.agegrp1_t1 <- crr(ftime = agegrp1$GA_LoS, fstatus = agegrp1$ga_transitions, cov1 = agegrp1$WaitingTime, failcode = 1)
    toc()
    tic("crr runtime trans 2")
    CI.agegrp1_t2 <- crr(ftime = agegrp1$GA_LoS, fstatus = agegrp1$ga_transitions, cov1 = agegrp1$WaitingTime, failcode = 2)
    toc()
    tic("crr runtime trans 3")
    CI.agegrp1_t3 <- crr(ftime = agegrp1$GA_LoS, fstatus = agegrp1$ga_transitions, cov1 = agegrp1$WaitingTime, failcode = 3)
    toc()
    
    print("end run")
    
    trans_age <- c()
    
    age_rows <- c((k*3) -2, (k*3) -1, (k*3))
    cif1 <- predict.crr(CI.agegrp1_t1, cov1 = 0)
    cif2 <- predict.crr(CI.agegrp1_t2, cov1 = 0)
    cif3 <- predict.crr(CI.agegrp1_t3, cov1 = 0)
    
    
    if(7 %in% cif1[,1]){
      seven_t1 <- cif1[cif1[,1] == 7, 2]
    }else if(nrow(cif1[cif1[,1] <7, ]) > 0){
      seven_t1 <- cif1[cif1[,1] < 7, 2]
      seven_t1 <- max(seven_t1)
    }else{
      seven_t1 <- 0
    }
    if(7 %in% cif2[,1]){
      seven_t2 <- cif2[cif2[,1] == 7, 2]
    }else if(nrow(cif2[cif2[,1] <7, ]) > 0){
      seven_t2 <- cif2[cif2[,1] < 7, 2]
      seven_t2 <- max(seven_t2)
    }else{
      seven_t2 <- 0
    }
    if(7 %in% cif3[,1]){
      seven_t3 <- cif3[cif3[,1] == 7, 2]
    }else if(nrow(cif3[cif3[,1] <7, ]) > 0){
      seven_t3 <- cif3[cif3[,1] < 7, 2]
      seven_t3 <- max(seven_t3)
    }else{
      seven_t3 <- 0
    }
    
    seven_seq <- seq(7, 77, by = 7)
    multi_7_t1  <- cif1[cif1[,1] %in% seven_seq, 2]
    multi_7_t2  <- cif2[cif2[,1] %in% seven_seq, 2]
    multi_7_t3  <- cif3[cif3[,1] %in% seven_seq, 2]
    
    multi_7_t1 <- c(0, multi_7_t1)
    multi_7_t2 <- c(0, multi_7_t2)
    multi_7_t3 <- c(0, multi_7_t3)
    
    multi_7_t1 <- multi_7_t1[2:length(multi_7_t1)] - multi_7_t1[1:(length(multi_7_t1) -1)]
    multi_7_t2 <- multi_7_t2[2:length(multi_7_t2)] - multi_7_t1[1:(length(multi_7_t2) -1)]
    multi_7_t3 <- multi_7_t3[2:length(multi_7_t3)] - multi_7_t1[1:(length(multi_7_t3) -1)]
    
    
    
    
    out_df$coeff[age_rows[1]] <- CI.agegrp1_t1$coef
    out_df$variance[age_rows[1]] <- CI.agegrp1_t1$var[1,1]
    out_df$day_7[age_rows[1]] <- seven_t1
    if(length(multi_7_t1) > 0)
      out_df$mean_7[age_rows[1]] <- mean(multi_7_t1)
    else
      out_df$mean_7[age_rows[1]] <- seven_t1
    if(length(multi_7_t1) > 0)
      out_df$median_7[age_rows[1]] <- median(multi_7_t1)
    else
      out_df$median_7[age_rows[1]] <- seven_t1
    
    
    
    out_df$coeff[age_rows[2]] <- CI.agegrp1_t2$coef
    out_df$variance[age_rows[2]] <- CI.agegrp1_t2$var[1,1]
    out_df$day_7[age_rows[2]] <- seven_t2
    if(length(multi_7_t2) > 0)
      out_df$mean_7[age_rows[2]] <- mean(multi_7_t2)
    else
      out_df$mean_7[age_rows[2]] <- seven_t2
    if(length(multi_7_t2) > 0)
      out_df$median_7[age_rows[2]] <- median(multi_7_t2)
    else
      out_df$median_7[age_rows[2]] <- seven_t2
    
    
    out_df$coeff[age_rows[3]] <- CI.agegrp1_t3$coef
    out_df$variance[age_rows[3]] <- CI.agegrp1_t3$var[1,1]
    out_df$day_7[age_rows[3]] <- seven_t3
    if(length(multi_7_t3) > 0)
      out_df$mean_7[age_rows[3]] <- mean(multi_7_t3)
    else
      out_df$mean_7[age_rows[3]] <- seven_t3
    if(length(multi_7_t3) > 0)
      out_df$median_7[age_rows[3]] <- median(multi_7_t3)
    else
      out_df$median_7[age_rows[3]] <- seven_t3
    
      
      
  }
  
  
  for(k in 1:3){
    print(paste("On CC age group:", k))
    
    agegrp1 <- subset(cohort1, subset = agegrp_v3 == k)
    tic("crr runtime")
    print(paste("start running crrs agegrp:", k))
    CI.agegrp1_t1 <- crr(ftime = agegrp1$cc_LoS, fstatus = agegrp1$cc_transitions, cov1 = agegrp1$WaitingTime, failcode = 1)
    toc()
    tic("crr runtime trans 2")
    CI.agegrp1_t2 <- crr(ftime = agegrp1$cc_LoS, fstatus = agegrp1$cc_transitions, cov1 = agegrp1$WaitingTime, failcode = 2)
    toc()
    tic("crr runtime trans 3")
    CI.agegrp1_t3 <- crr(ftime = agegrp1$cc_LoS, fstatus = agegrp1$cc_transitions, cov1 = agegrp1$WaitingTime, failcode = 3)
    toc()
    
    print("end run")
    
    trans_age <- c()
    
    age_rows <- c((k*3) -2, (k*3) -1, (k*3))
    
    cif1 <- predict.crr(CI.agegrp1_t1, cov1 = 0)
    cif2 <- predict.crr(CI.agegrp1_t2, cov1 = 0)
    cif3 <- predict.crr(CI.agegrp1_t3, cov1 = 0)
    
    
    if(7 %in% cif1[,1]){
      seven_t1 <- cif1[cif1[,1] == 7, 2]
    }else if(nrow(cif1[cif1[,1] <7, ]) > 0){
      seven_t1 <- cif1[cif1[,1] < 7, 2]
      seven_t1 <- max(seven_t1)
    }else{
      seven_t1 <- 0
    }
    if(7 %in% cif2[,1]){
      seven_t2 <- cif2[cif2[,1] == 7, 2]
    }else if(nrow(cif2[cif2[,1] <7, ]) > 0){
      seven_t2 <- cif2[cif2[,1] < 7, 2]
      seven_t2 <- max(seven_t2)
    }else{
      seven_t2 <- 0
    }
    if(7 %in% cif3[,1]){
      seven_t3 <- cif3[cif3[,1] == 7, 2]
    }else if(nrow(cif3[cif3[,1] <7, ]) > 0){
      seven_t3 <- cif3[cif3[,1] < 7, 2]
      seven_t3 <- max(seven_t3)
    }else{
      seven_t3 <- 0
    }
    
    seven_seq <- seq(7, 77, by = 7)
    multi_7_t1  <- cif1[cif1[,1] %in% seven_seq, 2]
    multi_7_t2  <- cif2[cif2[,1] %in% seven_seq, 2]
    multi_7_t3  <- cif3[cif3[,1] %in% seven_seq, 2]
    
    multi_7_t1 <- c(0, multi_7_t1)
    multi_7_t2 <- c(0, multi_7_t2)
    multi_7_t3 <- c(0, multi_7_t3)
    
    multi_7_t1 <- multi_7_t1[2:length(multi_7_t1)] - multi_7_t1[1:(length(multi_7_t1) -1)]
    multi_7_t2 <- multi_7_t2[2:length(multi_7_t2)] - multi_7_t1[1:(length(multi_7_t2) -1)]
    multi_7_t3 <- multi_7_t3[2:length(multi_7_t3)] - multi_7_t1[1:(length(multi_7_t3) -1)]
    
    
    
    
    cc_out_df$coeff[age_rows[1]] <- CI.agegrp1_t1$coef
    cc_out_df$variance[age_rows[1]] <- CI.agegrp1_t1$var[1,1]
    cc_out_df$day_7[age_rows[1]] <- seven_t1
    if(length(multi_7_t1) > 0)
      cc_out_df$mean_7[age_rows[1]] <- mean(multi_7_t1)
    else
      cc_out_df$mean_7[age_rows[1]] <- seven_t1
    if(length(multi_7_t1) > 0)
      cc_out_df$median_7[age_rows[1]] <- median(multi_7_t1)
    else
      cc_out_df$median_7[age_rows[1]] <- seven_t1
    
    
    
    cc_out_df$coeff[age_rows[2]] <- CI.agegrp1_t2$coef
    cc_out_df$variance[age_rows[2]] <- CI.agegrp1_t2$var[1,1]
    cc_out_df$day_7[age_rows[2]] <- seven_t2
    if(length(multi_7_t2) > 0)
      cc_out_df$mean_7[age_rows[2]] <- mean(multi_7_t2)
    else
      cc_out_df$mean_7[age_rows[2]] <- seven_t2
    if(length(multi_7_t2) > 0)
      cc_out_df$median_7[age_rows[2]] <- median(multi_7_t2)
    else
      cc_out_df$median_7[age_rows[2]] <- seven_t2
    
    
    cc_out_df$coeff[age_rows[3]] <- CI.agegrp1_t3$coef
    cc_out_df$variance[age_rows[3]] <- CI.agegrp1_t3$var[1,1]
    cc_out_df$day_7[age_rows[3]] <- seven_t3
    if(length(multi_7_t3) > 0)
      cc_out_df$mean_7[age_rows[3]] <- mean(multi_7_t3)
    else
      cc_out_df$mean_7[age_rows[3]] <- seven_t3
    if(length(multi_7_t3) > 0)
      cc_out_df$median_7[age_rows[3]] <- median(multi_7_t3)
    else
      cc_out_df$median_7[age_rows[3]] <- seven_t3
    
    
    
  }
  
  # run this one as an example because it takes less time than the above
  
  ## what we need from this output: coef, var
  
  
  return(list(out_df, cc_out_df))
}


survival_analysis_set_up <- function(cohorts_data, single_ICD = TRUE, base_dir, single_icd  ){
  require(tidyverse)
  require(survival)
  require(msSurv)
  require(plyr)
  require(timereg)
  require(ggplot2)
  require(ggpubr) ## need rtools installed for multi-plots
  require(haven)
  require(dplyr)
  require(tictoc)
  require(cmprsk)
  
  
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
    
    cohort12_failure_func <- failure_function(cohort23 = cohort12, csv_survial_name = failure_csv,
                                              pdf_to_export = failure_pdf)   
    cohort_3_transitions <- cohort_3_comp_risks(cohort3, ICD_group = single_icd, results_pdf = transitions_pdf)
    cohort1_transitions <- cohort_1_competing_risk(cohort1, current_ICD = single_icd)
    
    failure_res <- cohort12_failure_func
    cohort_3_ga <- cohort_3_transitions[[2]]
    cohort_3_cc <- cohort_3_transitions[[4]]
    cohort_1_res_ga <- cohort1_transitions[[1]]
    cohort_1_res_cc <- cohort1_transitions[[2]]
    
  
  }else{
    cols_needed <- c("WaitingTime","Elective2Emergency","GA_LoS","ga_transitions",
                     "agegrp_v3","cc_LoS","cc_transitions")
    cols_to_keep <- which(colnames(cohorts_data) %in% cols_needed)
    
    cohort1 <- subset(cohorts_data, subset = cohort == 1)
    cohort2 <- subset(cohorts_data, subset = cohort == 2)
    cohort3 <- subset(cohorts_data, subset = cohort == 3)
    cohort12 <- subset(cohorts_data, subset = cohort != 3)
    
    cohort_1_icds <- plyr::count(cohort1$ICD)
    cohort_3_icds <- plyr::count(cohort3$ICD)
    
    
    
    
    
    
  }
  
  return(list(failure_res, cohort_3_ga, cohort_3_cc, cohort_1_res_ga, cohort_1_res_cc))
}



