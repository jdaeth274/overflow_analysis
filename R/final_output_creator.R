###############################################################################
## Script to sum up results and output excel sheet ############################
###############################################################################

sum_up_function <- function(reg_res, time_series_data_res, time_series_forecasts, forecast_length, COVID_preds,
                            out_sheet){
  require(dplyr)  
  require(openxlsx)

  ## Forecasted_admissions 
  elec_admis <- time_series_forecasts[[4]]
  emerg_admis <- time_series_forecasts[[1]]
  
  
  icds_in_play <- union(elec_admis$icd_name, emerg_admis$icd_name)
  
  icd_cols <- rep(icds_in_play, each = 3)
  age_cols <- rep(c(1,2,3), length(icds_in_play))
  for(k in 1:length(icd_cols)){
    if(nchar(icd_cols[k]) == 1)
      icd_cols[k] <- paste("0",icd_cols[k], sep = "")
  }
  phi_colnames <- paste("ICD",icd_cols,"_AGE",age_cols, sep = "")
  phi_icds <- phi_colnames
  phi_colnames <- c("t","a",phi_colnames)
  phi_COVID <- c(phi_colnames, "COVID_AGE1", "COVID_AGE2","COVID_AGE3" )
  
  t_col <- rep(seq(0,length.out = forecast_length), 2)
  a_col <- rep(c("N","E"), each = forecast_length)

  phi_df <- data.frame(matrix(0, nrow = length(t_col), ncol = length(phi_COVID)))
  colnames(phi_df) <- phi_COVID
  phi_df$t <- t_col 
  phi_df$a <- a_col
  
  current_col <- 3
  for(icdee in icds_in_play){
    for(tiempo in c("1","2","3")){
      current_elec_df <- elec_admis[elec_admis$icd_name == as.character(icdee) &
                                      elec_admis$age == tiempo,]
      current_emerg_df <- emerg_admis[emerg_admis$icd_name == as.character(icdee) &
                                        emerg_admis$age == tiempo,]
      if(nrow(current_elec_df) != 0){
        phi_df[1:forecast_length,current_col] <- current_elec_df$median
      }
      if(nrow(current_emerg_df) != 0){
        if(length(current_emerg_df) == 208)
          browser()
        phi_df[(forecast_length + 1):(2*forecast_length), current_col] <- current_emerg_df$median
      }
      
      print(paste(icdee,tiempo, sep = "-"))
      print(phi_colnames[current_col])
      
      current_col <- current_col + 1
      
    }
  }
  
  
  phi_df$COVID_AGE1[(forecast_length + 1):(2*forecast_length)] <- COVID_preds[1:52,2]
  phi_df$COVID_AGE2[(forecast_length + 1):(2*forecast_length)] <- COVID_preds[1:52,3]
  phi_df$COVID_AGE3[(forecast_length + 1):(2*forecast_length)] <- COVID_preds[1:52,4]
  
  phi2_df <- phi_df
  phi2_df$COVID_AGE1[(forecast_length + 1):(2*forecast_length)] <- COVID_preds[1:52,5]
  phi2_df$COVID_AGE2[(forecast_length + 1):(2*forecast_length)] <- COVID_preds[1:52,6]
  phi2_df$COVID_AGE3[(forecast_length + 1):(2*forecast_length)] <- COVID_preds[1:52,7]
  
  phi3_df <- phi_df
  phi2_df$COVID_AGE1[(forecast_length + 1):(2*forecast_length)] <- COVID_preds[1:52,8]
  phi2_df$COVID_AGE2[(forecast_length + 1):(2*forecast_length)] <- COVID_preds[1:52,9]
  phi2_df$COVID_AGE3[(forecast_length + 1):(2*forecast_length)] <- COVID_preds[1:52,10]
  
  
  ## pi_x df for transitions ##
  
  failure_func_df <- reg_res[[3]]
  
  pi_x_df <- data.frame(matrix(0,nrow = 2*length(phi_icds),ncol = 3))
  colnames(pi_x_df) <- c("a","p","pi_x")
  pi_x_df$a <- rep(c("N","E"), each = length(phi_icds))
  pi_x_df$p <- rep(phi_icds, 2)
  
  current_row <- 1
  for(icdee in icds_in_play){
    for(age in c(1,2,3)){
      failure_row <- failure_func_df[failure_func_df$age == age &
                                       failure_func_df$ICD == icdee,]
      if(nrow(failure_row) == 1){
        pi_x_df$pi_x[current_row] <- failure_row$mean_7
      }
      
      current_row <- current_row + 1
    }
  }
  
  
  ## pi_y for transitions ##
  print("Working on the pi y transitions")
  elec_res <- reg_res[[1]][[1]]
  emerg_res <- reg_res[[2]][[1]]
  coef_df <- reg_res[[1]][[2]]
  
  ga_emerg <- grep("ga",x = emerg_res$patient_group)
  emerg_ga_sa <- emerg_res[ga_emerg,]
  emerg_cc_sa <- emerg_res[-ga_emerg,]
  
  ga_elec <- grep("ga", x = elec_res$patient_group)
  elec_ga_sa <- elec_res[ga_elec,]
  elec_cc_sa <- elec_res[-ga_elec,]
  
  ga_coef_rows <- grep("ga", x = coef_df$patient_group)
  ga_coef <- coef_df[ga_coef_rows,]
  cc_coef <- coef_df[-ga_coef_rows,]
  
  pi_y_a_row <- rep(c("N","E"), each = 4*length(pi_x_df$p))
  pi_y_p_row <- rep(pi_x_df$p, each = 8)
  pi_y_s_row <- rep(rep(c("G","C"), each = 4), length(pi_x_df$p))
  pi_y_sbar_row <- rep(c("H","C","D","G","H","G","D","C"), length(pi_x_df$p))
  
  
  pi_y_df <- data.frame(matrix(0, nrow = length(pi_y_p_row), ncol = 7))
  colnames(pi_y_df) <- c("a","p","s","sbar","pi_y","coeff","variance")
  pi_y_df$a <- pi_y_a_row
  pi_y_df$p <- pi_y_p_row
  pi_y_df$s <- pi_y_s_row
  pi_y_df$sbar <- pi_y_sbar_row
  
  current_row <- 1
  elec_emerg_break <- 3*length(icds_in_play)
  
  trans_icds_in_play <- rep(icds_in_play, 2)
  
  for(icdee in trans_icds_in_play){
    for(age in c(1,2,3)){
      current_rows <- seq(current_row*8 - 7, current_row*8) ## Looking at the rows for the current age & icd comboi
      if(current_row <= elec_emerg_break){ ## Checking if we're in the elective or emergency section of Results
        if(icdee %in% elec_ga_sa$ICD){ ## Check if comboi is present in GA and CC res separtely, just in case we decide to drop symmetry in the future
          icdee_rows <- elec_ga_sa[elec_ga_sa$ICD == icdee,]
          if(age %in% icdee_rows$age){
            age_rows_ga <- icdee_rows[icdee_rows$age == age,]
            seven_rows <- age_rows_ga[grep("seven",age_rows_ga$WT),]
            Discharge_trans <- seven_rows[seven_rows$variable == "Discharged","value"][1]
            CC_trans <- seven_rows[seven_rows$variable == "CC","value"][1]
            death_trans <- seven_rows[seven_rows$variable == "Dead","value"][1]
            GA_trans <- seven_rows[seven_rows$variable == "GA","value"][1]
            pi_y_df[current_rows[1:4],"pi_y"] <- c(Discharge_trans,CC_trans, death_trans, GA_trans)
            if(icdee %in% ga_coef$ICD){
              icdee_rows <- ga_coef[ga_coef$ICD == icdee,]
              if(age %in% icdee_rows$age){
                age_rows_ga <- icdee_rows[icdee_rows$age == age,]
                pi_y_df[current_rows[1:4],"coeff"] <- age_rows_ga$coef[c(3,2,4,1)]
              }
            }
            
          }
          
          
        }
        if(icdee %in% elec_cc_sa$ICD){
          icdee_rows_cc <- elec_cc_sa[elec_cc_sa$ICD == icdee,]
          if(age %in% icdee_rows_cc$age){
            age_rows_ga <- icdee_rows_cc[icdee_rows_cc$age == age,]
            seven_rows <- age_rows_ga[grep("seven",age_rows_ga$WT),]
            Discharge_trans <- seven_rows[seven_rows$variable == "Discharged","value"][1]
            CC_trans <- seven_rows[seven_rows$variable == "CC","value"][1]
            death_trans <- seven_rows[seven_rows$variable == "Dead","value"][1]
            GA_trans <- seven_rows[seven_rows$variable == "GA","value"][1]
            pi_y_df[current_rows[5:8],"pi_y"] <- c(Discharge_trans,GA_trans, death_trans, CC_trans)
            if(icdee %in% cc_coef$ICD){
              icdee_rows <- cc_coef[cc_coef$ICD == icdee,]
              if(age %in% icdee_rows$age){
                age_rows_cc <- icdee_rows[icdee_rows$age == age,]
                pi_y_df[current_rows[5:8],"coeff"] <- age_rows_cc$coef[c(3,2,4,1)]
              }
            }
            
          }
        }
      
      
      }else{
        if(icdee %in% emerg_ga_sa$ICD){
          icdee_rows_ga <- emerg_ga_sa[emerg_ga_sa$ICD == icdee,]
          if(age %in% icdee_rows_ga$age){
            age_rows_ga <- icdee_rows_ga[icdee_rows_ga$age == age,]
            seven_rows <- age_rows_ga[grep("mean",age_rows_ga$WT),]
            Discharge_trans <- seven_rows[seven_rows$variable == "Discharged","value"][1]
            CC_trans <- seven_rows[seven_rows$variable == "CC","value"][1]
            death_trans <- seven_rows[seven_rows$variable == "Dead","value"][1]
            GA_trans <- seven_rows[seven_rows$variable == "GA","value"][1]
            pi_y_df[current_rows[1:4],"pi_y"] <- c(Discharge_trans,CC_trans, death_trans, GA_trans)
            
          }
        }
        if(icdee %in% emerg_cc_sa$ICD){
          icdee_rows_cc <- emerg_cc_sa[emerg_cc_sa$ICD == icdee,]
          if(age %in% icdee_rows_cc$age){
            age_rows_cc <- icdee_rows_cc[icdee_rows_cc$age == age,]
            seven_rows <- age_rows_cc[grep("mean",age_rows_cc$WT),]
            Discharge_trans <- seven_rows[seven_rows$variable == "Discharged","value"][1]
            CC_trans <- seven_rows[seven_rows$variable == "CC","value"][1]
            death_trans <- seven_rows[seven_rows$variable == "Dead","value"][1]
            GA_trans <- seven_rows[seven_rows$variable == "GA","value"][1]
            pi_y_df[current_rows[5:8],"pi_y"] <- c(Discharge_trans,GA_trans, death_trans, CC_trans)
            
          }
        }
        
        
        
      }
      
      print(paste(icdee,age, pi_y_df[current_rows[1],"p"]))
      current_row <- current_row + 1
      
  }
  
  
  }

  
  ## pi z proportions ##
  print("On the pi z transitions")
  emerg_cc_props <- time_series_forecasts[[3]]
  elec_cc_props <- time_series_forecasts[[8]]
  
  t_col <- rep(seq(0,length.out = forecast_length), 2)
  a_col <- rep(c("N","E"), each = forecast_length)
  
  pi_z_df <- data.frame(matrix(0, nrow = length(t_col), ncol = length(phi_colnames)))
  colnames(pi_z_df) <- phi_colnames
  pi_z_df$t <- t_col
  pi_z_df$a <- a_col
  
  current_col <- 3
  for(icdee in icds_in_play){
    for(tiempo in c("1","2","3")){
      current_elec_df <- elec_cc_props[elec_cc_props$icd_name == as.character(icdee) &
                                         elec_cc_props$age == tiempo,]
      current_emerg_df <- emerg_cc_props[emerg_cc_props$icd_name == as.character(icdee) &
                                           emerg_cc_props$age == tiempo,]
      if(nrow(current_elec_df) != 0){
        pi_z_df[1:forecast_length,current_col] <- current_elec_df$median
      }
      if(nrow(current_emerg_df) != 0){
        pi_z_df[(forecast_length + 1):(2*forecast_length), current_col] <- current_emerg_df$median
      }
      
      print(paste(icdee,tiempo, sep = "-"))
      print(phi_colnames[current_col])
      
      current_col <- current_col + 1
      
    }
  }
  
  
  ## bundled icd proportions ##
  print("on the bundle now")
  
  elec_bundles <- time_series_forecasts[[9]]
  emerg_bundles <- time_series_forecasts[[10]]
  icds_in_play_bundle <- union(elec_bundles$icd_name, emerg_bundles$icd_name)
  
  icd_cols <- rep(icds_in_play_bundle, each = 3)
  age_cols <- rep(c(1,2,3), length(icds_in_play_bundle))
  for(k in 1:length(icd_cols)){
    if(nchar(icd_cols[k]) == 1)
      icd_cols[k] <- paste("0",icd_cols[k], sep = "")
  }
  prop_colnames <- paste("ICD",icd_cols,"_AGE",age_cols, sep = "")
  
  prop_colnames <- c("t","a",prop_colnames)
  
  
  t_col <- rep(seq(0,length.out = forecast_length), 2)
  a_col <- rep(c("N","E"), each = forecast_length)
  
  prop_df <- data.frame(matrix(0, nrow = length(t_col), ncol = length(prop_colnames)))
  colnames(prop_df) <- prop_colnames
  
  prop_df$t <- t_col
  prop_df$a <- a_col
  
  
  current_col <- 3
  for(icdee in icds_in_play_bundle){
    for(tiempo in c("1","2","3")){
      # if(icdee == 15 & tiempo == "3")
      #   browser()
      #   
      current_elec_df <- elec_bundles[elec_bundles$icd_name == as.character(icdee) &
                                        elec_bundles$age == tiempo,]
      current_emerg_df <- emerg_bundles[emerg_bundles$icd_name == as.character(icdee) &
                                        emerg_bundles$age == tiempo,]
      if(nrow(current_elec_df) != 0){
        prop_df[1:forecast_length,current_col] <- current_elec_df$median
      }
      if(nrow(current_emerg_df) != 0){
        if(length(current_emerg_df$median) == 208)
          browser()
        prop_df[(forecast_length + 1):(2*forecast_length), current_col] <- current_emerg_df$median
      }
      
      print(paste(icdee,tiempo, sep = "-"))
      print(prop_colnames[current_col])
      
      current_col <- current_col + 1
      
    }
  }
  
  ## x0 pool of patients 
  print("Waiting pool")
  waiting_pool <- time_series_data_res[[3]]
  x0_df <- data.frame(matrix(0,nrow = 2*length(phi_icds),ncol = 3))
  colnames(x0_df) <- c("a","p","x0")
  x0_df$a <- rep(c("N","E"), each = length(phi_icds))
  x0_df$p <- rep(phi_icds, 2)
  
  current_row <- 1
  for(icdee in icds_in_play){
    for(age in c("<25","25-64","65+")){
      waiting_row <- waiting_pool[waiting_pool$age == age &
                                    waiting_pool$icd == icdee,]
      if(nrow(waiting_row) == 1){
        x0_df$x0[current_row] <- waiting_row$pool
      }
      
      current_row <- current_row + 1
    }
  }
  
  ## y0 pool ##
  print("In hospital pool")
  
  in_hops_patients <- time_series_data_res[[6]]
  
  y0_df <- data.frame(matrix(0, nrow = 4*length(phi_icds), ncol = 4))
  colnames(y0_df) <- c("a","p","s","y0")
  y0_df$a <- rep(c("N","E"), each = 2*length(phi_icds))
  y0_df$p <- rep(phi_icds, 4)
  y0_df$s <- rep(c("G","C","G","C"), each = length(phi_icds))
  y0_df$combo <- paste(y0_df$a, y0_df$p, y0_df$s, sep = "-")
  
  in_hops_patients$combo <- paste(in_hops_patients$a, in_hops_patients$p, in_hops_patients$s, sep = "-")
  colnames(in_hops_patients) <- c("admit","patient_group","ward","one","combo")
  
  y0_df <- dplyr::left_join(y0_df, in_hops_patients, by = c("combo" = "combo"))
  y0_df <- y0_df[,c(1,2,3,9)]
  colnames(y0_df)[4] <- "y0"
  y0_df$y0[which(is.na(y0_df$y0))] <- 0
  
  
  
  df_list <- list("phi1" = phi_df, "phi2" = phi2_df, "phi3" = phi3_df,
                  "pi_x" = pi_x_df, "pi_y" = pi_y_df,
                  "pi_z" = pi_z_df, "ICDprop" = prop_df,
                  "x0" = x0_df, "y0" = y0_df)
  
  write.xlsx(df_list, file = out_sheet)
  
  
  return(pi_y_df)
  
  
}
