###############################################################################
## Quick crr for noeplasms ####################################################
###############################################################################

require(cmprsk)
require(fastcmprsk)
require(vroom)
require(tictoc)

transitions_data <- vroom("D:/Overflows/neo_trans_data_07_06_2020.csv",
                          delim = ",")


agegrp_3 <- transitions_data[transitions_data$agegrp_v3 == 3,]
transitions_3_grp <- agegrp_3[-which(is.na(agegrp_3$GA_LoS)),]
transitions_3_grp[which(is.na(agegrp_3$ga_transitions)),"ga_transitions"] <- 0

system.time(cif_2_ga <- fastCrr(Crisk(transitions_3_grp$GA_LoS, transitions_3_grp$ga_transitions, failcode = 1) ~ transitions_3_grp$WaitingTime,
                    variance = TRUE, returnDataFrame = TRUE))

agegrp_2 <- transitions_data[transitions_data$agegrp_v3 == 2,]
transitions_2_grp <- agegrp_2[-which(is.na(agegrp_2$GA_LoS)),]
transitions_2_grp[which(is.na(transitions_2_grp$ga_transitions)),"ga_transitions"] <- 0

system.time(cif_2_2_ga <- fastCrr(Crisk(transitions_2_grp$GA_LoS, transitions_2_grp$ga_transitions, failcode = 1) ~ transitions_2_grp$WaitingTime,
                                variance = TRUE, returnDataFrame = TRUE))

agegrp_1 <- transitions_data[transitions_data$agegrp_v3 == 1,]
transitions_1_grp <- agegrp_1[-which(is.na(agegrp_1$GA_LoS)),]
transitions_1_grp[which(is.na(transitions_1_grp$ga_transitions)),"ga_transitions"] <- 0

system.time(cif_1_1_ga <- fastCrr(Crisk(transitions_1_grp$GA_LoS, transitions_1_grp$ga_transitions, failcode = 1) ~ transitions_1_grp$WaitingTime,
                                  variance = TRUE, returnDataFrame = TRUE))

system.time(cif_1_2_ga <- fastCrr(Crisk(transitions_1_grp$GA_LoS, transitions_1_grp$ga_transitions, failcode = 2) ~ transitions_1_grp$WaitingTime,
                                  variance = TRUE, returnDataFrame = TRUE))

system.time(cif_1_2_ga_crr <- crr(ftime = transitions_1_grp$GA_LoS, fstatus = transitions_1_grp$ga_transitions,
                                  failcode = 2, cov1 = transitions_1_grp$WaitingTime))

## lets try recoding the ga_transitions for each type to run ##
agegrp_3_2_time <- agegrp_3$GA_LoS
na_tims <- which(is.na(agegrp_3_2_time))
agegrp_3_2_time <- agegrp_3_2_time[-na_tims]
agegrp_3_2_trans <- agegrp_3$ga_transitions
agegrp_3_2_trans <- agegrp_3_2_trans[-na_tims]
old_2ers <- which(agegrp_3_2_trans == 2)
old_1ers <- which(agegrp_3_2_trans == 1)
old_na <- which(is.na(agegrp_3_2_trans))
agegrp_3_2_trans[old_2ers] <- 1
agegrp_3_2_trans[old_1ers] <- 2
agegrp_3_2_trans[old_na] <- 0

system.time(cif_2_ga <- fastCrr(Crisk(agegrp_3_2_time, agegrp_3_2_trans, failcode = 1) ~ transitions_3_grp$WaitingTime,
                                variance = TRUE, returnDataFrame = TRUE))


cheeky_survival_func <- function(transitions_data, current_ICD = "2"){
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
  
  
  
  ages <- unique(transitions_data$agegrp_v3)
  
  for(k in 1:length(ages)){
    
    current_age_grp <- transitions_data[transitions_data$agegrp_v3 == ages[k],]
    transitions_1_grp <- current_age_grp[-which(is.na(current_age_grp$GA_LoS)),]
    transitions_1_grp[which(is.na(transitions_1_grp$ga_transitions)),"ga_transitions"] <- 0
    print(paste("Working on age group:", ages[k]))    
    tic("Working on trans 1 GA")
    
    cif_1_ga <- fastCrr(Crisk(transitions_1_grp$GA_LoS, transitions_1_grp$ga_transitions, failcode = 1) ~ transitions_1_grp$WaitingTime,
                        variance = TRUE, var.control = varianceControl(useMultipleCores = TRUE), returnDataFrame = TRUE)
    toc()
    tic("Transition 2 GA")
    agegrp_3_2_time <- current_age_grp$GA_LoS
    na_tims <- which(is.na(agegrp_3_2_time))
    agegrp_3_2_time <- agegrp_3_2_time[-na_tims]
    agegrp_3_2_trans <- current_age_grp$ga_transitions
    agegrp_3_2_trans <- agegrp_3_2_trans[-na_tims]
    old_2ers <- which(agegrp_3_2_trans == 2)
    old_1ers <- which(agegrp_3_2_trans == 1)
    old_na <- which(is.na(agegrp_3_2_trans))
    agegrp_3_2_trans[old_2ers] <- 1
    agegrp_3_2_trans[old_1ers] <- 2
    agegrp_3_2_trans[old_na] <- 0
    
    system.time(cif_2_ga <- fastCrr(Crisk(agegrp_3_2_time, agegrp_3_2_trans, failcode = 1) ~ transitions_3_grp$WaitingTime,
                                    variance = TRUE, returnDataFrame = TRUE))
    toc()
    tic("Transition 3 GA")
    agegrp_3_3_time <- current_age_grp$GA_LoS
    na_tims <- which(is.na(agegrp_3_3_time))
    agegrp_3_3_time <- agegrp_3_3_time[-na_tims]
    agegrp_3_3_trans <- current_age_grp$ga_transitions
    agegrp_3_3_trans <- agegrp_3_3_trans[-na_tims]
    old_3ers <- which(agegrp_3_3_trans == 3)
    old_1ers <- which(agegrp_3_3_trans == 1)
    old_na <- which(is.na(agegrp_3_3_trans))
    agegrp_3_3_trans[old_3ers] <- 1
    agegrp_3_3_trans[old_1ers] <- 2
    agegrp_3_3_trans[old_na] <- 0
    
    system.time(cif_3_ga <- fastCrr(Crisk(agegrp_3_3_time, agegrp_3_3_trans, failcode = 1) ~ transitions_3_grp$WaitingTime,
                                    variance = TRUE, returnDataFrame = TRUE))
    
    toc()
    tic("Transition 1 CC")
    
    cif_1_cc <- crr(ftime = cohort_data$ccLoS_cln, fstatus = cohort_data$cc_transitions,
                         cov1 = cohort_data$WaitingTime, failcode = 1)
    
    toc()
    tic("Transition 2 CC:")
    
    cif_2_cc <- crr(ftime = cohort_data$ccLoS_cln, fstatus = cohort_data$cc_transitions,
                         cov1 = cohort_data$WaitingTime, failcode = 2)
    
    toc()
    tic("Transition 3 CC:")
    cif_3_cc <- crr(ftime = cohort_data$ccLoS_cln, fstatus = cohort_data$cc_transitions,
                         cov1 = cohort_data$WaitingTime, failcode = 3)
    
    toc()
    
    age_rows <- c((k*3) -2, (k*3) -1, (k*3))
    tic("Predicting trans 1")
    cif1_pred <- predict(cif_1_ga, newdata = 0, tL = 1)
    toc()
    tic("Predicting trans 2")
    cif2_pred <- predict(cif_2_ga, newdata = 0, tL = 1)
    toc()
    tic("Predicting trans 3")
    cif3_pred <- predict(cif_3_ga, newdata = 0, tL = 1)
    toc()
    
    cif1 <- cbind(cif1_pred$ftime, cif1_pred$CIF)
    cif2 <- cbind(cif2_pred$ftime, cif2_pred$CIF)
    cif3 <- cbind(cif3_pred$ftime, cif3_pred$CIF)
    browser()    
    
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
    
    
    
    
    out_df$coeff[age_rows[1]] <- cif_1_ga$coef
    out_df$variance[age_rows[1]] <- cif_1_ga$var[1,1]
    out_df$day_7[age_rows[1]] <- seven_t1
    if(length(multi_7_t1) > 0)
      out_df$mean_7[age_rows[1]] <- mean(multi_7_t1)
    else
      out_df$mean_7[age_rows[1]] <- seven_t1
    if(length(multi_7_t1) > 0)
      out_df$median_7[age_rows[1]] <- median(multi_7_t1)
    else
      out_df$median_7[age_rows[1]] <- seven_t1
    
    
    
    out_df$coeff[age_rows[2]] <- cif_2_ga$coef
    out_df$variance[age_rows[2]] <- cif_2_ga$var[1,1]
    out_df$day_7[age_rows[2]] <- seven_t2
    if(length(multi_7_t2) > 0)
      out_df$mean_7[age_rows[2]] <- mean(multi_7_t2)
    else
      out_df$mean_7[age_rows[2]] <- seven_t2
    if(length(multi_7_t2) > 0)
      out_df$median_7[age_rows[2]] <- median(multi_7_t2)
    else
      out_df$median_7[age_rows[2]] <- seven_t2
    
    
    out_df$coeff[age_rows[3]] <- cif_3_ga$coef
    out_df$variance[age_rows[3]] <- cif_3_ga$var[1,1]
    out_df$day_7[age_rows[3]] <- seven_t3
    if(length(multi_7_t3) > 0)
      out_df$mean_7[age_rows[3]] <- mean(multi_7_t3)
    else
      out_df$mean_7[age_rows[3]] <- seven_t3
    if(length(multi_7_t3) > 0)
      out_df$median_7[age_rows[3]] <- median(multi_7_t3)
    else
      out_df$median_7[age_rows[3]] <- seven_t3
    
    
    age_rows <- c((k*3) -2, (k*3) -1, (k*3))
    
    tic("Predicting trans cc 1")
    cif1_cc <- predict.crr(cif_1_cc, cov1 = 0)
    toc()
    tic("Predicting trans cc 2")
    cif1_cc <- predict.crr(cif_2_cc, cov1 = 0)
    toc()
    tic("Predicting trans cc 3")
    cif1_cc <- predict.crr(cif_3_cc, cov1 = 0)
    toc()
    
    
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
    seven_seq <- seq(7, 77, by = 7)
    multi_7_t1_cc  <- cif1_cc[cif1_cc[,1] %in% seven_seq, 2]
    multi_7_t2_cc  <- cif2_cc[cif2_cc[,1] %in% seven_seq, 2]
    multi_7_t3_cc  <- cif3_cc[cif3_cc[,1] %in% seven_seq, 2]
    
    multi_7_t1_cc <- c(0, multi_7_t1_cc)
    multi_7_t2_cc <- c(0, multi_7_t2_cc)
    multi_7_t3_cc <- c(0, multi_7_t3_cc)
    
    multi_7_t1_cc <- multi_7_t1_cc[2:length(multi_7_t1_cc)] - multi_7_t1_cc[1:(length(multi_7_t1_cc) -1)]
    multi_7_t2_cc <- multi_7_t2_cc[2:length(multi_7_t2_cc)] - multi_7_t1_cc[1:(length(multi_7_t2_cc) -1)]
    multi_7_t3_cc <- multi_7_t3_cc[2:length(multi_7_t3_cc)] - multi_7_t1_cc[1:(length(multi_7_t3_cc) -1)]
    
    
    
    
    cc_out_df$coeff[age_rows[1]] <- cif_1_cc$coef
    cc_out_df$variance[age_rows[1]] <- cif_1_cc$var[1,1]
    cc_out_df$day_7[age_rows[1]] <- seven_t1_cc
    if(length(multi_7_t1_cc) > 0)
      cc_out_df$mean_7[age_rows[1]] <- mean(multi_7_t1_cc)
    else
      cc_out_df$mean_7[age_rows[1]] <- seven_t1_cc
    if(length(multi_7_t1_cc) > 0)
      cc_out_df$median_7[age_rows[1]] <- median(multi_7_t1_cc)
    else
      cc_out_df$median_7[age_rows[1]] <- seven_t1_cc
    
    
    
    cc_out_df$coeff[age_rows[2]] <- cif_2_cc$coef
    cc_out_df$variance[age_rows[2]] <- cif_2_cc$var[1,1]
    cc_out_df$day_7[age_rows[2]] <- seven_t2_cc
    if(length(multi_7_t2) > 0)
      cc_out_df$mean_7[age_rows[2]] <- mean(multi_7_t2_cc)
    else
      cc_out_df$mean_7[age_rows[2]] <- seven_t2_cc
    if(length(multi_7_t2) > 0)
      cc_out_df$median_7[age_rows[2]] <- median(multi_7_t2_cc)
    else
      cc_out_df$median_7[age_rows[2]] <- seven_t2_cc
    
    
    cc_out_df$coeff[age_rows[3]] <- cif_3_cc$coef
    cc_out_df$variance[age_rows[3]] <- cif_3_cc$var[1,1]
    cc_out_df$day_7[age_rows[3]] <- seven_t3_cc
    if(length(multi_7_t3_cc) > 0)
      cc_out_df$mean_7[age_rows[3]] <- mean(multi_7_t3_cc)
    else
      cc_out_df$mean_7[age_rows[3]] <- seven_t3_cc
    if(length(multi_7_t3_cc) > 0)
      cc_out_df$median_7[age_rows[3]] <- median(multi_7_t3_cc)
    else
      cc_out_df$median_7[age_rows[3]] <- seven_t3_cc
    
    
    
  }
  
  
  return(list(out_df, cc_out_df))
  
}

cheeky_survival_func(transitions_data)


