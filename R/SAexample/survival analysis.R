# Title: Survival Analysis for OverFlow Deaths Project
# Last updated: 30 May 2020

# install necessary packages
install.packages("tidyverse")
install.packages("survival")
install.packages("msSurv")
install.packages("ggpubr")
install.packages("plyr")
require(ggplot2)
require(ggpubr)

# set working directory
setwd("D:/Overflows/kl1215")

# import survival analysis data
library(haven)
data <- read_dta("survival analysis.dta")

# isolate cohorts
cohort1 <- subset(data, subset = cohort1 == 1)
cohort2 <- subset(data, subset = cohort2 == 1)
cohort3 <- subset(data, subset = cohort3 == 1)
cohort23 <- subset(data, subset = cohort1 !=1)

# 1. Failure function for Electives to Emergencies (Cohort = 2 & 3, event = Elective2Emergency == 1, time = WaitingTime)
library(survival)
library(dplyr)
surv_object <- Surv(time = cohort23$WaitingTime, event = cohort23$Elective2Emergency)
fit <- survfit(surv_object ~ agegrp_v3, data = cohort23, type = "kaplan-meier")
summary(fit)

# Plots survival function (CURVLAB NOT WORKING - JOSH, MAKE THIS FOR GGPLOT2 - BY TRANSITION TYPE)
survival_df <- cbind.data.frame(fit$time, fit$surv, rep(c("<24", "25-64", "65+"), each = length(fit$surv)/3),
                                fit$cumhaz)
colnames(survival_df) <- c("time","survival","age_group", "cumulative_hazard")

survival_plot_2 <- ggplot(data = survival_df, aes(x = time, y = survival, group = age_group, colour = age_group)) +
        geom_line() + xlab("Waiting time (days)") + ylab("Survival probability") + ylim(c(0,1)) + labs(colour = "Age group") +
        ggtitle("Survival function")
survival_plot_2

plot(fit, main = "survival function", 
     col = c("black","blue","green"),
             curvlab = c("<25", "25-64", "65+"),
             xlab = "Waiting Time (days)", 
             ylab = "Survival Probability")

# plots failure function

cumhaz_plot_2 <- ggplot(data = survival_df, aes(x = time, y = cumulative_hazard, group = age_group, colour = age_group)) +
        geom_line() + xlab("Waiting time (days)") + ylab("Failure probability") + ylim(c(0,1)) + labs(colour = "Age group") +
        ggtitle("Failure function (cumulative hazard)")
cumhaz_plot_2


plot(fit, fun="event", main = "Failure function (cumulative hazard)",
     col = c("black","blue","green"),
     curvlab = c("<25", "25-64", "65+"),
     xlab = "Waiting Time (days)", 
     ylab = "Failure Probability")

# puts failure function and standard errors in a dataframe
failurefunction <- data.frame(do.call(rbind, Map(cbind, fit[["time"]], fit[["cumhaz"]], fit[["lower"]], fit[["upper"]])))
failurefunction <- failurefunction %>% rename(Time = X1, FailureFunction = X2 , LowerCI = X3, UpperCI = X4)
failurefunction$agegrp <- c(rep(1:3, each = 366))
failurefunction$icd <- 2

# Use these probabilities to find cohort 2 from cohort 3
XXXXXXXX


# 2. Competing Risks (use cmprsk for cohort 1 with no WT covariate; use crr for cohort 3 with WT covariate)
# LOOPING: by ICD, agegrp, transition number, vector of GA_LoS & ga_transitions vs. vector of cc_LoS & cc_transitions)
# 2.1 Cohort 1 (True Emergencies) GA transitions (Cohort = 1, event = ga_transitions == 1, 2, or 3, time = GA_LoS)
library(cmprsk)
library(msSurv)
# give frequency of different GA transitions
table(cohort1$ga_transitions)                                            

# set survival datatset (grouped by age groups)
CI.byagegrp <- cuminc(ftime = cohort1$GA_LoS, fstatus = cohort1$ga_transitions, group = cohort1$agegrp_v3)

# plots cumulative hazard functions by transition type and age-group (JOSH - PUT THIS IN GGPLOT2 - BY TRANSITION TYPE)

cuminc_list_to_df <- function(CI_dataset, group_names = c("Recovered, <25", "Recovered, 25-64", "Recovered, 65+",
                                                          "CC, <25", "CC, 25-64", "CC, 65+","Died, <25",
                                                          "Died, 25-64","Died, 65+")){

        tests_pos <- grep("Tests", names(CI_dataset), ignore.case = TRUE)
        
        cuminc_df <- NULL
        
        for(k in 1:(tests_pos - 1)){
                current_title <- group_names[k]
                current_age <-stringr::str_split_fixed(current_title, ",", 2)[2]
                current_transition <- stringr::str_split_fixed(current_title, ",", 2)[1]
                current_list_out <- CI_dataset[[k]]
                new_rowz <- cbind.data.frame(current_list_out$time, current_list_out$est,
                                             rep(group_names[k], length(current_list_out$time)),
                                             rep(current_age, length(current_list_out$time)),
                                             rep(current_transition, length(current_list_out$time)))
                colnames(new_rowz) <- c("time","est",'age_tran', "age","transition")
                cuminc_df <- rbind.data.frame(cuminc_df, new_rowz)
                
        }
        
        return(cuminc_df)
}

cuminc_df <- cuminc_list_to_df(CI.byagegrp)

cuminc_plot <- ggplot(data = cuminc_df, aes(x = time, y = est, group = age_tran)) +
        geom_line(aes(linetype = transition, colour = age)) + xlab("Length of stay (days)") +
        ylab("Cumulative transition probability")
cuminc_plot

## Just recovered now ##

recovered_df <- cuminc_df[grep("recovered", cuminc_df$transition, ignore.case = TRUE),]

recovered_plot <- ggplot(data = recovered_df, aes(x = time, y = est, group = age_tran)) +
        geom_line(aes(colour = age)) + xlab("Length of stay (days)") +
        ylab("Cumulative Recovery probability") + ggtitle("Recovery probabilities")
recovered_plot

## Just CC transition ##

cc_df <- cuminc_df[grep("CC", cuminc_df$transition, ignore.case = TRUE),]

cc_plot <- ggplot(data = cc_df, aes(x = time, y = est, group = age_tran)) +
        geom_line(aes(colour = age)) + xlab("Length of stay (days)") +
        ylab("Cumulative move to CC probability") + ggtitle("G&A > CC probabilities")
cc_plot

## Just Died ##

died_df <- cuminc_df[grep("died", cuminc_df$transition, ignore.case = TRUE),]

died_plot <- ggplot(data = died_df, aes(x = time, y = est, group = age_tran)) +
        geom_line(aes(colour = age)) + xlab("Length of stay (days)") +
        ylab("Cumulative Death probability") + ggtitle("Death probabilities")
died_plot



plot(CI.byagegrp, lty = c(1,1,2,2,3,3),
     col = c("black","blue","black","blue","black","blue"),
     curvlab = c("Recovered, <25", "Recovered, 25-64", "Recovered, 65+", "CC, <25", "CC, 25-64", "CC, 65+",
                 "Died, <25", "Died, 25-64","Died, 65+"), xlab = "Length of Stay (days)")

# Tests for significant differences across something (COMPETING RISKS? AGE GROUPS? ASK PABLO)
CI.byagegrp$Tests

# (Age = <25, Event = Recovered)
CI_Transition1_agegrp1 <- data.frame(do.call(rbind, Map(cbind, CI.byagegrp[["1 1"]][["time"]],
                                            CI.byagegrp[["1 1"]][["est"]], CI.byagegrp[["1 1"]][["var"]])))
CI_Transition1_agegrp1 <- CI_Transition1_agegrp1 %>% rename(Time = X1 , CIF_T1 = X2, var_T1 = X3)
# drops duplicates (WHY SO MANY DUPLICATES? ASK PABLO)
CI_Transition1_agegrp1 <- CI_Transition1_agegrp1[!duplicated(CI_Transition1_agegrp1$CIF),]
# drop first row of zeros (WHY DOES THIS EXIST? ASK PABLO)
CI_Transition1_agegrp1 <- CI_Transition1_agegrp1[-1,]

# (Age = <25, Event = CC)
CI_Transition2_agegrp1 <- data.frame(do.call(rbind, Map(cbind, CI.byagegrp[["1 2"]][["time"]],
                                                        CI.byagegrp[["1 2"]][["est"]], CI.byagegrp[["1 2"]][["var"]])))
CI_Transition2_agegrp1 <- CI_Transition2_agegrp1 %>% rename(Time = X1 , CIF_T2 = X2, var_T2 = X3)
# drops duplicates (WHY SO MANY DUPLICATES? ASK PABLO)
CI_Transition2_agegrp1 <- CI_Transition2_agegrp1[!duplicated(CI_Transition2_agegrp1$CIF),]

# (Age = <25, Event = Died)
CI_Transition3_agegrp1 <- data.frame(do.call(rbind, Map(cbind, CI.byagegrp[["1 3"]][["time"]],
                                                        CI.byagegrp[["1 3"]][["est"]], CI.byagegrp[["1 3"]][["var"]])))
CI_Transition3_agegrp1 <- CI_Transition3_agegrp1 %>% rename(Time = X1 , CIF_T3 = X2, var_T3 = X3)
# drops duplicates (WHY SO MANY DUPLICATES? ASK PABLO)
CI_Transition3_agegrp1 <- CI_Transition3_agegrp1[!duplicated(CI_Transition3_agegrp1$CIF),]

# puts all Age = < 25 events in one df (NOT PICKING UP TRANSITION3 - JOSH FIND A WAY TO PUT THESE 3 DFS TOGETHER)
CI_agegrp1 <- merge(CI_Transition1_agegrp1,CI_Transition2_agegrp1, 
                    by.x = "Time",by.y = "Time", all = TRUE)
CI_agegrp1 <- merge(CI_agegrp1, CI_Transition3_agegrp1,
                    by.x = "Time",by.y = "Time", all = TRUE)

colnames(CI_agegrp1)

# DO THE SAME FOR AGE GROUPS 2 AND 3
# TAKE 7TH DAY LOS AND MEAN & MEDIAN OF 7TH DAY MULTIPLES LOS TRANSITION PROBABILITIES --> PUT IN NEW DF

multi_7 <- seq(7,max(CI_agegrp1$Time), by = 7)

los_df <- data.frame(matrix(ncol = 5, nrow = 3))
colnames(los_df) <- c("Transition","Day_7","mean_7_multiple","median_7_multiple", "age")
los_df$age <- c("<24","<24","<24")
los_df$Transition <- c("recovered","CC","die")

los_df$Day_7[1:3] <- as.numeric(CI_agegrp1[CI_agegrp1$Time == 7,c(2,4,6)])

average_df <- CI_agegrp1[CI_agegrp1$Time %in% multi_7, c(2,4,6)]

los_df$mean_7_multiple <- colMeans(average_df[sapply(average_df, is.numeric)], na.rm = TRUE)
los_df$median_7_multiple[1] <- median(average_df$CIF_T1, na.rm = TRUE)
los_df$median_7_multiple[2] <- median(average_df$CIF_T2, na.rm = TRUE)
los_df$median_7_multiple[3] <- median(average_df$CIF_T3, na.rm = TRUE)

los_df



# DO THE SAME FOR CC TRANSITIONS FOR ALL AGE-GROUPS

# Cohort 3 - Electives staying Electives
# give frequency of different GA transitions
table(cohort3$ga_transitions) 

# set survival datatset (grouped by age groups)
CI.byagegrp_cohort3 <- crr(ftime = cohort3$GA_LoS, fstatus = cohort3$ga_transitions, cov1 = cohort3$WaitingTime, cengroup = cohort3$agegrp_v3)



















