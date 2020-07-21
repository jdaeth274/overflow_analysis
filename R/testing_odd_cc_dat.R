###############################################################################
## 

require(fastcmprsk)
install.packages("cmprsk")

odd_cc_data <- read.csv(file = "~/Dropbox/Overflow/odd_icd3_age_3_cc_grouping.csv",
                        stringsAsFactors = FALSE)

odd_cc_data_no0 <- odd_cc_data[odd_cc_data$cc_LoS != 0,]

crr_res_1 <- crr(odd_cc_data_no0$cc_LoS, odd_cc_data_no0$cc_transitions, cov1 = odd_cc_data_no0$WaitingTime,
                 failcode = 1)
crr_res_2 <- crr(odd_cc_data_no0$cc_LoS, odd_cc_data_no0$cc_transitions, cov1 = odd_cc_data_no0$WaitingTime,
                 failcode = 2)
crr_res_3 <- crr(odd_cc_data_no0$cc_LoS, odd_cc_data_no0$cc_transitions, cov1 = odd_cc_data_no0$WaitingTime,
                 failcode = 3)


crr_pred_1 <- predict.crr(crr_res_1, cov1 = 0)
crr_pred_2 <- predict.crr(crr_res_2, cov1 = 0)
crr_pred_3 <- predict.crr(crr_res_3, cov1 = 0)



crr_res_1 <- crr(odd_cc_data$cc_LoS, odd_cc_data$cc_transitions, cov1 = odd_cc_data$WaitingTime,
               failcode = 1)
crr_res_2 <- crr(odd_cc_data$cc_LoS, odd_cc_data$cc_transitions, cov1 = odd_cc_data$WaitingTime,
                 failcode = 2)
crr_res_3 <- crr(odd_cc_data$cc_LoS, odd_cc_data$cc_transitions, cov1 = odd_cc_data$WaitingTime,
                 failcode = 3)


crr_pred_1 <- predict.crr(crr_res_1, cov1 = 0)
crr_pred_2 <- predict.crr(crr_res_2, cov1 = 0)
crr_pred_3 <- predict.crr(crr_res_3, cov1 = 0)

## cuminc 

CI.byagegrp <- cuminc(ftime = odd_cc_data$cc_LoS, fstatus = odd_cc_data$cc_transitions, group = odd_cc_data$agegrp_v3)

cuminc_df <- cuminc_list_to_df(CI.byagegrp, ICD = 3)

## alter the wait times for discharge from cc 

odd_cc_data_altered_dis_wait <- odd_cc_data
odd_cc_data_altered_dis_wait[which(odd_cc_data_altered_dis_wait$cc_transitions == 1),"WaitingTime"] <- odd_cc_data_altered_dis_wait[which(odd_cc_data_altered_dis_wait$cc_transitions == 1),"WaitingTime"] * 12

crr_res_1 <- crr(odd_cc_data_altered_dis_wait$cc_LoS, odd_cc_data_altered_dis_wait$cc_transitions,
                 cov1 = odd_cc_data_altered_dis_wait$WaitingTime,
                 failcode = 1)
crr_res_2 <- crr(odd_cc_data_altered_dis_wait$cc_LoS, odd_cc_data_altered_dis_wait$cc_transitions,
                 cov1 = odd_cc_data_altered_dis_wait$WaitingTime,
                 failcode = 2)
crr_res_3 <- crr(odd_cc_data_altered_dis_wait$cc_LoS, odd_cc_data_altered_dis_wait$cc_transitions,
                 cov1 = odd_cc_data_altered_dis_wait$WaitingTime,
                 failcode = 3)


crr_pred_1 <- predict.crr(crr_res_1, cov1 = 0)
crr_pred_2 <- predict.crr(crr_res_2, cov1 = 0)
crr_pred_3 <- predict.crr(crr_res_3, cov1 = 0)

crr_pred_1
crr_pred_2
crr_pred_3

crr_pred_1_WT <- predict.crr(crr_res_1, cov1 = seq(1,121))
crr_pred_2_WT <- predict.crr(crr_res_2, cov1 = seq(1,121))
crr_pred_3_WT <- predict.crr(crr_res_3, cov1 = seq(1,121))

crr_pred_wt_avg_1 <- rowMeans(crr_pred_1_WT[,2:ncol(crr_pred_1_WT)])
crr_pred_wt_avg_2 <- rowMeans(crr_pred_2_WT[,2:ncol(crr_pred_2_WT)])
crr_pred_wt_avg_3 <- rowMeans(crr_pred_3_WT[,2:ncol(crr_pred_3_WT)])

crr_pred_wt_avg_1
crr_pred_wt_avg_2
crr_pred_wt_avg_3


plot.predict.crr(crr_pred_1_WT)
plot.predict.crr(crr_pred_2_WT)


###############################################################################
## fastcmprsk #################################################################
###############################################################################



fail_1 <- fastCrr(Crisk(ftime = odd_cc_data$cc_LoS, fstatus = odd_cc_data$cc_transitions, failcode = 1) ~ odd_cc_data$WaitingTime,
        returnDataFrame = TRUE)
fail_1_pred <- predict(fail_1, newdata = odd_cc_data$WaitingTime,getBootstrapVariance = FALSE, tL = 1)
fail_1_pred$CIF


fail_2 <- fastCrr(Crisk(ftime = odd_cc_data$cc_LoS, fstatus = odd_cc_data$cc_transitions, failcode = 2) ~ odd_cc_data$WaitingTime,
                  returnDataFrame = TRUE)
fail_2_pred <- predict(fail_2, newdata = 0,getBootstrapVariance = FALSE)

fail_3 <- fastCrr(Crisk(ftime = odd_cc_data$cc_LoS, fstatus = odd_cc_data$cc_transitions, failcode = 3) ~ odd_cc_data$WaitingTime,
                  returnDataFrame = TRUE)
fail_3_pred <- predict(fail_3, newdata = 0,getBootstrapVariance = FALSE)


fail_1_pred$CIF
fail_2_pred$CIF
fail_3_pred$CIF

###############################################################################
## Lets have a quick go at this multinomial logisitic regression ##############
###############################################################################


require(foreign)
require(nnet)
require(ggplot2)
require(reshape2)


odd_cc_data$trans_2 <- odd_cc_data$cc_transitions
odd_cc_data$trans_2[odd_cc_data$cc_LoS > 7] <- 4
odd_cc_data$trans_2 <- as.factor(odd_cc_data$trans_2)

odd_cc_data$trans_2 <- relevel(odd_cc_data$trans_2, ref = 4)

test_multino <- multinom(trans_2 ~ WaitingTime, data = odd_cc_data)

probs_out <- fitted(test_multino)


## one age grp ##

odd_cc_data$trans_2 <- relevel(odd_cc_data$trans_2, ref = 4)

test_multino <- multinom(trans_2 ~ WaitingTime, data = odd_cc_data)

probs_out <- fitted(test_multino)




