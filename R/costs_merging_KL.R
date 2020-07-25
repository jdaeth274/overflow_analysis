###############################################################################
## Costing analysis ###########################################################
###############################################################################
require(stringr)
require(dplyr)
require(vroom)
require(dtplyr)
require(data.table)

## Load HES data
system.time(transitions_data_whole <- fread("D:/Overflows/data/HES_APC_CC_0913_transitions_all_ICD.csv"))
# keep variables you need
hes_data <- select(transitions_data_whole, c("procode3","admimeth_C","EpiEnd_FY","SUSHRGep", "ICD", "agegrp_v3", "admimeth_C"))

##############################################################################################################
## Load costs data
# Load organisational level unit costs
system.time(costs_dataset <- vroom("D:/Overflows/data/cleaned_1213_orgrcs.csv", delim = ","))

# Load previous years reference costs to match HES admissions that are left unmatched
#load national ref cost data for 1213 - for costs that are unmatched using organisational level cost data
national_rcs1213 <- read.csv("D:/Overflows/data/national_rcs1213.csv",fileEncoding = "UTF-8-BOM")

# load national ref cost data for 1112 - for costs that are unmatched after using national_rcs1213
national_rcs1112 <- read.csv("D:/Overflows/data/national_rcs1112.csv",fileEncoding = "UTF-8-BOM")

# load national ref cost data for 1011 - for costs that are unmatched after using national_rcs1112
national_rcs1011 <- read.csv("D:/Overflows/data/national_rcs1011.csv",fileEncoding = "UTF-8-BOM")

# load national ref cost data for 0910 - for costs that are unmatched after using national_rcs1011
national_rcs0910 <- read.csv("D:/Overflows/data/national_rcs0910.csv",fileEncoding = "UTF-8-BOM")

##############################################################################################################
# Notes about merging HES and costs data:
# ORG.CODE matched to procode3 
# Currency.code matched to HRG or susHRG (this is for post FY1213 onwards. will need to include a change in variable name for HF pre cleaning file)
# we need to keep the UNIT COST, MEAN COST, Actual Cost MAPPING POT and PERIOD (2009-10) EpiEnd_FY== "0910"
# mean cost = the national mean average unit cost
# unit cost = the average cost to the organisation of providing that activity
# actual cost = the organisation's activity (number of cases) x the organisation unit cost
# we are really only interested in matching the unit cost of each organisation's hrg
# Mapping pot is matched to admimeth_C
# Mapping to 01_EI and 02_NEI and 03_XS
# XS bedays with the bed.days variable

##############################################################################################################

# Pre-process HES data
# Only keep FY1213 in transitions_data_whole as it is the only year with SUSHRGep code (cost dataset has HRG as HRG4.0, imperial hes has hrg 3.5)
hes_data_1213 <- subset(hes_data, EpiEnd_FY == "1213")

# Now merge HES and organisational reference costs
mergedcosts_orglevel <- left_join(hes_data_1213, costs_dataset, 
                                  by = c("procode3" = "orgcode",
                                         "admimeth_C" = "admimeth_c",
                                         "EpiEnd_FY" = "EpiEnd_FY",
                                         "SUSHRGep" = "currencycode"))

# Not all HES episodes are matched with a cost, so we next move onto the national ref costs 1213
# subset MATCHED episodes @ org level
matched_orglevel <- mergedcosts_orglevel[!is.na(mergedcosts_orglevel$unitcost),]

# subset UNMATCHED episodes @ org level (procode3-SUSHRGep pair does not exist in costs data but does exist in hes data)
unmatched_orglevel <- mergedcosts_orglevel[is.na(mergedcosts_orglevel$unitcost),]

# remove UNIT.COST column from unmatched_orglevel
unmatched_orglevel <- select(unmatched_orglevel, -unitcost)

# match unmatched_orglevel with national_rcs1213
merged_rcs1213level <- left_join(unmatched_orglevel, national_rcs1213,
                                 by = c("admimeth_C" = "admimeth_C",
                                        "SUSHRGep" = "CURRENCY.CODE"))

# Not all HES episodes are matched with a cost, so we next move onto the national ref costs 1112
# subset MATCHED and UNMATCHED episodes
matched_rcs1213level <- merged_rcs1213level[!is.na(merged_rcs1213level$unitcost),]
unmatched_rcs1213level <- merged_rcs1213level[is.na(merged_rcs1213level$unitcost),]

# remove UNIT.COST column from unmatched_rcs1213level
unmatched_rcs1213level <- select(unmatched_rcs1213level, -unitcost)

# match unmatched_rcs1213level with national_rcs1112
merged_rcs1112level <- left_join(unmatched_rcs1213level, national_rcs1112,
                                 by = c("admimeth_C" = "admimeth_C",
                                        "SUSHRGep" = "CURRENCY.CODE"))

# Not all HES episodes are matched with a cost, so we next move onto the national ref costs 1011
# subset MATCHED and UNMATCHED episodes
matched_rcs1112level <- merged_rcs1112level[!is.na(merged_rcs1112level$unitcost),]
unmatched_rcs1112level <- merged_rcs1112level[is.na(merged_rcs1112level$unitcost),]

# remove UNIT.COST column
unmatched_rcs1112level <- select(unmatched_rcs1112level, -unitcost)

# match unmatched_rcs1112level with national_rcs1011
merged_rcs1011level <- left_join(unmatched_rcs1112level, national_rcs1011,
                                 by = c("admimeth_C" = "admimeth_C",
                                        "SUSHRGep" = "CURRENCY.CODE"))

# Not all HES episodes are matched with a cost, so we next move onto the national ref costs 1011
# subset MATCHED and UNMATCHED episodes
matched_rcs1011level <- merged_rcs1011level[!is.na(merged_rcs1011level$unitcost),]
unmatched_rcs1011level <- merged_rcs1011level[is.na(merged_rcs1011level$unitcost),]

# remove UNIT.COST column
unmatched_rcs1011level <- select(unmatched_rcs1011level, -unitcost)

# match unmatched_rcs1011level with national_rcs0910
merged_rcs0910level <- left_join(unmatched_rcs1011level, national_rcs0910,
                                 by = c("admimeth_C" = "admimeth_C",
                                        "SUSHRGep" = "CURRENCY.CODE"))

# append all matched dataframes to form complete hes-cost merge
full_mergedcosts <- rbind(matched_orglevel,
                          matched_rcs1213level,
                          matched_rcs1112level,
                          matched_rcs1011level,
                          merged_rcs0910level)

# figure out how many unmatched episodes we still have after this
print('Number of uncosted HES episodes')
unmatched <- sum(is.na(full_mergedcosts$unitcost))

print('Number of HES episodes')
total <- nrow(full_mergedcosts)

print('Percentage of HES episodes unmatched')
(unmatched/total)*100

# drop unmatched from full_mergedcosts$unitcost
full_mergedcosts <- full_mergedcosts[!is.na(full_mergedcosts$unitcost),]

# Now we move onto calculating unit costs
# first, calculate average unit cost per patient group (mean, sd, min, p25, p50, p75, max)
quantilecost <- full_mergedcosts %>% group_by(ICD, agegrp_v3, admimeth_C) %>%
  do(data.frame(t(quantile(.$unitcost, props = c(0.25, 0.50, 0.75)))))

avgcost <- do.call(data.frame, aggregate(full_mergedcosts$unitcost, 
                                         list(full_mergedcosts$ICD, full_mergedcosts$agegrp_v3, full_mergedcosts$admimeth_C),
                                         function(x) c(mean = mean(x), sd = sd(x))))

colnames(avgcost) <- c("ICD","agegrp","admimeth","avg_cost","sd_cost")

costs <- left_join(avgcost, quantilecost,
                   by = c("ICD" = "ICD",
                          "agegrp" = "agegrp_v3",
                          "admimeth" = "admimeth_C"))
costs <- rename(costs, min = X0., p25 = X25., median = X50., p75 = X75., max = X100.)

costs$disease <- "ICD"

costs$age <- paste("_AGE", costs$agegrp, sep = "")

costs$p <- paste0(costs$disease, costs$ICD,costs$age)
costs$s <- "G"
costs$a <- ifelse(costs$admimeth == 1, "N", "E")

costs_export <- select(costs, -c(ICD, agegrp, admimeth, disease, age))

# Write out costs
write.csv(costs_export,
          file = "D:/Overflows/output/costs.csv",
          row.names = FALSE,
          quote = FALSE)

###########################################################################################
# Now calculate average unit cost per patient group across waiting time quantiles for electives (eg. for p25 WT, what is the cost p25, p50, p75, mean, sd)
electives_merged <- subset(full_mergedcosts, admimeth_C == 1)

quantilewt <- electives_merged %>% group_by(ICD, agegrp_v3) %>%
  do(data.frame(t(quantile(.$WaitingTime, props = c(0.25, 0.50, 0.75)))))

quantilewt <- rename(quantilewt, min = X0., p25 = X25., median = X50., p75 = X75., max = X100.)

# at each quantile, what is the avg, sd, min, 25, 50, 75, max cost?
# identify which electives fall within each quantile
electives_wtquant <- left_join(electives_merged, quantilewt,
                               by = c("ICD","agegrp_v3"))

electives_wtquant$p025 <- ifelse(electives_wtquant$WaitingTime < electives_wtquant$p25, 1, 0)
electives_wtquant$p2550 <- ifelse((electives_wtquant$WaitingTime >= electives_wtquant$p25 & 
                                     electives_wtquant$WaitingTime < electives_wtquant$p50), 1, 0)
electives_wtquant$p5075 <- ifelse((electives_wtquant$WaitingTime >= electives_wtquant$p50 & 
                                     electives_wtquant$WaitingTime < electives_wtquant$p75), 1, 0)
electives_wtquant$p75100 <- ifelse((electives_wtquant$WaitingTime >= electives_wtquant$p75), 1, 0)

electives_p025 <- subset(electives_wtquant, p025 == 1)
electives_p2550 <- subset(electives_wtquant, p2550 == 1)
electives_p5075 <- subset(electives_wtquant, p5075 == 1)
electives_p75100 <- subset(electives_wtquant, p75100 == 1)

# quantiles of costs within each WT quantile
quant_costs_wtp025 <- electives_p025 %>% group_by(ICD, agegrp_v3, p025) %>%
  do(data.frame(t(quantile(.$UNIT.COST, props = c(0.25, 0.50, 0.75)))))

quant_costs_wtp2550 <- electives_p2550 %>% group_by(ICD, agegrp_v3, p2550) %>%
  do(data.frame(t(quantile(.$UNIT.COST, props = c(0.25, 0.50, 0.75)))))

quant_costs_wtp5075 <- electives_p5075 %>% group_by(ICD, agegrp_v3, p5075) %>%
  do(data.frame(t(quantile(.$UNIT.COST, props = c(0.25, 0.50, 0.75)))))

quant_costs_wtp75100 <- electives_p75100 %>% group_by(ICD, agegrp_v3, p75100) %>%
  do(data.frame(t(quantile(.$UNIT.COST, props = c(0.25, 0.50, 0.75)))))

# mean and sd cost within each WT quantile
avg_costs_wtp025 <- do.call(data.frame, aggregate(electives_p025$UNIT.COST, 
                                                  list(electives_p025$ICD, electives_p025$agegrp_v3),
                                                  function(x) c(mean = mean(x), sd = sd(x))))
colnames(avg_costs_wtp025) <- c("ICD","agegrp_v3","avg_cost","sd_cost")

avg_costs_wtp2550 <- do.call(data.frame, aggregate(electives_p2550$UNIT.COST, 
                                                   list(electives_p2550$ICD, electives_p2550$agegrp_v3),
                                                   function(x) c(mean = mean(x), sd = sd(x))))
colnames(avg_costs_wtp2550) <- c("ICD","agegrp_v3","avg_cost","sd_cost")

avg_costs_wtp5075 <- do.call(data.frame, aggregate(electives_p5075$UNIT.COST, 
                                                   list(electives_p5075$ICD, electives_p5075$agegrp_v3),
                                                   function(x) c(mean = mean(x), sd = sd(x))))
colnames(avg_costs_wtp5075) <- c("ICD","agegrp_v3","avg_cost","sd_cost")

avg_costs_wtp75100 <- do.call(data.frame, aggregate(electives_p75100$UNIT.COST, 
                                                    list(electives_p75100$ICD, electives_p75100$agegrp_v3),
                                                    function(x) c(mean = mean(x), sd = sd(x))))
colnames(avg_costs_wtp75100) <- c("ICD","agegrp_v3","avg_cost","sd_cost")

# combine mean, sd, and quantile costs at each WT quantile
costs_wtp025 <- left_join(avg_costs_wtp025, quant_costs_wtp025,
                          by = c("ICD" = "ICD","agegrp_v3" = "agegrp_v3"))
costs_wtp025 <- rename(costs_wtp025, min = X0., p25 = X25., median = X50., p75 = X75., max = X100.)
costs_wtp025 <- subset(costs_wtp025, select = -p025)

costs_wtp2550 <- left_join(avg_costs_wtp2550, quant_costs_wtp2550,
                           by = c("ICD" = "ICD","agegrp_v3" = "agegrp_v3"))
costs_wtp2550 <- rename(costs_wtp2550, min = X0., p25 = X25., median = X50., p75 = X75., max = X100.)
costs_wtp2550 <- subset(costs_wtp2550, select = -p2550)

costs_wtp5075 <- left_join(avg_costs_wtp5075, quant_costs_wtp5075,
                           by = c("ICD" = "ICD", "agegrp_v3" = "agegrp_v3"))
costs_wtp5075 <- rename(costs_wtp5075, min = X0., p25 = X25., median = X50., p75 = X75., max = X100.)
costs_wtp5075 <- subset(costs_wtp5075, select = -p5075)

costs_wtp75100 <- left_join(avg_costs_wtp75100, quant_costs_wtp75100,
                            by = c("ICD" = "ICD", "agegrp_v3" = "agegrp_v3"))
costs_wtp75100 <- rename(costs_wtp75100, min = X0., p25 = X25., median = X50., p75 = X75., max = X100.)
costs_wtp75100 <- subset(costs_wtp75100, select = -p75100)

# write out to csv
write.csv(costs_wtp025,
          file = "D:/Overflows/output/costs_p025.csv", 
          row.names = FALSE,
          quote = FALSE)

write.csv(costs_wtp2550,
          file = "D:/Overflows/output/costs_p2550.csv", 
          row.names = FALSE,
          quote = FALSE)

write.csv(costs_wtp5075,
          file = "D:/Overflows/output/costs_p5075.csv", 
          row.names = FALSE,
          quote = FALSE)

write.csv(costs_wtp75100,
          file = "D:/Overflows/output/costs_p75100.csv", 
          row.names = FALSE,
          quote = FALSE)
