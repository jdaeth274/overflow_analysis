
# Create main categories of ICD10 codes
# Requires a dataframe, a column with ICD 10 codes.
# Returns a modified dataframe.
categorise_ICD10 <- function(data, code, new_col){
  
  data <- data %>% 
    mutate(clean_code := toupper(!!rlang::sym(code)),
           clean_code = gsub("[[:punct:]]", "", clean_code),
           clean_code = gsub(" ", "", clean_code))
  
  detect_diagnoses <- function(x, codes){
    return(str_detect(replace_na(x, ""), str_c(codes, collapse = "|")))
  }
  
  data <- data %>%
    mutate(!!new_col := case_when(detect_diagnoses(clean_code, c("A.*", "B.*")) ~ 1,
                                  detect_diagnoses(clean_code, c("C.*", "D[1|2|3|4].*")) ~ 2,
                                  detect_diagnoses(clean_code, c("D[5|6|7|8].*"))~ 3,
                                  detect_diagnoses(clean_code, c("E.*")) ~ 4,
                                  detect_diagnoses(clean_code, c("F.*")) ~ 5,
                                  detect_diagnoses(clean_code, c("G.*")) ~ 6,
                                  detect_diagnoses(clean_code, c("H[1|2|3|4|5].*")) ~ 7, 
                                  detect_diagnoses(clean_code, c("H[6|7|8|9].*")) ~ 8,
                                  detect_diagnoses(clean_code, c("I.*")) ~ 9,
                                  detect_diagnoses(clean_code, c("J.*")) ~ 10,
                                  detect_diagnoses(clean_code, c("K.*")) ~ 11,
                                  detect_diagnoses(clean_code, c("L.*")) ~ 12,
                                  detect_diagnoses(clean_code, c("M.*")) ~ 13,
                                  detect_diagnoses(clean_code, c("N.*")) ~ 14,
                                  detect_diagnoses(clean_code, c("O.*")) ~ 15,
                                  detect_diagnoses(clean_code, c("P.*")) ~ 16,
                                  detect_diagnoses(clean_code, c("Q.*")) ~ 17,
                                  detect_diagnoses(clean_code, c("R.*")) ~ 18,
                                  detect_diagnoses(clean_code, c("S.*", "T.*")) ~ 19,
                                  detect_diagnoses(clean_code, c("V.*", "W.*", "X.*", "Y.*")) ~ 20,
                                  detect_diagnoses(clean_code, c("Z.*")) ~ 21,
                                  detect_diagnoses(clean_code, c("U.*")) ~ 22)) %>% 
    select(-clean_code) 
}


# Combined ICD10 groups based on main categories of ICD10 groups and admission method
combined_ICD10_admimeth <- function(data, ICD10_grp, admimeth, new_col){
  
  data <- data %>%
    mutate(!!new_col := case_when(!!rlang::sym(ICD10_grp) %in% c(1, 4, 5, 8, 15, 16, 17) & !!rlang::sym(admimeth) == 1 ~ 50,
                                  !!rlang::sym(ICD10_grp) %in% c(3, 7, 8, 16, 17, 21) & !!rlang::sym(admimeth) == 2 ~ 51,
                                  TRUE ~ !!rlang::sym(ICD10_grp)))
}


# Import one year of APC data from HES database
# Requires an open SQLite database connection and a string for the year to be extracted
# Returns a dataframe 
import_APC <- function(db, year_char){
  
  APC_db <- tbl(db, 'APC')
  
  # General and acute episodes
  GAEs <- APC_db %>% 
    filter(substr(EPISTART, 1, 4) == year_char &
             substr(PROCODE3, 1, 1) =='R' & 
             ADMIMETH %in% c(11, 12, 13, 21:28, 66, 67, 69)) %>% 
    select(hesid = ENCRYPTED_HESID, 
           diag_01 = DIAG_01,
           rttstart = RTTPERSTART,
           admidate_MDY = ADMIDATE_FILLED,
           admimeth = ADMIMETH,
           dismeth = DISMETH,
           disdate_MDY = DISDATE,
           Procode = PROCODE3,
           startage = STARTAGE_CALC,
           susrecid = SUSRECID,
           epistart_MDY = EPISTART, 
           epiend = EPIEND, 
           SUSHRG, SUSCOREHRG, SUSHRGVERS,
           IMPFRAILTY_DELIRIUM, IMPFRAILTY_DEMENTIA, IMPFRAILTY_DEPENDENCE, 
           IMPFRAILTY_FALLSFRAX, IMPFRAILTY_MOBILITY, IMPFRAILTY_ULCERS, IMPFRAILTY_SENILITY) %>% 
    collect()
  
  return(GAEs)
}


# Clean APC data and derive variables
# Requires a dataframe produced by the import_APC function
# Returns a dataframe
process_APC <- function(data, year_char){
  
  mhtrusts_codes <- c("R1A", "R1C", "R1E", "RAT", "RDY", "RGD", "RH5", "RJ8", "RKL", "RLY", "RMY", "RNK",
                      "RNN", "RNU", "RP1", "RP7", "RPG", "RQY", "RRE", "RRP", "RT1", "RT2", "RT5", "RT6",
                      "RTQ", "RTV", "RV3", "RV5", "RV9", "RVN", "RW1", "RW4", "RW5", "RWK", "RWN", "RWQ",
                      "RWV", "RWX", "RX2", "RX3", "RX4", "RXA", "RXG", "RXM", "RXT", "RXV", "RXX", "RXY",
                      "RYG", "RYK")
  
  data <- data %>% 
    filter(!is.na(disdate_MDY) &
             !is.element(dismeth, c(5, 8)) &
             !is.na(startage) &
             !is.element(Procode, mhtrusts_codes) &
             !is.na(hesid))
  
  
  data <- data %>% 
    mutate_at(vars(admidate_MDY, disdate_MDY, epistart_MDY, epiend, rttstart), ~as.Date(., "%Y-%m-%d")) 
  
  if(year_char == "2020"){
    data <- data %>% 
      filter(month(epistart_MDY) %in% c(1,2))
  }
  
  # the YYYY year variables and week variables refer to year and week assigned according
  # to ISO8601
  
  data <- data %>% 
    mutate(admimeth_C = case_when(admimeth %in% c(11, 12, 13) ~ 1, # elective
                                  admimeth %in% c(21:28, 66, 67, 69) ~ 2), # emergency
           WaitingTime = as.numeric(admidate_MDY - rttstart),
           agegrp_v3 = case_when(startage < 25 ~ 1,
                                 startage >= 25 & startage < 65 ~ 2,
                                 startage >= 65 ~ 3),
           Frail = ifelse(IMPFRAILTY_DELIRIUM == 1 | IMPFRAILTY_DEMENTIA == 1 |
                              IMPFRAILTY_DEPENDENCE == 1 | IMPFRAILTY_FALLSFRAX == 1 | 
                              IMPFRAILTY_MOBILITY == 1 | IMPFRAILTY_ULCERS == 1 | 
                              IMPFRAILTY_SENILITY == 1, 1, 0),
           admidate_YYYY = format(admidate_MDY, format = "%G"), 
           admidate_week = format(admidate_MDY, format = "%V"), 
           epistart_YYYY = format(epistart_MDY, format = "%G"), 
           epistart_week = format(epistart_MDY, format = "%V"), 
           rttstart_week = format(rttstart, format = "%V"), 
           rttstart_YYYY = format(rttstart, format = "%G"))  
  
  data <- data %>%
    categorise_ICD10("diag_01", new_col = "MainICD10Cat") %>% 
    combined_ICD10_admimeth(ICD10_grp = "MainICD10Cat", admimeth = "admimeth_C", new_col = "ICD") %>% 
    filter(!is.na(MainICD10Cat) & 
             !(MainICD10Cat == 15 & agegrp_v3 == 3) & 
             !(MainICD10Cat == 16 & agegrp_v3 > 1)  &
             epiend >= epistart_MDY &
             (is.na(WaitingTime) | WaitingTime >= 0))
  
  return(data)
  
}

# Import one year of CC data from HES databast
# Requires an open SQLite database connection and a string for the year to be extracted
# Returns a dataframe 
import_CC <- function(db, year_char){
  
  CC_db <- tbl(db, 'CC')
  
  data <- CC_db %>% 
    filter(substr(CCSTARTDATE, 1, 4) == year_char) %>% 
    select(ccstartdate = CCSTARTDATE,
           ccdisdate = CCDISDATE,
           ccdisdest = CCDISDEST,
           susrecid = SUSRECID, 
           ccunitfun = CCUNITFUN) %>% 
    collect()
  
  return(data)
}

# Clean CC data and derive variables
# Requires a dataframe produced by the import_CC function
# and the APC data that it will be added to later on
# Returns a dataframe
process_CC <- function(data, APC_data){
  
  data <- data %>% 
    mutate_at(vars(ccstartdate, ccdisdate), ~as.Date(., "%Y%m%d")) %>% 
    filter(ccdisdate >= ccstartdate & ccunitfun %in% c(1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 90, 91)) %>% 
    mutate(cc_LoS = as.numeric(ccdisdate - ccstartdate)) 
  
  # Link in APC epistart and epiend 
  
  data <- data %>% 
    left_join(APC_data[, c("susrecid", "epistart_MDY", "epiend")], by = "susrecid")
  
  # First figure out how the CC episode overlaps with the assigned APC episode
  # Exclude if they don't overlap at all
  data <- data %>% 
    filter(!is.na(epistart_MDY) & !is.na(epiend)) %>% 
    mutate(match_level = case_when(ccstartdate == epistart_MDY & ccdisdate == epiend ~ 1,
                                   ccstartdate > epistart_MDY & ccdisdate < epiend ~ 2,
                                   ccstartdate > epistart_MDY & ccdisdate == epiend ~ 3,
                                   ccstartdate == epistart_MDY & ccdisdate < epiend ~ 4,
                                   TRUE ~ NA_real_)) %>% 
    filter(!is.na(match_level))
  
  # Then figure out whether CC episodes overlap within each susrecid grouping 
  # or wether on episode is within another episode
  data <- data %>% 
    group_by(susrecid) %>% 
    arrange(susrecid, ccstartdate, desc(cc_LoS)) %>% 
    mutate(to_remove = ifelse(ccstartdate >= lag(ccstartdate) & ccdisdate <= lag(ccdisdate), 1, 0),
           overlap_with_prev = ifelse(ccstartdate < lag(ccdisdate), as.numeric(lag(ccdisdate)-ccstartdate), 0))
  
  # If there is an overlap, filter out the shorter episode
  data <- data %>% 
    ungroup() %>% 
    filter(is.na(to_remove) | to_remove == 0)
  
  # Summarise
  data <- data %>% 
    group_by(susrecid) %>% 
    arrange(ccstartdate, ccdisdate) %>% 
    summarise(ccstartdate = ccstartdate[row_number() == 1],
              ccdisdate = ccdisdate[row_number() == n()],
              ccdisdest = ccdisdest[row_number() == n()],
              cc_LoS = sum(cc_LoS) - sum(overlap_with_prev, na.rm = TRUE))
  
  return(data)
}

# Join APC and CC data and derive LoS Variables
# Requires a dataframe with APC data produced by the process_APC function
# and a dataframe with CC data produced by the process_CC function
# Returns a dataframe
join_APC_CC <- function(APC_data, CC_data){
  
  APC_data <- APC_data %>% 
    left_join(CC_data, by = c("susrecid"))
  

  APC_data <- APC_data %>% 
    mutate(cc_start_flag = ifelse(ccstartdate == epistart_MDY, 1, 0),
           cc_dis_flag = ifelse(ccdisdate == epiend, 1, 0),
           Total_LoS = as.numeric(epiend - epistart_MDY),
           GA_LOS = Total_LoS - cc_LoS,
           cc = ifelse(!is.na(ccstartdate), 1, 0)) 
  
  return(APC_data)
  
}


# Imports, processes and combined one year of HES APC and CC data
# Requires an open SQLite database connection, a string for the year to be extracted
# and a folder where the resulting Rds file will be saved
prepare_HES_year <- function(db, year_char, output_folder){
  
  APC <- import_APC(db, year_char) 
  APC <- process_APC(APC, year_char)

  CC <- import_CC(db, year_char)
  CC <- process_CC(data = CC, APC_data = APC)
  
  APC_CC <- join_APC_CC(APC_data = APC, CC_data = CC)
  
  saveRDS(APC_CC, here::here(output_folder, str_c("APC_CC_", year_char, ".Rds")))

}

  