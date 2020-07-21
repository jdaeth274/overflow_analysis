require(vroom)
require(tidyverse)
require(lubridate)
require(ggplot2)

data <- vroom("D:/Overflows/data/HES_APC_CC_0913_transitions_all_ICD.csv", delim = ",",
                  col_names = TRUE, num_threads = 8)

icd50 <- subset(data, subset = ICD == 50)
icd50$ICD01 <- ifelse(icd50$MainICD10Cat == 1, 1, 0)
icd50$ICD04 <- ifelse(icd50$MainICD10Cat == 4, 1, 0)
icd50$ICD05 <- ifelse(icd50$MainICD10Cat == 5, 1, 0)
icd50$ICD15 <- ifelse(icd50$MainICD10Cat == 15, 1, 0)
icd50$ICDBundle <- ifelse(icd50$MainICD10Cat == 8 | icd50$MainICD10Cat == 16 | icd50$MainICD10Cat == 17, 1, 0)

icd50$date <- as.POSIXct( paste(1, icd50$admidate_week, icd50$admidate_YYYY, sep = "-"), format="%u-%U-%Y", locale = "UK")
icd50 <- icd50[order(icd50$date),]

icd50 <- icd50[icd50$date >= "2009-05-01",]

bundleprop <- data.frame(matrix(vector(), 204, 1))
bundleprop <- bundleprop %>% mutate (date = unique(icd50$date))
bundleprop <- icd50 %>%
  group_by(date) %>%
  summarize(ICD01 = sum(ICD01),
            ICD04 = sum(ICD04),
            ICD05 = sum(ICD05),
            ICD15 = sum(ICD15),
            ICDBundle = sum(ICDBundle),
            totalbundle = sum(One))

bundleprop <- bundleprop %>%
  mutate(ICD01prop = ICD01/totalbundle,
         ICD04prop = ICD04/totalbundle,
         ICD05prop = ICD05/totalbundle,
         ICD15prop = ICD15/totalbundle,
         ICDBundleprop = ICDBundle/totalbundle)

setwd("D:/Overflows/output/")
pdf(file = "proportions_ICD_electives_bundle.pdf",paper='A4')

print(ggplot(bundleprop, aes(x = date)) +
        geom_line(aes(x = date, y = ICD01prop, color = 'red')) + 
        geom_line(aes(x = date, y = ICD04prop, color = 'blue')) + 
        geom_line(aes(x = date, y = ICD05prop, color = 'green')) +
        geom_line(aes(x = date, y = ICD15prop, color = 'purple')) +
        geom_line(aes(x = date, y = ICDBundleprop, color = 'orange')) + 
        scale_color_discrete(name = "Legend", labels = c("Proportion_ICD01", "Proportion_ICD04", "Proportion_ICD05",
                                                         "Proportion_ICD15", "Proportion_ICDBundle"))+
        ggtitle("Proportion of each ICD in Electives Bundle"))

dev.off()
