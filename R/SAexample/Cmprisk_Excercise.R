# Competing risk models excercise

library(cmprsk)
library(msSurv)
setwd("/Users/Monkeyface/Dropbox/00_COVID/Rcode")

Data <- read.csv("1471-2288-11-144-S4.csv", header=TRUE)

CI.overall <- cuminc(ftime=Data$LOS, fstatus=Data$CR.event)

plot(CI.overall, curvlab=c("Discharge", "Death"), xlab="Days")

CI.4vs5 <- cuminc(ftime=Data$LOS, fstatus=Data$CR.event, group=Data$RiskClass)

plot(CI.4vs5, lty=c(1,1,2,2),
     col=c("black", "blue", "black", "blue"),
     curvlab=c("Discharge, younger", "Discharge, older",
               "Death, younger", "Death, older"), xlab="Days")
CI.4vs5$Tests