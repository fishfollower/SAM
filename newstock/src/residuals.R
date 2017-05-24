library(stockassessment)
load("run/model.RData")
RES<-residuals(fit)
RESP<-procres(fit)
save(RES, RESP, file="run/residuals.RData")












