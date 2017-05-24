library(stockassessment)
load("run/model.RData")
LO<-leaveout(fit)
save(LO, file="run/leaveout.RData")
