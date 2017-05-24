library(stockassessment)
load("run/model.RData")
RETRO<-retro(fit, year=5)
save(RETRO, file="run/retro.RData")
