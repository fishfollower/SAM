library(stockassessment)
setwd("run")
load("data.RData")

conf <- defcon(dat)
saveConf(conf,"../conf/model.cfg")
conf<-loadConf(dat,"../conf/model.cfg")
par<-defpar(dat,conf)
fit<-sam.fit(dat,conf,par)

save(fit, file="model.RData")
