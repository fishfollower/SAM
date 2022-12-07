## Tests the combination of ALN likelihood with AR correlation structure
library(stockassessment)
data(nscodData)
data(nscodConf)

nscodConf$keyLogFmu <- nscodConf$keyLogFsta
nscodConf$keyLogFrho[1,] <- 0
nscodConf$obsLikelihoodFlag[3]<-"ALN"   
nscodConf$obsCorStruct[3]<-"AR"
nscodConf$keyCorObs[3,]<-c(0,1,1,-1,-1)

par<-defpar(nscodData,nscodConf)

fit<-sam.fit(nscodData,nscodConf,par)

cat(fit$opt$objective,"\n", file="res.out")

## Will warn about biopar
suppressWarnings(fv <- modelforecast(fit,rep(NA,100), nosim=100))
