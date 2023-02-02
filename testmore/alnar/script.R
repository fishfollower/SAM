## Tests the combination of ALN likelihood with AR correlation structure
library(stockassessment)
data(nscodData)
data(nscodConf)

nscodConf$obsLikelihoodFlag[3]<-"ALN"
nscodConf$obsCorStruct[3]<-"AR"
nscodConf$keyCorObs[3,]<-c(0,1,1,-1,-1)

par<-defpar(nscodData,nscodConf)

fit<-sam.fit(nscodData,nscodConf,par)

d <- deterministicReferencepoints(fit,"Max")

## Test that simulation does not give error
s <- simulate(fit,1)


cat(fit$opt$objective,"\n", file="res.out")
