## Tests the combination of ALN likelihood with ID correlation structure
library(stockassessment)
data(nscodData)
data(nscodConf)

nscodConf$obsLikelihoodFlag[3]<-"ALN"
nscodConf$obsCorStruct[3]<-"ID"

par<-defpar(nscodData,nscodConf)

fit<-sam.fit(nscodData,nscodConf,par)

cat(fit$opt$objective,"\n", file="res.out")
