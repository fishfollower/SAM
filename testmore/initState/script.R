## Tests the combination of ALN likelihood with ID correlation structure
library(stockassessment)
data(nscodData)
data(nscodConf)

nscodConf$obsLikelihoodFlag[3]<-"ALN"
nscodConf$obsCorStruct[3]<-"ID"

nscodConf$initState <- 0
par<-defpar(nscodData,nscodConf)
fit<-sam.fit(nscodData,nscodConf,par)


nscodConf$initState <- 1
par<-defpar(nscodData,nscodConf)
fit2<-sam.fit(nscodData,nscodConf,par)

AIC(fit)
AIC(fit2)

cat(logLik(fit),"\n", file = "res.out")
cat(logLik(fit2),"\n", file = "res.out", append=TRUE)

cat(AIC(fit),"\n", file = "res.out", append=TRUE)
cat(AIC(fit2),"\n", file = "res.out", append=TRUE)
