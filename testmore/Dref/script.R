library(stockassessment)

data(nscodData)
data(nscodConf)
data(nscodParameters)

fit <- sam.fit(nscodData, nscodConf, defpar(nscodData,nscodConf))

cat(class(try({deterministicReferencepoints(fit,c("MSY"))},silent=TRUE)),"\n", file = "res.out")


nscodConf$stockRecruitmentModelCode <- 69
fitA <- sam.fit(nscodData, nscodConf, defpar(nscodData,nscodConf))

cat(logLik(fitA),"\n", file = "res.out")

rpA <- deterministicReferencepoints(fitA,c("MSY"))

cat(logLik(fitA),"\n", file = "res.out")


nscodConf$
nscodConf$stockRecruitmentModelCode <- 69
fitB <- sam.fit(nscodData, nscodConf, defpar(nscodData,nscodConf))

system.time(rp2 <- deterministicReferencepoints(fitB,c("MSY")))


data(nscodData)
data(nscodConf)
data(nscodParameters)

nscodConf
nscodConf$stockRecruitmentModelCode <- 69
fitA <- sam.fit(nscodData, nscodConf, defpar(nscodData,nscodConf))

system.time(rp3 <- deterministicReferencepoints(fitC,c("MSY")))
