library(stockassessment)

data(nscodData)
data(nscodConf)

conf <- nscodConf
conf$stockRecruitmentModelCode <- 3
par <- defpar(nscodData,conf)
fit1 <- sam.fit(nscodData,conf,par)

conf$constRecBreaks <- c(1997.5)
par <- defpar(nscodData,conf)
fit2 <- sam.fit(nscodData,conf,par)

cat(c(fit1$opt$objective,fit2$opt$objective) ,"\n", file="res.out")
