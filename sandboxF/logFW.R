load("Fobs.RData")
matplot(Fobs$year, log(Fobs$Fobs), xlab="Year", ylab="logF", pch=colnames(Fobs$Fobs))

library(TMB)
compile("logFW.cpp")
dyn.load(dynlib("logFW"))

Fobs$cormode <- 0
 
par <- list()
par$logsdV <- -1
par$logsd <- -2

par$betaU = 0
par$betaV = 0
par$alphaV = 0
par$alphaU = rep(0,ncol(Fobs$Fobs)-1)
par$logsdU <- rep(-1,ncol(Fobs$Fobs)-1)
par$logF <- matrix(0,  nrow=nrow(Fobs$Fobs), ncol=ncol(Fobs$Fobs))
par$logsdW <- 0
par$logW <- matrix(0,  nrow=nrow(Fobs$Fobs), ncol=ncol(Fobs$Fobs))

map = list(logsdU=factor(rep(0,length(par$logsdU))), logsd=factor(NA))

obj <- MakeADFun(Fobs, par, random=c("logF", "logW"), DLL="logFW", map=map)
fit <- nlminb(obj$par, obj$fn, obj$gr)
sdr<-sdreport(obj)

Fest = as.list(sdr,"Est")$logF

matplot(Fobs$year, (Fest), type="l", add=TRUE)
2/(1 + exp(-2*fit$par["betaU"] )) - 1
2/(1 + exp(-2*fit$par["betaV"] )) - 1

## res
## ID: nll=441.9986
## CS: nll=278.1695
## AR: nll=230.7495
## parallel: nll ~ 485
## xsam nll 289.8875
