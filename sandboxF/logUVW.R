load("Fobs.RData")
matplot(Fobs$year, log(Fobs$Fobs), xlab="Year", ylab="logF", pch=colnames(Fobs$Fobs))

library(TMB)
compile("logUVW.cpp")
dyn.load(dynlib("logUVW"))

Fobs$cormode <- 0

par <- list()
par$logsdV <- -1
par$logsd <- -1
par$logsdObs <- -2

par$betaU = 0
par$betaV = 0
par$alphaV = 0
par$logU <- matrix(0,  nrow=ncol(Fobs$Fobs)-1 , ncol=nrow(Fobs$Fobs))
par$alphaU = rep(0,dim(par$logU)[1] )
par$logsdU <- rep(-1,dim(par$logU)[1])
par$logV <- rep(0,nrow(Fobs$Fobs))
par$logFUpper <- matrix(0,  nrow=nrow(Fobs$Fobs), ncol=ncol(Fobs$Fobs))


map = list(logsdU=factor(rep(0,length(par$logsdU))), logsdObs=factor(NA))

obj <- MakeADFun(Fobs, par, random=c("logU","logV","logFUpper"), DLL="logUVW", 
                 map=map)
fit <- nlminb(obj$par, obj$fn, obj$gr)
sdr<-sdreport(obj)

Fest = as.list(sdr,"Est")$logFUpper

matplot(Fobs$year, (Fest), type="l", add=TRUE)
2/(1 + exp(-2*fit$par["betaU"] )) - 1
2/(1 + exp(-2*fit$par["betaV"] )) - 1

## res
## ID: nll=441.9986
## CS: nll=278.1695
## AR: nll=230.7495
## parallel: nll ~ 485
