library(TMB)
compile("babysamUVW.cpp")
dyn.load(dynlib("babysamUVW"))

load("dat.RData")

dat$srmode <- 0
dat$fcormode <- 0

dat$keyF <- c(0,1,2,3,4,5,6,7,7)

dat$keyQ <- rbind(c(NA,NA,NA,NA,NA,NA,NA,NA,NA),
                  c(NA, 0, 1, 2, 3, 4, 5, 5,NA),
                  c( 6, 7, 8, 9,10,10,NA,NA,NA))

dat$keySd <- rbind(c( 0, 0, 0, 0, 0, 0, 0, 0, 0),
                   c(NA, 1, 1, 1, 1, 1, 1, 1,NA),
                   c( 2, 2, 2, 2, 2, 2,NA,NA,NA))

dat$covType <- c(0,1,1) # 0=ID, 1=IGAR, 2=US ...

#  0-1 1-2 2-3 3-4 4-5 5-6 6-7 7-8 
dat$keyIGAR <- rbind(c( -1, -1, -1, -1, -1, -1, -1, -1),
                     c( NA, 0,  1,  2,  3,  4,  5, NA),
                     c( 6,  7,  8,  9,  9, NA, NA, NA))

dat$noParUS <- sapply(1:length(dat$fleetTypes),
                      function(f){
                        A<-sum(!is.na(dat$keySd[f,]))
                        ifelse(dat$covType[f]==2, (A*A-A)/2, 0)
                      })

dat$ncolF = max(dat$keyF)+1
dat$nrowF = length(dat$year)

par <- list()
par$logsdR <-0
par$logsdS <- 0
par$rickerpar <- if(dat$srmode==1){c(1,1)}else{numeric(0)}
par$transRhoF <- if(dat$fcormode==0){numeric(0)}else{0.1}
par$bhpar <- if(dat$srmode==2){c(1,1)}else{numeric(0)}
par$logQ <- numeric(max(dat$keyQ, na.rm=TRUE)+1)
par$logsd <- numeric(max(dat$keySd, na.rm=TRUE)+1)
par$logIGARdist <- numeric(max(dat$keyIGAR, na.rm=TRUE)+1)
par$parUS <- numeric(sum(dat$noParUS))
par$logN <- matrix(0, nrow=length(dat$year), ncol=length(dat$age))
par$missing <- numeric(sum(is.na(dat$obs)))

par$betaU = -1
par$betaV = -1
par$alphaV = -1
par$alphaU = rep(-1,dim(par$logU)[2] )
par$logsdU <- rep(-1,dim(par$logU)[2])
par$logsdV <- -1
par$logsdW <- -1
par$logU <- matrix(0,  nrow=dat$nrowF , ncol=dat$ncolF-1)
par$logV <- rep(0,dat$nrowF)
par$logW <- matrix(0,  nrow=dat$nrowF, ncol=dat$ncolF)


obj <- MakeADFun(dat, par, random=c("logN", "logU","logV","logW", "missing"), DLL="babysamUVW", map= list(logsdU=factor(rep(0,length(par$logsdU)))))

fit <- nlminb(obj$par, obj$fn, obj$gr, control=list(eval.max=1000, iter.max=1000))
sdr<-sdreport(obj)
