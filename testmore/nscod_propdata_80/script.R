library(stockassessment)
cn<-read.ices("cn.dat")
cw<-read.ices("cw.dat")
dw<-read.ices("dw.dat")
lf<-read.ices("lf.dat")
lw<-read.ices("lw.dat")
mo<-read.ices("mo.dat")
nm<-read.ices("nm.dat")
pf<-read.ices("pf.dat")
pm<-read.ices("pm.dat")
sw<-read.ices("sw.dat")
surveys<-read.ices("survey.dat")

seasonTimes <- c(0,0.125,1)
p <- c(0.2,0.8)
## p <- rep(0.5,2)
set.seed(345)
v <- cbind(replicate(length(p)-1,as.numeric(arima.sim(list(ar=0.8),n=nrow(cn),sd=0.1))),0)
P <- plogis(qlogis(p)+t(v))


rdirichlet <- function(n,p,s){
    replicate(n,{xx <- rgamma(length(p),shape = p * s, rate = 1); xx/sum(xx)})
}

propC <- simplify2array(apply(P,2,function(pp){
    v <- rdirichlet(1,pp,1000)
    matrix(log(v[,1] / tail(v[,1],1)),ncol=1)
}, simplify=FALSE))[seq_len(length(p)-1),,,drop=FALSE]

aux <- as.matrix(data.frame(fleetType = 80,
                  year = as.numeric(rownames(cn))[slice.index(propC,3)],
                  age = -1, ##as.numeric(colnames(cn))[slice.index(propC,2)],
                  logobs = (as.vector(propC)),
                  catchFleet = 1,
                  seasonStart = seasonTimes[slice.index(propC,1)],
                  seasonEnd = seasonTimes[slice.index(propC,1)+1],
                  catchType = 1,
                  seasonNumber = as.vector(slice.index(propC,1)),
                  numberOfSeasons = dim(propC)[1]+1))

dat<-setup.sam.data(surveys=surveys,
                    residual.fleets=cn,
                    aux.fleet = aux,
                    prop.mature=mo, 
                    stock.mean.weight=sw, 
                    catch.mean.weight=cw, 
                    dis.mean.weight=dw, 
                    land.mean.weight=lw,
                    prop.f=pf, 
                    prop.m=pm, 
                    natural.mortality=nm, 
                    land.frac=lf)

conf<-defcon(dat)
## conf$keyLogFsta[1,] <- c(0, 1, 2, 3, 4, 5)
## conf$corFlag <- c(2)
## conf$keyVarF[1,] <- c(0, 1, 1, 1, 1, 1)
## conf$keyLogFpar <- matrix(
##     c(-1, -1, -1, -1, -1, -1,      
##       0,  1,  2,  3,  4, -1,
##       5,  6,  7,  8, -1, -1
##       -1, -1, -1, -1, -1, -1
##       ), nrow=4, byrow=TRUE)
## conf$keyVarObs <- matrix(
##     c( 0,  1,  2,  2,  2,  2,
##       3,  4,  4,  4,  4, -1,
##       5,  6,  6,  6, -1, -1,
##       7, 7, 7, 7, 7, 7
##       ), nrow=4, byrow=TRUE)
## conf$noScaledYears <- 13
## conf$keyScaledYears <- 1993:2005
## conf$keyParScaledYA <- row(matrix(NA, nrow=13, ncol=6))-1
## conf$fbarRange <- c(2,4)
#conf$keyLogFpar[3,] <- c(0,1,2,3,-1,-1)
conf$seasonTimes <- seasonTimes
conf$isFishingSeason <- rep(1,length(seasonTimes)-1)
conf$keyLogFseason[1,] <- 0
conf$obsLikelihoodFlag[4] <- "Dirichlet"

par<-defpar(dat,conf)
par$logSdLogObs[4] <- 5
fit<-sam.fit(dat,conf,par, newtonsteps=0, rm.unidentified=FALSE, run=TRUE)

resid <- residuals(fit)

cat(fit$opt$objective,"\n", file="res.out")
