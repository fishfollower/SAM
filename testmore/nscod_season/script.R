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

p <- c(0.8,0.1,0.1)
set.seed(345)
v <- cbind(replicate(length(p)-1,as.numeric(arima.sim(list(ar=0.8),n=nrow(cn),sd=0.1))),0)
cnS <- lapply(seq_along(p), function(pp) cn * rbinom(prod(dim(cn)),100,plogis(qlogis(p[pp]) + v[,pp])) / 100)
for(i in seq_along(cnS))
    attr(cnS[[i]],"time") <- c(1/3*(i-1),1/3*i)


dat<-setup.sam.data(surveys=surveys,
                    residual.fleet=cnS, 
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
conf$keyLogFsta[1,] <- c(0, 1, 2, 3, 4, 5)
conf$keyLogFsta[2,] <- c(0, 1, 2, 3, 4, 5)
conf$keyLogFsta[3,] <- c(0, 1, 2, 3, 4, 5)
conf$corFlag <- c(0,0,2)
conf$keyVarF[1,] <- c(0, 1, 1, 1, 1, 1)
conf$keyVarF[2,] <- c(0, 1, 1, 1, 1, 1)
conf$keyVarF[3,] <- c(0, 1, 1, 1, 1, 1)
conf$keyLogFpar <- matrix(
    c(-1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1,
      0,  1,  2,  3,  4, -1,
      5,  6,  7,  8, -1, -1
      ), nrow=5, byrow=TRUE)
conf$keyVarObs <- matrix(
    c( 0,  1,  2,  2,  2,  2,
      0,  1,  2,  2,  2,  2,
      0,  1,  2,  2,  2,  2,
      3,  4,  4,  4,  4, -1,
      5,  6,  6,  6, -1, -1
      ), nrow=5, byrow=TRUE)
conf$noScaledYears <- 13
conf$keyScaledYears <- 1993:2005
conf$keyParScaledYA <- row(matrix(NA, nrow=13, ncol=6))-1
conf$fbarRange <- c(2,4)
conf$seasonTimes <- c(0,1/3,2/3,1)
conf$isFishingSeason <- c(1,1,1)
conf$keyLogFseason[1,] <- 0
conf$keyLogFseason[2,] <- 0
conf$keyLogFseason[3,] <- 0

par<-defpar(dat,conf)
fit<-sam.fit(dat,conf,par)

cat(fit$opt$objective,"\n", file="res.out")
