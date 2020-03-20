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

dat<-setup.sam.data(surveys=surveys,
                    residual.fleet=cn, 
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
conf$corFlag <- 2
conf$keyLogFpar <- matrix(
         c(-1, -1, -1, -1, -1, -1,
            0,  1,  2,  3,  4, -1,
            5,  6,  7,  8, -1, -1
           ), nrow=3, byrow=TRUE)
conf$keyVarF[1,] <- c(0, 1, 1, 1, 1, 1)
conf$keyVarObs <- matrix(
         c( 0,  1,  2,  2,  2,  2,
            3,  4,  4,  4,  4, -1,
            5,  6,  6,  6, -1, -1
           ), nrow=3, byrow=TRUE)
conf$noScaledYears <- 13
conf$keyScaledYears <- 1993:2005
conf$keyParScaledYA <- row(matrix(NA, nrow=13, ncol=6))-1
conf$fbarRange <- c(2,4)

conf$predVarObsLink<-rbind(c(0,0,0,0,0,0),
                           c(1,1,1,1,1,NA),
                           c(2,2,2,2,NA,NA))

par<-defpar(dat,conf)
fit<-sam.fit(dat,conf,par, sim.condRE=FALSE)

set.seed(12345)
cc<-TMB:::checkConsistency(fit$obj, n=500)
summary(cc)
