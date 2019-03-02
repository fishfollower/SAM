library(stockassessment)
library(TMB)
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

conf$corFlag <- 3
par<-defpar(dat,conf)

lower = list();upper = list();
lower$sepFalpha = -5
upper$sepFalpha = 5

map = list()
map$sepFalpha = as.factor(c(0,1,2,3,4,NA))
par$sepFlogitRho = c(0,10)
map$sepFlogitRho = as.factor(c(0,NA))
fit<-sam.fit(dat,conf,par,lower = lower, upper = upper,map = map)


conf$corFlag <- 0
par<-defpar(dat,conf)
fit2<-sam.fit(dat,conf,par)

conf$corFlag <- 1
par<-defpar(dat,conf)
fit3<-sam.fit(dat,conf,par)

conf$corFlag <- 2
par<-defpar(dat,conf)
fit4<-sam.fit(dat,conf,par)


print(modeltable(c("separabel"=fit, "IID"=fit2, "cs"=fit3, "AR"=fit4)))

fselectivityplot(fit)#Strange selevtivityplot...
fselectivityplot(fit2)
