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

set.seed(123)
surveys$FIDX[,1]<-exp(log(surveys$FIDX[,1])+rnorm(26,sd=.1))

surveys<-surveys[1:2]

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

conf$corFlag <- 2
conf$keyVarF[1,] <- c(0, 1, 1, 1, 1, 1)
conf$noScaledYears <- 13
conf$keyScaledYears <- 1993:2005
conf$keyParScaledYA <- row(matrix(NA, nrow=13, ncol=6))-1
conf$fbarRange <- c(2,4)
conf$keyBiomassTreat[4]<-10
par<-defpar(dat,conf)
fit<-sam.fit(dat,conf,par)

cat(fit$opt$objective,"\n", file="res.out")
