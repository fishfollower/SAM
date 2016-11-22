stopifnot(require("testthat"),
		require("stockassessment"))

context("nsher example")
setwd("../nsher/")
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

conf$fbarRange<-c(2,6)
conf$corFlag<-1
conf$keyLogFpar<-matrix(c(
-1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
-1,    0,    1,    2,    3,    4,    5,    6,   -1,
-1,    7,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
 8,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1), nrow=4, byrow=TRUE)

par<-defpar(dat,conf)
par$logFpar<-rep(0,9)

fit<-sam.fit(dat,conf,par)
load("fit.expected.Rdata")#reads in fit.exp
compare.fits(fit, fit.exp)
