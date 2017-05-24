library(stockassessment)

oldwd<-setwd("data")

  cn<-read.ices("cn.dat")
  cw<-read.ices("cw.dat")
  dw<-read.ices("dw.dat")
  lw<-read.ices("lw.dat")
  mo<-read.ices("mo.dat")
  nm<-read.ices("nm.dat")
  pf<-read.ices("pf.dat")
  pm<-read.ices("pm.dat")
  sw<-read.ices("sw.dat")
  lf<-read.ices("lf.dat")
  surveys<-read.ices("survey.dat")

setwd(oldwd)

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


save(dat, file="run/data.RData")
