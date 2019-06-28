library(stockassessment)

cn<-read.ices("cn.dat")

tcn <- t(cn)
tcn[c(2, 4, 5, 6, 7, 8, 12, 14, 17, 19, 25)] <- NA
cn1 <- t(tcn)

cw1<-read.ices("cw.dat")
dw1<-read.ices("dw.dat")
lf1<-read.ices("lf.dat")
lw1<-read.ices("lw.dat")
mo1<-read.ices("mo.dat")
nm1<-read.ices("nm.dat")
pf1<-read.ices("pf.dat")
pm1<-read.ices("pm.dat")
sw1<-read.ices("sw.dat")
surveys<-read.ices("survey.dat")
dat<-setup.sam.data(surveys=surveys,
                    residual.fleet=cn1, 
                    prop.mature=mo1, 
                    stock.mean.weight=sw1, 
                    catch.mean.weight=cw1, 
                    dis.mean.weight=dw1, 
                    land.mean.weight=lw1,
                    prop.f=pf1, 
                    prop.m=pm1, 
                    natural.mortality=nm1, 
                    land.frac=lf1)

conf<-defcon(dat)
    
par<-defpar(dat,conf)
   
fit<-sam.fit(dat,conf,par)

fit2<-sam.fit(dat,conf,par)

cat(fit$opt$objective,fit2$opt$objective ,"\n", file="res.out")
