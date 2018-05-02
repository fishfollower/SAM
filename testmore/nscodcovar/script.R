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

## Estimate a correlation 

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
conf$fbarRange <- c(2,4)
conf$obsCorStruct[] <- c("ID","ID","AR")
conf$keyCorObs[3,1:3] <- 0

par<-defpar(dat,conf)
fit<-sam.fit(dat,conf,par)

## assign same correlation 

estcor<-cov2cor(fit$rep$obsCov[[3]])
attr(surveys[[2]], "cor") <- lapply(1:23, function(i)estcor)

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
conf$fbarRange <- c(2,4)
par<-defpar(dat,conf)
fit2<-sam.fit(dat,conf,par)

sink("res.out")
modeltable(c(ESTMATED=fit,ASSIGNED=fit2))
sink()
