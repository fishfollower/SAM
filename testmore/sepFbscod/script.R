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

#conf<-loadConf(dat,"model.cfg", patch=TRUE)
#par<-defpar(dat,conf)
#fit<-sam.fit(dat,conf,par)
#if(fit$opt$convergence!=0) stop("Model did not converge.")

conf<-loadConf(dat,"model.cfg", patch=TRUE)
conf$corFlag <- 3
par<-defpar(dat,conf)
fit2<-sam.fit(dat,conf,par)
if(fit2$opt$convergence!=0) stop("Model did not converge.")

conf<-loadConf(dat,"model.cfg", patch=TRUE)
conf$corFlag <- 4
par<-defpar(dat,conf)
fit3<-sam.fit(dat,conf,par)
if(fit3$opt$convergence!=0) stop("Model did not converge.")

#print(modeltable(c("AR"=fit, "separabel"=fit2,"separabelW"=fit3)))
cat(c(fit2$opt$objective,fit3$opt$objective), file="res.out")


