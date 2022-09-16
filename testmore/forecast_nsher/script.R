library(stockassessment)
cn<-read.ices("../nsher/cn.dat")
cw<-read.ices("../nsher/cw.dat")
dw<-read.ices("../nsher/dw.dat")
lf<-read.ices("../nsher/lf.dat")
lw<-read.ices("../nsher/lw.dat")
mo<-read.ices("../nsher/mo.dat")
nm<-read.ices("../nsher/nm.dat")
pf<-read.ices("../nsher/pf.dat")
pm<-read.ices("../nsher/pm.dat")
sw<-read.ices("../nsher/sw.dat")
surveys<-read.ices("../nsher/survey.dat")

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

tol <- 1e-6
## Test fscale
target <- c(0.9,0.9,0.9,1.1,1.1,1.1,1,1,1)
fv <- modelforecast(fit, fscale = target, biasCorrect = FALSE,ave.years = max(fit$data$years)+(-4:0))
cat("fscale tab:\t",paste(all.equal(tail(fbartable(fit)[,1],1) * cumprod(c(1,target)),
                     attr(fv,"tab")[,"fbar:median"], check.attributes = FALSE, tolerance = tol), collapse = ", "),"\n", file="res.out")
cat("fscale fulltab:\t",paste(all.equal(c(head(fbartable(fit)[,1],-1),tail(fbartable(fit)[,1],1) * cumprod(c(1,target))),
                                        attr(fv,"fulltab")[,"fbar:median"], check.attributes = FALSE, tolerance = tol), collapse = ", "),"\n", file="res.out", append = TRUE)


## Test fval
target <- c(0.1,0.2,0.3,0.1)
fv <- modelforecast(fit, fval = target, biasCorrect = FALSE,ave.years = max(fit$data$years)+(-4:0))
cat("fval tab:\t",paste(all.equal(c(tail(fbartable(fit)[,1],1), target),
                 attr(fv,"tab")[,"fbar:median"], check.attributes = FALSE, tolerance = tol)),"\n", file="res.out", append = TRUE) 
cat("fval fulltab:\t",paste(all.equal(c(head(fbartable(fit)[,1],-1),tail(fbartable(fit)[,1],1), target),
                     attr(fv,"fulltab")[,"fbar:median"], check.attributes = FALSE, tolerance = tol), collapse = ", "),"\n", file="res.out", append = TRUE)

## Test catchval
target <- round(c(0.5,0.3) * mean(tail(rep(catchtable(fit)[,1],4))))
fv <- modelforecast(fit, catchval = target, biasCorrect = FALSE,ave.years = max(fit$data$years)+(-4:0))
cat("catchval tab:\t",paste(all.equal(target,
                     round(attr(fv,"tab")[-1,"catch:median"]), check.attributes = FALSE, tolerance = tol), collapse = ", "),"\n", file="res.out", append = TRUE) 
cat("catchval fulltab:\t",paste(all.equal(head(catchtable(fit)[,1],-1),
                     attr(fv,"fulltab")[1:(nrow(catchtable(fit))-1),"catch:median"], check.attributes = FALSE, tolerance = tol), collapse = ", "),"\n", file="res.out", append = TRUE)

## Test landval
target <- round(c(0.1,0.1) * mean(tail(rep(catchtable(fit)[,1],4))))
fv <- modelforecast(fit, landval = target, biasCorrect = FALSE, splitLD = TRUE,ave.years = max(fit$data$years)+(-4:0))
cat("landval tab:\t",paste(all.equal(target,
                     round(attr(fv,"tab")[-1,"Land:median"]), check.attributes = FALSE, tolerance = tol), collapse = ", "),"\n", file="res.out", append = TRUE) 


## Test fmodel
target <- rep(NA,5)
fv <- modelforecast(fit, landval = target, biasCorrect = FALSE,ave.years = max(fit$data$years)+(-4:0))
cat("fmodel tab:\t",paste(all.equal(rep(tail(fbartable(fit)[,1],1),length(target)),
                     attr(fv,"tab")[-1,"fbar:median"], check.attributes = FALSE, tolerance = tol), collapse = ", "),"\n", file="res.out", append = TRUE) 
cat("fmodel tab:\t",paste(all.equal(c(fbartable(fit)[,1],rep(tail(fbartable(fit)[,1],1),length(target))),
                     attr(fv,"fulltab")[,"fbar:median"], check.attributes = FALSE, tolerance = tol), collapse = ", "),"\n", file="res.out", append = TRUE) 


