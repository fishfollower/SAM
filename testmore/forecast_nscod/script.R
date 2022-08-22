library(stockassessment)
cn<-read.ices("../nscod/cn.dat")
cw<-read.ices("../nscod/cw.dat")
dw<-read.ices("../nscod/dw.dat")
lf<-read.ices("../nscod/lf.dat")
lw<-read.ices("../nscod/lw.dat")
mo<-read.ices("../nscod/mo.dat")
nm<-read.ices("../nscod/nm.dat")
pf<-read.ices("../nscod/pf.dat")
pm<-read.ices("../nscod/pm.dat")
sw<-read.ices("../nscod/sw.dat")
surveys<-read.ices("../nscod/survey.dat")

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

par<-defpar(dat,conf)
fit<-sam.fit(dat,conf,par)

catchtable(fit)

tol <- 1e-6
## Test fscale
target <- c(0.9,0.9,0.9,1.1,1.1,1.1,1,1,1)
fv <- modelforecast(fit, fscale = target, biasCorrect = FALSE, ave.years = max(fit$data$years)+(-4:0))
cat("fscale tab:\t",paste(all.equal(tail(fbartable(fit)[,1],1) * cumprod(c(1,target)),
                     attr(fv,"tab")[,"fbar:Estimate"], check.attributes = FALSE, tolerance = tol), collapse = ", "),"\n", file="res.out")
cat("fscale fulltab:\t",paste(all.equal(c(head(fbartable(fit)[,1],-1),tail(fbartable(fit)[,1],1) * cumprod(c(1,target))),
                                        attr(fv,"fulltab")[,"fbar:Estimate"], check.attributes = FALSE, tolerance = tol), collapse = ", "),"\n", file="res.out", append = TRUE)


## Test fval
target <- c(0.1,0.2,0.3,0.1)
fv <- modelforecast(fit, fval = target, biasCorrect = FALSE,ave.years = max(fit$data$years)+(-4:0))
cat("fval tab:\t",paste(all.equal(c(tail(fbartable(fit)[,1],1), target),
                 attr(fv,"tab")[,"fbar:Estimate"], check.attributes = FALSE, tolerance = tol)),"\n", file="res.out", append = TRUE) 
cat("fval fulltab:\t",paste(all.equal(c(head(fbartable(fit)[,1],-1),tail(fbartable(fit)[,1],1), target),
                     attr(fv,"fulltab")[,"fbar:Estimate"], check.attributes = FALSE, tolerance = tol), collapse = ", "),"\n", file="res.out", append = TRUE)

## Test catchval
target <- round(c(0.5,0.3) * mean(tail(rep(catchtable(fit)[,1],4))))
fv <- modelforecast(fit, catchval = target, biasCorrect = FALSE,ave.years = max(fit$data$years)+(-4:0))
cat("catchval tab:\t",paste(all.equal(target,
                     round(attr(fv,"tab")[-1,"catch:Estimate"]), check.attributes = FALSE, tolerance = tol), collapse = ", "),"\n", file="res.out", append = TRUE) 
cat("catchval fulltab:\t",paste(all.equal(head(catchtable(fit)[,1],-1),
                     attr(fv,"fulltab")[1:(nrow(catchtable(fit))-1),"catch:Estimate"], check.attributes = FALSE, tolerance = tol), collapse = ", "),"\n", file="res.out", append = TRUE)

## Test landval
target <- round(c(0.1,0.1) * mean(tail(rep(catchtable(fit)[,1],4))))
fv <- modelforecast(fit, landval = target, biasCorrect = FALSE, splitLD = TRUE,ave.years = max(fit$data$years)+(-4:0))
cat("landval tab:\t",paste(all.equal(target,
                     round(attr(fv,"tab")[-1,"Land:Estimate"]), check.attributes = FALSE, tolerance = tol), collapse = ", "),"\n", file="res.out", append = TRUE) 


## Test fmodel
target <- rep(NA,5)
fv <- modelforecast(fit, landval = target, biasCorrect = FALSE,ave.years = max(fit$data$years)+(-4:0))
cat("fmodel tab:\t",paste(all.equal(rep(tail(fbartable(fit)[,1],1),length(target)),
                     attr(fv,"tab")[-1,"fbar:Estimate"], check.attributes = FALSE, tolerance = tol), collapse = ", "),"\n", file="res.out", append = TRUE) 
cat("fmodel tab:\t",paste(all.equal(c(fbartable(fit)[,1],rep(tail(fbartable(fit)[,1],1),length(target))),
                     attr(fv,"fulltab")[,"fbar:Estimate"], check.attributes = FALSE, tolerance = tol), collapse = ", "),"\n", file="res.out", append = TRUE) 


