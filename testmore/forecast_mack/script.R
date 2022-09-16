library(stockassessment)

cn<-read.ices("../mack/cn.dat")
cw<-read.ices("../mack/cw.dat")
dw<-read.ices("../mack/dw.dat")
lf<-read.ices("../mack/lf.dat")
lw<-read.ices("../mack/lw.dat")
mo<-read.ices("../mack/mo.dat")
nm<-read.ices("../mack/nm.dat")
pf<-read.ices("../mack/pf.dat")
pm<-read.ices("../mack/pm.dat")
sw<-read.ices("../mack/sw.dat")
surveys<-read.ices("../mack/survey.dat")

recap<-read.table("../mack/tag.dat", header=TRUE)
recap<-recap[recap$Type==1 & recap$RecaptureY<=2006,]
recap<-recap[recap[,1]>=min(as.numeric(rownames(sw))), ]

W<-matrix(NA,nrow=nrow(cn), ncol=ncol(cn))
W[as.numeric(rownames(cn))<2000]<-10
attr(cn,"weight")<-W

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
                    land.frac=lf,
                    recapture=recap)

conf<-defcon(dat)
conf$fbarRange <- c(4,8)
conf$corFlag <- 2
conf$fixVarToWeight<-1

par<-defpar(dat,conf)
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
target <- round(c(0.1,0.1)*mean(tail(rep(catchtable(fit)[,1],4))))
fv <- modelforecast(fit, landval=target, biasCorrect = FALSE, splitLD = TRUE,ave.years = max(fit$data$years)+(-4:0), rec.years = max(fit$data$years)+(-4:0))
cat("landval tab:\t",paste(all.equal(target,
                     round(attr(fv,"tab")[-1,"Land:median"]), check.attributes = FALSE, tolerance = tol), collapse = ", "),"\n", file="res.out", append = TRUE) 


## Test fmodel
target <- rep(NA,5)
fv <- modelforecast(fit, landval = target, biasCorrect = FALSE,ave.years = max(fit$data$years)+(-4:0))
cat("fmodel tab:\t",paste(all.equal(rep(tail(fbartable(fit)[,1],1),length(target)),
                     attr(fv,"tab")[-1,"fbar:median"], check.attributes = FALSE, tolerance = tol), collapse = ", "),"\n", file="res.out", append = TRUE) 
cat("fmodel tab:\t",paste(all.equal(c(fbartable(fit)[,1],rep(tail(fbartable(fit)[,1],1),length(target))),
                     attr(fv,"fulltab")[,"fbar:median"], check.attributes = FALSE, tolerance = tol), collapse = ", "),"\n", file="res.out", append = TRUE) 


