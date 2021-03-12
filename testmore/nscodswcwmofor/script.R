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
fit0<-sam.fit(dat,conf,par)

conf$stockWeightModel<-1
conf$keyStockWeightMean<-0:5
conf$keyStockWeightObsVar<-rep(0,6)

conf$catchWeightModel<-1
conf$keyCatchWeightMean<-0:5
conf$keyCatchWeightObsVar<-rep(0,6)

conf$matureModel<-1
conf$keyMatureMean<-0:5

par<-defpar(dat,conf)
fit1<-sam.fit(dat,conf,par)

#par(mfrow=c(3,1))
#matplot(fit1$data$stockMeanWeight)
#matplot(exp(fit1$pl$logSW), type="l", add=TRUE)
#matplot(exp(fit1$pl$logSW-2*fit1$plsd$logSW), type="l", add=TRUE)
#matplot(exp(fit1$pl$logSW+2*fit1$plsd$logSW), type="l", add=TRUE)
#matplot(fit1$data$catchMeanWeight)
#matplot(exp(fit1$pl$logCW), type="l", add=TRUE)
#matplot(exp(fit1$pl$logCW-2*fit1$plsd$logCW), type="l", add=TRUE)
#matplot(exp(fit1$pl$logCW+2*fit1$plsd$logCW), type="l", add=TRUE)
#matplot(fit1$data$propMat)
#matplot(plogis(fit1$pl$logitMO), type="l", add=TRUE)
#matplot(plogis(fit1$pl$logitMO-2*fit1$plsd$logitMO), type="l", add=TRUE)
#matplot(plogis(fit1$pl$logitMO+2*fit1$plsd$logitMO), type="l", add=TRUE)

set.seed(123);
fc000<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=FALSE, useMOmodel=FALSE, processNoiseF=FALSE)
set.seed(123);
fc001<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=FALSE, useMOmodel=FALSE, processNoiseF=TRUE)
set.seed(123);
fc010<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=FALSE, useMOmodel=TRUE, processNoiseF=FALSE)
set.seed(123);
fc100<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=TRUE, useMOmodel=FALSE, processNoiseF=FALSE)
set.seed(123);
fc110<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=TRUE, useMOmodel=TRUE, processNoiseF=FALSE)
set.seed(123);
fc101<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=TRUE, useMOmodel=FALSE, processNoiseF=TRUE)
set.seed(123);
fc011<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=FALSE, useMOmodel=TRUE, processNoiseF=TRUE)
set.seed(123);
fc111<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=TRUE, useMOmodel=TRUE, processNoiseF=TRUE)

cat(capture.output(attr(fc000,"shorttab")), file = 'res.out', sep = '\n')
cat(capture.output(attr(fc001,"shorttab")), file = 'res.out', sep = '\n', append=TRUE)
cat(capture.output(attr(fc010,"shorttab")), file = 'res.out', sep = '\n', append=TRUE)
cat(capture.output(attr(fc100,"shorttab")), file = 'res.out', sep = '\n', append=TRUE)
cat(capture.output(attr(fc110,"shorttab")), file = 'res.out', sep = '\n', append=TRUE)
cat(capture.output(attr(fc101,"shorttab")), file = 'res.out', sep = '\n', append=TRUE)
cat(capture.output(attr(fc011,"shorttab")), file = 'res.out', sep = '\n', append=TRUE)
cat(capture.output(attr(fc111,"shorttab")), file = 'res.out', sep = '\n', append=TRUE)

 
