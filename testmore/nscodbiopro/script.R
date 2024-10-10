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
conf$fbarRange <- c(2,4)

par<-defpar(dat,conf)
fit0<-sam.fit(dat,conf,par)

conf$stockWeightModel<-1
conf$keyStockWeightMean<-0:5
conf$keyStockWeightObsVar<-rep(0,6)

conf$catchWeightModel<-1
conf$keyCatchWeightMean[1,]<-0:5
conf$keyCatchWeightObsVar[1,]<-rep(0,6)

conf$matureModel<-1
conf$keyMatureMean<-0:5

conf$mortalityModel<-1
conf$keyMortalityMean<-0:5
conf$keyMortalityObsVar<-rep(0,6)

par<-defpar(dat,conf)
fit1<-sam.fit(dat,conf,par)

m0 <- modelforecast(fit0,rep("F=0.2",2),nosim=100,useModelLastN=TRUE,year.base=2018,resampleFirst=TRUE)

round(cov(m0[[2]]$sim),2)
round(fit0$sdrep$covY,2)

m0x <- modelforecast(fit0,rep("F=0.2",2),nosim=100,year.base=2019,resampleFirst=TRUE)
m1 <- modelforecast(fit1,rep("F=0.2",2),nosim=1000,useModelLastN=TRUE,year.base=2018)

matplot(fit1$dat$year[1] + seq_len(nrow(fit1$pl$logSW))-1,fit1$pl$logSW)
for(y in 1:3)
for(i in 1:6)
    points(m1[[y]]$year,mean(log(m1[[y]]$bio_stockMeanWeight[,i])),pch=5)
points(rep(2019,6),fit1$sdrep$estY[13:18])
points(rep(2018,6),fit1$sdrep$estYm1[13:18])

m1x <- modelforecast(fit1,rep("F=0.2",2),nosim=100,year.base=2019,resampleFirst=TRUE)

definit <- defpar(dat, conf)
if (!identical(par, relist(unlist(par), skeleton = definit))) 
        warning("Initial values are not consistent, so running with default init values from defpar()")


#par(mfrow=c(2,2))
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
#matplot(fit1$data$natMor)
#matplot(exp(fit1$pl$logNM), type="l", add=TRUE)
#matplot(exp(fit1$pl$logNM-2*fit1$plsd$logNM), type="l", add=TRUE)
#matplot(exp(fit1$pl$logNM+2*fit1$plsd$logNM), type="l", add=TRUE)

 

set.seed(123);
fc0000<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=FALSE, useCWmodel=FALSE, useMOmodel=FALSE, useNMmodel=FALSE, processNoiseF=FALSE)
set.seed(123);
fc0001<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=FALSE, useCWmodel=FALSE, useMOmodel=FALSE, useNMmodel=TRUE, processNoiseF=FALSE)
set.seed(123);
fc0010<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=FALSE, useCWmodel=FALSE, useMOmodel=TRUE, useNMmodel=FALSE, processNoiseF=FALSE)
set.seed(123);
fc0100<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=FALSE, useCWmodel=TRUE, useMOmodel=FALSE, useNMmodel=FALSE, processNoiseF=FALSE)
set.seed(123);
fc1000<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=TRUE, useCWmodel=FALSE, useMOmodel=FALSE, useNMmodel=FALSE, processNoiseF=FALSE)
set.seed(123);
fc1001<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=TRUE, useCWmodel=FALSE, useMOmodel=FALSE, useNMmodel=TRUE, processNoiseF=FALSE)
set.seed(123);
fc1010<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=TRUE, useCWmodel=FALSE, useMOmodel=TRUE, useNMmodel=FALSE, processNoiseF=FALSE)
set.seed(123);
fc1100<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=TRUE, useCWmodel=TRUE, useMOmodel=FALSE, useNMmodel=FALSE, processNoiseF=FALSE)
set.seed(123);
fc0101<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=FALSE, useCWmodel=TRUE, useMOmodel=FALSE, useNMmodel=TRUE, processNoiseF=FALSE)
set.seed(123);
fc0110<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=FALSE, useCWmodel=TRUE, useMOmodel=TRUE, useNMmodel=FALSE, processNoiseF=FALSE)
set.seed(123);
fc0011<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=FALSE, useCWmodel=FALSE, useMOmodel=TRUE, useNMmodel=TRUE, processNoiseF=FALSE)
set.seed(123);
fc1110<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=TRUE, useCWmodel=TRUE, useMOmodel=TRUE, useNMmodel=FALSE, processNoiseF=FALSE)
set.seed(123);
fc1101<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=TRUE, useCWmodel=TRUE, useMOmodel=FALSE, useNMmodel=TRUE, processNoiseF=FALSE)
set.seed(123);
fc1011<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=TRUE, useCWmodel=FALSE, useMOmodel=TRUE, useNMmodel=TRUE, processNoiseF=FALSE)
set.seed(123);
fc0111<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=FALSE, useCWmodel=TRUE, useMOmodel=TRUE, useNMmodel=TRUE, processNoiseF=FALSE)
set.seed(123);
fc1111<-forecast(fit1, fscale=c(1,1,1,1), useSWmodel=TRUE, useCWmodel=TRUE, useMOmodel=TRUE, useNMmodel=TRUE, processNoiseF=FALSE)


cat(capture.output(attr(fc0000,"shorttab")), file = 'res.out', sep = '\n')
cat(capture.output(attr(fc0001,"shorttab")), file = 'res.out', sep = '\n', append=TRUE)
cat(capture.output(attr(fc0010,"shorttab")), file = 'res.out', sep = '\n', append=TRUE)
cat(capture.output(attr(fc0100,"shorttab")), file = 'res.out', sep = '\n', append=TRUE)
cat(capture.output(attr(fc1000,"shorttab")), file = 'res.out', sep = '\n', append=TRUE)
cat(capture.output(attr(fc1001,"shorttab")), file = 'res.out', sep = '\n', append=TRUE)
cat(capture.output(attr(fc1010,"shorttab")), file = 'res.out', sep = '\n', append=TRUE)
cat(capture.output(attr(fc1100,"shorttab")), file = 'res.out', sep = '\n', append=TRUE)
cat(capture.output(attr(fc0101,"shorttab")), file = 'res.out', sep = '\n', append=TRUE)
cat(capture.output(attr(fc0110,"shorttab")), file = 'res.out', sep = '\n', append=TRUE)
cat(capture.output(attr(fc0011,"shorttab")), file = 'res.out', sep = '\n', append=TRUE)
cat(capture.output(attr(fc1110,"shorttab")), file = 'res.out', sep = '\n', append=TRUE)
cat(capture.output(attr(fc1101,"shorttab")), file = 'res.out', sep = '\n', append=TRUE)
cat(capture.output(attr(fc1011,"shorttab")), file = 'res.out', sep = '\n', append=TRUE)
cat(capture.output(attr(fc0111,"shorttab")), file = 'res.out', sep = '\n', append=TRUE)
cat(capture.output(attr(fc1111,"shorttab")), file = 'res.out', sep = '\n', append=TRUE)
