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

varC = as.matrix(read.table("varCommercial.txt", sep = " "))
attributes(cn)$weight = 1/varC
varS1 = as.matrix(read.table("varSurvey1.txt", sep = " "))
attributes(surveys[[1]])$weight = 1/varS1
varS2 = as.matrix(read.table("varSurvey2.txt", sep = " "))
attributes(surveys[[2]])$weight = 1/varS2
varS3 = as.matrix(read.table("varSurvey3.txt", sep = " "))
attributes(surveys[[3]])$weight = 1/varS3

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


conf = loadConf(dat,"modelConf.cfg")
par<-defpar(dat,conf)
par$logSdLogN = c(-0.35, -5)
map = list(logSdLogN = as.factor(c(0,NA)))
fit<-sam.fit(dat,conf,par,map =map)
cat(fit$opt$objective, file="res.out")
