library(stockassessment)
cn<-read.ices("bfte2014-CANUM-reported.txt")
dn<-read.ices("bfte2014-DISNUM.txt")
ln<-read.ices("bfte2014-LANUM.txt")
cw<-read.ices("bfte2014-WECA.txt")
dw<-read.ices("bfte2014-WEDIS.txt")
lw<-read.ices("bfte2014-WELAND.txt")
mo<-read.ices("bfte2014-MATPROP.txt")
nm<-read.ices("bfte2014-NATMOR.txt")
pf<-read.ices("bfte2014-FPROP.txt")
pm<-read.ices("bfte2014-MPROP.txt")
sw<-read.ices("bfte2014-WEST.txt")
lf<-ln/cn

surveys<-read.ices("index.txt")

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

conf <- defcon(dat)
conf$corFlag <- 2
conf$keyLogFsta[1,] <- c(0:8,8)
conf$keyLogFpar[5,] <- c(-1,   -1,   -1 ,   9,   10,   11,   12,   13,   14,    15)
conf$keyLogFpar[6,] <- c(-1,   -1,   -1,   -1,   16,   17,   -1,   -1,   -1,    -1)
conf$keyLogFpar[7,] <- c(-1,   18,   19,   -1,   -1,   -1,   -1,   -1,   -1,    -1)
conf$keyLogFpar[8,] <- c(-1,   -1,   20,   21,   22,   23,   -1,   -1,   -1,    -1)
conf$obsCorStruct[] <- c("ID", "AR", "AR", "ID", "ID", "ID", "ID", "ID")
conf$keyCorObs[2,6:9] <- c(0,1,1,1) 
conf$keyCorObs[3,6:9] <- c(2,3,3,3) 
conf$keyVarF[1,] <- c(0,1,1,1,1,1,1,1,1,1)
par <- defpar(dat,conf)

fit <- sam.fit(dat,conf,par)

cat(fit$opt$objective,"\n", file="res.out")
