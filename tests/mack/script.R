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

recap<-read.table("tag.dat", header=TRUE)

recap<-recap[recap[,1]>=min(as.numeric(rownames(sw))), ]

W<-matrix(1,nrow=nrow(cn), ncol=ncol(cn))
W[as.numeric(rownames(cn))<2000]<-100
attr(cn,"weight")<-W
#cn[1,1:11]<-NA
#cn[2:5,]<-NA


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

cat(fit$opt$objective,"\n", file="res.out") 
