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
recap<-recap[recap$Type==1 & recap$RecaptureY<=2006,]
recap<-recap[recap[,1]>=min(as.numeric(rownames(sw))), ]

recap2<-read.table("tag2.dat", header=TRUE)
recap2$r<-round(recap2$r)
recap2<-recap2[recap2$Nscan>0,]
recap<-rbind(recap,recap2)

age<- recap$ReleaseY - recap$Yearclass
recap  <- recap[age>1 & age<13,]

#### Remove the tags recapture in the release year
recap <- recap[recap$ReleaseY!=recap$RecaptureY,]     

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

#### model configuration

conf<-defcon(dat)
conf$keyLogFsta[1,] <- c(0,1,2,3,4,5,6,7,7,7,7,7,7)
conf$keyVarObs[4,]      <- c(-1,-1,-1, 3,4,4,4,4,4,4,4,4,-1)
conf$fbarRange <- c(4,8)
conf$corFlag <- 0
conf$fixVarToWeight<-1
conf$obsCorStruct[] <- c("ID", "ID", "ID", "AR", "ID")
conf$keyCorObs[4,4:11]<-0

par<-defpar(dat,conf)

# model fitting
fit <-sam.fit(dat,conf,par, newtonsteps=0)
