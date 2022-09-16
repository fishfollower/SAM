# to install the package from the multi branch 
#devtools::install_github("fishfollower/SAM/stockassessment", ref="multi")

library(stockassessment)
cn<-read.ices("old/cn.dat")
cw<-read.ices("old/cw.dat")
dw<-cw
lf<-cn; lf[]<-1
lw<-cw
mo<-read.ices("old/mo.dat")
nm<-read.ices("old/nm.dat")
pf<-read.ices("old/pf.dat")
pm<-read.ices("old/pm.dat")
sw<-read.ices("old/sw.dat")
surveys<-read.ices("old/survey.dat")

cnA<-read.ices("split/cn_A.dat")
cnC<-read.ices("split/cn_C.dat")
cnD<-read.ices("split/cn_D.dat")
cnF<-read.ices("split/cn_F.dat")

idx<-!rownames(cn)%in%rownames(cnA)


#sumACDF<-(cnA+cnC+cnD+cnF)[idx,]
sumACDF<-cn[idx,]
attr(sumACDF, "sumof")<-c(1,2,3,4)

# temp fix to convergence problems 
cnA[cnA==0]<-1; cnA[,1]<-0
cnC[cnC==0]<-1
cnD[cnD==0]<-1
cnF[cnF==0]<-1

dat<-setup.sam.data(surveys=surveys,
                    residual.fleets=list(cnA, cnC, cnD, cnF), # Notice list
                    sum.residual.fleets=sumACDF,
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
par<-defpar(dat,conf)
fit<-sam.fit(dat,conf,par)

cat(fit$opt$objective,"\n\n", file="res.out")
cat(capture.output(prmatrix(t(fit$pl$logF))), sep="\n", file="res.out", append=TRUE)
