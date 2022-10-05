library(stockassessment)

cn<-read.ices("dat/cn.dat")
cw<-read.ices("dat/cw.dat")
dw<-read.ices("dat/dw.dat")
lw<-read.ices("dat/lw.dat")
mo<-read.ices("dat/mo.dat")
nm<-read.ices("dat/nm.dat")
pf<-read.ices("dat/pf.dat")
pm<-read.ices("dat/pm.dat")
sw<-read.ices("dat/sw.dat")
lf<-read.ices("dat/lf.dat")
surveys<-read.ices("dat/survey.dat")
cnA<-read.ices("dat/cn_A.dat")
cnC<-read.ices("dat/cn_C.dat")
cnD<-read.ices("dat/cn_D.dat")
cnF<-read.ices("dat/cn_F.dat")
cwA<-read.ices("dat/cw_A.dat")
cwC<-read.ices("dat/cw_C.dat")
cwD<-read.ices("dat/cw_D.dat")
cwF<-read.ices("dat/cw_F.dat")

idx<-!rownames(cn)%in%rownames(cnA)
sumACDF<-cn[idx,]
attr(sumACDF, "sumof")<-c(1,2,3,4)
                                        # temp fix to convergence problems 
cnA[cnA==0]<-1; cnA[,1]<-0
                                        #cnC[cnC==0]<-1
cnD[cnD==0]<-1
cnF[cnF==0]<-1





                                        # patch back in time:
cwA<-rbind(t(t(cw[!rownames(cw)%in%rownames(cwA),])), cwA)
cwC<-rbind(t(t(cw[!rownames(cw)%in%rownames(cwC),])), cwC)
cwD<-rbind(t(t(cw[!rownames(cw)%in%rownames(cwD),])), cwD)
cwF<-rbind(t(t(cw[!rownames(cw)%in%rownames(cwF),])), cwF)


dat<-setup.sam.data(surveys=surveys,
                    residual.fleet=list(cnA, cnC, cnD, cnF), 
                    sum.residual.fleets=sumACDF,
                    prop.mature=mo, 
                    stock.mean.weight=sw, 
                    catch.mean.weight=list(cwA, cwC, cwD, cwF), 
                    dis.mean.weight=dw, 
                    land.mean.weight=lw,
                    prop.f=pf, 
                    prop.m=pm, 
                    natural.mortality=nm, 
                    land.frac=lf)

conf<-loadConf(dat,"dat/model.cfg", patch=TRUE)
par<-defpar(dat,conf)
par$itrans_rho <- c(5,1,1)
### Changes below compared to 2019 to solve convergence problem for logSdLogFsta fleet C in 2020
mymap <- list( itrans_rho = factor(c(NA,0,1)), logSdLogFsta=factor(c(1,1,1,1)))
fit<-sam.fit(dat,conf,par,map=mymap)


mymap <- list( itrans_rho = factor(c(NA,0,1)))
par2<-fit$pl
par2$missing<-NULL
attr(par2, "what")<-NULL
 
fit<-sam.fit(dat,conf,par2,map=mymap)

cat(fit$opt$objective,file = "res.out")
