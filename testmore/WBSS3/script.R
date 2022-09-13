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

# patch back in time:
cntot<-(cnA+cnC+cnD+cnF)
fracA<-colMeans(cnA/cntot)
fracC<-colMeans(cnC/cntot)
fracD<-colMeans(cnD/cntot)
fracF<-colMeans(cnF/cntot)
cnA<-rbind(t(t(cn[!rownames(cn)%in%rownames(cnA),])*fracA), cnA)
cnC<-rbind(t(t(cn[!rownames(cn)%in%rownames(cnC),])*fracC), cnC)
cnD<-rbind(t(t(cn[!rownames(cn)%in%rownames(cnD),])*fracD), cnD)
cnF<-rbind(t(t(cn[!rownames(cn)%in%rownames(cnF),])*fracF), cnF)

dat<-setup.sam.data(surveys=surveys,
                    residual.fleets=list(cnA, cnC, cnD, cnF), # Notice list 
                    prop.mature=mo, 
                    stock.mean.weight=sw, 
                    catch.mean.weight=list(cw,cw*2,cw*3,cw*4),
                    dis.mean.weight=dw, 
                    land.mean.weight=lw,
                    prop.f=list(pf,pf*2,pf*3,pf*4),
                    prop.m=pm, 
                    natural.mortality=nm, 
                    land.frac=lf)

conf<-defcon(dat)

conf$corFlag<-c(2,1,0,2)
par<-defpar(dat,conf)
fit<-sam.fit(dat,conf,par)
corPar <- fit$pl$itrans_rho
cor<-round(cov2cor(fit$obj$report(fit$obj$env$last.par.best)$fvar),3)
cat(fit$opt$objective,"\n\n", file="res.out")
cat(corPar, sep="\n", file="res.out", append=TRUE)
options(width=1000) 
cat(capture.output(prmatrix(cor)), sep="\n", file="res.out", append=TRUE)
cat(capture.output(prmatrix(fbartable(fit))), sep="\n", file="res.out", append=TRUE)
