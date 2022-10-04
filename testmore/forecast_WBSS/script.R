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
                    catch.mean.weight=cw, 
                    dis.mean.weight=dw, 
                    land.mean.weight=lw,
                    prop.f=pf, 
                    prop.m=pm, 
                    natural.mortality=nm, 
                    land.frac=lf)

conf<-defcon(dat)

conf$corFlag<-c(0,0,0,0)
par<-defpar(dat,conf)
fit<-sam.fit(dat,conf,par)

cat(fit$opt$objective,"\n", file="res.out")

checkNum <- function(val,test){
    test <- rep(test, length.out = length(val))
    isTRUE(all.equal.numeric(val,test, check.attributes=FALSE))
}

## Test of multi fleet constrained forecast
## Total F at ICES Fmsy: F = 0.31
## (A fleet #1) Unintended catch from North Sea at 6142t: C[1]=6142
## (F fleet #4) Directed catch in subdivisions 22-24 is half of total: C[4]=0.5*C
## (D fleet #3) Bycatch from other fisheries in subdivisions 20-24 is practically zero: C[3]=0.001
suppressWarnings(v <- modelforecast(fit,rep("F=0.31 & C[1]=6142 & C[4]=0.5*C & C[3]=1",5),
                   rec.years = tail(fit$data$years,5),
                   ave.years = tail(fit$data$years,5),
                   nosim = 0, newton_config=list(grad_tol=1e-7)))

cat(checkNum(tail(fbartable(v)[,1],-1),0.31),"\n", file="res.out", append=TRUE)
Ctab <- attr(v,"catchby")
CtabMed <- Ctab[rownames(Ctab) %in% "mostLikelyTrajectory",]
CtabRel <- CtabMed / rowSums(CtabMed)[row(CtabMed)]
cat(checkNum(CtabMed[-1,1],6142),"\n", file="res.out", append=TRUE)
cat(checkNum(CtabMed[-1,3],1),"\n", file="res.out", append=TRUE)
cat(checkNum(CtabRel[-1,4],0.5),"\n", file="res.out", append=TRUE)

suppressWarnings(v <- modelforecast(fit,rep("F=0.2",5),
                   rec.years = tail(fit$data$years,5),
                   ave.years = tail(fit$data$years,5),
                   nosim = 0, newton_config=list(grad_tol=1e-7)))
cat(checkNum(tail(fbartable(v)[,1],-1),0.2),"\n", file="res.out", append=TRUE)
Ftab <- attr(v,"fbarby")
FtabMed <- Ftab[rownames(Ftab) %in% "mostLikelyTrajectory",]
FtabRel <- FtabMed[-1,] / FtabMed[1,][col(FtabMed[-1,])]
cat(checkNum(as.vector(FtabRel[]),FtabRel[1]),"\n", file="res.out", append=TRUE)


v <- forecast(fit, fscale = c(1,rep(NA,5)), fval = c(NA,rep(0.31,5)))
cat(checkNum(tail(fbartable(v)[,1],-1),0.31),"\n", file="res.out", append=TRUE)
