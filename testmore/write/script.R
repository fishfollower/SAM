library(stockassessment)

toNA <- function(x){
    if(is.list(x)){
        return(lapply(x, function(y){
            y[y <= 0] <- NA
            y
        }))
    }
    x[x <= 0] <- NA
    x
}

tol <- 1e-5

dir <- tempdir()
cat("Read/Write","\n",file = "res.out",append=FALSE)

##########################################################################################
## Attribute test
##########################################################################################


m <- read.ices("test.dat")
write.ices(m,file.path(dir,"test.dat"))
m2 <- read.ices(file.path(dir,"test.dat"))

cat("Attributes","\n",file = "res.out",append=TRUE)
cat(all.equal(m,m2),"\n",file = "res.out",append=TRUE)
cat(all.equal(attr(m,"sumof"),1:4),"\n",file = "res.out",append=TRUE)
cat(all.equal(attr(m,"testAttrib"),"2-4"),"\n",file = "res.out",append=TRUE)

##########################################################################################
## Single fleet data
##########################################################################################

cn<-read.ices("sf/cn.dat")
cw<-read.ices("sf/cw.dat")
dw<-read.ices("sf/dw.dat")
lw<-read.ices("sf/lw.dat")
mo<-read.ices("sf/mo.dat")
nm<-read.ices("sf/nm.dat")
pf<-read.ices("sf/pf.dat")
pm<-read.ices("sf/pm.dat")
sw<-read.ices("sf/sw.dat")
lf<-read.ices("sf/lf.dat")
surveys<-read.ices("sf/survey.dat")
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

dir.create(file.path(dir,"sf"))
dir_sf <- file.path(dir,"sf")

write.data.files(dat,dir_sf)

cat("Single","\n",file = "res.out",append=TRUE)
cat(all.equal(toNA(cn),read.ices(file.path(dir_sf,"cn.dat")), tolerance= tol, check.attributes=FALSE),"\n",file = "res.out",append=TRUE)
cat(all.equal((cw),read.ices(file.path(dir_sf,"cw.dat")), tolerance= tol),"\n",file = "res.out", append=TRUE)
cat(all.equal((dw),read.ices(file.path(dir_sf,"dw.dat")), tolerance= tol), "\n", file = "res.out", append=TRUE)
cat(all.equal((lw),read.ices(file.path(dir_sf,"lw.dat")), tolerance= tol), "\n", file = "res.out", append=TRUE)
cat(all.equal(mo,read.ices(file.path(dir_sf,"mo.dat")), tolerance= tol), "\n", file = "res.out", append=TRUE)
cat(all.equal(nm,read.ices(file.path(dir_sf,"nm.dat")), tolerance= tol), "\n", file = "res.out", append=TRUE)
cat(all.equal(pf,read.ices(file.path(dir_sf,"pf.dat")), tolerance= tol), "\n", file = "res.out", append=TRUE)
cat(all.equal(pm,read.ices(file.path(dir_sf,"pm.dat")), tolerance= tol), "\n", file = "res.out", append=TRUE)
cat(all.equal((sw),read.ices(file.path(dir_sf,"sw.dat")), tolerance= tol), "\n", file = "res.out", append=TRUE)
cat(all.equal((lf),read.ices(file.path(dir_sf,"lf.dat")), tolerance= tol), "\n", file = "res.out", append=TRUE)
cat(all.equal(toNA(surveys),read.ices(file.path(dir_sf,"survey.dat")), tolerance= tol, check.attributes=FALSE), "\n", file = "res.out", append=TRUE)

##########################################################################################
## Multi fleet data w.o. sum fleet
##########################################################################################

cn<-read.ices("mf/cn.dat")
cw<-read.ices("mf/cw.dat")
dw<-cw
lf<-cn; lf[]<-1
lw<-cw
mo<-read.ices("mf/mo.dat")
nm<-read.ices("mf/nm.dat")
pf<-read.ices("mf/pf.dat")
pm<-read.ices("mf/pm.dat")
sw<-read.ices("mf/sw.dat")
surveys<-read.ices("mf/survey.dat")

cnA<-read.ices("mf/cn_A.dat")
cnC<-read.ices("mf/cn_C.dat")
cnD<-read.ices("mf/cn_D.dat")
cnF<-read.ices("mf/cn_F.dat")
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

dir <- tempdir()
dir.create(file.path(dir,"mf"))
dir_mf <- file.path(dir,"mf")

## Write fleet data to multiple files
write.data.files(dat,dir_mf, writeToOne = FALSE)

cat("Multi1","\n",file = "res.out",append=TRUE)
## CN
cat(all.equal(toNA(cnA),read.ices(file.path(dir_mf,"cn_00001.dat")), tolerance= tol, check.attributes=FALSE), "\n", file = "res.out", append = TRUE)
cat(all.equal(toNA(cnC),read.ices(file.path(dir_mf,"cn_00002.dat")), tolerance= tol, check.attributes=FALSE), "\n", file = "res.out", append = TRUE)
cat(all.equal(toNA(cnD),read.ices(file.path(dir_mf,"cn_00003.dat")), tolerance= tol, check.attributes=FALSE), "\n", file = "res.out", append = TRUE)
cat(all.equal(toNA(cnF),read.ices(file.path(dir_mf,"cn_00004.dat")), tolerance= tol, check.attributes=FALSE), "\n", file = "res.out", append = TRUE)
## CW
cat(all.equal((cw),read.ices(file.path(dir_mf,"cw_00001.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
cat(all.equal((cw),read.ices(file.path(dir_mf,"cw_00002.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
cat(all.equal((cw),read.ices(file.path(dir_mf,"cw_00003.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
cat(all.equal((cw),read.ices(file.path(dir_mf,"cw_00004.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)

## DW
cat(all.equal((dw),read.ices(file.path(dir_mf,"dw_00001.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
cat(all.equal((dw),read.ices(file.path(dir_mf,"dw_00002.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
cat(all.equal((dw),read.ices(file.path(dir_mf,"dw_00003.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
cat(all.equal((dw),read.ices(file.path(dir_mf,"dw_00004.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
## LW
cat(all.equal((lw),read.ices(file.path(dir_mf,"lw_00001.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
cat(all.equal((lw),read.ices(file.path(dir_mf,"lw_00002.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
cat(all.equal((lw),read.ices(file.path(dir_mf,"lw_00003.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
cat(all.equal((lw),read.ices(file.path(dir_mf,"lw_00004.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
## LF
cat(all.equal((lf),read.ices(file.path(dir_mf,"lf_00001.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
cat(all.equal((lf),read.ices(file.path(dir_mf,"lf_00002.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
cat(all.equal((lf),read.ices(file.path(dir_mf,"lf_00003.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
cat(all.equal((lf),read.ices(file.path(dir_mf,"lf_00004.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
## PF
cat(all.equal((pf),read.ices(file.path(dir_mf,"pf_00001.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
cat(all.equal((pf),read.ices(file.path(dir_mf,"pf_00002.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
cat(all.equal((pf),read.ices(file.path(dir_mf,"pf_00003.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
cat(all.equal((pf),read.ices(file.path(dir_mf,"pf_00004.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
## Stock data
cat(all.equal(mo,read.ices(file.path(dir_mf,"mo.dat")), tolerance= tol), "\n", file = "res.out", append=TRUE)
cat(all.equal(nm,read.ices(file.path(dir_mf,"nm.dat")), tolerance= tol), "\n", file = "res.out", append=TRUE)
cat(all.equal(pm,read.ices(file.path(dir_mf,"pm.dat")), tolerance= tol), "\n", file = "res.out", append=TRUE)
cat(all.equal((sw),read.ices(file.path(dir_mf,"sw.dat")), tolerance= tol), "\n", file = "res.out", append=TRUE)
cat(all.equal(toNA(surveys),read.ices(file.path(dir_mf,"survey.dat")), tolerance= tol, check.attributes=FALSE), "\n", file = "res.out", append=TRUE)


## Write fleet data to one files
dir.create(file.path(dir,"mf2"))
dir_mf2 <- file.path(dir,"mf2")
write.data.files(dat,dir_mf2, writeToOne = TRUE)

cat("Multi2","\n",file = "res.out",append=TRUE)
## CN
cat(all.equal(toNA(cnA),read.ices(file.path(dir_mf2,"cn_00001.dat")), tolerance= tol, check.attributes=FALSE), "\n", file = "res.out", append = TRUE)
cat(all.equal(toNA(cnC),read.ices(file.path(dir_mf2,"cn_00002.dat")), tolerance= tol, check.attributes=FALSE), "\n", file = "res.out", append = TRUE)
cat(all.equal(toNA(cnD),read.ices(file.path(dir_mf2,"cn_00003.dat")), tolerance= tol, check.attributes=FALSE), "\n", file = "res.out", append = TRUE)
cat(all.equal(toNA(cnF),read.ices(file.path(dir_mf2,"cn_00004.dat")), tolerance= tol, check.attributes=FALSE), "\n", file = "res.out", append = TRUE)
## CW
cat(all.equal((cw),read.ices(file.path(dir_mf2,"cw.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)

## DW
cat(all.equal((dw),read.ices(file.path(dir_mf2,"dw.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
## LW
cat(all.equal((lw),read.ices(file.path(dir_mf2,"lw.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
## LF
cat(all.equal((lf),read.ices(file.path(dir_mf2,"lf.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
## PF
cat(all.equal((pf),read.ices(file.path(dir_mf2,"pf.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
## Stock data
cat(all.equal(mo,read.ices(file.path(dir_mf2,"mo.dat")), tolerance= tol), "\n", file = "res.out", append=TRUE)
cat(all.equal(nm,read.ices(file.path(dir_mf2,"nm.dat")), tolerance= tol), "\n", file = "res.out", append=TRUE)
cat(all.equal(pm,read.ices(file.path(dir_mf2,"pm.dat")), tolerance= tol), "\n", file = "res.out", append=TRUE)
cat(all.equal((sw),read.ices(file.path(dir_mf2,"sw.dat")), tolerance= tol), "\n", file = "res.out", append=TRUE)
cat(all.equal(toNA(surveys),read.ices(file.path(dir_mf2,"survey.dat")), tolerance= tol, check.attributes=FALSE), "\n", file = "res.out", append=TRUE)



##########################################################################################
## Multi fleet data w. sum fleet
##########################################################################################


cn<-read.ices("mfs/cn.dat")
cw<-read.ices("mfs/cw.dat")
dw<-read.ices("mfs/dw.dat")
lw<-read.ices("mfs/lw.dat")
mo<-read.ices("mfs/mo.dat")
nm<-read.ices("mfs/nm.dat")
pf<-read.ices("mfs/pf.dat")
pm<-read.ices("mfs/pm.dat")
sw<-read.ices("mfs/sw.dat")
lf<-read.ices("mfs/lf.dat")
surveys<-read.ices("mfs/survey.dat")
cnA<-read.ices("mfs/cn_A.dat")
cnC<-read.ices("mfs/cn_C.dat")
cnD<-read.ices("mfs/cn_D.dat")
cnF<-read.ices("mfs/cn_F.dat")
cwA<-read.ices("mfs/cw_A.dat")
cwC<-read.ices("mfs/cw_C.dat")
cwD<-read.ices("mfs/cw_D.dat")
cwF<-read.ices("mfs/cw_F.dat")

idx<-!rownames(cn)%in%rownames(cnA)
sumACDF<-cn[idx,]
attr(sumACDF, "sumof")<-c(1,2,3,4)
                                        
cnA[cnA==0]<-1; cnA[,1]<-0
                                        
cnD[cnD==0]<-1
cnF[cnF==0]<-1

                                        
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


dir <- tempdir()
dir.create(file.path(dir,"mf3"))
dir_mf3 <- file.path(dir,"mf3")

## Write fleet data to one file
write.data.files(dat,dir_mf3, writeToOne = TRUE)

list.files(dir_mf3)


cat("Multi3","\n",file = "res.out",append=TRUE)
## CN
cat(all.equal(toNA(cnA),read.ices(file.path(dir_mf3,"cn_00001.dat")), tolerance= tol, check.attributes=FALSE), "\n", file = "res.out", append = TRUE)
cat(all.equal(toNA(cnC),read.ices(file.path(dir_mf3,"cn_00002.dat")), tolerance= tol, check.attributes=FALSE), "\n", file = "res.out", append = TRUE)
cat(all.equal(toNA(cnD),read.ices(file.path(dir_mf3,"cn_00003.dat")), tolerance= tol, check.attributes=FALSE), "\n", file = "res.out", append = TRUE)
cat(all.equal(toNA(cnF),read.ices(file.path(dir_mf3,"cn_00004.dat")), tolerance= tol, check.attributes=FALSE), "\n", file = "res.out", append = TRUE)
## CW
cat(all.equal((cwA),read.ices(file.path(dir_mf3,"cw_00001.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
cat(all.equal((cwC),read.ices(file.path(dir_mf3,"cw_00002.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
cat(all.equal((cwD),read.ices(file.path(dir_mf3,"cw_00003.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
cat(all.equal((cwF),read.ices(file.path(dir_mf3,"cw_00004.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
## CN sum
cat(all.equal(toNA(sumACDF),read.ices(file.path(dir_mf3,"cn_sum.dat")), tolerance= tol, check.attributes=FALSE), "\n", file = "res.out", append = TRUE)
## DW
cat(all.equal((dw),read.ices(file.path(dir_mf3,"dw.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
## LW
cat(all.equal((lw),read.ices(file.path(dir_mf3,"lw.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
## LF
cat(all.equal((lf),read.ices(file.path(dir_mf3,"lf.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
## PF
cat(all.equal((pf),read.ices(file.path(dir_mf3,"pf.dat")), tolerance= tol), "\n", file = "res.out", append = TRUE)
## Stock data
cat(all.equal(mo,read.ices(file.path(dir_mf3,"mo.dat")), tolerance= tol), "\n", file = "res.out", append=TRUE)
cat(all.equal(nm,read.ices(file.path(dir_mf3,"nm.dat")), tolerance= tol), "\n", file = "res.out", append=TRUE)
cat(all.equal(pm,read.ices(file.path(dir_mf3,"pm.dat")), tolerance= tol), "\n", file = "res.out", append=TRUE)
cat(all.equal((sw),read.ices(file.path(dir_mf3,"sw.dat")), tolerance= tol), "\n", file = "res.out", append=TRUE)
cat(all.equal(toNA(surveys),read.ices(file.path(dir_mf3,"survey.dat")), tolerance= tol, check.attributes=FALSE), "\n", file = "res.out", append=TRUE)
