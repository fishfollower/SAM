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

conf<-loadConf(dat,"model.cfg")
par<-defpar(dat,conf)
fit<-sam.fit(dat,conf,par)


conf2<-loadConf(dat,"model2.cfg")
par2<-defpar(dat,conf2)
fit2<-sam.fit(dat,conf2,par2)

#par2$logXtraSd<-log(2)
#fit2<-sam.fit(dat,conf2,par2, map=list(logXtraSd=factor(NA)))
#set.seed(123)

#cn3<-cn
#cn3[rownames(cn3)%in%1993:2005,]<-exp(log(cn3[rownames(cn3)%in%1993:2005,])+rnorm(length(cn3[rownames(cn3)%in%1993:2005,]),0,.2))
#dat3<-setup.sam.data(surveys=surveys,
#                    residual.fleet=cn3, 
#                    prop.mature=mo, 
#                    stock.mean.weight=sw, 
#                    catch.mean.weight=cw, 
#                    dis.mean.weight=dw, 
#                    land.mean.weight=lw,
#                    prop.f=pf, 
#                    prop.m=pm, 
#                    natural.mortality=nm, 
#                    land.frac=lf)
#
#conf3<-loadConf(dat3,"model2.cfg")
#par3<-defpar(dat3,conf3)
#fit3<-sam.fit(dat3,conf3,par3)

cat(capture.output(modeltable(c("Standard"=fit,"XtraSd"=fit2))), file = 'res.out', sep = '\n')
