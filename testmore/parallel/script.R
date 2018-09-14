## Tests runwithout using maps and bounds
library(stockassessment)

prefix <- "WHITEHAKE_"

cn <- read.ices(paste0(prefix,"cn.dat"))
cw <- read.ices(paste0(prefix,"cw.dat"))
dw <- read.ices(paste0(prefix,"dw.dat"))
lw <- read.ices(paste0(prefix,"lw.dat"))
mo <- read.ices(paste0(prefix,"mo.dat"))
nm <- read.ices(paste0(prefix,"nm.dat"))
pf <- read.ices(paste0(prefix,"pf.dat"))
pm <- read.ices(paste0(prefix,"pm.dat"))
sw <- read.ices(paste0(prefix,"sw.dat"))
lf <- read.ices(paste0(prefix,"lf.dat"))
surveys <- read.ices(paste0(prefix,"survey.dat"))

dat <- setup.sam.data(surveys=surveys,
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

conf <- loadConf(dat,"model.cfg")

par <- defpar(dat,conf)

rhoF = 0.99999
par$itrans_rho = -log(2/(1+rhoF)-1)/2
lower = list(logSdLogN = c(-4,-4))
fit <- sam.fit(dat,conf,par, map = list(itrans_rho = as.factor(NA)),lower=lower)

RETRO <- retro(fit, year=2, newtonsteps=0) 

res <- c(  RETRO[[1]]$pl$itrans_rho - par$itrans_rho, ## should be zero if map is working 
           RETRO[[1]]$pl$logSdLogN) ## should not be lower than -4
cat(round(res,4),"\n", file="res.out")



