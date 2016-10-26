# SAM
Git page for SAM R-package

Can be installed by typing: 

```R
devtools::install_github("fishfollower/SAM/stockassessment")
```

The following is a quick tour og the package, but let's first download a few needed files (make sure you are not overwriting existing files). 

```R
setwd(tempdir())
filestoget <- c("cn.dat", "cw.dat", "dw.dat", "lf.dat", "lw.dat", 
                "mo.dat", "nm.dat", "pf.dat", "pm.dat", "sw.dat", 
                "survey.dat")
url<-"https://raw.githubusercontent.com/fishfollower/SAM/master/tests/nsher/"
d<-lapply(filestoget, function(f)download.file(paste(url,f,sep=""), f))
```
Now the files should be downloaded to our current working directory. Next we read in the data files and create a collected data object: 

```R
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
```

From this defined data object it is possible to generate a default/minimalistic model configuration.

```R
conf<-defcon(dat)
```

This configuration can be changed by modifying the elements in the list. Here we set the fbar-range to 2-4, allow correlated F processes, and modify the catchability couplings. 

```R
conf$fbarRange<-c(2,6)
conf$corFlag<-1
conf$keyLogFpar<-matrix(c(
-1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
-1,    0,    1,    2,    3,    4,    5,    6,   -1,
-1,    7,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
 8,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1), nrow=4, byrow=TRUE)
``` 

 