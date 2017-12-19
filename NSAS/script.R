# to install the package from the multi branch 
#devtools::install_github("fishfollower/SAM/stockassessment", ref="multi")
library(TMB)
setwd("D:/Repository/SAM/SAMMULTILAI/NSAS/")
#compile("../stockassessment/src/stockassessment.cpp")
library(stockassessment)
library(reshape2)
cn<-read.ices("old/CANUM.txt")
cw<-read.ices("old/WECA.txt")
dw<-cw
lf<-cn; lf[]<-1
lw<-cw
mo<-read.ices("old/MATPROP.txt")
nm<-read.ices("old/NATMOR.txt")
pf<-read.ices("old/FPROP.txt")
pm<-read.ices("old/MPROP.txt")
sw<-read.ices("old/WEST.txt")
surveys<-read.ices("old/FLEET.txt")[1:3]
surveyLAI <- read.table("old/lai.txt",stringsAsFactors=F,header=T)
ORSH <- subset(surveyLAI,Area == "Or/Sh")
CNS  <- subset(surveyLAI,Area == "CNS")
BUN  <- subset(surveyLAI,Area == "Buchan")
SNS  <- subset(surveyLAI,Area == "SNS")

formatLAI <- function(x){
              x <- dcast(x[,c("Year","LAIUnit","L..9")],Year ~ LAIUnit,value.var="L..9")
              rownames(x) <- x$Year+1900
              x <- x[,-grep("Year",colnames(x))]
              colnames(x) <- 0:(ncol(x)-1)
              x <- as.matrix(x)
              attr(x,"time") <- c(0.67,0.67)
             return(x)}

ORSH <- formatLAI(ORSH); attr(ORSH,"part") <- 1
CNS <- formatLAI(CNS);   attr(CNS, "part") <- 2
BUN <- formatLAI(BUN);   attr(BUN, "part") <- 3
SNS <- formatLAI(SNS);   attr(SNS, "part") <- 4

surveys[[1]] <- surveys[[1]][,-ncol(surveys[[1]])]; attr(surveys[[1]],"time") <- c(0.54,0.56)
surveys[[2]] <- surveys[[2]][-nrow(surveys[[2]]),1:2]; attr(surveys[[2]],"time") <- c(0.08,0.17)
surveys[[3]] <- matrix(surveys[[3]][-nrow(surveys[[3]]),],ncol=1,dimnames=list(1992:2016,0)); attr(surveys[[3]],"time") <- c(0.08,0.17)
surveys[[4]] <- ORSH
surveys[[5]] <- CNS
surveys[[6]] <- BUN
surveys[[7]] <- SNS
names(surveys)[4:7] <- c("LAI_ORSH","LAI_CNS","LAI_BUN","LAI_SNS")

source("../stockassessment/R/reading.R")
source("../stockassessment/R/conf.R")
source("../stockassessment/R/init.R")
source("../stockassessment/R/run.R")



dat<-setup.sam.data(surveys=surveys,
                    residual.fleets=list(cn), # Notice list
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
save.image("./alldat.RData")
load("./alldat.RData")

dyn.load("../stockassessment/src/stockassessment.dll")

  data <- dat
  parameters <- par
  data<-clean.void.catches(dat,conf)
  tmball <- c(data, conf, simFlag=T)
tmball$resFlag <- 0
  nmissing <- sum(is.na(dat$logobs))
  parameters$missing <- numeric(nmissing)
  ran <- c("logN", "logF","logP", "missing")
  random <- ran
  data <- tmball
  DLL <- "stockassessment"
#  obj <- MakeADFun(tmball, parameters, random=ran, DLL="stockassessment")


#fit<-sam.fit(dat,conf,par)

#cat(fit$opt$objective,"\n\n", file="res.out")
#cat(capture.output(prmatrix(t(fit$pl$logF))), sep="\n", file="res.out", append=TRUE)

