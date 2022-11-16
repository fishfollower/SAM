library(stockassessment)

cn<-read.ices("data/CANUM.txt")
cw<-read.ices("data/WECA.txt")
dw<-cw
lf<-cn; lf[]<-1
lw<-cw
mo<-read.ices("data/MATPROP.txt")
nm<-read.ices("data/NATMOR.txt")
pf<-read.ices("data/FPROP.txt")
pm<-read.ices("data/MPROP.txt")
sw<-read.ices("data/WEST.txt")
surveys<-read.ices("data/FLEET.txt")[1:3]

surveyLAI <- read.table("data/lai.txt",stringsAsFactors=F,header=T)
ORSH <- subset(surveyLAI,Area == "Or/Sh")
CNS  <- subset(surveyLAI,Area == "CNS")
BUN  <- subset(surveyLAI,Area == "Buchan")
SNS  <- subset(surveyLAI,Area == "SNS")

formatLAI <- function(x){
    x <- reshape(x[,c("Year","LAIUnit","L..9")], timevar="LAIUnit", idvar="Year",direction="wide")
    rownames(x) <- x$Year+1900
    x <- x[,-grep("Year",colnames(x))]
    colnames(x) <- 0:(ncol(x)-1)
    x <- as.matrix(x)
    attr(x,"time") <- c(0.67,0.67)
    return(x)}

ORSH <- formatLAI(ORSH); attr(ORSH,"part") <- 2
CNS <- formatLAI(CNS);   attr(CNS, "part") <- 1
BUN <- formatLAI(BUN);   attr(BUN, "part") <- 3
SNS <- formatLAI(SNS);   attr(SNS, "part") <- 4

surveys[[1]] <- surveys[[1]][,-ncol(surveys[[1]])]; attr(surveys[[1]],"time") <- c(0.54,0.56)
surveys[[2]] <- surveys[[2]][-nrow(surveys[[2]]),1:2]; attr(surveys[[2]],"time") <- c(0.08,0.17)
surveys[[3]] <- matrix(surveys[[3]][-nrow(surveys[[3]]),],ncol=1,dimnames=list(1992:2016,0)); attr(surveys[[3]],"time") <- c(0.08,0.17)
surveys[[4]] <- CNS
surveys[[5]] <- ORSH
surveys[[6]] <- BUN
surveys[[7]] <- SNS
names(surveys)[4:7] <- c("LAI_CNS","LAI_ORSH","LAI_BUN","LAI_SNS")

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
conf$keyLogFpar[2,-1]  <- c(1,1,2,2,2,3,3,3)-1
conf$keyLogFpar[3,2:3] <- c(4,4)-1
conf$keyLogFpar[4,1]   <- 5-1
conf$keyLogFpar[5:8,1] <- 6-1
#conf$keyVarObs[4:8,1]  <- 3
#conf$keyVarLogP <- rep(0,3)
par<-defpar(dat,conf)

par$logFpar<-c(0, 0.3, 0.4, -9, -12.5, -4.1)

suppressWarnings(fit<-sam.fit(dat,conf,par))

cat(fit$opt$objective,"\n\n", file="res.out")

cat(capture.output(prmatrix(fit$rep$comps)), sep="\n", file="res.out", append=TRUE)
