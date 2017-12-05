# to install the package from the multi branch 
#devtools::install_github("fishfollower/SAM/stockassessment", ref="multi")

library(stockassessment)
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
surveys<-read.ices("old/FLEET.txt")


cnA<-read.ices("split/cn_A.dat");
cnB<-read.ices("split/cn_B.dat");
cnC<-read.ices("split/cn_C.dat");
cnD<-read.ices("split/cn_D.dat");
cwA<-read.ices("split/cw_A.dat");
cwB<-read.ices("split/cw_B.dat");
cwC<-read.ices("split/cw_C.dat");
cwD<-read.ices("split/cw_D.dat");

cwA[,as.character(8)] <- rowSums(cnA[,as.character(8:9)] * cwA[,as.character(8:9)])/ rowSums(cnA[,as.character(8:9)])
cwB[,as.character(8)] <- rowSums(cnB[,as.character(8:9)] * cwB[,as.character(8:9)])/ rowSums(cnB[,as.character(8:9)])
cwC[,as.character(8)] <- rowSums(cnC[,as.character(8:9)] * cwC[,as.character(8:9)])/ rowSums(cnC[,as.character(8:9)])
cwD[,as.character(8)] <- rowSums(cnD[,as.character(8:9)] * cwD[,as.character(8:9)])/ rowSums(cnD[,as.character(8:9)])

cnA[,as.character(8)] <- rowSums(cnA[,as.character(8:9)]);
cnB[,as.character(8)] <- rowSums(cnB[,as.character(8:9)]);
cnC[,as.character(8)] <- rowSums(cnC[,as.character(8:9)]);
cnD[,as.character(8)] <- rowSums(cnD[,as.character(8:9)]);

cnA <- cnA[,-ncol(cnA)]
cnB <- cnB[,-ncol(cnB)]
cnC <- cnC[,-ncol(cnC)]
cnD <- cnD[,-ncol(cnD)]

cwA <- cwA[,-ncol(cwA)]
cwB <- cwB[,-ncol(cwB)]
cwC <- cwC[,-ncol(cwC)]
cwD <- cwD[,-ncol(cwD)]

idx<-!rownames(cn)%in%rownames(cnA)

cwAF <- cwBF <- cwCF <- cwDF <- cw
cwAF[-which(idx),] <- cwA; cwAF[is.nan(cwAF)] <- 0
cwBF[-which(idx),] <- cwB; cwBF[is.nan(cwBF)] <- 0
cwCF[-which(idx),] <- cwC; cwCF[is.nan(cwCF)] <- 0
cwDF[-which(idx),] <- cwD; cwDF[is.nan(cwDF)] <- 0

#sumACDF<-(cnA+cnC+cnD+cnF)[idx,]
sumABCD<-cn[idx,]

attr(sumABCD, "sumof")<-c(1,2,3,4)

# temp fix to convergence problems 
cnA[cnA==0]<-1
cnB[cnB==0]<-1
cnC[cnC==0]<-1
cnD[cnD==0]<-1

surveys[[1]] <- surveys[[1]][,-ncol(surveys[[1]])]; attr(surveys[[1]],"time") <- c(0.54,0.56)
surveys[[2]] <- surveys[[2]][-nrow(surveys[[2]]),1:2]; attr(surveys[[2]],"time") <- c(0.08,0.17)
surveys[[3]] <- matrix(surveys[[3]][-nrow(surveys[[3]]),],ncol=1,dimnames=list(1992:2016,0)); attr(surveys[[3]],"time") <- c(0.08,0.17)

dat<-setup.sam.data(surveys=surveys[1:3],
                    residual.fleets=list(cnA, cnB, cnC, cnD), # Notice list
                    sum.residual.fleets=sumABCD,
                    prop.mature=mo, 
                    stock.mean.weight=sw, 
                    catch.mean.weight=cw,#list(cwAF,cwBF,cwCF,cwDF),
                    dis.mean.weight=dw, 
                    land.mean.weight=lw,#list(cwAF,cwBF,cwCF,cwDF),
                    prop.f=pf, 
                    prop.m=pm, 
                    natural.mortality=nm, 
                    land.frac=lf)
conf<-defcon(dat)
conf$keyLogFsta[2,] <- c(8,9,10,rep(11,6))
conf$keyLogFsta[3,] <- c(12,13,14,rep(15,6))
conf$keyLogFsta[4,] <- c(16,17,18,rep(19,6))

par<-defpar(dat,conf)
fit<-sam.fit(dat,conf,par)

cat(fit$opt$objective,"\n\n", file="res.out")
cat(capture.output(prmatrix(t(fit$pl$logF))), sep="\n", file="res.out", append=TRUE)
