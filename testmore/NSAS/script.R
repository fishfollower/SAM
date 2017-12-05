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

#- Original split data
cnA<-read.ices("split/cn_A.dat");
cnB<-read.ices("split/cn_B.dat");
cnC<-read.ices("split/cn_C.dat");
cnD<-read.ices("split/cn_D.dat");
cwA<-read.ices("split/cw_A.dat");
cwB<-read.ices("split/cw_B.dat");
cwC<-read.ices("split/cw_C.dat");
cwD<-read.ices("split/cw_D.dat");

#- Set plusgroup to 8+ (data in 9+)
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

#- Fake split
cnAFS <- rbind(colMeans(cnA / (cnA + cnB + cnC + cnD)) * cn[1:50,],cnA)
cnBFS <- rbind(colMeans(cnB / (cnA + cnB + cnC + cnD)) * cn[1:50,],cnB)
cnCFS <- rbind(colMeans(cnC / (cnA + cnB + cnC + cnD)) * cn[1:50,],cnC)
cnDFS <- rbind(colMeans(cnD / (cnA + cnB + cnC + cnD)) * cn[1:50,],cnD)
cnAFS[cnAFS==0]<-1
cnBFS[cnBFS==0]<-1
cnCFS[cnCFS==0]<-1
cnDFS[cnDFS==0]<-1

idx<-!rownames(cn)%in%rownames(cnA)

#- get catch weight by fleet full length
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
sumABCD[sumABCD==0] <- 1

cwAF[which(cwAF==0)] <- 0.01
cwBF[which(cwBF==0)] <- 0.01
cwCF[which(cwCF==0)] <- 0.01
cwDF[which(cwDF==0)] <- 0.01

surveys[[1]] <- surveys[[1]][,-ncol(surveys[[1]])]; attr(surveys[[1]],"time") <- c(0.54,0.56)
surveys[[2]] <- surveys[[2]][-nrow(surveys[[2]]),1:2]; attr(surveys[[2]],"time") <- c(0.08,0.17)
surveys[[3]] <- matrix(surveys[[3]][-nrow(surveys[[3]]),],ncol=1,dimnames=list(1992:2016,0)); attr(surveys[[3]],"time") <- c(0.08,0.17)

iSY         <- 1
dat<-setup.sam.data(surveys=surveys[1:3],
                    #residual.fleets=list(cn), # Notice list
                    #residual.fleets=list(cnAFS, cnBFS, cnCFS, cnDFS), # Notice list
                    residual.fleets=list(cnA, cnB, cnC, cnD), # Notice list
                    sum.residual.fleets=sumABCD[iSY:50,],
                    prop.mature=mo[iSY:70,],
                    stock.mean.weight=sw[iSY:70,],
                    catch.mean.weight=cw[iSY:70,],#list(cwAF,cwBF,cwCF,cwDF),
                    dis.mean.weight=dw[iSY:70,],
                    land.mean.weight=lw[iSY:70,],#list(cwAF,cwBF,cwCF,cwDF),
                    prop.f=pf[iSY:70,],
                    prop.m=pm[iSY:70,],
                    natural.mortality=nm[iSY:70,],
                    land.frac=lf[iSY:70,])
conf<-defcon(dat)
# Potential improvements to get the model to fit
conf$keyLogFsta[1,] <- c(0:5,6,6,6)
conf$keyLogFsta[2,] <- c(8,9,10,rep(11,6))-1
conf$keyLogFsta[3,] <- c(12,13,14,15,16,rep(17,4))-1
conf$keyLogFsta[4,] <- c(18,19,20,rep(21,6))-1

conf$keyVarObs[1,] <- c(rep(0,4),rep(1,5))
conf$keyVarObs[2,] <- c(rep(2,4),rep(3,5))
conf$keyVarObs[3,] <- c(rep(4,4),rep(5,5))
conf$keyVarObs[4,] <- c(rep(6,4),rep(7,5))
conf$keyVarObs[5,-1] <- conf$keyVarObs[5,-1]+4
conf$keyVarObs[6,2:3] <- conf$keyVarObs[6,2:3]+4
conf$keyVarObs[7,1] <- conf$keyVarObs[7,1]+4

conf$keyVarF[1,] <- c(rep(0,4),rep(1,5))
conf$keyVarF[2,] <- c(rep(2,3),rep(3,6))
conf$keyVarF[3,] <- c(rep(4,3),rep(5,6))
conf$keyVarF[4,] <- c(rep(6,3),rep(7,6))

conf$corFlag <- c(2,0,0,0)

#
conf$fbarRange   <- c(2,6)
par<-defpar(dat,conf)
fit<-sam.fit(dat,conf,par)

save(fit,file="./truesplit.RData")
save(fit,file="./nosplit.RData")
save(fit,file="./fakesplit.RData")

catchplot(fit)
ssbplot(fit)
fbarplot(fit)

cat(fit$opt$objective,"\n\n", file="res.out")
cat(capture.output(prmatrix(t(fit$pl$logF))), sep="\n", file="res.out", append=TRUE)
