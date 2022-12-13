library(stockassessment)

## Only for illustration

obs<-read.table("data/obs.dat")
obs$Year <- floor(obs$t1)
obs$Q <- (obs$t1-obs$Year)*4+1

CN <- xtabs(obs ~ Year + ageFrom + Q, data = subset(obs, fleet == 1))
CN[dim(CN)[1],,4] <- NA
CN <- aperm(apply(CN, c(1,3), function(x) c(x[1:3],x[4]+x[5])),c(2,1,3))
CNq <- list(Q1 = local({x <- CN[,,1]; attr(x,"time") <- c(0,0.25); x[x==0] <- NA; x[,1] <- 0; x}),
            Q2 = local({x <- CN[,,2]; attr(x,"time") <- c(0.25,0.5); x[x==0] <- NA; x[,1] <- 0; x}),
            Q3 = local({x <- CN[,,3]; attr(x,"time") <- c(0.5,0.75); x[x==0] <- NA; x}),
            Q4 = local({x <- CN[,,4]; attr(x,"time") <- c(0.75,1); x[x==0] <- NA; x[nrow(x),] <- NA; x})
            )
CNyear <- apply(CN,c(1:2),sum)
CNp <- lapply(CNq, function(x){
    x[is.na(x)] <- 0
    x / ifelse(CNyear==0,1,CNyear)
})

obsS <- subset(obs, fleet > 2 & ageFrom < 4)
surveys <- lapply(split(obsS,obsS$fleet),function(dd){
    v <- as.matrix(xtabs(obs ~ Year + ageFrom, data = dd))
    v <- v[rowSums(v) > 0, colSums(v) > 0, drop = FALSE]
    attr(v,"time") <- c((dd$Q[1]-1),(dd$Q[1]))*0.25
    v
})


aux<-read.table("data/aux.dat")
aux$Year <- floor(aux$t1)
aux$Q <- (aux$t1-aux$Year)*4+1

xsa=read.table("data/xsa.dat")

## All the same
NM <- apply(xtabs(M ~ Year + age + Q, data = aux),1:2,mean)[,1:4]

## All the same
PM <- xtabs(PM ~ Year + age + Q, data = aux)
PM["2022",,4] <- PM["2022",,3]
PM <- apply(PM, 1:2, mean)[,1:4]

## Spawning beginning of Q3
SW <- xtabs(SW ~ Year + age + Q, data = aux)[,1:4,3]

CW <- xtabs(CW ~ Year + age + Q, data = aux)[,1:4,]

DW <- xtabs(DW ~ Year + age + Q, data = aux)[,1:4,]

LF <- xtabs(LF ~ Year + age + Q, data = aux)[,1:4,]
class(LF) <- "array"

propM <- SW*0 + 0.5
propF<- SW*0 + 0.5

## Fit yearly data

datY<-setup.sam.data(surveys=surveys,
                    residual.fleet=CNyear, 
                    prop.mature=PM, 
                    stock.mean.weight=SW, 
                    catch.mean.weight=Reduce("+",list(CW[,,1],CW[,,2],CW[,,3],CW[,,4]))/4, 
                    dis.mean.weight=Reduce("+",list(DW[,,1],DW[,,2],DW[,,3],DW[,,4]))/4, 
                    land.mean.weight=Reduce("+",list(CW[,,1],CW[,,2],CW[,,3],CW[,,4]))/4,
                    prop.f=propM, 
                    prop.m=propF, 
                    natural.mortality=NM, 
                    land.frac=Reduce("+",list(LF[,,1],LF[,,2],LF[,,3],LF[,,4]))/4)

confY <- defcon(datY)
confY$keyLogFmu <- confY$keyLogFsta
confY$keyLogFmu[1,] <- c(0,1,1,1)
confY$keyLogFrho[1,] <- 0
confY$corFlag[] <- 2
confY$keyVarObs[1,] <- c(0,1,1,2)
confY$keyVarObs[2,] <- c(-1,3,3,4)
confY$keyVarObs[3,] <- c(5,6,-1,-1)
confY$keyVarObs[4,] <- c(7,8,-1,-1)
confY$keyVarObs[5,] <- c(-1,-1,9,10)
confY$keyLogFpar[2,] <- c(-1,0,0,0)
confY$keyLogFpar[3,] <- c(1,2,-1,-1)
confY$keyLogFpar[4,] <- c(3,4,-1,-1)
confY$keyLogFpar[5,] <- c(-1,-1,5,5)


parY <- defpar(datY,confY)

fitY <- sam.fit(datY,confY,parY)
cat(AIC(fitY),"\n", file = "res.out")
## Fit quarterly fleets

datQ<-setup.sam.data(surveys=surveys,
                    residual.fleet=CNq, 
                    prop.mature=PM, 
                    stock.mean.weight=SW, 
                    catch.mean.weight=list(CW[,,1],CW[,,2],CW[,,3],CW[,,4]), 
                    dis.mean.weight=list(DW[,,1],DW[,,2],DW[,,3],DW[,,4]), 
                    land.mean.weight=list(CW[,,1],CW[,,2],CW[,,3],CW[,,4]),
                    prop.f=list(propF*0+1,propF*0+1,propF*0,propF*0), 
                    prop.m=propM, 
                    natural.mortality=NM, 
                    land.frac=list(LF[,,1],LF[,,2],LF[,,3],LF[,,4]))

confQ <- defcon(datQ)
## confQ$keyLogFmu <- confQ$keyLogFsta
## confQ$keyLogFmu[1,] <- c(-1,0,0,0)
## confQ$keyLogFmu[2,] <- c(-1,1,1,1)
## confQ$keyLogFmu[3,] <- c(2,3,3,3)
## confQ$keyLogFmu[4,] <- c(4,5,5,5)
## confQ$keyLogFrho[1,] <- 0
## confQ$keyLogFrho[2,] <- 1
## confQ$keyLogFrho[3,] <- 2
## confQ$keyLogFrho[4,] <- 3
confQ$keyVarF[1:4,] <- 0
#confQ$corFlag[] <- c(4,4,4,4)
confQ$keyVarObs[1,] <- c(-1,0,0,1)
confQ$keyVarObs[2,] <- c(-1,0,0,1)
confQ$keyVarObs[3,] <- c(2,0,0,3)
confQ$keyVarObs[4,] <- c(4,0,0,5)
confQ$keyVarObs[5,] <- c(-1,6,6,6)
confQ$keyVarObs[6,] <- c(7,8,-1,-1)
confQ$keyVarObs[7,] <- c(9,10,-1,-1)
confQ$keyVarObs[8,] <- c(-1,-1,11,12)
confQ$keyLogFpar[5,] <- c(-1,0,0,0)
confQ$keyLogFpar[6,] <- c(1,0,-1,-1)
confQ$keyLogFpar[7,] <- c(2,0,-1,-1)
confQ$keyLogFpar[8,] <- c(-1,-1,0,0)

parQ <- defpar(datQ,confQ)

fitQ <- sam.fit(datQ,confQ,parQ, run = TRUE,
               lower = list(itrans_rho = -5),
               upper = list(itrans_rho = 5))

cat(AIC(fitQ),"\n", file = "res.out", append=TRUE)

