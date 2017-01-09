##' small helper function
##' @param min from 
##' @param max to 
##' @details ...
setSeq<-function(min,max){
  if(min==max){
    ret <- 1
  }else{
    ret <- c(1:(max-min),max-min)
  }
  return(ret)
}

##' Setup basic minimal configuration for sam assessment
##' @param dat sam data object 
##' @return a list containing 
##' \item{minAge}{}
##' \item{maxAge}{}
##' \item{maxAgePlusGroup}{}
##' \item{keyLogFsta}{}
##' \item{keyLogFpar}{}
##' \item{keyQpow}{}
##' \item{keyVarF}{}
##' \item{keyVarLogN}{}
##' \item{keyVarObs}{}
##' \item{obsCorStruct}{}
##' \item{keyCorObs}{}
##' \item{corFlag}{set to zero as a placeholder here.}
##' \item{stockRecruitmentModelCode}{set to zero as a placeholder here.}
##' \item{noScaledYears}{set to zero as a placeholder here.}
##' \item{keyScaledYears}{a scalar set to zero as a placeholder here.}
##' \item{keyParScaledYA}{an array set to zero as a placeholder here.}
##' \item{fbarRange}{}
##' \item{obsLikelihoodFlag}{A vector with an element for each fleet. "LN" means, and "ALN" means}
##' @details ...
##' @export
defcon<-function(dat){
  fleetTypes <- dat$fleetTypes
  ages <- do.call(rbind,tapply(dat$aux[,3], INDEX=dat$aux[,2], FUN=range))
  ages[fleetTypes%in%c(3,5),] <- NA
  minAge <- min(ages, na.rm=TRUE)
  maxAge <- max(ages, na.rm=TRUE)
  ages[is.na(ages)] <- minAge
  nAges <- maxAge-minAge+1
  nFleets <- nrow(ages)
  ret <- list()
  ret$minAge <- minAge
  ret$maxAge <- maxAge  
  ret$maxAgePlusGroup <- 1
  x <- matrix(0, nrow=nFleets, ncol=nAges)
  lastMax <- 0
  for(i in 1:nrow(x)){
    if(fleetTypes[i]==0){
      x[i,(ages[i,1]-minAge+1):(ages[i,2]-minAge+1)] <- setSeq(ages[i,1],ages[i,2])+lastMax
      lastMax <- max(x)
    }
  }  
  ret$keyLogFsta <- x - 1
  ret$corFlag <- 0
  x <- matrix(0, nrow=nFleets, ncol=nAges)
  lastMax <- 0
  for(i in 1:nrow(x)){
    if(fleetTypes[i]%in%c(1,2,3)){
      x[i,(ages[i,1]-minAge+1):(ages[i,2]-minAge+1)] <- setSeq(ages[i,1],ages[i,2])+lastMax
      lastMax <- max(x)
    }
  }  
  ret$keyLogFpar <- x - 1 
  ret$keyQpow <- matrix(-1, nrow=nFleets, ncol=nAges)
  x<-matrix(0, nrow=nFleets, ncol=nAges)
  lastMax <- 0
  for(i in 1:nrow(x)){
    if(fleetTypes[i]==0){
      x[i,(ages[i,1]-minAge+1):(ages[i,2]-minAge+1)]<-lastMax+1
      lastMax <- max(x)
    }
  }  
  ret$keyVarF <- x - 1
  ret$keyVarLogN <- c(1,rep(2,nAges-1)) - 1 
  x <- matrix(0, nrow=nFleets, ncol=nAges)
  lastMax <- 0
  for(i in 1:nrow(x)){
    if(fleetTypes[i]%in%c(0,1,2,3)){
      x[i,(ages[i,1]-minAge+1):(ages[i,2]-minAge+1)] <- lastMax+1
      lastMax <- max(x)
    }
  }  
  ret$keyVarObs <- x - 1
  ret$obsCorStruct <- factor(rep("ID",nFleets),levels=c("ID","AR","US"))
  ret$keyCorObs <- matrix(-1, nrow=nFleets, ncol=nAges-1)
  colnames(ret$keyCorObs)<-paste(minAge:(maxAge-1),(minAge+1):maxAge,sep="-")
  for(i in 1:nrow(x)){
      ret$keyCorObs[i,(ages[i,1]-minAge+1):(ages[i,2]-minAge)]<-NA
  }
    
  ret$stockRecruitmentModelCode <- 0
  ret$noScaledYears <- 0
  ret$keyScaledYears <- numeric(0)
  ret$keyParScaledYA <- array(0,c(0,0))
  ret$fbarRange <- ages[fleetTypes==0,]

  ret$obsLikelihoodFlag <- factor(rep("LN",nFleets),levels=c("LN","ALN"))
  ret$cutReleaseSurvival <- integer(0)
  ret$fixVarToWeight <- 0
  return(ret) 
}
