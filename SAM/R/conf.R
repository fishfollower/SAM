##' Setup basic minimal configuration for sam assessment
##' @param dat sam data object 
##' @details ...
##' @export
defcon<-function(dat){
  setSeq<-function(min,max){
    if(min==max){
      ret<-1
    }else{
      ret<-c(1:(max-min),max-min)
    }
    return(ret)
  }
  fleetTypes<-dat$fleetTypes
  ages<-do.call(rbind,tapply(dat$obs[,3], INDEX=dat$obs[,2], FUN=range))
  minAge<-min(ages)
  maxAge<-max(ages)
  nAges<-maxAge-minAge+1
  nFleets<-nrow(ages)
  ret<-list()
  ret$minAge <- minAge
  ret$maxAge <- maxAge  
  ret$maxAgePlusGroup <- 1
  x<-matrix(0, nrow=nFleets, ncol=nAges)
  lastMax<-0
  for(i in 1:nrow(x)){
    if(fleetTypes[i]==0){
      x[i,(ages[i,1]-minAge+1):(ages[i,2]-minAge+1)]<-setSeq(ages[i,1],ages[i,2])+lastMax
      lastMax<-max(x)
    }
  }  
  ret$keyLogFsta <- x
  ret$corFlag <- 0
  x<-matrix(0, nrow=nFleets, ncol=nAges)
  lastMax<-0
  for(i in 1:nrow(x)){
    if(fleetTypes[i]%in%c(1,2)){
      x[i,(ages[i,1]-minAge+1):(ages[i,2]-minAge+1)]<-setSeq(ages[i,1],ages[i,2])+lastMax
      lastMax<-max(x)
    }
  }  
  ret$keyLogFpar <- x 
  ret$keyQpow <- matrix(0, nrow=nFleets, ncol=nAges)
  x<-matrix(0, nrow=nFleets, ncol=nAges)
  lastMax<-0
  for(i in 1:nrow(x)){
    if(fleetTypes[i]==0){
      x[i,(ages[i,1]-minAge+1):(ages[i,2]-minAge+1)]<-lastMax+1
      lastMax<-max(x)
    }
  }  
  ret$keyVarF <- x
  ret$keyVarLogN <- c(1,rep(2,nAges-1))
  x<-matrix(0, nrow=nFleets, ncol=nAges)
  lastMax<-0
  for(i in 1:nrow(x)){
    if(fleetTypes[i]%in%c(0,1,2)){
      x[i,(ages[i,1]-minAge+1):(ages[i,2]-minAge+1)]<-lastMax+1
      lastMax<-max(x)
    }
  }  
  ret$keyVarObs <- x
  ret$stockRecruitmentModelCode <- 0
  ret$noScaledYears <- 0
  ret$keyScaledYears <- numeric(0)
  ret$keyParScaledYA <- numeric(0)
  ret$fbarRange <- ages[fleetTypes==0,]

  return(ret)
}
