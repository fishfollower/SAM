##' Setup minimal initial parameters
##' @param dat sam data object
##' @param conf configuration list  
##' @details ...
##' @export
defpar <- function(dat,conf){
  ret<-list()
  ret$logFpar=numeric(max(conf$keyLogFpar)+1)-5 
  ret$logQpow=numeric(max(conf$keyQpow)+1)
  ret$logSdLogFsta=numeric(max(conf$keyVarF)+1)-.7
  ret$logSdLogN=numeric(max(conf$keyVarLogN)+1)-.35
  ret$logSdLogObs=numeric(max(conf$keyVarObs)+1)-.35
  ret$logSdLogTotalObs=numeric(sum(conf$obsLikelihoodFlag %in% c("ALN")))
  ret$transfIRARdist=if(all(is.na(conf$keyCorObs)))numeric(0) else numeric(max(conf$keyCorObs,na.rm=TRUE)+1)+0.05
  nbyfleet = (conf$obsCorStruct=="US")*(dat$maxAgePerFleet-dat$minAgePerFleet+1)
  ret$sigmaObsParUS=numeric(sum(nbyfleet*(nbyfleet-1)/2))
  ret$rec_loga=if(conf$stockRecruitmentModelCode==0){numeric(0)}else{numeric(1)}
  ret$rec_logb=if(conf$stockRecruitmentModelCode==0){numeric(0)}else{numeric(1)} 
  ret$itrans_rho=if(conf$corFlag==0){numeric(0)}else{numeric(1)+.5}
  ret$logScale=numeric(conf$noScaledYears)
  ret$logScaleSSB=if(any(dat$fleetTypes%in%c(3,4))){numeric(1)}else{numeric(0)}
  ret$logPowSSB=if(any(dat$fleetTypes==4)){numeric(1)}else{numeric(0)}
  ret$logSdSSB=if(any(dat$fleetTypes%in%c(3,4))){numeric(1)}else{numeric(0)}
  ret$logF=matrix(0, nrow=max(conf$keyLogFsta)+1,ncol=dat$noYears)
  ret$logN=matrix(0, nrow=conf$maxAge-conf$minAge+1, ncol=dat$noYears)
  return(ret)
}
