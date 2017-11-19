##' Setup minimal initial parameters
##' @param dat sam data object
##' @param conf configuration list  
##' @details ...
##' @return a list containing the following
##' \item{logFpar}{}
##' \item{logQpow}{}
##' \item{logSdLogFsta}{}
##' \item{logSdLogN}{}
##' \item{logSdLogObs}{}
##' \item{logSdLogTotalObs}{}
##' \item{transfIRARdist}{}
##' \item{sigmaObsParUS}{}
##' \item{rec_loga}{}
##' \item{rec_logb}{}
##' \item{itrans_rho}{}
##' \item{logScale}{}
##' \item{logScaleSSB}{}
##' \item{logPowSSB}{}
##' \item{logSdSSB}{}
##' \item{logitReleaseSurvival}{}
##' \item{logitRecapturePhi}{}
##' \item{logF}{}
##' \item{logN}{}
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
  nbyfleet = (conf$obsCorStruct=="US")*(dat$maxAgePerFleet-dat$minAgePerFleet+1-(conf$obsLikelihoodFlag=="ALN"))
  ret$sigmaObsParUS=numeric(sum(nbyfleet*(nbyfleet-1)/2))
  ret$rec_loga=if(conf$stockRecruitmentModelCode==0){numeric(0)}else{numeric(1)}
  ret$rec_logb=if(conf$stockRecruitmentModelCode==0){numeric(0)}else{numeric(1)} 
  ret$itrans_rho=unlist(lapply(as.list(conf$corFlag),function(x){if(x==0){ ret <- numeric()} else { ret <- numeric(1)+.5}; return(ret)}))
  ret$logScale=if(conf$noScaledYears==0){numeric(0)}else{numeric(max(conf$keyParScaledYA)+1)}
  ret$logitReleaseSurvival=if(any(dat$fleetTypes==5)){numeric(length(unique(dat$aux[!is.na(dat$aux[,8]),8])))
                           }else{numeric(0)}
  ret$logitRecapturePhi=if(any(dat$fleetTypes==5)){numeric(length(ret$logitReleaseSurvival))
                        }else{numeric(0)}
  ret$logF=matrix(0, nrow=max(conf$keyLogFsta)+1,ncol=dat$noYears)
  ret$logN=matrix(0, nrow=conf$maxAge-conf$minAge+1, ncol=dat$noYears)
  return(ret)
}
