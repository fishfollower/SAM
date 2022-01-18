##' Setup initial values for all model parameters and random effects.
##' @param dat sam data object as returned from the function \code{setup.sam.data}
##' @param conf sam configuration list, which could be read from a configuration file via the \code{loadConf} function. A default/dummy configuration can be generated via the \code{defcon} function.
##' @param spinoutyear Technical setting only used for biological parameter process models to insure equilibrium distribution in final edge year 
##' @details The model parameters and random effects are not initialized in any clever way - most are simply set to zero. If convergence problems occour different initial values can be tested, but it is more likely a problem with the model configuration. 
##' @return a list containing initial values for all model parameters and random effects in the model.
##' @export
defpar <- function(dat,conf,spinoutyear=10){
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



  M <- max(dat$natMor)
  F <- M
  C <- unname(quantile(dat$logobs[dat$aux[,"fleet"] %in% which(dat$fleetTypes == 0)],0.95, na.rm = TRUE))
  if(conf$stockRecruitmentModelCode==0){ # Random walk
      ret$rec_pars <- numeric(0)
  }else if(conf$stockRecruitmentModelCode==3){ # Constant mean
      ret$rec_pars <- numeric(length(unique(conf$constRecBreaks))+1)      
  }else if(conf$stockRecruitmentModelCode==50){ # depensatory logistic Hockey stick
      ret$rec_pars <- c(log(1),C + log(F+M) - log(F) -log(1 - exp(-F -M)),0,0)
      if(!is.na(conf$hockeyStickCurve))
          ret$rec_pars[3] <- log(conf$hockeyStickCurve)   
      ret$rec_pars[4] <- C
  }else if(conf$stockRecruitmentModelCode %in% c(51,52)){ # Shepherd, Deriso
      ret$rec_pars <- numeric(3)
      ret$rec_pars[3] <- C
  }else if(conf$stockRecruitmentModelCode==60){ # logistic Hockey stick
      ret$rec_pars <- c(log(1),C + log(F+M) - log(F) -log(1 - exp(-F -M)),0)
      if(!is.na(conf$hockeyStickCurve))
          ret$rec_pars[3] <- log(conf$hockeyStickCurve)
  }else if(conf$stockRecruitmentModelCode==61){ # Hockey stick
      ret$rec_pars <- c(C + log(F+M) - log(F) -log(1 - exp(-F -M))  + log(2),C + log(F+M) - log(F) -log(1 - exp(-F -M)))
  }else if(conf$stockRecruitmentModelCode==63){ # Bent hypoerbola / Hockey-stick-like
      ret$rec_pars <- c(C + log(F+M) - log(F) -log(1 - exp(-F -M)) ,log(0.5),3)
      if(!is.na(conf$hockeyStickCurve))
          ret$rec_pars[3] <- log(conf$hockeyStickCurve)
  }else if(conf$stockRecruitmentModelCode %in% c(66,67,68,69)){ # Shepherd, Deriso
      ret$rec_pars <- numeric(3)
  }else if(conf$stockRecruitmentModelCode %in% c(90,91,92)){
      ret$rec_pars <- numeric(length(conf$constRecBreaks) + as.integer(conf$stockRecruitmentModelCode==90))
      if(conf$stockRecruitmentModelCode==90)
          ret$rec_pars[length(ret$rec_pars)] <- 5
  }else{ # The rest
      ret$rec_pars <- numeric(2)
  }
  
  ret$itrans_rho=if(conf$corFlag==0){numeric(0)}else{numeric(1)+.5}
  ret$logScale=if(conf$noScaledYears==0){numeric(0)}else{numeric(max(conf$keyParScaledYA)+1)}
  ret$logitReleaseSurvival=if(any(dat$fleetTypes==5)){numeric(length(unique(dat$aux[!is.na(dat$aux[,8]),8])))
                           }else{numeric(0)}
  ret$logitRecapturePhi=if(any(dat$fleetTypes==5)){numeric(length(ret$logitReleaseSurvival))
                        }else{numeric(0)}
  
  if(conf$corFlag ==3 ){
    ret$sepFalpha=rep(0,max(conf$keyLogFsta)+1)
    ret$sepFlogitRho = rep(1,2)
    ret$sepFlogSd = rep(-1,2)
    ret$logSdLogFsta = numeric(0)
    ret$itrans_rho = numeric(0)

    if(conf$corFlag==3){
      ret$logSdLogFsta = numeric(0)
    }
  }else{
    ret$sepFalpha=numeric(0)
    ret$sepFlogitRho = numeric(0)
    ret$sepFlogSd = numeric(0)
    
  }
  
  ret$predVarObs =if(length(conf$predVarObsLink)!=0){
    if(max(conf$predVarObsLink,na.rm=TRUE)>(-1)){
      rep(log(1),(max(conf$predVarObsLink,na.rm=TRUE) + 1))
    }else{numeric(0)}
  }else{numeric(0)}

  ret$logPhiSW=if(conf$stockWeightModel==0){numeric(0)}else{numeric(2)}
  ret$logSdProcLogSW=if(conf$stockWeightModel==0){numeric(0)}else{numeric(1)}
  ret$meanLogSW=if(conf$stockWeightModel==0){numeric(0)}else{numeric(max(conf$keyStockWeightMean,na.rm=TRUE)+1)}
  ret$logSdLogSW=if(conf$stockWeightModel==0){numeric(0)}else{numeric(max(conf$keyStockWeightObsVar,na.rm=TRUE)+1)}

  ret$logPhiCW=if(conf$catchWeightModel==0){numeric(0)}else{numeric(2)}
  ret$logSdProcLogCW=if(conf$catchWeightModel==0){numeric(0)}else{numeric(1)}
  ret$meanLogCW=if(conf$catchWeightModel==0){numeric(0)}else{numeric(max(conf$keyCatchWeightMean,na.rm=TRUE)+1)}
  ret$logSdLogCW=if(conf$catchWeightModel==0){numeric(0)}else{numeric(max(conf$keyCatchWeightObsVar,na.rm=TRUE)+1)}

  ret$logPhiMO=if(conf$matureModel==0){numeric(0)}else{numeric(2)}
  ret$logSdProcLogitMO=if(conf$matureModel==0){numeric(0)}else{numeric(1)}
  ret$meanLogitMO=if(conf$matureModel==0){numeric(0)}else{numeric(max(conf$keyMatureMean,na.rm=TRUE)+1)}
  ret$logSdMO=if(conf$matureModel==0){numeric(0)}else{numeric(1)}

  ret$logPhiNM=if(conf$mortalityModel==0){numeric(0)}else{numeric(2)}
  ret$logSdProcLogNM=if(conf$mortalityModel==0){numeric(0)}else{numeric(1)}
  ret$meanLogNM=if(conf$mortalityModel==0){numeric(0)}else{numeric(max(conf$keyMortalityMean,na.rm=TRUE)+1)}
  ret$logSdLogNM=if(conf$mortalityModel==0){numeric(0)}else{numeric(max(conf$keyMortalityObsVar,na.rm=TRUE)+1)}
  ret$logXtraSd=if(nrow(conf$keyXtraSd)==0){numeric(0)}else{numeric(length(unique(conf$keyXtraSd[,4])))}
  

  ## Reference points
  ret$logFScaleMSY <- 0
  ret$implicitFunctionDelta <- 0
  ret$logScaleFmsy <- 0
  ret$logScaleFmypyl <- 0
  ret$logScaleFmdy <- 0
  ret$logScaleFmax <- 0
  ret$logScaleFxdYPR <- numeric(0)
  ret$logScaleFxB0 <- numeric(0)
  ret$logScaleFcrash <- 0
  ret$logScaleFext <- 0
  ret$logScaleFxPercent <- numeric(0)
  ret$logScaleFlim <- 0
  ret$logScaleFmsyRange <- matrix(0,2,0)
  ret$splinePenalty <- 0
  
  ## Latent variables
  ret$logF=matrix(0, nrow=max(conf$keyLogFsta)+1,ncol=dat$noYears)
  ret$logN=matrix(0, nrow=conf$maxAge-conf$minAge+1, ncol=dat$noYears)
  ret$logSW=if(conf$stockWeightModel==0){matrix(0, nrow=0, ncol=0)}else{matrix(0, ncol=ncol(dat$stockMeanWeight), nrow=nrow(dat$stockMeanWeight)+spinoutyear)}
  ret$logCW=if(conf$catchWeightModel==0){matrix(0, nrow=0, ncol=0)}else{matrix(0, ncol=ncol(dat$catchMeanWeight), nrow=nrow(dat$catchMeanWeight)+spinoutyear)}  
  ret$logitMO=if(conf$matureModel==0){matrix(0, nrow=0, ncol=0)}else{matrix(0, ncol=ncol(dat$propMat), nrow=nrow(dat$propMat)+spinoutyear)}
  ret$logNM=if(conf$mortalityModel==0){matrix(0, nrow=0, ncol=0)}else{matrix(0, ncol=ncol(dat$natMor), nrow=nrow(dat$natMor)+spinoutyear)}
  return(ret)
}
