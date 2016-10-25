##' runwithout helper function
##' @param fit a fittes model objest as returned from sam.fit
##' @param year a vector of years to be excluded.  When both fleet and year are supplied they need to be of same length as only the pairs are excluded
##' @param fleet a vector of fleets to be excluded.  When both fleet and year are supplied they need to be of same length as only the pairs are excluded
##' @param ... extra arguments to sam.fit
##' @details ...
##' @export
runwithout <- function(fit, year=NULL, fleet=NULL, ...){
  data <- fit$data
  if(length(year)==0 & length(fleet)==0){
    idx <- 1:nrow(data$obs)
  }
  if(length(year)>0 & length(fleet)==0){
    idx <- which(!data$obs[,"year"]%in%year)
  }
  if(length(fleet)>0 & length(year)==0){
    idx <- which(!data$obs[,"fleet"]%in%fleet)
  }
  if(length(fleet)>0 & length(year)>0){
    if(length(fleet)!=length(year))stop("length of fleet not equal to length of year: When both fleet and year are supplied they need to be of same length as only the pairs are excluded")
    idx<-!apply(apply(cbind(year,fleet),1,function(s)data$obs[,"year"]%in%s[1] & data$obs[,"fleet"]%in%s[2]),1,any) # could this be done simpler???
  }
  data$obs <- data$obs[idx,]
  data$logobs <- data$logobs[idx]
  suf <- sort(unique(data$obs[,"fleet"])) # sort-unique-fleet
  data$noFleets <- length(suf)
  data$fleetTypes <- data$fleetTypes[suf]
  data$sampleTimes <- data$sampleTimes[suf]
  data$years <- min(as.numeric(data$obs[,"year"])):max(as.numeric(data$obs[,"year"]))
  data$noYears <- length(data$years)
  data$idx1 <- sapply(data$years, function(y){idx<-which(as.numeric(data$obs[,"year"])==y);ifelse(length(idx)==0,-1,min(idx))}) 
  data$idx2 <- sapply(data$years, function(y){idx<-which(as.numeric(data$obs[,"year"])==y);ifelse(length(idx)==0,-1,max(idx))})
  data$nobs <- length(data$logobs[idx])
  
  data$propMat <- data$propMat[rownames(data$propMat)%in%data$years, ]
  data$stockMeanWeight <- data$stockMeanWeight[rownames(data$stockMeanWeight)%in%data$years, ]
  data$natMor <- data$natMor[rownames(data$natMor)%in%data$years, ]
  data$propF <- data$propF[rownames(data$propF)%in%data$years, ]
  data$propM <- data$propM[rownames(data$propM)%in%data$years, ]
  # there may be a later issue with the year range of the last four
  data$landFrac <- data$landFrac[rownames(data$landFrac)%in%data$years, ]  
  data$catchMeanWeight <- data$catchMeanWeight[rownames(data$catchMeanWeight)%in%data$years, ]
  data$disMeanWeight <- data$disMeanWeight[rownames(data$disMeanWeight)%in%data$years, ]
  data$landMeanWeight <- data$landMeanWeight[rownames(data$landMeanWeight)%in%data$years, ]
  data$obs[,"fleet"] <- match(data$obs[,"fleet"],suf)

  .reidx <- function(x){
    if(any(x >= (-0.5))){
      xx<-x[x >= (-.5)]
      x[x >= (-.5)] <- match(xx,sort(unique(xx)))-1
    };
    x
  }
                                  
  conf <- fit$conf
  conf$keyLogFsta <- .reidx(conf$keyLogFsta[suf,,drop=FALSE])
  conf$keyLogFpar <- .reidx(conf$keyLogFpar[suf,,drop=FALSE])
  conf$keyQpow <- .reidx(conf$keyQpow[suf,,drop=FALSE])
  conf$keyVarF <- .reidx(conf$keyVarF[suf,,drop=FALSE])
  conf$keyVarObs <- .reidx(conf$keyVarObs[suf,,drop=FALSE])
  yidx <- conf$keyScaledYears%in%data$obs[,'year']
  conf$noScaledYears[1,1] <- sum(yidx)
  conf$keyScaledYears <- conf$keyScaledYears[,yidx,drop=FALSE]
  conf$keyParScaledYA <- .reidx(conf$keyParScaledYA[yidx,,drop=FALSE])
  par<-defpar(data,conf)
  
  ret <- sam.fit(data,conf,par)
  return(ret)
}
