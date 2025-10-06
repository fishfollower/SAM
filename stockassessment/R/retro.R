##' reduce helper function to reduce data
##' @param data a data object as returned by the function setup.sam.data
##' @param year a vector of years to be excluded.  
##' @param fleet a vector of fleets to be excluded.
##' @param age a vector of ages fleets to be excluded.
##' @param conf an optional corresponding configuration to be modified along with the data change. Modified is returned as attribute "conf"  
##' @details When more than one vector is supplied they need to be of same length, as only the pairs are excluded
##' @export
reduce<-function(data, year=NULL, fleet=NULL, age=NULL, conf=NULL, onlyobs=FALSE){
    yOrig <- min(as.numeric(data$aux[,"year"])):max(as.numeric(data$aux[,"year"]))
    aOrig <- min(as.numeric(data$aux[data$aux[,"age"]>=0,"age"])):max(as.numeric(data$aux[,"age"]))
    nam<-c("year", "fleet", "age")[c(length(year)>0,length(fleet)>0,length(age)>0)]
    if((length(year)==0) & (length(fleet)==0) & (length(age)==0)){
        idx <- rep(TRUE,nrow(data$aux))
    }else{
        idx <- !do.call(paste, as.data.frame(data$aux[,nam,drop=FALSE])) %in% do.call(paste, as.data.frame(cbind(year=year, fleet=fleet, age=age)))
    }
    data$aux <- data$aux[idx,]
    data$auxData <- data$auxData[idx,]
    data$logobs <- data$logobs[idx]
    data$weight <- data$weight[idx]
    suf <- sort(unique(data$aux[,"fleet"])) # sort-unique-fleet
    data$noFleets <- length(suf)
    data$fleetTypes <- data$fleetTypes[suf]
    data$fleetCovarianceSize <- data$fleetCovarianceSize[suf]
    data$sampleTimesStart <- data$sampleTimesStart[suf]
    data$sampleTimesEnd <- data$sampleTimesEnd[suf]
    oldYears <- data$years
    data$years <- min(as.numeric(data$aux[,"year"])):max(as.numeric(data$aux[,"year"]))
    ages <- min(as.numeric(data$aux[data$aux[,"age"]>=0,"age"])):max(as.numeric(data$aux[,"age"]))
    data$noYears <- length(data$years)
    mmfun<-function(f,y, ff){idx<-which(data$aux[,"year"]==y & data$aux[,"fleet"]==f); ifelse(length(idx)==0, NA, ff(idx)-1)}
    data$idx1 <- outer(suf, data$years, Vectorize(mmfun,c("f","y")), ff=min)
    data$idx2 <- outer(suf, data$years, Vectorize(mmfun,c("f","y")), ff=max)
    data$idxCor <- data$idxCor[suf,match(data$years,oldYears)]
    data$nobs <- length(data$logobs[idx])  
    ##data$propMat <- data$propMat[rownames(data$propMat)%in%data$years, colnames(data$propMat)%in%ages,drop=FALSE]
    if(!onlyobs){
        data$propMat <- data$propMat[match(data$years,yOrig), match(ages,aOrig),drop=FALSE]
        data$stockMeanWeight <- data$stockMeanWeight[match(data$years,yOrig), match(ages,aOrig),drop=FALSE]
        data$natMor <- data$natMor[match(data$years,yOrig), match(ages,aOrig),drop=FALSE]
        if(length(dim(data$propF))==3){
            data$propF <- data$propF[match(data$years,yOrig), match(ages,aOrig),,drop=FALSE]
        }else{
            data$propF <- data$propF[match(data$years,yOrig), match(ages,aOrig),drop=FALSE]
        }
        data$propM <- data$propM[match(data$years,yOrig), match(ages,aOrig),drop=FALSE]
        if(length(dim(data$landFrac))==3){
            data$landFrac <- data$landFrac[na.omit(match(data$years,yOrig[1:nrow(data$landFrac)])), match(ages,aOrig),,drop=FALSE]
        }else{
            data$landFrac <- data$landFrac[na.omit(match(data$years,yOrig[1:nrow(data$landFrac)])), match(ages,aOrig),drop=FALSE]
        }
        if(length(dim(data$catchMeanWeight))==3){
            data$catchMeanWeight <- data$catchMeanWeight[na.omit(match(data$years,yOrig[1:nrow(data$landFrac)])), match(ages,aOrig),,drop=FALSE]
        }else{
            data$catchMeanWeight <- data$catchMeanWeight[na.omit(match(data$years,yOrig[1:nrow(data$landFrac)])), match(ages,aOrig),drop=FALSE]
        }
        if(length(dim(data$disMeanWeight))==3){
            data$disMeanWeight <- data$disMeanWeight[na.omit(match(data$years,yOrig[1:nrow(data$landFrac)])), match(ages,aOrig),,drop=FALSE]
        }else{
            data$disMeanWeight <- data$disMeanWeight[na.omit(match(data$years,yOrig[1:nrow(data$landFrac)])), match(ages,aOrig),drop=FALSE]
        }
        if(length(dim(data$landMeanWeight))==3){
            data$landMeanWeight <- data$landMeanWeight[na.omit(match(data$years,yOrig[1:nrow(data$landFrac)])), match(ages,aOrig),,drop=FALSE]
        }else{
            data$landMeanWeight <- data$landMeanWeight[na.omit(match(data$years,yOrig[1:nrow(data$landFrac)])), match(ages,aOrig),drop=FALSE]
        }
        data$TAC <- data$TAC[match(data$years,yOrig),,drop=FALSE]
        data$RecruitClimate <- data$RecruitClimate[match(data$years,yOrig),,,drop=FALSE]
    }
    data$aux[,"fleet"] <- match(data$aux[,"fleet"],suf) 
    data$minAgePerFleet <- tapply(as.integer(data$aux[,"age"]), INDEX=data$aux[,"fleet"], FUN=min)
    data$maxAgePerFleet <- tapply(as.integer(data$aux[,"age"]), INDEX=data$aux[,"fleet"], FUN=max)
    if(any(data$fleetTypes >= 80)){
        xtrMin <- as.integer(tapply(data$auxData[,5], INDEX=data$aux[,"fleet"], FUN=min))
        xtrMax <- as.integer(tapply(data$auxData[,5], INDEX=data$aux[,"fleet"], FUN=max))
        data$fleetCovarianceSize <- ifelse(data$fleetTypes >= 80,xtrMax-xtrMin+1,data$maxAgePerFleet - data$minAgePerFleet + 1)
    }else{
        data$fleetCovarianceSize <- data$maxAgePerFleet - data$minAgePerFleet + 1
    }
    attr(data,"fleetNames") <- attr(data,"fleetNames")[suf]
    if(!missing(conf)){
        .reidx <- function(x){
            if(any(x >= (-0.5) & !is.na(x))){
                xx <- x[x >= (-.5) & !is.na(x)]
                x[x >= (-.5) &!is.na(x)] <- match(xx,sort(unique(xx)))-1
            }
            x
        }
        conf$keyLogFsta <- .reidx(conf$keyLogFsta[suf,,drop=FALSE])
        conf$keyLogFpar <- .reidx(conf$keyLogFpar[suf,,drop=FALSE])
        conf$keyLogFmu <- .reidx(conf$keyLogFmu[suf,,drop=FALSE])
        conf$keyLogFrho <- .reidx(conf$keyLogFrho[suf,,drop=FALSE])
        conf$keyLogFbound_kappa <- .reidx(conf$keyLogFbound_kappa[suf,,drop=FALSE])
        conf$keyLogFbound_alpha <- .reidx(conf$keyLogFbound_alpha[suf,,drop=FALSE])
        conf$keyLogFbound_tau <- .reidx(conf$keyLogFbound_tau[suf,,drop=FALSE])
        conf$keyLogFboundTAC_kappa <- .reidx(conf$keyLogFboundTAC_kappa[suf,,drop=FALSE])
        conf$keyLogFboundTAC_alpha <- .reidx(conf$keyLogFboundTAC_alpha[suf,,drop=FALSE])
        conf$keyQpow <- .reidx(conf$keyQpow[suf,,drop=FALSE])
        conf$keyVarF <- .reidx(conf$keyVarF[suf,,drop=FALSE])
        conf$keyVarObs <- .reidx(conf$keyVarObs[suf,,drop=FALSE])
        conf$predVarObsLink <- .reidx(conf$predVarObsLink[suf,,drop=FALSE])
        conf$obsCorStruct <- conf$obsCorStruct[suf]
        conf$fixVarToWeight <- conf$fixVarToWeight[suf]
        conf$maxAgePlusGroup <- conf$maxAgePlusGroup[suf]
        conf$keyCorObs <- .reidx(conf$keyCorObs[suf,,drop=FALSE])
        conf$keyLogFseason <- .reidx(conf$keyLogFseason[suf,,drop=FALSE])
        yidx <- conf$keyScaledYears%in%data$aux[data$aux[,'fleet']==1,'year']
        if(length(conf$keyScaledYears)>0){
            conf$noScaledYears <- sum(yidx)
            conf$keyScaledYears <- as.vector(conf$keyScaledYears)[yidx]
            conf$keyParScaledYA <- .reidx(conf$keyParScaledYA[yidx,,drop=FALSE])
        }
        conf$keyBiomassTreat <- conf$keyBiomassTreat[suf]
        conf$obsLikelihoodFlag <- conf$obsLikelihoodFlag[suf]
        conf$fracMixObs <- conf$fracMixObs[suf]
        attr(data, "conf") <- conf
    }
    data
}


##' runwithout helper function
##' @param fit a fitted model object as returned from sam.fit
##' @param year a vector of years to be excluded.  When both fleet and year are supplied they need to be of same length, as only the pairs are excluded
##' @param fleet a vector of fleets to be excluded.  When both fleet and year are supplied they need to be of same length, as only the pairs are excluded
##' @param map map from original fit 
##' @param ... extra arguments to sam.fit
##' @details ...
##' @export
runwithout <- function(fit, year, fleet, ...){
    UseMethod("runwithout")
}

##' @param map map to use
##' @rdname runwithout
##' @method runwithout sam
##' @export
runwithout.sam <- function(fit, year=NULL, fleet=NULL, map=fit$obj$env$map,rm.unidentified=TRUE,initialize.par=TRUE,initialize.re=FALSE, ...){
  data <- reduce(fit$data, year=year, fleet=fleet, conf=fit$conf)      
  conf <- attr(data, "conf")
  fakefile <- file()
  sink(fakefile); on.exit(if(sink.number()>0){sink()})
  saveConf(conf, file="")
  sink()
  conf <- loadConf(data, fakefile, patch=TRUE)
  close(fakefile)
  par <- defpar(data,conf)
  if(initialize.par)
      par[!names(par)%in%c("logN", "logF", "logSW", "logCW", "logitMO", "logNM","logP","logitFseason")]<-fit$pl[!names(fit$pl)%in%c("missing", "logN", "logF", "logSW", "logCW", "logitMO", "logNM","logP","logitFseason")]
  if(initialize.re){
      par$logN <- fit$pl$logN[,match(data$years,fit$data$years)]
      par$logF <- fit$pl$logF[,match(data$years,fit$data$years)]
  }
  ret <- sam.fit(data, conf, par, rm.unidentified=rm.unidentified, map=map, lower=fit$low, upper=fit$hig, ...)
  return(ret)
}

##' retro run 
##' @param fit a fitted model object as returned from sam.fit
##' @param year either 1) a single integer n in which case runs where all fleets are reduced by 1, 2, ..., n are returned, 2) a vector of years in which case runs where years from and later are excluded for all fleets, and 3 a matrix of years were each column is a fleet and each column corresponds to a run where the years and later are excluded.    
##' @param ncores the number of cores to attempt to use
##' @param ... extra arguments to \code{\link{sam.fit}}
##' @details ...
##' @export
retro <- function(fit, year=NULL, ncores=detectCores(), ...){
    UseMethod("retro")
}

##' @rdname retro
##' @method retro sam
##' @importFrom parallel detectCores makeCluster clusterExport parLapply stopCluster clusterEvalQ
##' @export
retro.sam <- function(fit, year=NULL, ncores=detectCores(), ...){
  data <- fit$data
  y <- fit$data$aux[,"year"]
  f <- fit$data$aux[,"fleet"]
  suf <- sort(unique(f))
  maxy <- sapply(suf, function(ff)max(y[f==ff]))
  if(length(year)==1){
    mat <- sapply(suf,function(ff){my<-maxy[ff];my:(my-year+1)})
    if(year==1)mat <- matrix(mat, nrow=1)
  }
  if(is.vector(year) & length(year)>1){
    mat <- sapply(suf,function(ff)year)
  }
  if(is.matrix(year)){
    mat <- year
  }

  if(nrow(mat)>length(unique(y)))stop("The number of retro runs exceeds number of years")
  if(ncol(mat)!=length(suf))stop("Number of retro fleets does not match")

  setup <- lapply(1:nrow(mat),function(i)do.call(rbind,lapply(suf,function(ff)if(mat[i,ff]<=maxy[ff]){cbind(mat[i,ff]:maxy[ff], ff)})))
  if(ncores>1){
    cl <- makeCluster(ncores) #set up nodes
    on.exit(stopCluster(cl)) #shut it down
    clusterExport(cl, varlist="fit", envir=environment())
    lib.ver <- dirname(path.package("stockassessment"))
    clusterExport(cl, varlist="lib.ver", envir=environment())
    clusterEvalQ(cl, {library(stockassessment, lib.loc=lib.ver)})
    runs <- parLapply(cl, setup, function(s)runwithout(fit, year=s[,1], fleet=s[,2], ...))
  } else {
    runs <- lapply( setup, function(s)runwithout(fit, year=s[,1], fleet=s[,2], ...))
  }
  converg <- unlist(lapply(runs, function(x)x$opt$conv))
  if(any(converg!=0)) warning(paste0("retro run(s) ", paste0(which(converg!=0),collapse=",")," did not converge."))
  attr(runs, "fit") <- fit
  class(runs)<-"samset"
  runs
}

##' leaveout run 
##' @param fit a fitted model object as returned from sam.fit
##' @param fleet a list of vectors. Each element in the list specifies a run where the fleets mentioned are omitted 
##' @param ncores the number of cores to attemp to use
##' @param ... extra arguments to \code{\link{sam.fit}}
##' @details ...
##' @importFrom parallel detectCores makeCluster clusterExport parLapply stopCluster clusterEvalQ
##' @export
leaveout <- function(fit, fleet=as.list(2:fit$data$noFleets), ncores=detectCores(), ...){
  if(ncores>1){
    cl <- makeCluster(ncores) #set up nodes
    on.exit(stopCluster(cl)) #shut it down
    clusterExport(cl, varlist="fit", envir=environment())
    lib.ver <- dirname(path.package("stockassessment"))
    clusterExport(cl, varlist="lib.ver", envir=environment())
    clusterEvalQ(cl, {library(stockassessment, lib.loc=lib.ver)})
    runs <- parLapply(cl, fleet, function(f)runwithout(fit, fleet=f, ...))
  } else {
    runs <- lapply(fleet, function(f)runwithout(fit, fleet=f, ...))
  }
  converg <- unlist(lapply(runs, function(x)x$opt$conv))
  if(any(converg!=0)) warning(paste0("leavout run(s) ", paste0(which(converg!=0),collapse=",")," did not converge."))
  names(runs) <- paste0("w.o. ", lapply(fleet, function(x)paste(attr(fit$data,"fleetNames")[x], collapse=" and ")))
  attr(runs, "fit") <- fit
  class(runs)<-"samset"
  runs
}


##' Mohn's rho calculation
##' @param fits a samset object as returned from the retro function.
##' @param what a function computing the quantity to calculate Mohn's rho for (default NULL computes Fbar, SSB, and R).
##' @param lag lag applied to fits and reference fit.
##' @param ... extra arguments
##' @details ...
##' @export
mohn <- function(fits, what=NULL, lag=0, ...){
    UseMethod("mohn")
}

##' @rdname mohn
##' @method mohn samset
##' @export
mohn.samset <- function(fits, what=NULL, lag=0, ...){
  if(is.null(what)){
    what <- function(fit){
      ret <- cbind(rectable(fit,...)[,1], ssbtable(fit)[,1], fbartable(fit)[,1])
      add <- 0
      dots <- list(...)
      if(!is.null(dots$lagR)){
        if(dots$lagR == TRUE){
          add <- 1
        }
      }
      colnames(ret) <- c(paste("R(age ", fit$conf$minAge + add, ")", sep = ""), "SSB",
                         paste("Fbar(", fit$conf$fbarRange[1], "-", fit$conf$fbarRange[2], ")", sep = ""))
      ret
    }
  }
  ref <- what(attr(fits,"fit"))
  ret <- lapply(fits, what)
  bias <- lapply(ret, function(x){y<-rownames(x)[nrow(x)-lag]; (x[rownames(x)==y,]-ref[rownames(ref)==y,])/ref[rownames(ref)==y,]})
  colMeans(do.call(rbind,bias))
}
