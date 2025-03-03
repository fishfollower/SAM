##' Table helper 
##' @param fit returned object from sam.fit
##' @param what quoted name of what to extract
##' @param x rownames of table
##' @param trans function to be applied
##' @param ... extra arguments not currently used
##' @details ...
tableit <-function (fit, what, x=fit$data$years, trans=function(x)x,...){
    UseMethod("tableit")
}
##' @rdname tableit
##' @method tableit sam
##' @export
tableit.sam <- function (fit, what, x=fit$data$years, trans=function(x)x,...){
   idx<-names(fit$sdrep$value)==what
   y<-fit$sdrep$value[idx]
   ci<-y+fit$sdrep$sd[idx]%o%c(-2,2)
   ret<-trans(cbind(y,ci))
   rownames(ret)<-x
   colnames(ret)<-c("Estimate","Low","High")
   return(ret)
}

##' @rdname tableit
##' @method tableit samforecast
##' @export
tableit.samforecast <- function (fit, what, x=fit$data$years, trans=function(x)x,...){
    tab <- attr(fit,"tab")
    if(what == "logssb"){
        ret <- tab[,7:9,drop=FALSE]
    }else if(what == "logfbar"){
        ret <- tab[,1:3,drop=FALSE]
    }else if(what == "logCatch"){
        ret <- tab[,10:12,drop=FALSE]
    }else if(what == "logR"){
        ret <- tab[,4:6,drop=FALSE]
    }else if(what == "logtsb"){
        if(any(grepl("^tsb:",colnames(tab)))){
            ret <- tab[,grepl("^tsb:",colnames(tab)), drop=FALSE]
        }else{
            stop("The forecast was not made with addTSB=TRUE")
        }
    }else{
        stop("what not implemented for samforecast")
    }
    ret <- trans(log(ret))
    colnames(ret)<-c("Estimate","Low","High")
    return(ret)
}


##' Life expectancy table 
##' @param  fit ...
##' @param atRecruit If true, show life expectancy given survival until minAge, otherwise show life expectancy at birth
##' @param ... extra arguments not currently used
##' @details ...
##' @export
lifeexpectancytable <- function(fit, atRecruit=TRUE, ...){
     UseMethod("lifeexpectancytable")
}

##' @rdname lifeexpectancytable
##' @method lifeexpectancytable default
##' @export
lifeexpectancytable.default <- function(fit,atRecruit=TRUE, ...){
    if(!.checkFullDerived(fit)) stop("This function needs a fit with all derived values. Fit with `fullDerived=TRUE` or update with `getAllDerivedValues`.")
    ret<-tableit(fit, ifelse(atRecruit,"logLifeExpectancyRec","logLifeExpectancy"), trans=exp,...)
    return(ret)
}

##' Years Lost table 
##' @param  fit ...
##' @param cause Fisning, Other, or LifeExpectancy
##' @param ... extra arguments not currently used
##' @details ...
##' @export
yearslosttable <- function(fit, cause, ...){
    UseMethod("yearslosttable")
}

##' @rdname yearslosttable
##' @method yearslosttable default
##' @export
yearslosttable.default <- function(fit,cause = c("Fishing","Other","LifeExpectancy"), ...){
    if(!.checkFullDerived(fit)) stop("This function needs a fit with all derived values. Fit with `fullDerived=TRUE` or update with `getAllDerivedValues`.")
    cv <- match.arg(cause)
    what <- switch(cv,
                   "Fishing"="logYLTF",
                   "Other"="logYLTM",
                   "logYNL")
    ret<-tableit(fit, what=what, trans=exp,...)
    return(ret)
}

##' SSB table 
##' @param  fit ...
##' @param ... extra arguments not currently used
##' @details ...
##' @export
ssbtable<-function(fit,...){
    UseMethod("ssbtable")
}
##' @rdname ssbtable
##' @method ssbtable default
##' @export
ssbtable.default <- function(fit,...){
   ret<-tableit(fit, "logssb", trans=exp,...)
   return(ret)
}

##' TSB table 
##' @param  fit ... 
##' @param ... extra arguments not currently used
##' @details ...
##' @export
tsbtable<-function(fit,...){
    UseMethod("tsbtable")
}
##' @rdname tsbtable
##' @method tsbtable default
##' @export
tsbtable.default <- function(fit,...){
   ret<-tableit(fit, "logtsb", trans=exp,...)
   return(ret)
}

##' Fbar table 
##' @param  fit ... 
##' @param ... extra arguments not currently used
##' @details ...
##' @export
fbartable<-function(fit,...){
    UseMethod("fbartable")
}
##' @rdname fbartable
##' @method fbartable default
##' @export
fbartable.default <- function(fit,...){
   ret<-tableit(fit, "logfbar", trans=exp)
   return(ret)
}

##' Recruit table 
##' @param  fit ...
##' @param lagR use the age after the youngest as R 
##' @param ... extra arguments not currently used
##' @details ...
##' @export
rectable<-function(fit, lagR=FALSE, ...){
    UseMethod("rectable")
}
##' @rdname rectable
##' @method rectable default
##' @export
rectable.default <- function(fit, lagR=FALSE,...){ 
   ret<-tableit(fit, ifelse(lagR,"logLagR","logR"), trans=exp)
   return(ret)
}

##' Catch table 
##' @param  fit object returned from sam.fit
##' @param obs.show logical add a column with catch sum of product rowsums(C*W)
##' @param ... extra arguments not currently used
##' @details ...
##' @export
catchtable<-function(fit, obs.show=FALSE,...){
    UseMethod("catchtable")
}
##' @rdname catchtable
##' @method catchtable sam
##' @export
catchtable.sam <- function(fit, obs.show=FALSE,...){
   CW <- fit$data$catchMeanWeight
   ## CW <- CW[apply(!is.na(CW),1,all),] 
   ## xx <- as.integer(rownames(CW))
   ## ret <- tableit(fit, x=xx, "logCatch", trans=exp)
   if(length(dim(CW))==2)CW<-array(CW,dim=c(dim(CW),1), dimnames = c(dimnames(CW),list("Residual catch")))
   nf <- dim(CW)[3]
   yy <- as.integer(rownames(CW[,,1]))
   aa <- as.integer(colnames(CW[,,1]))
   ret <- tableit(fit, x=yy, "logCatch", trans=exp)
   idx <- fit$data$fleetTypes[fit$data$aux[,'fleet']]%in%c(0)
   o <- exp(fit$data$logobs[idx])
   f <- fit$data$aux[idx,"fleet"]
   a <- match(fit$data$aux[idx,"age"],aa)
   y <- match(fit$data$aux[idx,"year"],yy)
   w <- CW[cbind(y,a,f)]
   cw <- o*w
   sop <- tapply(cw, INDEX=y, FUN=sum, na.rm=TRUE)
   if(obs.show){
      sop.catch <- rep(NA, nrow(ret))
      sop.catch[as.integer(names(sop))] <- sop
      ret<-cbind(ret,sop.catch=sop.catch)
   }
   return(ret)
}

##' CatchByFleet table 
##' @param  fit object returned from sam.fit
##' @param obs.show logical add a column with catch sum of product rowsums(C*W)
##' @details ...
##' @export
catchbyfleettable<-function(fit, obs.show=FALSE){
   CW <- fit$data$catchMeanWeight
   if(length(dim(CW))==2)CW<-array(CW,dim=c(dim(CW),1))
   nf <- dim(CW)[3]
   yy <- as.integer(rownames(CW[,,1]))
   aa <- as.integer(colnames(CW[,,1]))
   ret <- tableit(fit, x=rep(yy, fit$data$noFleets), "logCatchByFleet", trans=exp)
   ret <- cbind(ret, fleet=rep(1:fit$data$noFleets, each=length(yy)))
   allzero<-which(tapply(ret[,"Estimate"], INDEX=ret[,"fleet"], FUN=sum)==0)
   ret<-ret[!ret[,"fleet"]%in%allzero,]
   idx <- fit$data$fleetTypes[fit$data$aux[,'fleet']]%in%c(0)  
   o <- exp(fit$data$logobs[idx])
   f <- fit$data$aux[idx,"fleet"]
   a <- match(fit$data$aux[idx,"age"],aa)
   y <- match(fit$data$aux[idx,"year"],yy)
   w <- CW[cbind(y,a,f)]
   cw <- o*w
   fnam <- attr(fit$data, "fleetNames")[unique(f)]
   est <- matrix(ret[,"Estimate"], ncol=nf)
   colnames(est)<-paste0("Catch(",fnam,")")
   low <- matrix(ret[,"Low"], ncol=nf)
   colnames(low)<-paste0("Low(",fnam,")")
   hig <- matrix(ret[,"High"], ncol=nf)
   colnames(hig)<-paste0("High(",fnam,")")
   sop <- xtabs(cw~y+f)
   ret <- cbind(est,low,hig)
   rownames(ret)<-yy
   if(obs.show){
      sop.catch <- matrix(NA, nrow=nrow(est), ncol=ncol(est))
      sop.catch[as.integer(rownames(sop)),] <- sop
      colnames(sop.catch)<-paste0("Obs(",fnam,")")
      ret<-cbind(ret,sop.catch=sop.catch)
   }
   return(ret)
}

##' N table 
##' @param  fit ... 
##' @param ... extra arguments not currently used
##' @details ...
##' @export
ntable <- function(fit,...){
    UseMethod("ntable")
}
##' @rdname ntable
##' @method ntable sam
##' @export
ntable.sam <- function(fit,...){
   ret <- exp(t(fit$pl$logN))
   colnames(ret) <- fit$conf$minAge:fit$conf$maxAge
   rownames(ret) <- fit$data$years
   return(ret)
}

##' F-at-age table 
##' @param  fit a fitted object of class 'sam' as returned from sam.fit
##' @param ... extra arguments not currently used
##' @details ...
##' @export
faytable <- function(fit,...){
    UseMethod("faytable")
}
##' @rdname faytable
##' @param fleet the fleet number(s) to return F summed for (default is to return the sum of all residual fleets). 
##' @method faytable sam
##' @export
faytable.sam <- function(fit, fleet=which(fit$data$fleetTypes==0), ...){
   getfleet <- function(f){
     idx <- fit$conf$keyLogFsta[f,]+2    
     ret <- cbind(NA,exp(t(fit$pl$logF)))[,idx]
     ret[is.na(ret)] <- 0
     ret
   }
   ret <- Reduce("+",lapply(fleet,getfleet))
   colnames(ret) <- fit$conf$minAge:fit$conf$maxAge
   rownames(ret) <- fit$data$years
   return(ret)
}

##' Catch-at-age in numbers table 
##' @param fit a fitted object of class 'sam' as returned from sam.fit
##' @param fleet the fleet number(s) to return catch summed for (default is to return the sum of all residual fleets).  
##' @details ...
##' @export
caytable <- function(fit, fleet=which(fit$data$fleetTypes==0)){
   getfleet <- function(f){
     idx <- fit$conf$keyLogFsta[f,]+2    
     F <- cbind(NA,exp(t(fit$pl$logF)))[,idx]
     F[is.na(F)] <- 0
     M <- fit$data$natMor
     N <- exp(t(fit$pl$logN))
     F/(F+M)*N*(1-exp(-F-M))
   }
   ret <- Reduce("+",lapply(fleet,getfleet)) 
   colnames(ret) <- fit$conf$minAge:fit$conf$maxAge
   rownames(ret) <- fit$data$years
   return(ret)
}

##' parameter table 
##' @param  fit ... 
##' @param ... extra arguments not currently used
##' @details ...
##' @export
partable <- function(fit,...){
    UseMethod("partable")
}
##' @rdname partable
##' @method partable sam
##' @export
partable.sam <- function(fit,...){
  param <- coef(fit)
  nam <- names(param)
  dup <- duplicated(nam)
  namadd <- rep(0, length(nam))
  for (i in 2:length(dup)) {
    if(dup[i])namadd[i] <- namadd[i - 1] + 1
  }
  nam <- paste(nam, namadd, sep = "_")
  ret<-cbind(param, attr(param,"sd"))
  ex<-exp(ret[,1])
  lo<-exp(ret[,1]-2*ret[,2])
  hi<-exp(ret[,1]+2*ret[,2])
  ret<-cbind(ret,ex,lo,hi)
  colnames(ret)<-c("par", "sd(par)", "exp(par)", "Low", "High")
  rownames(ret)<-nam
  return(ret)
}

##' model table 
##' @param fits A sam fit as returned from the sam.fit function, or a collection c(fit1, fit2, ...) of such fits  
##' @param ... extra arguments not currently used
##' @details ...
##' @importFrom stats AIC pchisq
##' @export
modeltable <- function(fits,...){
    UseMethod("modeltable")
}
##' @rdname modeltable
##' @method modeltable sam
##' @export
modeltable.sam <- function(fits,...){
    modeltable(c(fits))
}
##' @rdname modeltable
##' @method modeltable samset
##' @export
modeltable.samset <- function(fits,...){
    if(!is.null(attr(fits,"fit"))){
      fits[[length(fits)+1]] <- attr(fits,"fit")
      fits <- fits[c(length(fits),1:(length(fits)-1))]
    }
    fits <- fits[!sapply(fits, is.null)]
    if(is.null(names(fits))){
        nam <- paste("M", 1:length(fits), sep="")
    }else{
        nam <-ifelse(names(fits)=="",paste("M", 1:length(fits), sep=""), names(fits))
    }
    logL <- sapply(fits, logLik)
    npar <- sapply(fits, function(f)attr(logLik(f),"df"))
    aic <- sapply(fits, AIC)
    res <- cbind("log(L)"=logL, "#par"=npar, "AIC"=aic)
    rownames(res) <- nam
    o <- 1:length(fits)
    if(length(fits)==2){
        o <- order(npar, decreasing=TRUE)
        if(npar[o[1]]>npar[o[2]]){
            df <- npar[o[1]]-npar[o[2]]
            D <- 2*(logL[o[1]]-logL[o[2]])
            P <- 1-pchisq(D,df)
            cnam <- paste0("Pval( ",nam[o[1]]," -> ",nam[o[2]], " )")
            res <- cbind(res, c(NA, P)[o])
            colnames(res)[ncol(res)] <- cnam
        }
    }
    return(res[o,,drop=FALSE])
}

##' table of survey catchabilities
##' @param  fit ...
##' @param ... extra arguments not currently used
##' @details ...
##' @export
qtable<-function(fit,...){
    UseMethod("qtable")
}

##' table of survey catchabilities
##' @param fit A sam fit as returned from the sam.fit function
##' @param ... extra arguments not currently used
##' @export
qtable.sam<-function(fit,...){
  key <- fit$conf$keyLogFpar[-fit$data$fleetTypes!=0, ] + 1
  key[key==0]<-NA
  cf<-coef(fit)
  cfsd <- attr(cf,"sd")
  qt <- matrix( cf[ names(cf)=="logFpar" ][key],nrow=nrow(key),ncol=ncol(key) )
  rownames(qt) <- attr(fit$data, "fleetNames")[-fit$data$fleetTypes!=0]
  colnames(qt)<-fit$conf$minAge:fit$conf$maxAge
  sds <- qt
  sds[] <- cfsd[ names(cf)=="logFpar" ][key]
  attr(qt,"sd") <- sds
  class(qt) <- "samqtable"
  qt
}

##' Yield per recruit calculation
##' @param fit the object returned from sam.fit
##' @param Flimit Upper limit for Fbar
##' @param Fdelta increments on the Fbar axis 
##' @param aveYears Number of years back to use when calculating averages (selection, weights, ...)
##' @param sprProp Proportion of SPR at F=0, for example 0.35 if F0.35SPR
##' @param ... extra arguments not currently used
##' @export
ypr<-function(fit, Flimit=2, Fdelta=0.01, aveYears=min(15,length(fit$data$years)), sprProp=0.35, ...){
    UseMethod("ypr")
}
##' @rdname ypr
##' @method ypr sam
##' @export
ypr.sam <- function(fit, Flimit=2, Fdelta=0.01, aveYears=min(15,length(fit$data$years)), sprProp=0.35, ...){
    if(fit$conf$stockWeightModel==1){
      stockMeanWeight <- exp(fit$pl$logSW)[1:nrow(fit$data$stockMeanWeight),]
      attr(stockMeanWeight, "dimnames")<-attr(fit$data$stockMeanWeight, "dimnames")  
    }else{
      stockMeanWeight <- fit$data$stockMeanWeight
    }
    if(fit$conf$catchWeightModel==1){
      catchMeanWeight <- exp(fit$pl$logCW)[1:dim(fit$data$catchMeanWeight)[1],,,drop=FALSE]
      attr(catchMeanWeight, "dimnames")<-attr(fit$data$catchMeanWeight, "dimnames")
    }else{
      catchMeanWeight <- fit$data$catchMeanWeight
    }
    if(fit$conf$mortalityModel==1){
      natMor <- exp(fit$pl$logNM)[1:nrow(fit$data$natMor),]
      attr(natMor, "dimnames")<-attr(fit$data$natMor, "dimnames")  
    }else{
      natMor <- fit$data$natMor
    }
    if(fit$conf$matureModel==1){
      propMat <- plogis(fit$pl$logitMO)[1:nrow(fit$data$propmat),]
      attr(propMat, "dimnames")<-attr(fit$data$propMat, "dimnames")  
    }else{
      propMat <- fit$data$propMat
    }

    barAges <- do.call(":", as.list(fit$conf$fbarRange)) + (1 - fit$conf$minAge)
    last.year.used = max(fit$data$years)
    idxno <- which(fit$data$years == last.year.used)
    F <- t(faytable(fit))
    F[is.na(F)] <- 0
    sel <- function() {
        Sa <- rep(0, nrow(F))
        K <- 0
        for (i in 0:(aveYears - 1)) {
            thisF <- F[, idxno - i]
            Sa <- Sa + thisF
            K <- K + fbartable(fit)[idxno - i]
        }
        return(Sa/K)
    }
    aveByCatch <- function(X) {
        fun <- function(f) {
            cay <- caytable(fit, f)
            XX <- X[, , f]
            cay <- cay[rownames(cay) %in% rownames(XX), ]
            cay.safe <- cay
            cay.safe[cay == 0] <- 1
            cay * XX/cay.safe
        }
        Reduce("+", lapply(which(fit$data$fleetTypes == 0), fun))
    }
    ave.sl <- sel()
    ave.sw <- colMeans(stockMeanWeight[(idxno - aveYears + 1):idxno, , drop = FALSE])
    ave.cw <- colMeans(aveByCatch(catchMeanWeight)[(idxno - aveYears + 1):(idxno - 1), , drop = FALSE])
    ave.pm <- colMeans(propMat[(idxno - aveYears + 1):idxno, , drop = FALSE])
    ave.nm <- colMeans(natMor[(idxno - aveYears + 1):idxno, , drop = FALSE])
    ave.lf <- colMeans(aveByCatch(fit$data$landFrac)[(idxno - aveYears + 1):(idxno - 1), , drop = FALSE])
    ave.cw.land <- colMeans(aveByCatch(fit$data$landMeanWeight)[(idxno - aveYears + 1):(idxno - 1), , drop = FALSE])
    N <- numeric(length(ave.nm))
    deltafirst <- 1e-05
    delta <- Fdelta
    scales <- c(0, deltafirst, seq(0.01, Flimit, by = delta))
    yields <- numeric(length(scales))
    ssbs <- numeric(length(scales))
    for (i in 1:length(scales)) {
        scale <- scales[i]
        F <- ave.sl * scale
        Z <- ave.nm + F
        N <- exp(-cumsum(c(0,Z[-length(Z)])))
        N[length(N)]<-N[length(N)]/(1-exp(-Z[length(N)]))
        C <- F/Z * (1 - exp(-Z)) * N * ave.lf
        Y <- sum(C * ave.cw)
        yields[i] <- Y
        ssbs[i] <- sum(N * ave.pm * ave.sw)
    }
    fmaxidx <- which.max(yields)
    fmax <- scales[fmaxidx]
    deltaY <- diff(yields)
    f01idx <- which.min((deltaY/delta - 0.1 * deltaY[1]/deltafirst)^2) + 1
    f01 <- scales[f01idx]
    f35spridx <- which.min((ssbs - sprProp * ssbs[1])^2) + 1
    f35 <- scales[f35spridx]
    fbarlab <- substitute(bar(F)[X - Y], list(X = fit$conf$fbarRange[1], 
        Y = fit$conf$fbarRange[2]))
    ret <- list(fbar = scales, ssb = ssbs, yield = yields, fbarlab = fbarlab, 
        fsprProp = f35, f01 = f01, fmax = fmax, fsprPropIdx = f35spridx, 
        f01Idx = f01idx, fmaxIdx = fmaxidx, sprProp = sprProp)
    class(ret) <- "samypr"
    ret
}






##' rmax table 
##' @param  fit ...
##' @param ... extra arguments not currently used
##' @details ...
##' @export
rmaxtable<-function(fit,...){
    UseMethod("rmaxtable")
}
##' @rdname rmaxtable
##' @method rmaxtable default
##' @export
rmaxtable.default <- function(fit,...){
    if(!.checkFullDerived(fit)) stop("This function needs a fit with all derived values. Fit with `fullDerived=TRUE` or update with `getAllDerivedValues`.")
    ret<-tableit(fit, "logrmax", trans=exp,...)
   return(ret)
}



##' Generation length table 
##' @param  fit ...
##' @param ... extra arguments not currently used
##' @details ...
##' @export
generationlengthtable<-function(fit,...){
    UseMethod("generationlengthtable")
}
##' @rdname generationlengthtable
##' @method generationlengthtable default
##' @export
generationlengthtable.default <- function(fit,...){
    if(!.checkFullDerived(fit)) stop("This function needs a fit with all derived values. Fit with `fullDerived=TRUE` or update with `getAllDerivedValues`.")
    ret<-tableit(fit, "logGenerationLength", trans=exp,...)
    return(ret)
}


##' YPR table 
##' @param  fit ...
##' @param ... extra arguments not currently used
##' @details ...
##' @export
yprtable<-function(fit,...){
    UseMethod("yprtable")
}
##' @rdname yprtable
##' @method yprtable default
##' @export
yprtable.default <- function(fit,...){
    if(!.checkFullDerived(fit)) stop("This function needs a fit with all derived values. Fit with `fullDerived=TRUE` or update with `getAllDerivedValues`.")
    ret<-tableit(fit, "logYPR", trans=exp,x=as.numeric(rownames(fit$data$catchMeanWeight)),...)
    return(ret)
}


##' SPR table 
##' @param  fit ...
##' @param ... extra arguments not currently used
##' @details ...
##' @export
sprtable<-function(fit,...){
    UseMethod("sprtable")
}
##' @rdname sprtable
##' @method sprtable default
##' @export
sprtable.default <- function(fit,...){
    if(!.checkFullDerived(fit)) stop("This function needs a fit with all derived values. Fit with `fullDerived=TRUE` or update with `getAllDerivedValues`.")
    ret<-tableit(fit, "logSPR", trans=exp,x=as.numeric(rownames(fit$data$catchMeanWeight)),...)
   return(ret)
}


##' equilibrium biomass table 
##' @param  fit ...
##' @param ... extra arguments not currently used
##' @details ...
##' @export
equilibriumbiomasstable<-function(fit,...){
    UseMethod("equilibriumbiomasstable")
}
##' @rdname equilibriumbiomasstable
##' @method equilibriumbiomasstable default
##' @export
equilibriumbiomasstable.default <- function(fit,...){
    if(!.checkFullDerived(fit)) stop("This function needs a fit with all derived values. Fit with `fullDerived=TRUE` or update with `getAllDerivedValues`.")
    ret<-tableit(fit, "logSe", trans=exp,x=as.numeric(rownames(fit$data$catchMeanWeight)),...)
    return(ret)
}



##' B0 biomass table 
##' @param  fit ...
##' @param ... extra arguments not currently used
##' @details ...
##' @export
b0table<-function(fit,...){
    UseMethod("b0table")
}
##' @rdname b0table
##' @method b0table default
##' @export
b0table.default <- function(fit,...){
    if(!.checkFullDerived(fit)) stop("This function needs a fit with all derived values. Fit with `fullDerived=TRUE` or update with `getAllDerivedValues`.")
    ret<-tableit(fit, "logB0", trans=exp,x=as.numeric(rownames(fit$data$catchMeanWeight)),...)
    return(ret)
}
