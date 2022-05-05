
.refpointEnum <- c(None = -99,
                   FixedF = -1,
                   StatusQuo = 0,
                   MSY = 1,
                   MSYRange = 2,
                   Max = 3,
                   xdYPR = 4,
                   xSPR = 5,
                   xB0 = 6,
                   MYPYLdiv = 7,
                   MYPYL = 8,
                   MDY = 9,
                   Crash = 10,
                   Ext = 11,
                   Lim = 12)


.logFtoSel <- function(logF, selYears, conf){
    logSel <- log(rowSums(exp(logF[,selYears + 1,drop=FALSE])))
    frng <- conf$fbarRange[1]:conf$fbarRange[2]
    logFsum <- log(sum(colMeans(exp(logF)[conf$keyLogFsta[1,frng-conf$minAge+1]+1,selYears + 1,drop=FALSE])))
    logSel - logFsum
}

.perRecruitR <- function(logf, fit, nYears, aveYears, selYears, pl = fit$pl, ct = 0, logCustomSel = numeric(0)){
    if(length(logCustomSel) > 0){
        sel <- exp(logCustomSel)
    }else{
        sel <- exp(.logFtoSel(pl$logF, selYears, fit$conf))
    }
    as.data.frame(.Call("perRecruitR",
                        logFbar = logf,
                        tmbdat = fit$obj$env$data,
                        pl = pl,
                        sel = sel,
                        aveYears = aveYears,
                        nYears = ifelse(nYears==0,150,nYears),
                        CT = ct
                        ))
}


.refpointParser <- function(x, ...){
    if(length(x) > 1)
        stop("Only one at a time!")
    ## Allow shorthand
    if(grepl("^StatusQuo$",x)){
        x <- "StatusQuo-0"
    }else if(grepl("^F=\\.$",x)){
        x <- "F=0.0"
    }else if(grepl("^F=\\.[[:digit:]]+$",x)){
        x <- gsub("=\\.","=0.",x)
    }else if(grepl("^F=[[:digit:]]+\\.$",x)){
        x <- gsub("\\.$",".0",x)
    }else if(grepl("^MYPYL$",x)){
        x <- "1.0MYPYL"
    }else if(grepl("^[[:digit:]]+MYPYL$",x)){
        x <- gsub("MYPYL",".0MYPYL",x)
    }else if(grepl("^\\.[[:digit:]]+MYPYL$",x)){
        x <- gsub("\\.","0.",x)
    }
    ## NOTE: enum, typePatterns, and xvalPatterns are assumed to have the same order and list elements
  
    typePatterns <- list(None = "^$",       # E.g., ""
                         FixedF = "^F=([[:digit:]]+\\.)?([[:digit:]]+)$", # E.g., "F=0"
                         StatusQuo = "^StatusQuo(-[[:digit:]]+)?$", # E.g., "StatusQuo" or "StatusQuo-2"
                         MSY = "^MSY$", # E.g., "MSY"
                         MSYRange = "^0\\.[[:digit:]]+MSYRange$", # E.g., "0.95MSYRange"
                         Max = "^Max$",
                         xdYPR = "^0\\.[[:digit:]]+dYPR$", # E.g., "0.1dYPR"
                         xSPR = "^0\\.[[:digit:]]+SPR$", # E.g., "0.35SPR"
                         xB0 =  "^0\\.[[:digit:]]+B0$", # E.g., "0.2B0"
                         MYPYLdiv = "^MYPYLdiv$",
                         MYPYL = "^[[:digit:]]+\\.[[:digit:]]+MYPYL$",
                         MDY = "^MDY$",
                         Crash = "^Crash$",
                         Ext = "^Ext$",
                         Lim = "^Lim$")
    xvalPatterns <- list(None = NA,
                         FixedF = c("(F=)(([[:digit:]]+\\.)?([[:digit:]]+))","\\2"),
                         StatusQuo = c("(StatusQuo)(-?)(([[:digit:]]+)?)","\\3"),
                         MSY = NA,
                         MSYRange = c("(0\\.[[:digit:]]+)(MSYRange)","\\1"),
                         Max = NA,
                         xdYPR = c("(0\\.[[:digit:]]+)(dYPR)","\\1"),
                         xSPR = c("(0\\.[[:digit:]]+)(SPR)","\\1"),
                         xB0 =  c("(0\\.[[:digit:]]+)(B0)","\\1"),
                         MYPYLdiv = NA,
                         MYPYL = c("([[:digit:]]+\\.[[:digit:]])(MYPYL)","\\1"),
                         MDY = NA,
                         Crash = NA,
                         Ext = NA,
                         Lim = NA)
    rpIndx <- which(sapply(typePatterns, function(p) grepl(p, x)))
    if(length(rpIndx) == 0)
        stop(sprintf("Error in reference point specification. %s not recognized.",x))
    xVal <- numeric(0)
    if(!is.na(xvalPatterns[[rpIndx]][1])){
        pp <- xvalPatterns[[rpIndx]]
        xVal <- as.numeric(gsub(pp[1], pp[2], x))
    }
    c(list(rpType = .refpointEnum[rpIndx],
           xVal = xVal),
           list(...))    
}

.refpointMerger <- function(x,y, keep.order = FALSE, keep.all = FALSE){
    ## If y is the first, just return
    if(length(x) == 0)
        return(list(y))
    ## Check if all but xVal is the same
    nms <- setdiff(names(y),"xVal")
    isMatch <- sapply(x, function(xx) isTRUE(all.equal(xx[nms], y[nms])))
    ## If no match, add y to list
    if(!any(isMatch))
        return(c(x,list(y)))
    ## If match, add xVal
    .unique <- if(keep.all){ identity }else{ unique }
    .sort <- if(keep.order){ identity }else{ sort }
    .process <- function(v) sort(unique(v))
    x[isMatch] <- lapply(x[isMatch], function(r){r$xVal <- .sort(.unique(c(r$xVal,y$xVal))); r})
    return(x)
}


.refpointStartingValue <- function(rparg, fit, Fsequence, fay = faytable(fit), fbar = fbartable(fit)[,1]){

    ## Check validity of SRR-refpoint combination
    srmc <- fit$conf$stockRecruitmentModelCode
    if(srmc == 0){
        invalidRP <- c(-1,1,2,6,7,8,9,10,11,12)
        srName <- "random walk (no equilibrium biomass)"
    }else{
        invalidRP <- c(-Inf)
        srName <- "the selected"
    }
    
    i <- match(rparg$rpType, invalidRP)
    if(length(i) > 0 && !is.na(i))
        stop(sprintf("%s cannot be estimated for %s recruitment",names(.refpointEnum)[match(invalidRP[i],.refpointEnum)],srName))
    
    if(rparg$rpType == -99){            # None
        rparg$logF0 <- numeric(0)
        return(rparg)
    }else if(rparg$rpType == -1){       # FixedF
        rparg$logF0 <- log(rparg$xVal)
        return(rparg)
    }else if(rparg$rpType == 0){        # Status quo (is not using logF0)
        rparg$logF0 <- numeric(0)
        return(rparg)        
    }
    ## Wrapper for perRecruit calculations    
    dd <- fit$obj$env$data #lapply(tmbargs$data,dataSanitize)
    pl <- fit$pl #lapply(tmbargs$parameters,parameterSanitize)
    getPR <- function(logf){ .perRecruitR(logf,fit,rparg$nYears, rparg$aveYears, rparg$selYears, fit$pl, rparg$catchType, rparg$logCustomSel) }
    ## Apply to Fsequence for grid search
    .na2 <- function(x) ifelse(is.na(x) | is.nan(x), -Inf, x)
    logF <- log(Fsequence)
    prv <- do.call("rbind",lapply(logF, getPR))    
    logYPR <- .na2(prv$logYPR)
    logSPR <- .na2(prv$logSPR)
    logSe <- .na2(prv$logSe)
    logRe <- .na2(prv$logRe)
    logYe <- .na2(prv$logYe)
    logLifeExpectancy <- .na2(prv$logLifeExpectancy)
    logYearsLost <- .na2(prv$logYearsLost)
    logDiscYe <- .na2(prv$logDiscYe)
    dSR0 <- exp(stockRecruitmentModelR(-20, fit$pl$rec_pars, fit$conf$stockRecruitmentModelCode, fit$conf$constRecBreaks)$Gradient_logssb) / exp(-20)

    
    rpType <- rparg$rpType
    if(rpType == 1){ ## MSY
        logF0 <- logF[which.max(logYe)]
        if(logF0 == max(logF))
            stop("The stock does not appear to have a well-defined F~MSY~. Increase the upper bound of Fsequence to try again, or remove MSY from the list.")
    }else if(rpType == 2){ ## MSYRange

    }else if(rpType == 3){ ## Max
        logF0 <- logF[which.max(logSPR)]
        if(logF0 == max(logF))
            stop("The stock does not appear to have a well-defined F~Max~. Increase the upper bound of Fsequence to try again, or remove Max from the list.")
    }else if(rpType == 4){ ## xdYPR
        if(any(rparg$xVal <= 0 | rparg$xVal >= 1))
            stop("Proportion must be strictly between 0 and 1")
        ypr <- exp(logYPR)
        ff <- exp(logF)
        logF0 <- sapply(rparg$xVal, function(x) log(ff)[which.min((diff(ypr)/diff(ff) - x * diff(ypr)[1] / diff(ff)[1])^2)])
    }else if(rpType == 5){ ## xSPR
        if(any(rparg$xVal <= 0 | rparg$xVal >= 1))
            stop("Proportion must be strictly between 0 and 1")
        spr <- exp(logSPR)
        logF0 <- sapply(rparg$xVal, function(x) logF[which.min((spr - x * spr[1])^2)])
    }else if(rpType == 6){ ## xB0
        if(any(rparg$xVal <= 0 | rparg$xVal >= 1))
            stop("Proportion must be strictly between 0 and 1")
        Se <- exp(logSe)
        logF0 <- sapply(rparg$xVal, function(x) logF[which.min((Se - x * Se[1])^2)])  
    }else if(rpType == 7){ ## MYPYLdiv
        Arng <- fit$conf$maxAge - fit$conf$minAge + 1
        v <- logYe - log(1.0 + exp(logYearsLost - log(Arng)))        
        logF0 <- logF[which.max(v)]
        if(logF0 == max(logF))
            stop("The stock does not appear to have a well-defined F~MYPYLdiv~. Increase the upper bound of Fsequence to try again, or remove MYPYLdiv from the list.")
    }else if(rpType == 8){ ## MYPYLprod
        if(any(rparg$xVal < 0))
            stop("MYPYL power must be non-negative.")
        Arng <- fit$conf$maxAge - fit$conf$minAge + 1
        logF0 <- sapply(rparg$xVal, function(x){
            v <- logYe + log(1.0 - exp(x*(logYearsLost - log(Arng))))
            logF[which.max(v)]
        })
        if(any(logF0 == max(logF)))
            stop("The stock does not appear to have a well-defined F~MYPYL~. Increase the upper bound of Fsequence to try again, or remove MYPYL from the list.")
    }else if(rpType == 9){ ## MDY
        logF0 <- logF[which.max(logDiscYe)]
        if(logF0 == max(logF))
            stop("The stock does not appear to have a well-defined F~MDY~. Increase the upper bound of Fsequence to try again, or remove MDY from the list.")
    }else if(rpType == 10){ ## Crash (dSR(0) = 1 / SPR)
        logF0 <- logF[which.min((log(dSR0) - (-logSPR))^2)]
    }else if(rpType == 11){ ## Ext
        logF0 <- logF[which.min(logSe^2)]
    }else{
        stop("Reference point type not implemented yet")
    }    
    rparg$logF0 <- logF0
    return(rparg)
}

.refpointNames <- function(rpType, xVal){
    baseName <- names(.refpointEnum)[match(rpType,.refpointEnum)]
    if(length(xVal) == 0 | rpType %in% c(1,3,7,9,10,11,12))
        return(baseName)
    if(rpType == -1){                   #F=x
        return(sapply(xVal, function(xx) sprintf("F=%s",xx)))
    }else if(rpType == 0){              #baseName - x
        return(gsub("-0$","",sapply(xVal, function(xx) sprintf("%s-%d", baseName,xx))))
    }else if(rpType == 2){ # MSYRange
        return(sapply(xVal, function(xx) sapply(c("Lower","Upper"), function(y) sprintf("%s%s %s",xx, gsub("^x","",baseName),y))))
    }else{  ## Value before baseName
        return(sapply(xVal, function(xx) sprintf("%s%s",xx, gsub("^x","",baseName))))
    }
    return("Unknown")
}

.refpointOutput <- function(ssdr, rpArgs, fit, biasCorrect, aveYearsIn, selYearsIn, Fsequence, referencepoints){
    rwnms <- unlist(lapply(rpArgs,function(x) .refpointNames(x$rpType, x$xVal)))
    rpRename <- sapply(lapply(referencepoints,.refpointParser),function(args)do.call(.refpointNames,args))
    outputOrder <- match(rwnms, rpRename)
    toCI <- function(pattern){
        indx <- 1:2
        if(biasCorrect)
            indx <- 3:4
        tmp <- ssdr[grepl(pattern,rownames(ssdr)),indx,drop=FALSE]
        CI <- exp(tmp %*% cbind(CIL=c(1,-2),CIH=c(1,2)))
        Est <- exp(tmp[,1])
        CI[is.na(Est),] <- NA
        cbind(Estimate = Est, CI)
    }

    ## Tables assume that the first reference point is used for Fsequence
    Ftab <- toCI("referencepoint_[1-9][[:digit:]]*_.+_logF$")[outputOrder,,drop=FALSE]
    Btab <- toCI("referencepoint_[1-9][[:digit:]]*_.+_logSe$")[outputOrder,,drop=FALSE]
    Rtab <- toCI("referencepoint_[1-9][[:digit:]]*_.+_logRe$")[outputOrder,,drop=FALSE]
    Ytab <- toCI("referencepoint_[1-9][[:digit:]]*_.+_logYe$")[outputOrder,,drop=FALSE]
    SPRtab <- toCI("referencepoint_[1-9][[:digit:]]*_.+_logSPR$")[outputOrder,,drop=FALSE]
    YPRtab <- toCI("referencepoint_[1-9][[:digit:]]*_.+_logYPR$")[outputOrder,,drop=FALSE]
    YLtab <- toCI("referencepoint_[1-9][[:digit:]]*_.+_logYearsLost$")[outputOrder,,drop=FALSE]
    LEtab <- toCI("referencepoint_[1-9][[:digit:]]*_.+_logLifeExpectancy$")[outputOrder,,drop=FALSE]
    colnames(Ftab) <- colnames(Btab) <- colnames(Rtab) <- colnames(Ytab) <- colnames(SPRtab) <- colnames(YPRtab) <- colnames(YLtab) <- colnames(LEtab) <- c("Estimate","Low","High")
    rownames(Ftab) <- rownames(Btab) <- rownames(Rtab) <- rownames(Ytab) <- rownames(SPRtab) <- rownames(YPRtab) <- rownames(YLtab) <- rownames(LEtab) <- referencepoints[outputOrder]

    ## First reference point is used for Fsequence
    YPRseq <- toCI("referencepoint_0_FixedF_logYPR$")
    SPRseq <- toCI("referencepoint_0_FixedF_logSPR$")
    Yieldseq <- toCI("referencepoint_0_FixedF_logYe$")
    Bseq <- toCI("referencepoint_0_FixedF_logSe$")
    Rseq <- toCI("referencepoint_0_FixedF_logRe$")
    LLseq <- toCI("referencepoint_0_FixedF_logYearsLost$")
    LEseq <- toCI("referencepoint_0_FixedF_logLifeExpectancy$")

    
    rownames(YPRseq) <- rownames(SPRseq) <- rownames(Yieldseq) <- rownames(Bseq) <- rownames(Rseq) <- rownames(LLseq) <- rownames(LEseq) <- Fsequence
    
    res <- list(tables = list(F = Ftab,
                              Yield = Ytab,
                              YieldPerRecruit = YPRtab,
                              SpawnersPerRecruit = SPRtab,
                              Biomass = Btab,
                              Recruitment = Rtab,
                              LifeExpectancy = LEtab,
                              LifeYearsLost = YLtab
                              ),
                graphs = list(F = Fsequence,
                              Yield = Yieldseq,
                              YieldPerRecruit = YPRseq,
                              SpawnersPerRecruit = SPRseq,
                              Biomass = Bseq,
                              Recruitment = Rseq,
                              YearsLost = LLseq,
                              LifeExpectancy = LEseq),
                ## opt = NA,
                ## ssdr = sdr,
                fbarlabel = substitute(bar(F)[X - Y], list(X = fit$conf$fbarRange[1], Y = fit$conf$fbarRange[2])),
                stochastic = FALSE
                )
    
    attr(res,"aveYears") <-  aveYearsIn
    attr(res,"selYears") <- selYearsIn
    
    attr(res,"fit") <- fit
    class(res) <- "sam_referencepoints"
    res
}


##' Function to calculate reference points for the embedded deterministic model of a SAM fit
##'
##' @param fit A fitted SAM model
##' @param referencepoints list of reference points to calculate (See details)
##' @param ... other arguments not used
##' @return 
##' @author Christoffer Moesgaard Albertsen
##' @export
deterministicReferencepoints <- function(fit,
                                         referencepoints,
                                         ...){
    UseMethod("deterministicReferencepoints")
}

##' Function to calculate reference points for the embedded deterministic model of a SAM fit
##'
##' @param fit A fitted SAM model
##' @param referencepoints list of reference points to calculate (See details)
##' @param catchType landing, catch, or discard
##' @param nYears Number of years in per-recruit calculations
##' @param Fsequence Sequence of F values for plotting and starting values
##' @param aveYears Years to average over for biological input
##' @param selYears Years to average over for selectivity
##' @param biasCorrect Should bias correction be used in \link[TMB]{sdreport}?
##' @param newton.control Control arguments passed to the newton optimizer (See \link[TMB]{newton})
##' @param ... other arguments not used
##' @return 
##' @author Christoffer Moesgaard Albertsen
##' @export
deterministicReferencepoints.sam <- function(fit,
                                             referencepoints,
                                             catchType = "catch",
                                             nYears = 300,
                                             Fsequence = seq(0,2, len = 50),
                                             aveYears = max(fit$data$years)+(-9:0),
                                             selYears = max(fit$data$years),
                                             biasCorrect = FALSE,
                                             newton.control = list(),
                                             ...){
    if(!all(diff(Fsequence) > 0) || !all(Fsequence >= 0))
        stop("Values of Fsequence must be positive and increasing.")
    if(!isTRUE(all.equal(Fsequence[1],0, check.attributes = FALSE, use.names = FALSE)))
        warning("The first value of Fsequence should be 0.")

    catchType <- pmatch(catchType,c("catch","landing","discard"))
    if(is.na(catchType))
        stop("Invalid catch type")

     aveYearsIn <- aveYears
    aveYears <- match(aveYears, fit$data$years) - 1
    if(any(is.na(aveYears)))
        stop("aveYears has years without data.")

    selYearsIn <- selYears
    selYears <- match(selYears, fit$data$years) - 1
    if(any(is.na(selYears)))
        stop("selYears has years without data.")

    ## Parse input reference points
    rpArgs <- Reduce(.refpointMerger,
                     lapply(referencepoints, .refpointParser, nYears = nYears, aveYears = aveYears, selYears = selYears, logCustomSel = numeric(0), catchType = catchType - 1),
                     list())

    ## Add starting values    
    rpArgs <- lapply(rpArgs, .refpointStartingValue, fit = fit, Fsequence = Fsequence)
    ## Add Fsequence for plotting
    rp0 <- list(rpType = -1,
                xVal = log(Fsequence),
                nYears = nYears,
                aveYears = aveYears,
                selYears = selYears,
                logCustomSel = numeric(0),
                catchType = catchType - 1,
                logF0 = log(Fsequence))

    ## Make list for TMB
    obj0 <- fit$obj
    argsIn <- as.list(obj0$env)[methods::formalArgs(TMB::MakeADFun)[methods::formalArgs(TMB::MakeADFun) != "..."]]
    argsIn$silent <- fit$obj$env$silent
    argsIn$parameters <- fit$pl
    argsIn$random <- unique(names(obj0$env$par[obj0$env$random]))

    argsIn$data$referencepoints <- c(list(rp0), rpArgs)
    attr(argsIn$data$referencepoints,"newton_config") <- newton.control
    args <- argsIn

     objSDR <- do.call(TMB::MakeADFun, args)
    objSDR$fn(fit$opt$par)
    sdr <- TMB::sdreport(objSDR, objSDR$par, fit$opt$he,
                         bias.correct= biasCorrect,
                         skip.delta.method = biasCorrect,
                         bias.correct.control = list(sd = TRUE,
                                                     split = objSDR$env$ADreportIndex()[grepl("referencepoint_[[:digit:]]+_.+",names(objSDR$env$ADreportIndex()))]
                                                     ))
    ssdr <- summary(sdr)
    ## Make tables
    res <- .refpointOutput(ssdr,rpArgs, fit, biasCorrect, aveYearsIn, selYearsIn, Fsequence, referencepoints)
    return(res)
}

