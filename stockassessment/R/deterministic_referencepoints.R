

recruitmentProperties <- function(fit){
    r <- .Call(C_recruitmentProperties,fit$obj$env$data,fit$pl)
    r$ARpars <- numeric(0)
    if(length(fit$pl$rec_transphi) > 0)
        r$ARpars <- .Call(C_logitroots2ARpar, fit$pl$rec_transphi)
    r
}





.refpointEnum <- c(None = -99,
                   FixedSSB = -2,
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



.logTotalF <- function(logF, conf){
    t(log(apply(conf$keyLogFsta+1,2,function(i) colSums(exp(logF)[i[i>0],,drop=FALSE]))))
}
.logFbar <- function(logF, selYears, conf){
    frng <- conf$fbarRange[1]:conf$fbarRange[2]
    totalF <- exp(.logTotalF(logF,conf))[,selYears+1,drop=FALSE]
    log(sum(colMeans(totalF[frng-conf$minAge+1,,drop=FALSE])))
}

.logFtoSel <- function(logF, selYears, conf){
    logSel <- log(rowMeans(exp(logF[,selYears + 1,drop=FALSE])))
    logFbar <- .logFbar(logF, selYears, conf)
    logSel - logFbar
}

.perRecruitR <- function(logf, fit, nYears, aveYears, selYears, pl = fit$pl, ct = 0, logCustomSel = numeric(0)){
    if(length(logCustomSel) > 0){
        sel <- exp(logCustomSel)
    }else{
        sel <- exp(.logFtoSel(pl$logF, selYears, fit$conf))
    }
    as.data.frame(.Call(C_perRecruitR,
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
    }else if(grepl("^B=\\.[[:digit:]]+$",x)){
        x <- gsub("=\\.","=0.",x)
    }else if(grepl("^B=[[:digit:]]+\\.$",x)){
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
                         FixedSSB = "^SSB=([[:digit:]]+\\.)?([[:digit:]]+)$", # E.g., "B=100000"
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
                         FixedSSB = c("(SSB=)(([[:digit:]]+\\.)?([[:digit:]]+))","\\2"),
                         FixedF = c("(F=)(([[:digit:]]+\\.)?([[:digit:]]+))","\\2"),
                         StatusQuo = c("(StatusQuo)(-?)(([[:digit:]]+)?)","\\3"),
                         MSY = NA,
                         MSYRange = c("(0\\.[[:digit:]]+)(MSYRange)","\\1"),
                         Max = NA,
                         xdYPR = c("(0\\.[[:digit:]]+)(dYPR)","\\1"),
                         xSPR = c("(0\\.[[:digit:]]+)(SPR)","\\1"),
                         xB0 =  c("(0\\.[[:digit:]]+)(B0)","\\1"),
                         MYPYLdiv = NA,
                         MYPYL = c("([[:digit:]]+\\.[[:digit:]]+)(MYPYL)","\\1"),
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


.refpointCheckRecruitment <- function(rparg,fit){
    srmc <- fit$conf$stockRecruitmentModelCode
    recProp <- recruitmentProperties(fit)

    ## Check validity of SRR-refpoint combination
    srmc <- fit$conf$stockRecruitmentModelCode
    recProp <- recruitmentProperties(fit)

    ## FixedSSB = -2          # hasEquilibrium
    ## FixedF = -1          # hasEquilibrium
    ## MSY = 1              # hasEquilibrium
    ## MSYRange = 2         # hasEquilibrium
    ## Max = 3              # all
    ## xdYPR = 4            # all
    ## xSPR = 5,            # all
    ## xB0 = 6,             # hasEquilibrium
    ## MYPYLdiv = 7,        # hasEquilibrium    
    ## MYPYL = 8,           # hasEquilibrium
    ## MDY = 9,             # hasEquilibrium
    ## Crash = 10,          # isCompensatory
    ## Ext = 11             # hasFinitedMaxGradient
    ## Lim = 12             # none

    invalidRP <- c(12)
    if(!recProp$hasEquilibrium)
        invalidRP <- sort(unique(c(invalidRP,-2,-1, 1,2,6,7,8,9)))
    if(!recProp$isCompensatory)
        invalidRP <- sort(unique(c(invalidRP,10)))
    if(!recProp$hasFiniteMaxGradient)
        invalidRP <- sort(unique(c(invalidRP,10,11)))
    srName <- recProp$name
    
    i <- match(rparg$rpType, invalidRP)
    if(length(i) > 0 && !is.na(i))
        stop(sprintf("%s cannot be estimated for %s recruitment",names(.refpointEnum)[match(invalidRP[i],.refpointEnum)],srName))
    return(NULL)
    
}

.refpointStartingValue <- function(rparg, fit, Fsequence, fay = faytable(fit), fbar = fbartable(fit)[,1], checkValidity = TRUE, stochasticType = 0, q = 0){
    .refpointCheckRecruitment(rparg,fit)
    
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
    if(stochasticType == 0){ # Deterministic
        getPR <- function(logf){ .perRecruitR(logf,fit,rparg$nYears, rparg$aveYears, rparg$selYears, fit$pl, rparg$catchType, rparg$logCustomSel) }
    }else if(stochasticType == 1){ # Stochastic Median
        getPR <- function(logf){
            v <- .perRecruitSR_Calc(logf,fit,rparg$nYears, rparg$aveYears, rparg$selYears, fit$pl, rparg$catchType, rparg$logCustomSel)
            list(logFbar = v$E_logFbar,
                 logYPR = v$E_logYPR,
                 logSPR = v$E_logSPR,
                 logSe = v$E_logSe,
                 logRe = v$E_logRe,
                 dSR0 = NA_real_,
                 logLifeExpectancy = v$E_logLifeExpectancy,
                 logYearsLost = v$E_logYearsLost,
                 logDiscYPR = NA_real_,
                 logDiscYe = NA_real_)
        }
    }else if(stochasticType == 2){ # Stochastic Mean
        getPR <- function(logf){
            v <- .perRecruitSR_Calc(logf,fit,rparg$nYears, rparg$aveYears, rparg$selYears, fit$pl, rparg$catchType, rparg$logCustomSel)
            list(logFbar = v$E_logFbar + 0.5 * v$V_logFbar,
                 logYPR = v$E_logYPR + 0.5 * v$V_logYPR,
                 logSPR = v$E_logSPR + 0.5 * v$V_logSPR,
                 logSe = v$E_logSe + 0.5 * v$V_logSe,
                 logRe = v$E_logRe + 0.5 * v$V_logRe,
                 dSR0 = NA_real_,
                 logLifeExpectancy = v$E_logLifeExpectancy + 0.5 * v$V_logLifeExpectancy,
                 logYearsLost = v$E_logYearsLost + 0.5 * v$V_logYearsLost,
                 logDiscYPR = NA_real_,
                 logDiscYe = NA_real_)
        }
    }else if(stochasticType == 3){ # Stochastic Mode
        getPR <- function(logf){
            v <- .perRecruitSR_Calc(logf,fit,rparg$nYears, rparg$aveYears, rparg$selYears, fit$pl, rparg$catchType, rparg$logCustomSel)
            list(logFbar = v$E_logFbar - v$V_logFbar,
                 logYPR = v$E_logYPR - v$V_logYPR,
                 logSPR = v$E_logSPR - v$V_logSPR,
                 logSe = v$E_logSe - v$V_logSe,
                 logRe = v$E_logRe - v$V_logRe,
                 dSR0 = NA_real_,
                 logLifeExpectancy = v$E_logLifeExpectancy - v$V_logLifeExpectancy,
                 logYearsLost = v$E_logYearsLost - v$V_logYearsLost,
                 logDiscYPR = NA_real_,
                 logDiscYe = NA_real_)
        }
    }else if(stochasticType == 4){ # Stochastic Quantile
        getPR <- function(logf){
            v <- .perRecruitSR_Calc(logf,fit,rparg$nYears, rparg$aveYears, rparg$selYears, fit$pl, rparg$catchType, rparg$logCustomSel)
            list(logFbar = qnorm(q, v$E_logFbar, v$V_logFbar),
                 logYPR = qnorm(q, v$E_logYPR, v$V_logYPR),
                 logSPR = qnorm(q, v$E_logSPR, v$V_logSPR),
                 logSe = qnorm(q, v$E_logSe, v$V_logSe),
                 logRe = qnorm(q, v$E_logRe, v$V_logRe),
                 dSR0 = NA_real_,
                 logLifeExpectancy = qnorm(q, v$E_logLifeExpectancy, v$V_logLifeExpectancy),
                 logYearsLost = qnorm(q, v$E_logYearsLost, v$V_logYearsLost),
                 logDiscYPR = NA_real_,
                 logDiscYe = NA_real_)
        }
    }
    
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
    if(rpType == -2){ ## FixedSSB
        if(any(log(rparg$xVal) > max(logSe) | log(rparg$xVal) < min(logSe)))
            stop("The stock does not appear to have a well-defined F corresponding to the requested B. Increase the upper bound of Fsequence to try again, or remove SSB=. from the list.")
        logFx <- sapply(rparg$xVal, function(x){
            logF[which.min((logSe - log(x))^2)]
        })
        logF0 <- as.vector(logFx)
    }else if(rpType == 1){ ## MSY
        logF0 <- logF[which.max(logYe)]
        if(checkValidity && logF0 == max(logF))
            stop("The stock does not appear to have a well-defined F~MSY~. Increase the upper bound of Fsequence to try again, or remove MSY from the list.")
    }else if(rpType == 2){ ## MSYRange
        logFmsy <- logF[which.max(logYe)]
        YeMSY <- max(logYe)
        logFx <- sapply(rparg$xVal, function(x){
            ff <- logF[logF < logFmsy]
            yy <- logYe[logF < logFmsy]
            lfLow <- ff[which.min((yy - (YeMSY+log(x)))^2)]
            ff <- logF[logF > logFmsy]
            yy <- logYe[logF > logFmsy]
            if(length(yy) == 0){
                lfUp <- logFmsy + 0.1
            }else{
                lfUp <- ff[which.min((yy - (YeMSY+log(x)))^2)]
            }
            c(-log(logFmsy-lfLow), log(lfUp-logFmsy))
        })
        logF0 <- c(logFmsy, as.vector(logFx))
    }else if(rpType == 3){ ## Max
        logF0 <- logF[which.max(logYPR)]
        if(checkValidity && logF0 == max(logF))
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
        if(checkValidity && logF0 == max(logF))
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
            stop(checkValidity && "The stock does not appear to have a well-defined F~MYPYL~. Increase the upper bound of Fsequence to try again, or remove MYPYL from the list.")
    }else if(rpType == 9){ ## MDY
        logF0 <- logF[which.max(logDiscYe)]
        if(logF0 == max(logF))
            stop("The stock does not appear to have a well-defined F~MDY~. Increase the upper bound of Fsequence to try again, or remove MDY from the list.")
    }else if(rpType == 10){ ## Crash (dSR(0) = 1 / SPR)
        logF0 <- logF[which.min((log(dSR0) - (-logSPR))^2)]
    }else if(rpType == 11){ ## Ext
        ## If Hard crash, we want a value before
        lf2 <- logF[logSe > -5]
        ls2 <- logSe[logSe > -5]
        logF0 <- lf2[which.min(ls2^2)] - 0.01
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
        return(as.vector(sapply(xVal, function(xx) as.vector(sapply(c("(Lower)","(Upper)"), function(y) sprintf("%s%s %s",xx, gsub("^x","",baseName),y))))))
    }else{  ## Value before baseName
        return(sapply(xVal, function(xx) sprintf("%s%s",xx, gsub("^x","",baseName))))
    }
    return("Unknown")
}



.refpointObjective <- function(rpArgs,fit){
    makeOne <- function(rp){
        if(rp$stochasticType == 0){     #Deterministic
            getPR <- function(logF, pl = fit$pl){
                as.list(.perRecruitR(logF, fit, rp$nYears, rp$aveYears, rp$selYears, pl, ct = rp$catchType))
            }
        }else if(rp$stochasticType ==  1){ #Stochastic MEDIAN
            getPR <- function(logF, pl = fit$pl){
                if(is.null(rp$logN0) || length(rp$logN0) == 0)
                    rp$logN0 <- pl$logN[,ncol(pl$logN)]
                r <- .perRecruitSR_Calc(logF, fit, rp$nYears, rp$aveYears, rp$selYears, pl, ct = rp$catchType, logNinit = rp$logN0, DT = ifelse(is.null(rp$DT) || is.na(rp$DT),0,rp$DT))
                list(logF = r$E_logFbar,
                     logYPR = r$E_logYPR,
                     logSPR = r$E_logSPR,
                     logSe = r$E_logSe,
                     logRe = r$E_logRe,
                     logYe = r$E_logYe,
                     dSR0 = NA_real_,
                     logLifeExpectancy = r$E_logLifeExpectancy,
                     logYearsLost = r$E_logYearsLost,
                     logDiscYe = NA_real_,
                     logDiscYPR = NA_real_
                     )
            }
        }else if(rp$stochasticType ==  2){ #Stochastic MEAN
            getPR <- function(logF, pl = fit$pl){
                if(is.null(rp$logN0) || length(rp$logN0) == 0)
                    rp$logN0 <- pl$logN[,ncol(pl$logN)]
                r <- .perRecruitSR_Calc(logF, fit, rp$nYears, rp$aveYears, rp$selYears, pl, ct = rp$catchType, logNinit = rp$logN0, DT = ifelse(is.null(rp$DT) || is.na(rp$DT),0,rp$DT))
                list(logF = r$E_logFbar + 0.5 * r$V_logFbar,
                     logYPR = r$E_logYPR + 0.5 * r$V_logYPR,
                     logSPR = r$E_logSPR + 0.5 * r$V_logSPR,
                     logSe = r$E_logSe + 0.5 * r$V_logSe,
                     logRe = r$E_logRe + 0.5 * r$V_logRe,
                     logYe = r$E_logYe + 0.5 * r$V_logYe,
                     dSR0 = NA_real_,
                     logLifeExpectancy = r$E_logLifeExpectancy + 0.5 * r$V_logLifeExpectancy,
                     logYearsLost = r$E_logYearsLost + 0.5 * r$V_logYearsLost,
                     logDiscYe = NA_real_,
                     logDiscYPR = NA_real_
                     )
            }
        }else if(rp$stochasticType ==  3){ #Stochastic MODE
            getPR <- function(logF, pl = fit$pl){
                if(is.null(rp$logN0) || length(rp$logN0) == 0)
                    rp$logN0 <- pl$logN[,ncol(pl$logN)]
                r <- .perRecruitSR_Calc(logF, fit, rp$nYears, rp$aveYears, rp$selYears, pl, ct = rp$catchType, logNinit = rp$logN0, DT = ifelse(is.null(rp$DT) || is.na(rp$DT),0,rp$DT))
                list(logF = r$E_logFbar - r$V_logFbar,
                     logYPR = r$E_logYPR - r$V_logYPR,
                     logSPR = r$E_logSPR - r$V_logSPR,
                     logSe = r$E_logSe - r$V_logSe,
                     logRe = r$E_logRe - r$V_logRe,
                     logYe = r$E_logYe - r$V_logYe,
                     dSR0 = NA_real_,
                     logLifeExpectancy = r$E_logLifeExpectancy - r$V_logLifeExpectancy,
                     logYearsLost = r$E_logYearsLost - r$V_logYearsLost,
                     logDiscYe = NA_real_,
                     logDiscYPR = NA_real_
                     )
            }
        }else if(rp$stochasticType ==  4){ #Stochastic QUANTILE
            getPR <- function(logF, pl = fit$pl){
                if(is.null(rp$logN0) || length(rp$logN0) == 0)
                    rp$logN0 <- pl$logN[,ncol(pl$logN)]
                r <- .perRecruitSR_Calc(logF, fit, rp$nYears, rp$aveYears, rp$selYears, pl, ct = rp$catchType, logNinit = rp$logN0, DT = ifelse(is.null(rp$DT) || is.na(rp$DT),0,rp$DT))
                q <- ifelse(is.null(rp$q), 0.5, ifelse(is.na(rp$q), 0.5, rp$q))
                list(logF = qnorm(q, r$E_logFbar, r$V_logFbar),
                     logYPR = qnorm(q, r$E_logYPR, r$V_logYPR),
                     logSPR = qnorm(q, r$E_logSPR, r$V_logSPR),
                     logSe = qnorm(q, r$E_logSe, r$V_logSe),
                     logRe = qnorm(q, r$E_logRe, r$V_logRe),
                     logYe = qnorm(q, r$E_logYe, r$V_logYe),
                     dSR0 = NA_real_,
                     logLifeExpectancy = qnorm(q, r$E_logLifeExpectancy, r$V_logLifeExpectancy),
                     logYearsLost = qnorm(q, r$E_logYearsLost, r$V_logYearsLost),
                     logDiscYe = NA_real_,
                     logDiscYPR = NA_real_
                     )
            }
        }else{
            stop("Unknown equilibrium type")
        }
        par2logF <- identity
        if(rp$rpType == -99){ ## None
            fn <- NA
        }else if(rp$rpType == -2){ ## FixedSSB
            fn <- function(x, pl = fit$pl){
                logF <- x[1]
                r <- getPR(logF,pl)
                return((r$logSe - log(rp$xVal))^2)
            }
        }else if(rp$rpType == -1){ ## FixedF (Known F)
            fn <- NA
            attr(fn,"logF") <- function(pl = fit$pl) pmax(rp$logF0,-20)
        }else if(rp$rpType == 0){ ## StatusQuo (Known F)
            fn <- NA
            attr(fn,"logF") <- function(pl = fit$pl) sapply(rp$xVal, function(y).logFbar(pl$logF,ncol(pl$logF)-y-1,fit$conf))
        }else if(rp$rpType == 1){ ## MSY
            fn <- function(x, pl = fit$pl){
                logF <- x[1]
                r <- getPR(logF, pl)
                return(-r$logYe)
            }
        }else if(rp$rpType == 2){ ## MSYRange
            par2logF <- function(x){
                if(!(length(x)%%2 == 1))
                    stop("In reference point MSYrange, length of F must be odd.")
                r <- numeric(length(x)-1)
                for(i in 2:length(x)){
                    if(((i-1)%%2)==1){
                        r[i-1] <- x[1] - exp(x[i])
                    }else{
                        r[i-1] <- x[1] + exp(x[i])
                    }
                }
                r
            }
            fn <- function(x, pl = fit$pl){
                logF <- par2logF(x)
                ## MSY
                r <- getPR(x[1], pl)
                k1 <- -r$logYe
                ## Ranges
                k2 <- sapply(1:length(logF), function(i){
                    xvi <- (i-1)%/%2+1
                    r2 <- getPR(logF[i], pl)
                    tmp <- r2$logYe - (log(rp$xVal[xvi]) + r$logYe)
                    tmp^2
                })
                return(k1 + sum(k2))
            }
        }else if(rp$rpType == 3){ ## Max
            fn <- function(x, pl = fit$pl){
                logF <- x[1]
                r <- getPR(logF, pl)
                return(-r$logYPR)
            }
        }else if(rp$rpType == 4){ ## xdYPR
            fn <- function(x, pl = fit$pl){
                if(!(length(x) == length(rp$xVal)))
                    stop("In reference point xdYPR, length of F does not match length of fractions.")
                ## Gradient at 0
                ypr <- function(logf) getPR(logf, pl)$logYPR
                dYPR <- function(logf){
                    h <- 0.001
                    u <- logf
                    v <- -ypr(u + 2 * h) + 8 * ypr(u + h) - 8 * ypr(u - h) + ypr(u - 2 * h)
                    g <- v / (12 * h)
                    g / exp(u)
                }
                dYPR0 <- dYPR(-20)
                logF <- x                
                k <- sapply(seq_along(logF), function(i){
                    v <- dYPR(logF[i])
                    tmp <- v - rp$xVal[i] * dYPR0
                    tmp ^ 2
                })
                return(sum(k))
            }            
        }else if(rp$rpType == 5){ ## xSPR
             fn <- function(x, pl = fit$pl){
                if(!(length(x) == length(rp$xVal)))
                    stop("In reference point xSPR, length of F does not match length of fractions.")
                ## Gradient at 0
                spr <- function(logf) getPR(logf, pl)$logSPR
                SPR0 <- spr(-20)
                logF <- x                
                k <- sapply(seq_along(logF), function(i){
                    v <- spr(logF[i])
                    tmp <- v - rp$xVal[i] * SPR0
                    tmp ^ 2
                })
                return(sum(k))
            }            
        }else if(rp$rpType == 6){ ## xB0
            fn <- function(x, pl = fit$pl){
                if(!(length(x) == length(rp$xVal)))
                    stop("In reference point xB0, length of F does not match length of fractions.")
                ## Gradient at 0
                EquiB <- function(logf) getPR(logf, pl)$logSe
                B0 <- EquiB(-20)
                logF <- x                
                k <- sapply(seq_along(logF), function(i){
                    v <- EquiB(logF[i])
                    tmp <- v - rp$xVal[i] * B0
                    tmp ^ 2
                })
                return(sum(k))
            }     
        }else if(rp$rpType == 7){ ## MYPYLdiv
            fn <- function(x, pl = fit$pl){
                logF <- x[1]
                logAgeRange <- log(fit$conf$maxAge - fit$conf$minAge + 1)
                r <- getPR(logF, pl)
                tmp <- r$logYe - log(1 + exp(r$logYearsLost - logAgeRange))
                return(-tmp)
            }
        }else if(rp$rpType == 8){ ## MYPYLprod
            fn <- function(x, pl = fit$pl){
                if(!(length(x) == length(rp$xVal)))
                    stop("In reference point xB0, length of F does not match length of fractions.")
                logF <- x
                logAgeRange <- log(fit$conf$maxAge - fit$conf$minAge + 1)
                k <- sapply(seq_along(logF), function(i){
                    r <- getPR(logF[i], pl)
                    tmp <- r$logYe + log(1 - exp(rp$xVal[i] * (r$logYearsLost - logAgeRange)))
                    -tmp
                })
                return(sum(k))
            }
        }else if(rp$rpType == 9){ ## MDY
            fn <- function(x, pl = fit$pl){
                logF <- x[1]
                r <- getPR(logF, pl)
                return(-r$logDiscYe)
            }
        }else if(rp$rpType == 10){ ## Crash
            fn <- function(x, pl = fit$pl){
                logF <- x[1]
                r <- getPR(logF, pl)
                logdSR0 <- log(r$dSR0)
                tmp <- logdSR0 - (-r$logSPR)
                return(tmp * tmp)
            }
        }else if(rp$rpType == 11){ ## Ext
            fn <- function(x, pl = fit$pl){
                logF <- x[1]
                r <- getPR(logF, pl)
                logdSR0 <- log(r$dSR0)
                tmp <- logdSR0 - (-r$logSPR)
                return(tmp * tmp)
            }
        }else if(rp$rpType == 12){ ## Lim
            stop("Not implemented yet")
        }else{
            stop("Unknown reference point type")
        }
        list(rpType = rp$rpType, xVal = rp$xVal,
             fn = fn, getPR = getPR, par2logF = par2logF, p0 = rp$logF0)
    }
    lapply(rpArgs, makeOne)
}

    .getJointCovariance <- function(fit){
        ## Adapted from TMB::sdreport
        obj <- fit$obj
        par <- fit$obj$env$last.par.best
        r <- fit$obj$env$random
        hessian.fixed <- fit$opt$he
        hessian.random <- fit$obj$env$spHess(par,random=TRUE)
        f <- obj$env$f
        w <- rep(0, length(par))
        nonr <- setdiff(seq_along(par), r)
        tmp <- f(par, order = 1, type = "ADGrad", keepx=nonr, keepy=r) ## TMBad only !!!
        if(!is.matrix(tmp)) ## Happens if length(r)==1
            tmp <- matrix(tmp, ncol=length(nonr) )
        A <- solve(hessian.random, tmp)
        ##diag.term2 <- rowSums((A %*% Vtheta)*A)
        G <- hessian.random %*% A
        G <- as.matrix(G) ## Avoid Matrix::cbind2('dsCMatrix','dgeMatrix')
        M1 <- cbind2(hessian.random,G)
        M2 <- cbind2(t(G), as.matrix(t(A)%*%G)+hessian.fixed )
        M <- rbind2(M1,M2)
        M <- Matrix::forceSymmetric(M,uplo="L")
        dn <- c(names(par)[r],names(par[-r]))
        dimnames(M) <- list(dn,dn)
        p <- Matrix::invPerm(c(r,(1:length(par))[-r]))
        InvSigma <- M[p,p]
        Matrix::solve(InvSigma)
    }
    
.refpointOptimizer <- function(fit, rpArgs, nsim, incpb = function(x){}){
    objs <- .refpointObjective(rpArgs, fit)
    fitOne <- function(i, pl = fit$pl){
        obj <- objs[[i]]
        if(!is.function(obj$fn)){
            logF <- attr(obj$fn,"logF")(pl)
        }else{
            opt <- nlminb(obj$p0, obj$fn, pl = pl)
            logF <- obj$par2logF(opt$par)
        }
        prv <- as.list(do.call("rbind",lapply(lapply(logF, obj$getPR, pl=pl), as.data.frame)))
        rtab <- unlist(prv)
        names(rtab) <- paste0("referencepoint_",i-1,"_XX_",rep(names(prv),each = length(logF)))
        return(rtab)
    }
    ## Point estimates
    pointEst <- unlist(lapply(seq_along(objs), fitOne))
    incpb("Point estimate")
    ## Uncertainty
    if(nsim > 0){
        Sigma <- .getJointCovariance(fit)
        C <- Matrix::t(Matrix::chol(Sigma))
#### Simulate new parlists
        simPar <- fit$obj$env$last.par.best + C %*% matrix(rnorm(ncol(C) * nsim),ncol=nsim)
        simParList <- apply(simPar,2, function(p)fit$obj$env$parList(par=p), simplify=FALSE)
#### Optimize
        uncEst <- lapply(simParList, function(pl){
            v <- unlist(lapply(seq_along(objs), fitOne, pl = pl))
            incpb("Uncertainty")
            v
        })
        uncEst <- do.call("cbind",uncEst)
    }else{
        uncEst <- NULL
    }
    ## Output
    r <- list(pointEst = pointEst,
         uncertainty = uncEst 
         )
    class(r) <- "refR"
    r
}

.refpointOutput <- function(ssdr, rpArgs, fit, biasCorrect, aveYearsIn, selYearsIn, Fsequence, referencepoints, stochastic = FALSE){
    rwnms <- unlist(lapply(rpArgs,function(x) .refpointNames(x$rpType, x$xVal)))
    rpRename <- unlist(sapply(lapply(referencepoints,.refpointParser),function(args)do.call(.refpointNames,args)))
    rpLabels <- unlist(lapply(referencepoints,function(x){
        v <- .refpointParser(x)
        if(v$rpType == 2)
            return(c(paste(x,c("(Lower)","(Upper)"))))
        return(x)
    }))
    outputOrder <- match(rpRename, rwnms)
    if(is(ssdr,"refR")){
        toCI <- function(pattern){
            ii <- grepl(pattern,names(ssdr$pointEst))
            Est <- exp(ssdr$pointEst[ii])
            if(!is.null(ssdr$uncertainty)){
                CI <- t(exp(apply(ssdr$uncertainty[ii,,drop=FALSE],1,quantile,prob=c(0.025,0.975))))
            }else{
                CI <- matrix(NA_real_,nrow = length(Est), ncol = 2)
            }
            cbind(Estimate = Est,CI)
        }
    }else{
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
    }
    if(length(Fsequence) > 0){
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
        rownames(Ftab) <- rownames(Btab) <- rownames(Rtab) <- rownames(Ytab) <- rownames(SPRtab) <- rownames(YPRtab) <- rownames(YLtab) <- rownames(LEtab) <- rpLabels #[outputOrder]

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
                stochastic = stochastic
                )
    }else{
        ## No Fsequence
        Ftab <- toCI("referencepoint_[[:digit:]]*_.+_logF$")[outputOrder,,drop=FALSE]
        Btab <- toCI("referencepoint_[[:digit:]]*_.+_logSe$")[outputOrder,,drop=FALSE]
        Rtab <- toCI("referencepoint_[[:digit:]]*_.+_logRe$")[outputOrder,,drop=FALSE]
        Ytab <- toCI("referencepoint_[[:digit:]]*_.+_logYe$")[outputOrder,,drop=FALSE]
        SPRtab <- toCI("referencepoint_[[:digit:]]*_.+_logSPR$")[outputOrder,,drop=FALSE]
        YPRtab <- toCI("referencepoint_[[:digit:]]*_.+_logYPR$")[outputOrder,,drop=FALSE]
        YLtab <- toCI("referencepoint_[[:digit:]]*_.+_logYearsLost$")[outputOrder,,drop=FALSE]
        LEtab <- toCI("referencepoint_[[:digit:]]*_.+_logLifeExpectancy$")[outputOrder,,drop=FALSE]
        
        colnames(Ftab) <- colnames(Btab) <- colnames(Rtab) <- colnames(Ytab) <- colnames(SPRtab) <- colnames(YPRtab) <- colnames(YLtab) <- colnames(LEtab) <- c("Estimate","Low","High")
        rownames(Ftab) <- rownames(Btab) <- rownames(Rtab) <- rownames(Ytab) <- rownames(SPRtab) <- rownames(YPRtab) <- rownames(YLtab) <- rownames(LEtab) <- rpLabels #[outputOrder]

  res <- list(tables = list(F = Ftab,
                              Yield = Ytab,
                              YieldPerRecruit = YPRtab,
                              SpawnersPerRecruit = SPRtab,
                              Biomass = Btab,
                              Recruitment = Rtab,
                              LifeExpectancy = LEtab,
                              LifeYearsLost = YLtab
                              ),
                graphs = NULL,
                ## opt = NA,
                ## ssdr = sdr,
                fbarlabel = substitute(bar(F)[X - Y], list(X = fit$conf$fbarRange[1], Y = fit$conf$fbarRange[2])),
                stochastic = stochastic
                )
        
    }


    attr(res,"stochasticType") <- sapply(rpArgs,function(x) x$stochasticType)
    attr(res,"stochasticQ") <- sapply(rpArgs,function(x) x$q)
    attr(res,"nYears") <- sapply(rpArgs,function(x) x$nYears)
    
    attr(res,"aveYears") <-  aveYearsIn
    attr(res,"selYears") <- selYearsIn
    
    attr(res,"fit") <- fit
    class(res) <- "sam_referencepoints"
    res
}


##' Function to calculate reference points for the embedded deterministic model of a SAM fit
##'
##' The function estimates reference points based on deterministic per-recruit calculations with no process variance.
##' The following reference points are implemented:
##' \describe{
##'    \item{F=x}{F fixed to x, e.g., \code{"F=0.3"}}
##'    \item{StatusQuo}{F in the last year of the assessment}
##'    \item{StatusQuo-y}{F in the y years before the last in the assessment, e.g., \code{"StatusQuo-1"}}
##'    \item{MSY}{F that maximizes yield}
##'    \item{0.xMSY}{Fs that gives 0.x*100\% of MSY, e.g., \code{"0.95MSY"}}
##'    \item{Max}{F that maximizes yield per recruit}
##'    \item{0.xdYPR}{F such that the derivative of yield per recruit is 0.x times the derivative at F=0, e.g., \code{"0.1dYPR"}}
##'    \item{0.xSPR}{F such that spawners per recruit is 0.x times spawners per recruit at F=0, e.g., \code{"0.35SPR"}}
##'    \item{0.xB0}{F such that biomass is 0.x times the biomass at F=0, e.g., \code{"0.2B0"}}
##' }
##'
##' @examples
##' \dontrun{
##'   deterministicReferencepoints(fit, c("MSY","0.95MSY","Max","0.35SPR","0.1dYPR","StatusQuo-3"))
##' }
##' @param fit A fitted SAM model
##' @param referencepoints list of reference points to calculate (See details)
##' @param ... other arguments not used
##' @return List of estimated reference points
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
##' @param catchType Type of yield to optimize: landing, catch, or discard
##' @param nYears Number of years in per-recruit calculations
##' @param Fsequence Sequence of F values for plotting and starting values
##' @param aveYears Years to average over for biological input
##' @param selYears Years to average over for selectivity
##' @param biasCorrect Should bias correction be used in \link[TMB]{sdreport}?
##' @param newton.control Control arguments passed to the newton optimizer (See \link[TMB]{newton})
##' @param run Run estimation? If false, a list of arguments to MakeADFun is returned.
##' @param ... other arguments not used
##' @return List of estimated reference points
##' @method deterministicReferencepoints sam
##' @rdname deterministicReferencepoints
##' @export
deterministicReferencepoints.sam <- function(fit,
                                             referencepoints,
                                             catchType = "catch",
                                             nYears = 100,
                                             Fsequence = seq(0,2, len = 50),
                                             aveYears = max(fit$data$years)+(-9:0),
                                             selYears = max(fit$data$years),
                                             biasCorrect = FALSE,
                                             newton.control = list(),
                                             run = TRUE,
                                             equilibriumMethod = c("AD","EC"),
                                             nosim_ci = 200,
                                             ...){

    equilibriumMethod <- match.arg(equilibriumMethod)
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
                     lapply(referencepoints, .refpointParser, nYears = nYears, aveYears = aveYears, selYears = selYears, logCustomSel = numeric(0), catchType = catchType - 1,logN0=numeric(0),stochasticType=0,q=NA_real_),
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
                logF0 = log(Fsequence),
                logN0=numeric(0),
                stochasticType=0,
                q=NA_real_)
    ##rpArgs <- c(list(rp0), rpArgs)

    if(equilibriumMethod == "AD"){
        ## Make list for TMB
        obj0 <- fit$obj
        argsIn <- as.list(obj0$env)[setdiff(methods::formalArgs(TMB::MakeADFun),"...")]
        argsIn$silent <- obj0$env$silent
        argsIn$parameters <- fit$pl
        argsIn$random <- unique(names(obj0$env$par[obj0$env$random]))
        argsIn$data$reportingLevel <- -1

        argsIn$data$referencepoints <- rpArgs
        attr(argsIn$data$referencepoints,"newton_config") <- newton.control
        args <- argsIn

        if(!run) return(args)
        
        objSDR <- do.call(TMB::MakeADFun, args)
        objSDR$fn(fit$opt$par)
        sdr <- TMB::sdreport(objSDR, objSDR$par, fit$opt$he,
                             bias.correct= biasCorrect,
                             skip.delta.method = biasCorrect,
                             bias.correct.control = list(sd = TRUE,
                                                         split = objSDR$env$ADreportIndex()[grepl("referencepoint_[[:digit:]]+_.+",names(objSDR$env$ADreportIndex()))]
                                                         ))
        ssdr <- summary(sdr)
    }else if(equilibriumMethod == "EC"){
        pb <- .SAMpb(min = 0, max = nosim_ci + 1, label="Point estimate")
        incpb <- function(label="") .SAM_setPB(pb, pb$getVal()+1,label)
        ssdr <- .refpointOptimizer(fit, rpArgs, nosim_ci, incpb)
    }
    ## Make tables        
    res <- .refpointOutput(ssdr,rpArgs, fit, biasCorrect, aveYearsIn, selYearsIn,
                           c(), #Fsequence,
                           referencepoints)
    attr(res,"equilibriumMethod") <- equilibriumMethod
       attr(res,"aveYears") <-  aveYearsIn
        attr(res,"selYears") <- selYearsIn
        
        attr(res,"fit") <- fit
    return(res)
}

