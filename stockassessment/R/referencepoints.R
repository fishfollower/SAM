##' Calculate jacobian of a function
##'
##' @param func function
##' @param x parameter values
##' @param ... passed to func
##' @return jacobian matrix
##' @author Christoffer Moesgaard Albertsen
jacobian <- function(func, x, ...){
         r <- .Call("jacobian",
                   function(x)func(x,...),
                   x,
                   globalenv(),
                   30L,
                   0.1 * 10^floor(log10(abs(x))) + 1e-4,##abs(1e-4 * x) + 1e-4 * (abs(x) < 1e-8),
                   1e-12)
        do.call("cbind",r[-1])
}
##' Calculate gradient of a function
##'
##' @param func function
##' @param x parameter values
##' @param ... passed to func
##' @return gradient vector
##' @author Christoffer Moesgaard Albertsen
grad <- function(func, x, ...){
         r <- .Call("jacobian",
                   function(x)func(x,...),
                   x,
                   globalenv(),
                   30L,
                   0.1 * 10^floor(log10(abs(x))) + 1e-4,##abs(1e-4 * x) + 1e-4 * (abs(x) < 1e-8),
                   1e-12)
         v <- do.call("cbind",r[-1])
         if(nrow(v) == 1)
             return(as.vector(v))
         return(diag(v))         
}



##' Add stock-recruitment curve to srplot
##'
##' @param fit Object to show SR-curve for
##' @param CI Add confidence intervals?
##' @param col Color of fitted line
##' @param cicol Color of confidence intervals
##' @param plot Add the curve to a plot?
##' @param ... not used
##' @seealso srplot
##' @author Christoffer Moesgaard Albertsen
##' @export
addRecruitmentCurve <- function(fit,
                                CI = TRUE,
                                col = rgb(0.6,0,0),
                                cicol = rgb(0.6,0,0,0.3),
                                plot = TRUE,
                                ...){
    UseMethod("addRecruitmentCurve")
} 

##' @rdname addRecruitmentCurve
##' @method addRecruitmentCurve sam
##' @export
addRecruitmentCurve.sam <- function(fit,
                    CI = TRUE,
                    col = rgb(0.6,0,0),
                    cicol = rgb(0.6,0,0,0.3),
                    plot = TRUE,
                    ...){
       X <- summary(fit)
       R <- X[, 1]
       S <- X[, 4]
       cf <- fit$sdrep$cov.fixed
       covEst <- cf[rownames(cf) %in% c("rec_pars"), colnames(cf) %in% c("rec_pars"), drop = FALSE]
       m <- fit$obj$env$map$rec_pars
       if(is.null(m)){
           covar <- covEst
       }else{
           covar <- covEst[m,m, drop = FALSE]
           covar[is.na(covar)] <- 0
       }
  
       srfit <- function(ssb){
           v <- .Call("stockRecruitmentModelR",
                      ssb,
                      fit$pl$rec_pars,
                      fit$conf$stockRecruitmentModelCode)
           res <- v$Recruits
           g <- matrix(v$Gradient, 1)
           attr(res,"sd") <- as.vector(sqrt(g %*% covar %*% t(g)))
           return(res)
       }

       ssb <- seq(0, max(S), len = 2000)
       tab <- sapply(ssb, function(x) {
           tmp <- srfit(x)
           sd <- attr(tmp, "sd")
           c(Estimate = as.vector(tmp),
             CIlow = as.vector(tmp - 2 * sd),
             CIhigh = as.vector(tmp + 2 * sd))
       })
       if(plot){
           if(CI)
               polygon(c(ssb, rev(ssb)),
                       c(tab["CIlow",],
                         rev(tab["CIhigh",])),
                       col = cicol,
                       border = NA)
           lines(ssb,tab["Estimate",], col = col, lwd = 3)
       }
       invisible(srfit)
} 


##' Estimating Fmsy
##'
##' Work in progress - do not use
##' @param fit a SAM fit
##' @param nYears Number of years to forecast
##' @param nlminb.control list of control variables for nlminb
##' @param rec.years Numeric vector of years to use (to calculate mean and standard deviation) for recruitment. An empty vector will use the recruitment model.
##' @param processNoiseF Should random walk process noise be used for F?
##' @param ... other arguments
##' @author Christoffer Moesgaard Albertsen
##' @seealso forecast
##' @export
forecastMSY <- function(fit,
                    nYears = 100,
                    nlminb.control = list(eval.max = 2000, iter.max = 2000),
                    rec.years = c(),
                    processNoiseF = FALSE,
                    ...){
    UseMethod("forecastMSY")
}


##' @rdname forecastMSY
##' @method forecastMSY sam
##' @export
forecastMSY.sam <- function(fit,
                    nYears = 500,
                    nlminb.control = list(eval.max = 2000, iter.max = 2000, trace = 1),
                    rec.years = c(),
                    processNoiseF = FALSE,
                    ...){

    argsIn <- forecast(fit,
                       findMSY = rep(1,nYears),
                       rec.years = rec.years,
                       processNoiseF = processNoiseF,
                       ...)
    argsIn$DLL <- "stockassessment"

    args <- argsIn
    args$map$logFScaleMSY <- factor(NA)
    objForecast <- do.call(TMB::MakeADFun, args)
    objForecast$fn()

    ## Get joint precision
    jointPrecision <- TMB::sdreport(objForecast,
                                    fit$opt$par,
                                    solve(fit$sdrep$cov.fixed),
                                    getJointPrecision = TRUE)$jointPrecision


    rp <- c("logFScaleMSY")
    
    ## Find MSY value
    args <- argsIn
    args$parameters$logFScaleMSY <- 0
    args$parameters$implicitFunctionDelta <- 0
    map0 <- args$map
    fix <- setdiff(names(args$parameters), args$random)
    args$map <- lapply(args$parameters[fix], function(x)factor(x*NA))
    args$map <- args$map[names(args$map) != rp]
    args$map$implicitFunctionDelta <- NULL


    objOptim <- do.call(TMB::MakeADFun, args)

    objOptim$fn(objOptim$par)
    ## obj$gr(obj$par)

    fn <- Vectorize(function(x){
        objOptim$fn(c(x,0))
        theta <- objOptim$env$last.par
        objOptim$env$f(theta, order = 1)[1,which(names(objOptim$env$last.par) == "implicitFunctionDelta")]
    })


    opt <- nlminb(objOptim$par[names(objOptim$par) != "implicitFunctionDelta"], fn, control = nlminb.control)

    ## Object to do Delta method (no map, no random, delta = 1)
    args <- argsIn
    args$parameters <- objOptim$env$parList(x = opt$par)
    args$parameters$implicitFunctionDelta <- 1
    args$map$logFScaleMSY <- NULL
    args$map$implicitFunctionDelta <- factor(NA)
    args$random <- NULL
    
    objDelta <- do.call(TMB::MakeADFun, args)

    ## Get Jacobian
    gridx <- which(names(objDelta$par) %in% rp)
    JacAll <- jacobian(function(x) objDelta$gr(x)[gridx], objDelta$par)

    ## Implicit function gradient
    dCdTheta <- -solve(JacAll[,gridx,drop=FALSE]) %*% JacAll[,-gridx,drop=FALSE]
    rownames(dCdTheta) <- rp
    colnames(dCdTheta) <- names(objForecast$env$last.par)


    ## Do delta method
    xtra <- diag(1,length(objForecast$env$last.par.best))
    diag(xtra)[objForecast$env$random] <- 0
    dG <- rbind(xtra[diag(xtra) != 0,,drop = FALSE],dCdTheta)
    covAll <- dG %*% solve(jointPrecision) %*% t(dG)
    covAllOld <- covAll
    i <- 21
    tv <- ((10^(-i))*10^floor(log10(diag(covAll)[gridx])))
    while(tryCatch({solve(covAll);FALSE},error=function(e)TRUE)){
        i <- i-1
        covAll <- covAllOld
        tv <- ((10^(-i))*10^floor(log10(diag(covAll)[gridx])))
        diag(covAll)[gridx] <- diag(covAll)[gridx] + tv
    }

    ## Object to do sdreport (delta = 0)
    args <- argsIn
    ## Remvoe rp from map
    args$map <- args$map[!(names(args$map) %in% rp)]
    ## Use optimized parameters
    args$parameters <- objOptim$env$parList(x = opt$par)
    args$parameters$implicitFunctionDelta <- 0

    objSDR <- do.call(TMB::MakeADFun, args)

    
    sdr <- TMB::sdreport(objSDR, objSDR$par, solve(covAll))
    ssdr <- summary(sdr)

    toCI <- function(what){
        exp(ssdr[rownames(ssdr) %in% what,] %*% cbind(Estimate=c(1,0),Low=c(1,-2),High=c(1,2)))
    }
    
    Fmsy <- toCI("logFMSY")
    ssb <- toCI("logssb")
    yield <- toCI("logCatch")
    rec <- toCI("logR")
    fbar <- toCI("logfbar")

    tab <- rbind(Fmsy, ssb[nrow(ssb),], yield[nrow(yield),], rec[nrow(rec),])
    rownames(tab) <- c("Fmsy", "Bmsy", "Yield", "Rmsy")

    rownames(ssb) <- rownames(yield) <- rownames(rec) <- rownames(fbar) <- min(fit$data$years) + 0:(nrow(ssb)-1)

    return(list(table = tab, timeseries = list(SSB = ssb, Catch = yield, Recruitment = rec, Fbar = fbar), opt = opt, sdr = sdr))

}


##' Estimate reference points
##'
##' Work in progress - do not use
##' @param fit an object to calculate reference points for
##' @param nYears Number of years to use in per-recruit calculations
##' @param Fsequence Sequence of F values used to report per-recruit and equilibrium values
##' @param aveYears Vector of year indices used to calculate average natural mortality, weights, etc. (starting at 0)
##' @param selYears Vector of year indices used to calculate selectivity (starting at 0)
##' @param SPRpercent Vector of x values for F[x * 100\%] reference points. Default is 0.35.
##' @param catchType Catch type used: (total) catch, landings, discard.
##' @param ... not used
##' @return a sam_referencepoints fit
##' @author Christoffer Moesgaard Albertsen
##' @export
referencepoints <- function(fit,
                            nYears,
                            Fsequence,
                            aveYears,
                            selYears,
                            SPRpercent,
                            catchType,
                            ...){
    UseMethod("referencepoints")
}

##' @rdname referencepoints
##' @method referencepoints sam
#' @export
referencepoints.sam <- function(fit,
                                nYears = 100,
                            Fsequence = seq(1e-5,4, len = 200),
                            aveYears = max(fit$data$years)+(-4:0),
                            selYears = max(fit$data$years),
                            SPRpercent = c(0.35),
                            catchType = "catch",
                            ...){

    ## Get joint precision
    jointPrecision <- TMB::sdreport(fit$obj,
                               fit$opt$par,
                               solve(fit$sdrep$cov.fixed),
                               getJointPrecision = TRUE)$jointPrecision

    ## Prepare arguments to calculate reference points (fix parameters and latent variables, delta = 1)
    obj0 <- fit$obj
    argsIn <- as.list(obj0$env)[methods::formalArgs(TMB::MakeADFun)[methods::formalArgs(TMB::MakeADFun) != "..."]]
    argsIn$parameters <- fit$pl
    argsIn$random <- unique(names(obj0$env$par[obj0$env$random]))
    ## Add referencepointSet
    catchType <- pmatch(catchType,c("catch","landing","discard"))
    if(is.na(catchType))
        stop("Invalid catch type")

    aveYears <- match(aveYears, fit$data$years) - 1
    if(any(is.na(aveYears)))
        stop("aveYears has years without data.")

    selYears <- match(selYears, fit$data$years) - 1
    if(any(is.na(selYears)))
        stop("selYears has years without data.")

    
    argsIn$data$referencepoint <- list(nYears = nYears,
                                       aveYears = aveYears,
                                       selYears = selYears,
                                       Fsequence = Fsequence,
                                       xPercent = SPRpercent,
                                       catchType = catchType-1
                                       )

    args <- argsIn
    ## Remove random
    args$random <- NULL
    ## Remove referencepoint parameters from map
    map0 <- args$map
    fix <- names(args$parameters) ## setdiff(names(args$parameters), args$random)
    args$map <- lapply(args$parameters[fix], function(x)factor(x*NA))

    if(fit$conf$stockRecruitmentModelCode %in% c(0,3)){ # RW, constant mean
        rp <- c("logScaleFmax",
                "logScaleF01",
                "logScaleFxPercent")
    }else if(fit$conf$stockRecruitmentModelCode %in% c(61,63)){ # Hockey-sticks
        rp <- c("logScaleFmsy",
                "logScaleFmax",
                "logScaleF01",
                "logScaleFcrash",
                ## "logScaleFext",
                "logScaleFxPercent",
                "logScaleFlim")
    }else if(fit$conf$stockRecruitmentModelCode %in% c(62)){ # AR
        rp <- c("logScaleFmsy",
                "logScaleFmax",
                "logScaleF01",
                "logScaleFxPercent")
    }else if(fit$conf$stockRecruitmentModelCode %in% c(64)){ # Pow CMP
        rp <- c("logScaleFmsy",
                "logScaleFmax",
                "logScaleF01",
                "logScaleFxPercent"
                )
    }else if(fit$conf$stockRecruitmentModelCode %in% c(65)){ # Pow Non-CMP
        rp <- c("logScaleFmsy",
                "logScaleFmax",
                "logScaleF01",
                "logScaleFxPercent"
                )
    }else if(fit$conf$stockRecruitmentModelCode %in% c(68,69) && fit$pl$rec_par[3] > 0){ ## depensatory recruitment; Fcrash does not work.
        rp <- c("logScaleFmsy",
                "logScaleFmax",
                "logScaleF01",
                ## "logScaleFext",
                "logScaleFxPercent"
                )
    }else{
        rp <- c("logScaleFmsy",
                "logScaleFmax",
                "logScaleF01",
                "logScaleFcrash",
                ## "logScaleFext",
                "logScaleFxPercent"
                )
    }
    ## Referencepoints to estimate
    args$map <- args$map[-which(names(args$map) %in% rp)]

    args$parameters$logScaleFmsy <- -2
    args$parameters$logScaleF01 <- -2
    args$parameters$logScaleFmax <- -2
    args$parameters$logScaleFcrash <- -2
    args$parameters$logScaleFext <- -2
    args$parameters$logScaleFxPercent <- rep(-2, length(SPRpercent))
    args$parameters$logScaleFlim <- -2
    args$parameters$implicitFunctionDelta <- 1

    objOptim <- do.call(TMB::MakeADFun, args)

    ## Take inital look at YPR / SPR to determine if Fmax makes sense
    rep <- objOptim$report()
    tryAgain <- FALSE
    if(which.max(rep$logYPR[is.finite(rep$logYPR)]) == length(rep$logYPR[is.finite(rep$logYPR)]) && any(rp %in% "logScaleFmax")){
        warning("The stock does not appear to have a well-defined Fmax. Fmax will not be estimated. Increase the upper bound of Fsequence to try again.")
        rp <- rp[-which(rp %in% "logScaleFmax")]
        args$map$logScaleFmax <- factor(NA)
        tryAgain <- TRUE
    }

      if(min(rep$logSe[is.finite(rep$logSe)]) > -10 && any(rp %in% "logScaleFcrash")){
        warning("The stock does not appear to have a well-defined Fcrash. Fmax will not be estimated. Increase the upper bound of Fsequence to try again.")
        rp <- rp[-which(rp %in% "logScaleFcrash")]
        args$map$logScaleFcrash <- factor(NA)
        tryAgain <- TRUE
    }

    if(which.max(rep$logYe[is.finite(rep$logYe)]) == length(rep$logYe[is.finite(rep$logYe)]) && any(rp %in% "logScaleFmsy")){
        warning("The stock does not appear to have a well-defined Fmsy. Fmsy will not be estimated. Increase the upper bound of Fsequence to try again.")
        rp <- rp[-which(rp %in% "logScaleFmsy")]
        args$map$logScaleFmsy <- factor(NA)
        tryAgain <- TRUE
    }


    if(tryAgain)
        objOptim <- do.call(TMB::MakeADFun, args)

    opt <- nlminb(objOptim$par, objOptim$fn, objOptim$gr)#, objOptim$he)
    ii <- 0
    while(max(abs(objOptim$gr(opt$par))) > 1e-4 && ii < 20){
        g <- as.numeric( objOptim$gr(opt$par) )
        h <- objOptim$he(opt$par)
        opt$par <- opt$par - solve(h, g)
        opt$objective <- objOptim$fn(opt$par)
        ii <- ii + 1
    }
    ## Object to do Delta method (nothing mapped (that's not mapped in fit$obj, nothing random, delta = 1)
    args <- argsIn
    ## Remove random
    args$random <- NULL
    ## Remvoe rp from map
    args$map <- args$map[!(names(args$map) %in% rp)]
    ## Use optimized parameters
    args$parameters <- objOptim$env$parList(x = opt$par)

    objDelta <- do.call(TMB::MakeADFun, args)

    ## Get Jacobian
    gridx <- which(names(objDelta$par) %in% rp)
    JacAll <- jacobian(function(x) objDelta$gr(x)[gridx], objDelta$par)

    ## Implicit function gradient
    dCdTheta <- -solve(JacAll[,gridx,drop=FALSE]) %*% JacAll[,-gridx,drop=FALSE]
    rownames(dCdTheta) <- gsub("^logScale","",names(objDelta$par)[gridx])
    colnames(dCdTheta) <- names(fit$obj$env$last.par)

    ## Do delta method
    xtra <- diag(1,length(fit$obj$env$last.par.best))
    diag(xtra)[fit$obj$env$random] <- 0
    dG <- rbind(xtra[diag(xtra) != 0,,drop = FALSE],dCdTheta)
    covAll <- dG %*% solve(jointPrecision) %*% t(dG)
    covAllOld <- covAll
    i <- 21
    tv <- ((10^(-i))*10^floor(log10(diag(covAll)[gridx])))
    while(tryCatch({solve(covAll);FALSE},error=function(e)TRUE)){
        i <- i-1
        covAll <- covAllOld
        tv <- ((10^(-i))*10^floor(log10(diag(covAll)[gridx])))
        diag(covAll)[gridx] <- diag(covAll)[gridx] + tv
    }

    ## Object to do sdreport (delta = 0)
    args <- argsIn
    ## Remvoe rp from map
    args$map <- args$map[!(names(args$map) %in% rp)]
    ## Use optimized parameters
    args$parameters <- objOptim$env$parList(x = opt$par)
    args$parameters$implicitFunctionDelta <- 0

    objSDR <- do.call(TMB::MakeADFun, args)

    
    sdr <- TMB::sdreport(objSDR, objSDR$par, solve(covAll))
    ssdr <- summary(sdr)

    toCI <- function(what){
        exp(ssdr[rownames(ssdr) == what,,drop=FALSE] %*% cbind(Estimate=c(1,0),CIL=c(1,-2),CIH=c(1,2)))
    }

    Ftab <- do.call("rbind",sapply(unique(rownames(ssdr)[grepl("referencepoint.logF",rownames(ssdr))]),
                     toCI, simplify = FALSE))
    Btab <- do.call("rbind",sapply(unique(rownames(ssdr)[grepl("referencepoint.logB",rownames(ssdr))]),
                                   toCI, simplify = FALSE))
    Rtab <- do.call("rbind",sapply(unique(rownames(ssdr)[grepl("referencepoint.logR",rownames(ssdr))]),
                                   toCI, simplify = FALSE))
    Ytab <- do.call("rbind",sapply(unique(rownames(ssdr)[grepl("referencepoint.logY",rownames(ssdr)) & !grepl("referencepoint.logYPR",rownames(ssdr))]),
                     toCI, simplify = FALSE))
    SPRtab <- do.call("rbind",sapply(unique(rownames(ssdr)[grepl("referencepoint.logSPR",rownames(ssdr))]),
                                   toCI, simplify = FALSE))
    YPRtab <- do.call("rbind",sapply(unique(rownames(ssdr)[grepl("referencepoint.logYPR",rownames(ssdr))]),
                                   toCI, simplify = FALSE))
    colnames(Ftab) <- colnames(Btab) <- colnames(Ytab) <- colnames(SPRtab) <- colnames(YPRtab) <- c("Estimate","Low","High")

    toRowNames <- Vectorize(function(x){
        switch(gsub("^referencepoint.logF","",x),
               "sq"="Status quo",
               "0"="Zero catch",
               "msy"="MSY",
               "max"="Max",
               "01"="0.1",
               "crash"="Crash",
               "ext"="Ext",
               "xPercent"=NA,
               "lim"="lim",
               x
        )               
    })
    rn <- toRowNames(rownames(Ftab))
    rn[is.na(rn)] <- sapply(SPRpercent,function(x)sprintf("%s%%",x * 100))    
    rownames(Ftab) <- rownames(Btab) <- rownames(Rtab) <- rownames(Ytab) <- rownames(SPRtab) <- rownames(YPRtab) <- rn

    ## Ftab["Ext",c("Low","High")] <- NA
    ## Btab["Ext",c("Low","High")] <- NA
    ## Rtab["Ext",c("Low","High")] <- NA
    ## Ytab["Ext",c("Low","High")] <- NA
    ## SPRtab["Ext",c("Low","High")] <- NA
    ## YPRtab["Ext",c("Low","High")] <- NA
    
    YPRseq <- toCI("logYPR")
    SPRseq <- toCI("logSPR")
    Yieldseq <- toCI("logYe")
    Bseq <- toCI("logSe")
    Rseq <- toCI("logRe")

    rownames(YPRseq) <- rownames(SPRseq) <- rownames(Yieldseq) <- rownames(Bseq) <- rownames(Rseq) <- argsIn$data$referencepoint$Fsequence
    
    res <- list(tables = list(F = Ftab,
                              B = Btab,
                              Yield = Ytab,
                              Recruitment = Rtab,
                              SpawnersPerRecruit = SPRtab,
                              YieldPerRecruit = YPRtab),
                graphs = list(F = argsIn$data$referencepoint$Fsequence,
                              Yield = Yieldseq,
                              YieldPerRecruit = YPRseq,
                              SpawnersPerRecruit = SPRseq,
                              Biomass = Bseq,
                              Recruitment = Rseq),
                opt = opt,
                sdr = sdr,
                diagonalCorrection = tv
                )
                              


    class(res) <- "sam_referencepoints"
    return(res)
}

   
