##' Calculate jacobian of a function
##'
##' @param func function
##' @param x parameter values
##' @param h step size
##' @param ... passed to func
##' @return jacobian matrix
##' @author Christoffer Moesgaard Albertsen
jacobian <- function(func, x,
                     h = abs(1e-04 * x) + 1e-04 * (abs(x) < sqrt(.Machine$double.eps/7e-07)),
                     ## 0.1 * 10^floor(log10(abs(x))) + 1e-4,
                     ...){
         r <- .Call("jacobian",
                   function(x)func(x,...),
                   x,
                   globalenv(),
                   30L,
                   h,##abs(1e-4 * x) + 1e-4 * (abs(x) < 1e-8),
                   1e-12)
        do.call("cbind",r[-1])
}
##' Calculate gradient of a function
##'
##' @param func function
##' @param x parameter values
##' @param h step size
##' @param ... passed to func
##' @return gradient vector
##' @author Christoffer Moesgaard Albertsen
grad <- function(func, x,
                 h = abs(1e-04 * x) + 1e-04 * (abs(x) < sqrt(.Machine$double.eps/7e-07)),
                 ##0.1 * 10^floor(log10(abs(x))) + 1e-4,
                 ...){
         r <- .Call("jacobian",
                   function(x)func(x,...),
                   x,
                   globalenv(),
                   30L,
                   h,##abs(1e-4 * x) + 1e-4 * (abs(x) < 1e-8),
                   1e-12)
         v <- do.call("cbind",r[-1])
         if(nrow(v) == 1)
             return(as.vector(v))
         return(diag(v))         
}

svd_solve <- function(x){
    ss <- svd(x)
    ss$v %*% diag(1/ss$d, length(ss$d), length(ss$d)) %*% t(ss$u)
}


##' Estimating Fmsy
##'
##' @param fit a SAM fit
##' @param nYears Number of years to forecast
##' @param nlminb.control list of control variables for nlminb
##' @param rec.years Numeric vector of years to use (to calculate mean and standard deviation) for recruitment. An empty vector will use the recruitment model.
##' @param ave.years vector of years to average for weights, maturity, M and such. Following ICES guidelines, the default is the last 10 years.
##' @param processNoiseF Should random walk process noise be used for F?
##' @param ... other arguments passed to forecast
##' @author Christoffer Moesgaard Albertsen
##' @seealso \link{forecast} \link{referencepoints}
##' @references
##' Albertsen, C. M. and Trijoulet, V. (2020) Model-based estimates of reference points in an age-based state-space stock assessment model. Fisheries Research, 230, 105618. doi: 10.1016/j.fishres.2020.105618
##' @export
forecastMSY <- function(fit,
                    nYears = 100,
                    nlminb.control = list(eval.max = 2000, iter.max = 2000),
                    rec.years = c(),
                    ave.years = max(fit$data$years)+(-9:0),
                    processNoiseF = FALSE,
                    ...){
    UseMethod("forecastMSY")
}


##' @rdname forecastMSY
##' @method forecastMSY sam
##' @param jacobianHScale Scale step size in jacobian calculation
##' @param nCatchAverageYears Number of years to average catch over for finding MSY
##' @export
forecastMSY.sam <- function(fit,
                            nYears = 100,                            
                            nlminb.control = list(eval.max = 2000, iter.max = 2000, trace = 1),
                            rec.years = c(),
                            ave.years = max(fit$data$years)+(-9:0),
                            processNoiseF = FALSE,
                            jacobianHScale = 0.5,
                            nCatchAverageYears = 20,
                            ...){

    argsIn <- modelforecast(fit,
                       findMSY = rep(1,nYears),
                       rec.years = rec.years,
                       ave.years = ave.years,
                       processNoiseF = processNoiseF,
                       nCatchAverageYears = nCatchAverageYears,
                       ...)
    argsIn$DLL <- "stockassessment"
    argsIn$silent <- fit$obj$env$silent
    
    args <- argsIn
    args$map$logFScaleMSY <- factor(NA)
    objForecast <- do.call(TMB::MakeADFun, args)
    objForecast$fn()

    ## Get joint precision
    jointPrecision <- TMB::sdreport(objForecast,
                                    fit$opt$par,
                                    svd_solve(fit$sdrep$cov.fixed),
                                    getJointPrecision = TRUE)$jointPrecision

    rp <- c("logFScaleMSY")
    
    ## Find MSY value
    args <- argsIn
    args$parameters$logFScaleMSY <- -1
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

    ## Try different values??
    testStart <- c(-2, -1, -0.5, -0.1, 0, 0.1, 0.25, 0.3, 0.5, 1, 2)
    fnTestStart <- sapply(testStart, fn)
    opt <- nlminb(testStart[which(order(fnTestStart) == 1)], fn, control = nlminb.control)
    ## p0 <- objOptim$par[names(objOptim$par) != "implicitFunctionDelta"]
    ## opt <- nlminb(p0, fn, control = nlminb.control)

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
    JacAll <- jacobian(function(x) objDelta$gr(x)[gridx], objDelta$par, h = jacobianHScale * abs(1e-04 * objDelta$par) + 1e-04 * (abs(objDelta$par) < sqrt(.Machine$double.eps/7e-07)))

    ## Implicit function gradient
    dCdTheta <- -svd_solve(JacAll[,gridx,drop=FALSE]) %*% JacAll[,-gridx,drop=FALSE]
    rownames(dCdTheta) <- rp
    colnames(dCdTheta) <- names(objForecast$env$last.par)


    ## Do delta method
    xtra <- diag(1,length(objForecast$env$last.par.best))
    diag(xtra)[objForecast$env$random] <- 0
    dG <- rbind(xtra[diag(xtra) != 0,,drop = FALSE],dCdTheta)
    covAll <- dG %*% svd_solve(jointPrecision) %*% t(dG)
    covAllOld <- covAll
    ## i <- 21
    ## tv <- ((10^(-i))*10^floor(log10(diag(covAll)[gridx])))
    ## while(tryCatch({solve(covAll);FALSE},error=function(e)TRUE)){
    ##     i <- i-1
    ##     covAll <- covAllOld
    ##     tv <- ((10^(-i))*10^floor(log10(diag(covAll)[gridx])))
    ##     diag(covAll)[gridx] <- diag(covAll)[gridx] + tv
    ## }

    ## Object to do sdreport (delta = 0)
    args <- argsIn
    ## Remvoe rp from map
    args$map <- args$map[!(names(args$map) %in% rp)]
    ## Use optimized parameters
    args$parameters <- objOptim$env$parList(x = opt$par)
    args$parameters$implicitFunctionDelta <- 0

    objSDR <- do.call(TMB::MakeADFun, args)

    
    sdr <- TMB::sdreport(objSDR, objSDR$par, svd_solve(covAll))
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

    res <- list(table = tab,
                timeseries = list(SSB = ssb,
                                  Catch = yield,
                                  Recruitment = rec,
                                  Fbar = fbar),
                opt = opt,
                sdr = sdr,
                fbarlabel = substitute(bar(F)[X - Y], list(X = fit$conf$fbarRange[1], Y = fit$conf$fbarRange[2])))
    class(res) <- "sam_forecast_msy"
    res
}


na2 <- function(x,a=0) ifelse(is.na(x) | is.nan(x),a,x)

##' Estimate reference points
##'
##' @param fit an object to calculate reference points for
##' @param nYears Number of years to use in per-recruit calculations
##' @param Fsequence Sequence of F values used to report per-recruit and equilibrium values
##' @param aveYears Vector of year indices used to calculate average natural mortality, weights, etc. Following ICES guidelines, the default is the last 10 years (starting at 0)
##' @param selYears Vector of year indices used to calculate selectivity (starting at 0)
##' @param SPRpercent Vector of x values for F[x * 100\%] reference points. Default is 0.35.
##' @param catchType Catch type used: (total) catch, landings, discard.
##' @param MSYreduction Vector of proportions for MSY ranges. Default is 0.05 giving an MSY range corresponding to no more than a 5\% yield reduction.
##' @param newtonSteps Number of additional Newton steps at the end of the reference point optimization.
##' @param optN N used for numerical optimizers to find equilibrium biomass
##' @param jacobianHScale Scale step size in jacobian calculation
##' @param ... not used
##' @return a sam_referencepoints fit
##' @author Christoffer Moesgaard Albertsen
##' @seealso \link{forecastMSY}
##' @references
##' Albertsen, C. M. and Trijoulet, V. (2020) Model-based estimates of reference points in an age-based state-space stock assessment model. Fisheries Research, 230, 105618. doi: 10.1016/j.fishres.2020.105618
##' @export
referencepoints <- function(fit,
                            nYears,
                            Fsequence,
                            aveYears,
                            selYears,
                            SPRpercent,
                            catchType,
                            MSYreduction,
                            newtonSteps = 3,
                            optN = 100,
                            jacobianHScale = 0.5,
                            ...){
    UseMethod("referencepoints")
}

##' @rdname referencepoints
##' @method referencepoints sam
##' @export
referencepoints.sam <- function(fit,
                                nYears = 100,
                                Fsequence = seq(0,4, len = 200),
                                aveYears = max(fit$data$years)+(-9:0),
                                selYears = max(fit$data$years),
                                SPRpercent = c(0.35),
                                dYPRpercent = c(0.1),
                                B0percent = c(0.2),
                                catchType = "catch",
                                MSYreduction = c(0.05),
                                newtonSteps = 3,
                                optN = 20,
                                jacobianHScale = 0.5,
                                fixRP = c(),
                                RecCorrection = "median",
                                biasCorrect = FALSE,
                                ...){
    if(!all(diff(Fsequence) > 0) || !all(Fsequence >= 0))
        stop("Values of Fsequence must be positive and increasing.")
    if(!isTRUE(all.equal(Fsequence[1],0, check.attributes = FALSE, use.names = FALSE)))
        warning("The first value of Fsequence should be 0.")

    if(!(all(MSYreduction > 0) && all(MSYreduction < 1)))
        stop("MSYreduction must be between 0 and 1.")
    MSYfraction <- 1 - MSYreduction
    
    ## Prepare arguments to calculate reference points (fix parameters and latent variables, delta = 1)
    obj0 <- fit$obj
    argsIn <- as.list(obj0$env)[methods::formalArgs(TMB::MakeADFun)[methods::formalArgs(TMB::MakeADFun) != "..."]]
    argsIn$silent <- fit$obj$env$silent
    argsIn$parameters <- fit$pl
    argsIn$random <- unique(names(obj0$env$par[obj0$env$random]))
    ## Add referencepointSet
    catchType <- pmatch(catchType,c("catch","landing","discard"))
    if(is.na(catchType))
        stop("Invalid catch type")

    RecCorrection <- pmatch(RecCorrection,c("mean","median","mode"))
    if(is.na(RecCorrection))
        stop("Invalid recruitment correction")


    aveYearsIn <- aveYears
    aveYears <- match(aveYears, fit$data$years) - 1
    if(any(is.na(aveYears)))
        stop("aveYears has years without data.")

    selYearsIn <- selYears
    selYears <- match(selYears, fit$data$years) - 1
    if(any(is.na(selYears)))
        stop("selYears has years without data.")

    
    argsIn$data$referencepoint <- list(nYears = nYears,
                                       aveYears = aveYears,
                                       selYears = selYears,
                                       Fsequence = Fsequence,
                                       xPercent = SPRpercent,
                                       xB0 = B0percent,
                                       xdYPR = dYPRpercent,
                                       MSYRange = MSYfraction,
                                       catchType = catchType-1,
                                       RecCorrection = RecCorrection-1,
                                       optN = optN
                                       )

    args <- argsIn
    ## Remove random
    args$random <- NULL
    ## Remove referencepoint parameters from map
    map0 <- args$map
    fix <- names(args$parameters) ## setdiff(names(args$parameters), args$random)
    args$map <- lapply(args$parameters[fix], function(x)factor(x*NA))

    if(fit$conf$stockRecruitmentModelCode %in% c(0)){ # RW
        rp <- c("logScaleFmax",
                "logScaleFxdYPR",
                "logScaleFxPercent")
        MSYfraction <- MSYreduction <- c()
    }else if(fit$conf$stockRecruitmentModelCode %in% c(61,63)){ # Hockey-sticks
        rp <- c("logScaleFmsy",
                "logScaleFmypyl",
                "logScaleFmdy",
                "logScaleFmax",
                "logScaleFxdYPR",
                "logScaleFcrash",
                "logScaleFext",
                "logScaleFxPercent",
                "logScaleFlim",
                "logScaleFmsyRange")
    }else if(fit$conf$stockRecruitmentModelCode %in% c(3,62)){ # constant mean, AR
        rp <- c("logScaleFmsy",
                "logScaleFmypyl",
                "logScaleFmdy",
                "logScaleFmax",
                "logScaleFxdYPR",
                "logScaleFxPercent",
                "logScaleFmsyRange")
    }else if(fit$conf$stockRecruitmentModelCode %in% c(64)){ # Pow CMP
        rp <- c("logScaleFmsy",
                "logScaleFmypyl",
                "logScaleFmdy",
                "logScaleFmax",
                "logScaleFxdYPR",
                "logScaleFxPercent",
                "logScaleFmsyRange"
                )
    }else if(fit$conf$stockRecruitmentModelCode %in% c(65)){ # Pow Non-CMP
        rp <- c("logScaleFmsy",
                "logScaleFmypyl",
                "logScaleFmdy",
                "logScaleFmax",
                "logScaleFxdYPR",
                "logScaleFxPercent",
                "logScaleFmsyRange"
                )
    }else if(fit$conf$stockRecruitmentModelCode %in% c(68,69) && fit$pl$rec_par[3] > 0){ ## depensatory recruitment; Fcrash does not work.
        rp <- c("logScaleFmsy",
                "logScaleFmypyl",
                "logScaleFmdy",
                "logScaleFmax",
                "logScaleFxdYPR",
                "logScaleFext",
                "logScaleFxPercent",
                "logScaleFmsyRange"
                )
    }else if(fit$conf$stockRecruitmentModelCode %in% c(90)){
        rp <- c("logScaleFmsy",
                "logScaleFmypyl",
                "logScaleFmdy",
                "logScaleFmax",
                "logScaleFxdYPR",
                "logScaleFcrash",
                "logScaleFext",
                "logScaleFxPercent",
                "logScaleFmsyRange"
                )
       }else if(fit$conf$stockRecruitmentModelCode %in% c(91,92)){
        rp <- c("logScaleFmsy",
                "logScaleFmypyl",
                "logScaleFmdy",
                "logScaleFmax",
                "logScaleFxdYPR",
                ##"logScaleFcrash",
                "logScaleFext",
                "logScaleFxPercent",
                "logScaleFmsyRange"
                )
    }else{
        rp <- c("logScaleFmsy",
                "logScaleFmypyl",
                "logScaleFmdy",
                "logScaleFmax",
                "logScaleFxdYPR",
                "logScaleFcrash",
                "logScaleFext",
                "logScaleFxPercent",
                "logScaleFmsyRange"
                )
    }

    rp <- setdiff(rp, fixRP)
    
    ## Referencepoints to estimate
    args$map <- args$map[-which(names(args$map) %in% rp)]

    args$parameters$logScaleFmsy <- -1
    args$parameters$logScaleFmypyl <- -1
    args$parameters$logScaleFmdy <- -1
    args$parameters$logScaleFxdYPR <- rep(-1, length(dYPRpercent))
    args$parameters$logScaleFxB0 <- rep(-1, length(B0percent))
    args$parameters$logScaleFmax <- -1
    args$parameters$logScaleFcrash <- -1
    args$parameters$logScaleFext <- -1
    args$parameters$logScaleFxPercent <- rep(-1, length(SPRpercent))
    if(length(MSYfraction) > 0){
        args$parameters$logScaleFmsyRange <- rbind(log(-log(MSYfraction)),
                                                   log(log(1+(1-MSYfraction))))
    }
    args$parameters$logScaleFlim <- -1
    args$parameters$implicitFunctionDelta <- 1

    objOptim <- do.call(TMB::MakeADFun, args)
    pStart <- list()
    
    rep <- objOptim$report()
    tryAgain <- FALSE

    ## Fmax
    if(which.max(rep$refpointseq_logYPR[is.finite(rep$refpointseq_logYPR)]) == length(rep$refpointseq_logYPR[is.finite(rep$refpointseq_logYPR)]) && any(rp %in% "logScaleFmax")){
        warning("The stock does not appear to have a well-defined Fmax. Fmax will not be estimated. Increase the upper bound of Fsequence to try again.")
        rp <- rp[-which(rp %in% "logScaleFmax")]
        args$map$logScaleFmax <- factor(NA)
        tryAgain <- TRUE
    }else if(any(rp %in% "logScaleFmax")){
        indx <- which(is.finite(rep$refpointseq_logYPR) & Fsequence > 0)
        ypr <- rep$refpointseq_logYPR[indx]
        ff <- Fsequence[indx]
        pStart$logScaleFmax <- log(ff[which.max(ypr)]) - log(tail(fbartable(fit)[,"Estimate"],1))
    }

    ## FCrash
    if(any(rp %in% "logScaleFcrash") && min(rep$refpointseq_logSe[!is.na(rep$refpointseq_logSe)], na.rm = TRUE) > -4){
        warning("The stock does not appear to have a well-defined Fcrash. Fcrash will not be estimated. Increase the upper bound of Fsequence to try again.")
        rp <- rp[-which(rp %in% "logScaleFcrash")]
        args$map$logScaleFcrash <- factor(NA)
        tryAgain <- TRUE
    }else if(any(rp %in% "logScaleFcrash")){
        if(fit$conf$stockRecruitmentModelCode == 90){
            indx <- which(!is.na(rep$refpointseq_logSe))
            vi <- rep$refpointseq_logSe[indx]
            ff <- Fsequence[indx]
            pStart$logScaleFcrash <- log(ff[which.min(vi)]) - log(tail(fbartable(fit)[,"Estimate"],1))
        }else{
            vSR <- .Call("stockRecruitmentModelR",
                         exp(-10),
                         fit$pl$rec_pars,
                         fit$conf$stockRecruitmentModelCode)
            v <- (tail(vSR$Gradient,1) - 1 / rep$refpointseq_logSPR)^2
            indx <- which(!is.na(v))
            vi <- v[indx]
            ff <- Fsequence[indx]
            pStart$logScaleFcrash <- log(ff[which.min(vi)]) - log(tail(fbartable(fit)[,"Estimate"],1))
        }
    }

    ## Fext
    if(any(rp %in% "logScaleFext") && min(rep$refpointseq_logSe[!is.na(rep$refpointseq_logSe)], na.rm = TRUE) > -4){
        warning("The stock does not appear to have a well-defined Fext. Fext will not be estimated. Increase the upper bound of Fsequence to try again.")
        rp <- rp[-which(rp %in% "logScaleFext")]
        args$map$logScaleFext <- factor(NA)
        tryAgain <- TRUE
    }else if(any(rp %in% "logScaleFext")){
        indx <- which(!is.na(rep$refpointseq_logSe))
        vi <- rep$refpointseq_logSe[indx]^2 + 10 * Fsequence[indx]
        ff <- Fsequence[indx]
        pStart$logScaleFext <- log(ff[which.min(vi)]) - log(tail(fbartable(fit)[,"Estimate"],1))
    }

    
   ## Fmsy
    if(any(rp %in% "logScaleFmsy") && which.max(rep$refpointseq_logYe[is.finite(rep$refpointseq_logYe)]) == length(rep$refpointseq_logYe[is.finite(rep$refpointseq_logYe)])){
        warning("The stock does not appear to have a well-defined Fmsy. Fmsy will not be estimated. Increase the upper bound of Fsequence to try again.")
        rp <- rp[-which(rp %in% "logScaleFmsy")]
        rp <- rp[-which(rp %in% "logScaleFmsyRange")]
        args$map$logScaleFmsy <- factor(NA)
        MSYfraction <- MSYreduction <- numeric(0)
        args$parameters$logScaleFmsyRange <- matrix(0,2,0)
        tryAgain <- TRUE
    }else if(any(rp %in% "logScaleFmsy")){
        indx <- which(is.finite(rep$refpointseq_logYe) & Fsequence > 0)
        ye <- rep$refpointseq_logYe[indx]
        ff <- Fsequence[indx]
        pStart$logScaleFmsy <- log(ff[which.max(ye)]) - log(tail(fbartable(fit)[,"Estimate"],1))
        if(any(rp %in% "logScaleFmsyRange")){
            indx2 <- which(is.finite(rep$refpointseq_logYe) & Fsequence > 0 & Fsequence < ff[which.max(ye)])
            ye2 <- rep$refpointseq_logYe[indx2]
            ff2 <- Fsequence[indx2]
            fmsy <- ff[which.max(ye)]
            FmsyRangeLower <- sapply(MSYfraction, function(x){
                fL <- ff2[which.min((ye2 - x * max(ye))^2)]
                -log(-(log(fL) - log(fmsy)))
            })
            indx3 <- which(is.finite(rep$refpointseq_logYe) & Fsequence > 0 & Fsequence > ff[which.max(ye)])
            ye3 <- rep$refpointseq_logYe[indx3]
            ff3 <- Fsequence[indx3]
            fmsy <- ff[which.max(ye)]
            FmsyRangeUpper <- sapply(MSYfraction, function(x){
                fL <- ff3[which.min((ye3 - x * max(ye))^2)]
                  log(log(fL) - log(fmsy))
            })
            pStart$logScaleFmsyRange <- rbind(FmsyRangeLower,
                                              FmsyRangeUpper)
        }
    }

    ## Fmypyl
    if(any(rp %in% "logScaleFmypyl")){
        v <- rep$refpointseq_logYe - rep$refpointseq_logYearsLost
        indx <- which(is.finite(v) & Fsequence > 0)
        vi <- v[indx]
        ff <- Fsequence[indx]
        pStart$logScaleFmypyl <- log(ff[which.max(vi)]) - log(tail(fbartable(fit)[,"Estimate"],1))
    }

    ## Fmdy
    if(any(rp %in% "logScaleFmdy")){
        v <- rep$refpointseq_logDiscYe
        indx <- which(is.finite(v) & Fsequence > 0)
        vi <- v[indx]
        ff <- Fsequence[indx]
        pStart$logScaleFmdy <- log(ff[which.max(vi)]) - log(tail(fbartable(fit)[,"Estimate"],1))
    }

    
    ## FxdYPR
    if(any(rp %in% "logScaleFxdYPR")){
        indx <- which(is.finite(rep$refpointseq_logYPR) & Fsequence > 0)
        ypr <- rep$refpointseq_logYPR[indx]
        ff <- Fsequence[indx]
        pStart$logScaleFxdYPR <- sapply(dYPRpercent,function(x){
            log(ff[which.min((diff(ypr)/diff(ff) - x * diff(ypr)[1] / diff(ff)[1])^2)]) - log(tail(fbartable(fit)[,"Estimate"],1))
        })
    }

    ## FxB0
    if(any(rp %in% "logScaleFxB0")){
        indx <- which(is.finite(rep$refpointseq_logSe))
        Se <- rep$refpointseq_logSe[indx]
        ff <- Fsequence[indx]
        pStart$logScaleFxdYPR <- sapply(B0percent,function(x){
            log(ff[which.min((Se - x * Se[1])^2)]) - log(tail(fbartable(fit)[,"Estimate"],1))
        })
    }

    ## Fx%
    if(any(rp %in% "logScaleFxPercent")){
        indx <- which(is.finite(rep$refpointseq_logSPR) & Fsequence > 0)
        spr <- rep$refpointseq_logSPR[indx]
        ff <- Fsequence[indx]
        pStart$logScaleFxPercent <- sapply(SPRpercent,function(x){
            log(ff[which.min((spr - x * spr[1])^2)]) - log(tail(fbartable(fit)[,"Estimate"],1))
        })
    }

    if(tryAgain)                        
        objOptim <- do.call(TMB::MakeADFun, args)

    p0 <- objOptim$par
    for(ii in names(pStart))
        p0[names(p0) %in% ii] <- pStart[[match(ii,names(pStart))]]

    opt <- nlminb(p0, objOptim$fn, objOptim$gr, objOptim$he)
    h0 <- objOptim$he(opt$par)
    ee <- eigen(h0)
    if(any(ee$values < sqrt(.Machine$double.eps))){
        indx <- which(apply(ee$vectors[,ee$values < sqrt(.Machine$double.eps), drop = FALSE], 1, function(x) any(abs(x)>sqrt(.Machine$double.eps))))
        warning(sprintf("Singular convergence. Consider not estimating one or more of: %s",paste(names(opt$par)[indx], collapse = ", ")))
    }else{
        for(ii in seq_len(newtonSteps)){
            g <- as.numeric( objOptim$gr(opt$par) )
            h <- objOptim$he(opt$par)
            opt$par <- opt$par - solve(h, g)
            opt$objective <- objOptim$fn(opt$par)
        }
    }
    opt$gr <- as.numeric( objOptim$gr(opt$par) )
    opt$he <- objOptim$he(opt$par)
    
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
    JacAll <- jacobian(function(x) objDelta$gr(x)[gridx], objDelta$par, h = jacobianHScale * abs(1e-04 * objDelta$par) + 1e-04 * (abs(objDelta$par) < sqrt(.Machine$double.eps/7e-07)))

    ## Implicit function gradient
    dCdTheta <- -svd_solve(JacAll[,gridx,drop=FALSE]) %*% JacAll[,-gridx,drop=FALSE]
    rownames(dCdTheta) <- gsub("^logScale","",names(objDelta$par)[gridx])
    colnames(dCdTheta) <- names(fit$obj$env$last.par)
    dCdTheta2 <- dCdTheta
    aa <- as.character(seq(fit$conf$minAge,fit$conf$maxAge,1))
    yy <- as.character(fit$data$years)
    colnames(dCdTheta2)[colnames(dCdTheta2) == "logF"] <- paste0("logF_","",row(objDelta$env$parameters$logF),"Y",yy[col(objDelta$env$parameters$logF)])
    colnames(dCdTheta2)[colnames(dCdTheta2) == "logN"] <- paste0("logN_","A",aa[row(objDelta$env$parameters$logN)],"Y",yy[col(objDelta$env$parameters$logN)])
    ## Gradient
    gradient <- dCdTheta2[,apply(dCdTheta2,2,function(x)sum(abs(x))) > 0]

    ## Get joint precision
    jointPrecision <- TMB::sdreport(fit$obj,
                                    fit$opt$par,
                                    svd_solve(fit$sdrep$cov.fixed),
                                    getJointPrecision = TRUE)$jointPrecision

    ## Do delta method
    xtra <- diag(1,length(fit$obj$env$last.par.best))
    diag(xtra)[fit$obj$env$random] <- 0
    dG <- rbind(xtra[diag(xtra) != 0,,drop = FALSE],dCdTheta)
    covAll <- dG %*% svd_solve(jointPrecision) %*% t(dG)
    covAllOld <- covAll
    i <- 21    

    if(any(na2(abs(covAll[lower.tri(covAll)]) > 1-1e-8,FALSE))){
        warning("Some estimates are highly correlated. The covariance matrix was modified slighly to be invertible.")
        ii <- unique(as.vector(which(abs(covAll) > 1-1e-8 & row(covAll) < col(covAll), arr.ind = TRUE)))
        CX <- local({x <- diag(0, nrow(covAll)); diag(x)[ii] <- 1e-5; x})
        covAll <- covAll + CX
    }
            
    ## Object to do sdreport (delta = 0)
    args <- argsIn
    ## Remvoe rp from map
    args$map <- args$map[!(names(args$map) %in% rp)]
    ## Use optimized parameters
    args$parameters <- objOptim$env$parList(x = opt$par)
    args$parameters$implicitFunctionDelta <- 0

    objSDR <- do.call(TMB::MakeADFun, args)
    ## objSDR$fn(objSDR$par)
    
    sdr <- TMB::sdreport(objSDR, objSDR$par, svd_solve(covAll),
                         bias.correct= biasCorrect,
                         skip.delta.method = biasCorrect,
                         bias.correct.control = list(sd = TRUE,
                                                     split = objSDR$env$ADreportIndex()[grepl("referencepoint.log(F|B|R|Y|S)",names(objSDR$env$ADreportIndex()))]
                                                     ))
    sdr$covAll <- covAll
    ssdr <- summary(sdr)

    toCI <- function(what){
        indx <- 1:2
        if(biasCorrect)
            indx <- 3:4
        tmp <- ssdr[rownames(ssdr) == what,indx,drop=FALSE]
        CI <- exp(tmp %*% cbind(CIL=c(1,-2),CIH=c(1,2)))
        cbind(Estimate = exp(tmp[,1]), CI)
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
    colnames(Ftab) <- colnames(Btab) <- colnames(Rtab) <- colnames(Ytab) <- colnames(SPRtab) <- colnames(YPRtab) <- c("Estimate","Low","High")

    toRowNames <- Vectorize(function(x){
        switch(gsub("^referencepoint.logF","",x),
               "sq"="Status quo",
               "0"="Zero catch",
               "msy"="MSY",
               "msyRange"="xMR",
               "max"="Max",
               "01"="0.1",
               "crash"="Crash",
               "ext"="Ext",
               "xPercent"="xP",
               "lim"="lim",
               "mypyl"="Max Yield per Year Lost",
               "mdy"="Max Discounted Yield per Year Lost",
               x
        )               
    })
    rn <- toRowNames(rownames(Ftab))
    rn[which(rn == "xP")] <- sapply(SPRpercent,function(x)sprintf("%s%%",x * 100))
    rn[which(rn == "xdYPR")] <- sapply(dYPRpercent,function(x)sprintf("%s%%",x * 100))
    rn[which(rn == "xB0")] <- sapply(B0percent,function(x)sprintf("%s%%",x * 100))
    rn[which(rn == "xMR")] <- sapply(MSYreduction,function(x){
        c(sprintf("MSY %s%% reduction range (Lower)",x * 100),
          sprintf("MSY %s%% reduction range (Upper)",x * 100))
    })
    rownames(Ftab) <- rownames(Btab) <- rownames(Rtab) <- rownames(Ytab) <- rownames(SPRtab) <- rownames(YPRtab) <- rn

    ## Ftab["Ext",c("Low","High")] <- NA
    ## Btab["Ext",c("Low","High")] <- NA
    ## Rtab["Ext",c("Low","High")] <- NA
    ## Ytab["Ext",c("Low","High")] <- NA
    ## SPRtab["Ext",c("Low","High")] <- NA
    ## YPRtab["Ext",c("Low","High")] <- NA
    
    YPRseq <- toCI("refpointseq_logYPR")
    SPRseq <- toCI("refpointseq_logSPR")
    Yieldseq <- toCI("refpointseq_logYe")
    Bseq <- toCI("refpointseq_logSe")
    Rseq <- toCI("refpointseq_logRe")
    LLseq <- toCI("refpointseq_logYearsLost")
    LEseq <- toCI("refpointseq_logLifeExpectancy")
    
    rownames(YPRseq) <- rownames(SPRseq) <- rownames(Yieldseq) <- rownames(Bseq) <- rownames(Rseq) <- rownames(LLseq) <- rownames(LEseq) <- argsIn$data$referencepoint$Fsequence
    
    res <- list(tables = list(F = Ftab,
                              Yield = Ytab,
                              YieldPerRecruit = YPRtab,
                              SpawnersPerRecruit = SPRtab,
                              Biomass = Btab,
                              Recruitment = Rtab
                             ),
                graphs = list(F = argsIn$data$referencepoint$Fsequence,
                              Yield = Yieldseq,
                              YieldPerRecruit = YPRseq,
                              SpawnersPerRecruit = SPRseq,
                              Biomass = Bseq,
                              Recruitment = Rseq,
                              YearsLost = LLseq,
                              LifeExpectancy = LEseq),
                opt = opt,
                sdr = sdr,
                fbarlabel = substitute(bar(F)[X - Y], list(X = fit$conf$fbarRange[1], Y = fit$conf$fbarRange[2]))## ,
                ## diagonalCorrection = tv
                )
                              
    attr(res,"aveYears") <-  aveYearsIn
    attr(res,"selYears") <- selYearsIn
                              
    attr(res,"fit") <- fit
    class(res) <- "sam_referencepoints"
    return(res)
}


##' @importFrom grDevices col2rgb devAskNewPage
##' @importFrom graphics legend segments polygon
##' @importFrom methods is
##' @importFrom stats na.omit
##' @importFrom utils head
##' @export   
plot.sam_referencepoints <- function(x,
                                     show = c(1L:3L,5L),
                                     estimates = c("Status quo", "MSY"),
                                     ask = TRUE,
                                     legend.args = list(x = "top", ncol = length(estimates)),
                                     zoomToCurve = TRUE,
                                     ...){

    estimates <- match.arg(estimates,
                           choices = rownames(x$tables$F),
                           several.ok = TRUE)

    if(ask){
        oask <- grDevices::devAskNewPage(TRUE)
        on.exit(grDevices::devAskNewPage(oask))
    }


      toReportNames <- Vectorize(function(x,type = "F"){
        paste0("referencepoint.log",type,switch(x,
                                                "Status quo"="sq",
                                                "Zero catch"="0",
                                                "MSY" = "msy",
                                                "Max" = "max",
                                                "0.1"="01",
                                                "Crash"="crash",
                                                "Ext"="ext",
                                                "lim"="lim",
                                                "xPercent"
                                                )
               )
        })
 
    

    doPlot <- function(y){
        getACol <- function(col, a){
            do.call("rgb",c(as.list(grDevices::col2rgb(col,TRUE))[-4],
                            list(alpha = a * 255),
                            list(maxColorValue = 255)))
        }
        yp <- pmatch(y,names(x$graphs))
        yt <- pmatch(y,names(x$tables))
        if(is.na(yp) || is.na(yt))
            stop("Error message")
        tmpf <- x$tables$F["Crash",1]
        tmpval <- x$graphs[[yp]]
        if(!is.na(tmpf))
            tmpval <- x$graphs[[yp]][abs(x$graphs$F -tmpf) > 0.1,]
        if(zoomToCurve)
            tmpval <- tmpval[,1]
        plot(x$graphs$F, x$graphs[[yp]][,1], type = "n", ylim = range(tmpval, finite = TRUE, na.rm = TRUE),
             xlab = x$fbarlabel, ylab = gsub("([a-z])([A-Z])", "\\1-\\2",names(x$graphs)[yp]),
             main = paste("Equilibrium",tolower(gsub("([a-z])([A-Z])", "\\1 \\2",names(x$graphs)[yp]))))
        polygon(c(x$graphs$F, rev(x$graphs$F)),
                c(x$graphs[[yp]][,2], rev(x$graphs[[yp]][,3])),
                col = getACol(1,0.3),
                border = NA)
        lines(x$graphs$F, x$graphs[[yp]][,1], col = 1, lwd = 3)
        for(j in seq_along(estimates)){
            r <- estimates[j]
            r <- pmatch(r, rownames(x$tables$F))
            if(is.na(r))
                stop("Error message 2")
            if(!is.na(x$tables$F[r,"Estimate"])){
                acol <- getACol(j+1, 0.2)
                usr <- par("usr")
                xvals <- x$tables$F[r,]
                yvals <- x$tables[[yt]][r,]
                graphics::polygon(c(usr[1], xvals[2], xvals[2], xvals[3], xvals[3], usr[1]),
                        c(yvals[2], yvals[2], usr[3], usr[3], yvals[3], yvals[3]),
                        border = NA,
                        col = acol)
                graphics::segments(xvals[1], usr[3], xvals[1], yvals[1],
                         col = j+1, lwd = 3, lty = 1)
                graphics::segments(usr[1], yvals[1], xvals[1], yvals[1],
                         col = j+1, lwd = 3, lty = 1)
                
            }
        }

        if(!is.null(legend.args)){
            defArgs <- list(legend = estimates,
                            lwd = 3,
                            fill = sapply(seq_along(estimates)+1, getACol, a = 0.2),
                            border = NA,
                            merge = TRUE,
                            title = "Reference points",
                            col = seq_along(estimates)+1)
            do.call(graphics::legend,
                    c(legend.args,defArgs[which(!names(defArgs) %in% names(legend.args))])
                    )
        }
    }
    
    bshow <- rep(FALSE, 5)
    if(is.character(show)){
        show <- stats::na.omit(pmatch(show, c("YieldPerRecruit", "SpawnersPerRecruit", "Yield", "Biomass", "Recruitment")))
    }
        
    bshow[show[show %in% (1:length(bshow))]] <- TRUE    
    ## YPR
    if(bshow[1L]){
        doPlot("YieldPerRecruit")
    }
    ## SPR
    if(bshow[2L]){ 
        doPlot("SpawnersPerRecruit")
    }
    ## Yield
    if(bshow[3L]){                       
        doPlot("Yield")
    }
    ## Biomass
    if(bshow[4L]){                       
        doPlot("Biomass")
    }
    ## Recruitment
    if(bshow[5L]){                       
        doPlot("Recruitment")
    }
}

##' Function to estimate ICES values
##'
##' The function calculates MSYBtrigger and Bpa from estimated reference points through a simulation forecast. If a sam object is used, other reference points (e.g., Fmsy) are estimated first.
##' 
##' @param x object to calculate for
##' @param nosim Number of simulations
##' @param nyears Number of years to forecast before equilibrium
##' @param ntail Number of years to use for calculations
##' @param ... Other parameters passed to \link{referencepoints}
##' @return ICES values
##' @author Christoffer Moesgaard Albertsen
##' @export
ICESvalues <- function(x,
                        nosim,
                       nyears,
                       ntail,
                        ...){
    UseMethod("ICESvalues")
}


##' @rdname ICESvalues
##' @method ICESvalues sam
##' @export
ICESvalues.sam <- function(x,
                            nosim,
                            nyears,
                            ...){
    rp <- referencepoints(x,...)
    ICESvalues(rp, nosim, nyears)
}

stat_mode <- function(x, ...){
    d <- density(x, n = 100000, ...)
    d$x[which.max(d$y)]
}


getFp05 <- function(fit, rp, Btrigger, nosim, nyears, ntail){
    
    getP <- function(Fval){
        h1 <- hcr(fit, Ftarget = Fval, Btrigger = as.numeric(Btrigger), nosim = nosim, nYears = nyears)
        mean(unlist(lapply(tail(h1$forecast, ntail),function(x)x$ssb)) < rp$tables$Biomass["lim","Estimate"], na.rm = TRUE)
    }
    getObj <- function(F){
        getP(F) - 0.05
    }

    F1 <- rp$tables$F["MSY","Estimate"]
    p1 <- getP(F1)
    while(p1 > 0.05){
        F2 <- F1
        F1 <- 0.5 * F1
        p1 <- getP(F1)
    }
    F2 <- rp$tables$F["lim","Estimate"]
    p2 <- getP(F2)
    while(p2 < 0.05){
        F1 <- F2
        F2 <- 1.5 * F2
        p2 <- getP(F2)
    }    
    opt <- uniroot(getObj, c(F1,F2))
    h <- hcr(fit, Ftarget = opt$root, Btrigger = as.numeric(Btrigger), nosim = nosim, nYears = nyears)
    list(Fp05 = opt$root,
         HCR = h,
         opt <- opt)
}


##' @rdname ICESvalues
##' @method ICESvalues sam_referencepoints
##' @export
ICESvalues.sam_referencepoints <- function(x,
                                           nosim = 1000,
                                           nyears = 100,
                                           ntail = 10,
                                           quantile_buffer_Fmsy = 0.5,
                                           calculate_Fp05 = TRUE,
                                            ...){
    fit <- attr(x,"fit")
    Fmsy <- summary(x$sdr)["referencepoint.logFmsy",1:2]
    FmsyUse <- Fmsy[1]
    if(is.finite(Fmsy[2]) && quantile_buffer_Fmsy > 0 && quantile_buffer_Fmsy < 1)
        FmsyUse <- qnorm(quantile_buffer_Fmsy, Fmsy[1], Fmsy[2])
    cat("\nCalculating MSYBtrigger\n")
    forecast_MSYBtrigger = modelforecast(fit,
                                         fval = exp(rep(FmsyUse, nyears)),
                                         nosim = nosim,
                                         processNoiseF=FALSE,
                                         deterministicF=TRUE,
                                         resampleFirst=FALSE,                                     
                                         rec.years=c(),
                                         ave.years=attr(x,"aveYears"),
                                         overwriteSelYears = attr(x,"selYears"),
                                         biasCorrect=FALSE)
    Blim <- x$tables$Biomass["lim","Estimate"]
    SigmaB <- diff(log(tail(ssbtable(fit),1))[c(1,3)])/2
    Bpa <- Blim * exp(1.645 * SigmaB)
    Flim <- x$tables$F["lim","Estimate"]
    SigmaF <- diff(log(tail(fbartable(fit),1))[c(1,3)])/2
    Fpa2017 <- Flim * exp(1.645 * SigmaF)
    SSB_5pct <- unname(quantile(unlist(lapply(tail(forecast_MSYBtrigger,ntail),
                                              function(x) x$ssb)), 0.05, na.rm = TRUE))
    SSB_mean <- unname(mean(unlist(lapply(tail(forecast_MSYBtrigger,ntail),
                                                 function(x) x$ssb)), na.rm = TRUE))
    SSB_median <- unname(quantile(unlist(lapply(tail(forecast_MSYBtrigger,ntail),
                                                 function(x) x$ssb)), 0.5, na.rm = TRUE))
    SSB_mode <- unname(stat_mode(unlist(lapply(tail(forecast_MSYBtrigger,ntail),
                                               function(x) x$ssb)), na.rm = TRUE))    
    MSYBt <- pmax(Bpa, SSB_5pct, na.rm = TRUE)
    attr(MSYBt,"SigmaB") <- SigmaB
    attr(MSYBt,"SSB_5pct") <- SSB_5pct
    attr(MSYBt,"Bpa") <- Bpa
    class(MSYBt) <- "sam_icesval_MSYBtrigger"
    
    r <- list(MSYBtrigger = MSYBt,
              SSB_summary = c(Mean = SSB_mean,
                              Mode = SSB_mode,
                              Median = SSB_median),
              Bpa = Bpa,
              Fpa2017 = Fpa2017)

    class(r) <- "sam_icesval"    
    attr(r,"fit") <- fit
    attr(r,"referencepoints") <- x
    attr(r,"ssb_samples") <- unlist(lapply(tail(forecast_MSYBtrigger,ntail),
                                           function(x) x$ssb))
    attr(r,"catch_samples") <- unlist(lapply(tail(forecast_MSYBtrigger,ntail),
                                             function(x) x$catch))
    if(calculate_Fp05){
        cat("\nCalculating Fp.05\n")
        Fp050 <- getFp05(fit, x, as.numeric(MSYBt), nosim, nyears, ntail)
        r$Fp05 <- Fp050$F
        attr(r,"HCR_Fp05") <- Fp050$F$HCR
        r$Fpa2021 <- Fp050$F        
    }else{
        r$Fp05 <- NA
        attr(r,"HCR_Fp05") <- NULL
        r$Fpa2021 <- NA        
    }
    r
}

#' @export
print.sam_icesval <- function(x, ...){
    cat("\nEstimated ICES values\n\n")
    v1 <- c("MSYBtrigger:","Bpa:")
    v2 <- sapply(c(x$MSYBtrigger,x$Bpa),function(x)formatC(x,format="f",...))
    nc <- pmax(40, max(nchar(v1))+max(nchar(v2))+1)
    for(i in seq_along(v1))
        cat(v1[i], sprintf(sprintf("%%%ds",30-nchar(v1[i])),v2[i]),"\n")
}
