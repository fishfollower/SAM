
##' Calculate jacobian of a function
##'
##' @param func function
##' @param x parameter values
##' @param h step size
##' @param ... passed to func
##' @return jacobian matrix
jacobian <- function(func, x,
                     h = abs(1e-04 * x) + 1e-04 * (abs(x) < sqrt(.Machine$double.eps/7e-07)),
                     ## 0.1 * 10^floor(log10(abs(x))) + 1e-4,
                     ...){
         r <- .Call(C_jacobian,
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
grad <- function(func, x,
                 h = abs(1e-04 * x) + 1e-04 * (abs(x) < sqrt(.Machine$double.eps/7e-07)),
                 ##0.1 * 10^floor(log10(abs(x))) + 1e-4,
                 ...){
         r <- .Call(C_jacobian,
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
    argsIn$data$reportingLevel <- 0
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



##' @rdname referencepoints
##' @title Estimate reference points
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
    .Defunct("deterministicReferencepoints","stockassessment")
    ##UseMethod("referencepoints")
}


##' @rdname referencepoints
##' @method referencepoints sam
##' @param dYPRpercent Defunct
##' @param B0percent Defunct 
##' @param fixRP Defunct
##' @param RecCorrection Defunct
##' @param biasCorrect Defunct
##' @param nlminb.control Defunct
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
                                nlminb.control = list(eval.max = 1000, iter.max = 1000),
                                ...){
    .Defunct("deterministicReferencepoints","stockassessment")
}

#' @rdname stockassessment-deprecated
#' @name stockassessment-deprecated
#' @title Deprecated and defunct functions
#' @section \code{referencepoints}:
#' For \code{referencepoints}, use \code{\link{deterministicReferencepoints}}.
NULL

