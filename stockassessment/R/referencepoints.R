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

#' @rdname stockassessment-deprecated
#' @section \code{referencepoints}:
#' For \code{referencepoints}, use \code{\link{deterministicReferencepoints}}.

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
    .Defunct("deterministicReferencepoints","stockassessment")
    ##UseMethod("referencepoints")
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
                                nlminb.control = list(eval.max = 1000, iter.max = 1000),
                                ...){
    .Defunct("deterministicReferencepoints","stockassessment")
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
