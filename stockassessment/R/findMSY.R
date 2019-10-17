

##' Add stock-recruitment curve to srplot
##'
##' @param fit Object to show SR-curve for
##' @param CI Add confidence intervals?
##' @param col Color of fitted line
##' @param cicol Color of confidence intervals
##' @param ... not used
##' @seealso srplot
##' @author Christoffer Moesgaard Albertsen
##' @export
addRecruitmentCurve <- function(fit,
                                CI = TRUE,
                                col = rgb(0.6,0,0),
                                cicol = rgb(0.6,0,0,0.3),
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
                    ...){
       X <- summary(fit)
       R <- X[, 1]
       S <- X[, 4]
       cf <- fit$sdrep$cov.fixed
       covEst <- cf[rownames(cf) == "rec_pars", colnames(cf) == "rec_pars"]
       m <- fit$obj$env$map$rec_pars
       if(is.null(m)){
           covar <- covEst
       }else{
           covar <- covEst[m,m]
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
       polygon(c(ssb, rev(ssb)),
               c(tab["CIlow",],
                 rev(tab["CIhigh",])),
               col = cicol,
               border = NA)
       lines(ssb,tab["Estimate",], col = col, lwd = 3)
       invisible(srfit)
} 


##' Estimating Fmsy
##'
##' Work in progress - do not use
##' @param fit a SAM fit
##' @param nYears Number of years to forecast
##' @param nlminb.control list of control variables for nlminb
##' @param ... other arguments
##' @author Christoffer Moesgaard Albertsen
##' @export
forecastMSY <- function(fit,
                    nYears = 100,
                    nlminb.control = list(eval.max = 100, iter.max = 100),
                    rec.years = c(),
                    processNoiseF = FALSE,
                    ...){
    UseMethod("forecastMSY")
}


##' @rdname forecastMSY
##' @method forecastMSY sam
##' @export
forecastMSY.sam <- function(fit,
                    nYears = 100,
                    nlminb.control = list(eval.max = 100, iter.max = 100),
                    rec.years = c(),
                    processNoiseF = FALSE,
                    ...){

    argsIn <- forecast(fit,
                       findMSY = rep(1,nYears),
                       rec.years = rec.years,
                       processNoiseF = processNoiseF,
                     ...)

    ## Find MSY value
    args <- argsIn
    map0 <- args$map
    fix <- setdiff(names(args$parameters), args$random)
    args$map <- lapply(args$parameters[fix], function(x)factor(x*NA))
    args$map$logFScaleMSY <- NULL
   
    obj <- do.call(TMB::MakeADFun, args)
    
    obj$fn(obj$par)
    obj$gr(obj$par)

    opt <- nlminb(obj$par, obj$fn, obj$gr, control = nlminb.control)
    ## sdr1 <- TMB::sdreport(obj, opt$par)

    ## Get standard errors
    args <- argsIn
    args$map$logFScaleMSY <- NULL
    args$parameters$logFScaleMSY <- opt$par
    obj2 <- do.call(TMB::MakeADFun, args)

    Sigma <- solve(fit$opt$he)

    env <- environment()
    Jacobian <- function(x){
        r <- .Call("jacobian",
                   function(x)tail(obj2$gr(x),1),
                   x,
                   env,
                   30,
                   0.1 * 10^floor(log10(abs(x))) + 1e-4, #abs(1e-4 * x) + 1e-4 * (abs(x) < 1e-5),
                   1e-12)
        do.call("cbind",r[-1])
    }

    JacAll <- Jacobian(obj2$par)

    dCdTheta <- solve(JacAll[,length(obj2$par),drop=FALSE]) %*% JacAll[,-length(obj2$par),drop=FALSE]


    ## Reuse old fit 
    dG <- rbind(diag(1,length(fit$opt$par)),dCdTheta)
    covAll <- dG %*% solve(fit$opt$he) %*% t(dG)
    covAll[nrow(covAll), ncol(covAll)] <- covAll[nrow(covAll), ncol(covAll)] +  1e-8

    sdr2 <- TMB::sdreport(obj2, obj2$par, solve(covAll))
    
 
    return(list(opt = opt, sdr = sdr2))

}

##' Estimate reference points
##'
##' Work in progress - do not use
##' @param fit an object to calculate reference points for
##' @param ... not used
##' @return a sam_referencepoints fit
##' @author Christoffer Moesgaard Albertsen
##' @export
referencepoints <- function(fit,
                            nYears,
                            Fsequence,
                            aveYears,
                            selYears,
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
                            aveYears = as.numeric(c(length(fit$data$years)-2:1)),
                            selYears = as.numeric(c(length(fit$data$years)-1)),
                            catchType = "catch",
                            ...){

    obj0 <- fit$obj
    argsIn <- as.list(obj0$env)[methods::formalArgs(TMB::MakeADFun)[methods::formalArgs(TMB::MakeADFun) != "..."]]
    argsIn$parameters <- fit$pl
    argsIn$random <- unique(names(obj0$env$par[obj0$env$random]))
    ## Add referencepointSet
    catchType <- pmatch(catchType,c("catch","landing","discard"))
    if(is.na(catchType))
        stop("Invalid catch type")
    argsIn$data$referencepoint <- list(nYears = nYears,
                                       aveYears = aveYears,
                                       selYears = selYears,
                                       Fsequence = Fsequence,
                                       catchType = catchType-1
                                       )

    args <- argsIn
    ## Remove referencepoint parameters from map
    map0 <- args$map
    fix <- setdiff(names(args$parameters), args$random)
    args$map <- lapply(args$parameters[fix], function(x)factor(x*NA))

    if(fit$conf$stockRecruitmentModelCode %in% c(0,3)){ # RW, constant mean
        rp <- c("logScaleFmax",
                "logScaleF01",
                "logScaleF35")
    }else if(fit$conf$stockRecruitmentModelCode %in% c(61,63)){ # Hockey-sticks
        rp <- c("logScaleFmsy",
                "logScaleFmax",
                "logScaleF01",
                "logScaleFcrash",
                "logScaleF35",
                "logScaleFlim")
    }else if(fit$conf$stockRecruitmentModelCode %in% c(62)){ # AR
        rp <- c("logScaleFmsy",
                "logScaleFmax",
                "logScaleF01",
                "logScaleF35")
    }else if(fit$conf$stockRecruitmentModelCode %in% c(64)){ # Pow CMP
        rp <- c("logScaleFmsy",
                "logScaleFmax",
                "logScaleF01",
                "logScaleF35"
                )
    }else if(fit$conf$stockRecruitmentModelCode %in% c(65)){ # Pow Non-CMP
        rp <- c("logScaleFmsy",
                "logScaleFmax",
                "logScaleF01",
                "logScaleF35"
                )
    }else{
        rp <- c("logScaleFmsy",
                "logScaleFmax",
                "logScaleF01",
                "logScaleFcrash",
                "logScaleF35"
                )
    }
    ## Referencepoints to estimate
    args$map <- args$map[-which(names(args$map) %in% rp)]

    args$parameters$logScaleFmsy <- -2
    args$parameters$logScaleF01 <- -2
    args$parameters$logScaleFmax <- -2
    args$parameters$logScaleFcrash <- -2
    args$parameters$logScaleF35 <- -2
    args$parameters$logScaleFlim <- -2

    
    obj <- do.call(TMB::MakeADFun, args)

    ## Take inital look at YPR / SPR to determine if Fmax makes sense
    rep <- obj$report()
    tryAgain <- FALSE
    if(which.max(rep$logYPR) == length(rep$logYPR) && any(rp %in% "logScaleFmax")){
        warning("The stock does not appear to have a well-defined Fmax. Fmax will not be estimated. Increase the upper bound of Fsequence to try again.")
        rp <- rp[-which(rp %in% "logScaleFmax")]
        args$map$logScaleFmax <- factor(NA)
        tryAgain <- TRUE
    }

      if(min(rep$logSe) > -10 && any(rp %in% "logScaleFcrash")){
        warning("The stock does not appear to have a well-defined Fcrash. Fmax will not be estimated. Increase the upper bound of Fsequence to try again.")
        rp <- rp[-which(rp %in% "logScaleFcrash")]
        args$map$logScaleFcrash <- factor(NA)
        tryAgain <- TRUE
    }

    if(which.max(rep$logYe) == length(rep$logYe) && any(rp %in% "logScaleFmsy")){
        warning("The stock does not appear to have a well-defined Fmsy. Fmsy will not be estimated. Increase the upper bound of Fsequence to try again.")
        rp <- rp[-which(rp %in% "logScaleFmsy")]
        args$map$logScaleFmsy <- factor(NA)
        tryAgain <- TRUE
    }


    if(tryAgain)
        obj <- do.call(TMB::MakeADFun, args)

    
    opt <- nlminb(obj$par, obj$fn, obj$gr)

    ## Get standard errors
    args <- argsIn
    args$map <- args$map[-which(names(args$map) %in% rp)]

    opl <- obj$env$parList(par = opt$par)
    
    args$parameters$logScaleFmsy <- opl$logScaleFmsy
    args$parameters$logScaleF01 <- opl$logScaleF01
    args$parameters$logScaleFmax <- opl$logScaleFmax
    args$parameters$logScaleFcrash <- opl$logScaleFcrash
    args$parameters$logScaleF35 <- opl$logScaleF35
    if(fit$conf$stockRecruitmentModelCode %in% c(61,63))
        args$parameters$logScaleFlim <- opl$logScaleFlim

    obj2 <- do.call(TMB::MakeADFun, args)

    Sigma <- solve(fit$opt$he)

      
    gridx <- which(names(obj2$par) %in% rp)
    env <- environment()
    Jacobian <- function(x){
        r <- .Call("jacobian",
                   function(x) obj2$gr(x)[gridx],
                   x,
                   env,
                   30,
                   0.1 * 10^floor(log10(abs(x))) + 1e-4, #abs(1e-4 * x) + 1e-4 * (abs(x) < 1e-5),
                   1e-12)
        do.call("cbind",r[-1])
    }

    JacAll <- Jacobian(obj2$par)

    dCdTheta <- solve(JacAll[,gridx,drop=FALSE]) %*% JacAll[,-gridx,drop=FALSE]

    varLogRefIn <- dCdTheta %*% solve(fit$opt$he) %*% t(dCdTheta)

    ## Reuse old fit 
    dG <- rbind(diag(1,length(fit$opt$par)),dCdTheta)
    covAll <- dG %*% solve(fit$opt$he) %*% t(dG)
    diag(covAll)[gridx] <- diag(covAll)[gridx] + 1e-8

    sdr2 <- TMB::sdreport(obj2, obj2$par, solve(covAll))
    ssdr <- summary(sdr2)

    toCI <- function(what){
        exp(ssdr[rownames(ssdr) == what,,drop=FALSE] %*% cbind(Estimate=c(1,0),CIL=c(1,-2),CIH=c(1,2)))
    }

    Ftab <- t(sapply(rownames(ssdr)[grepl("referencepoint.logF",rownames(ssdr))],
                     toCI))
    Btab <- t(sapply(rownames(ssdr)[grepl("referencepoint.logB",rownames(ssdr))],
                     toCI))
    Ytab <-  t(sapply(rownames(ssdr)[grepl("referencepoint.logY",rownames(ssdr)) & !grepl("referencepoint.logYPR",rownames(ssdr))],
                      toCI))
    SPRtab <-  t(sapply(rownames(ssdr)[grepl("referencepoint.logSPR",rownames(ssdr))],
                      toCI))
    YPRtab <-  t(sapply(rownames(ssdr)[grepl("referencepoint.logYPR",rownames(ssdr))],
                      toCI))
    colnames(Ftab) <- colnames(Btab) <- colnames(Ytab) <- colnames(SPRtab) <- colnames(YPRtab) <- c("Estimate","Low","High")

    toRowNames <- Vectorize(function(x){
        switch(gsub("^referencepoint.logF","",x),
               "sq"="Status quo",
               "0"="Zero catch",
               "msy"="MSY",
               "max"="Max",
               "01"="0.1",
               "crash"="Crash",
               "35"="35%",
               "lim"="lim"
        )               
    })
    rownames(Ftab) <- rownames(Btab) <- rownames(Ytab) <- rownames(SPRtab) <- rownames(YPRtab) <- toRowNames(rownames(Ftab))

    YPRseq <- toCI("logYPR")
    SPRseq <- toCI("logSPR")
    Yieldseq <- toCI("logYe")
    Bseq <- toCI("logSe")
    Rseq <- toCI("logRe")

    res <- list(tables = list(F = Ftab,
                              B = Btab,
                              Yield = Ytab,
                              ## Recruitment = Rtab,
                              SpawnersPerRecruit = SPRtab,
                              YieldPerRecruit = YPRtab),
                graphs = list(F = argsIn$data$referencepoint$Fsequence,
                              Yield = Yieldseq,
                              YieldPerRecruit = YPRseq,
                              SpawnersPerRecruit = SPRseq,
                              Biomass = Bseq,
                              Recruitment = Rseq),
                opt = opt,
                sdr = sdr2
                )
                              


    class(res) <- "sam_referencepoints"
    return(res)
}

   
