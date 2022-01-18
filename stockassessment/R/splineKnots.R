getSplineRecBreaks <- function(dat,
                               conf,
                               par = defpar(dat,conf),
                               srmc = 90,
                               map = list(),
                               knotRange = 3:5,
                               returnFit = TRUE,
                               ...){
    if(!match(srmc, c(90,91,92), FALSE)){
        stop("Not a spline recruitment model code")
    }
    constRecBreaks <- numeric(0)
    conf$stockRecruitmentModelCode <- 0
    if(!is.null(attr(par,"what")) && attr(par,"what") == "Estimate"){
        par$missing <- NULL
        attr(par,"what") <- NULL
    }
    fitRW <- sam.fit(dat,conf,par, map = map, silent = TRUE)
    logit <- Vectorize(function(x) log(x/(1-x)))
    findknot <- lapply(knotRange, function(i){
        cat("\n\n\n\t",i,"knots\n\n\n")
        cnfF0 <- fitRW$conf
        cnfF0$stockRecruitmentModelCode <- srmc
        ssb <- ssbtable(fitRW)[,1]
        n <- length(ssb)
        alph <- 0
        ## Rule of
        ## Stone, C. (1986). [Generalized Additive Models]: Comment. Statistical Science, 1(3), 312-314
        prb <- tail(head(plogis(seq(logit(1 / (n+1)), logit(n / (n+1)), len = i)),-1),-1)
        F0 <- local({xx <- list(); class(xx) <- "try-error";xx})
        while(alph < 0.25*n && is(F0, "try-error")){
            kn <- c(quantile(log(ssb),probs = c(0+alph/n,prb,1-alph/n))) ##c(0+5 * alph,prb,1 - 5 * alph)))
            cnfF0$constRecBreaks <- kn
            dp0 <- defpar(fitRW$data,cnfF0)
            for(nn in names(dp0))
                if(class(dp0[[nn]]) == class(fitRW$pl[[nn]]) &&
                   length(dp0[[nn]]) == length(fitRW$pl[[nn]]))
                    dp0[[nn]][] <- fitRW$pl[[nn]][]          
            F0 <- try({sam.fit(fitRW$data, cnfF0, dp0, map = map, ...)})
            alph <- alph+1
        }
        list(AIC = ifelse(is(F0, "try-error"),NA,AIC(F0)),
             kn = kn,
             fit = F0)
    })
    kk <- which.min(sapply(findknot,function(x)x$AIC))
    if(returnFit)
        return(findknot[[kk]]$fit)
    findknot[[kk]]$kn
}


fitPenalizedSpline <- function(dat,
                               conf,
                               par = defpar(dat,conf),
                               srmc = 90,
                               map = list(),
                               maxKnots = 20,
                               ...){
    constRecBreaks <- numeric(0)
    conf$stockRecruitmentModelCode <- 0
    if(!is.null(attr(par,"what")) && attr(par,"what") == "Estimate"){
        par$missing <- NULL
        attr(par,"what") <- NULL
    }
    fitRW <- sam.fit(dat,conf,par, map = map, silent = TRUE)
    cnfF0 <- fitRW$conf
    cnfF0$stockRecruitmentModelCode <- srmc
    ssb <- ssbtable(fitRW)[,1]
    F0 <- local({xx <- list(); class(xx) <- "try-error";xx})
    nKnots <- maxKnots
    while(nKnots >= 3 && is(F0, "try-error")){
        kn <- log(seq(min(ssb), max(ssb), len = nKnots))
        cnfF0$constRecBreaks <- kn
        dp0 <- defpar(fitRW$data,cnfF0)
        for(nn in names(dp0))
            if(class(dp0[[nn]]) == class(fitRW$pl[[nn]]) &&
               length(dp0[[nn]]) == length(fitRW$pl[[nn]]))
                dp0[[nn]][] <- fitRW$pl[[nn]][]          
        F0 <- try({sam.fit(fitRW$data, cnfF0, dp0, map = map, penalizeSpline = TRUE, ...)})
        nKnots <- nKnots - 1
    }
    return(F0)
}


functionalSR_startingvalue <- function(fit, stockRecruitmentModelCode){
    srAll <- c(1,2,50,51,52,60,61,63,64,66,67,68,69)
    if(missing(stockRecruitmentModelCode)){
        srUse <- srAll
    }else{
        ii <- match(stockRecruitmentModelCode,srAll)
        if(any(is.na(ii))){
            if(all(is.na(ii)))
                stop("Wrong model codes")
            warning("Unused stockRecruitmentModelCodes")
            ii <- na.omit(ii)
        }
        srUse <- srAll[ii]
    }
    X <- summary(fit)
    n <- nrow(X)
    lag <- fit$conf$minAge
    idxR <- (lag + 1):n
    idxS <- 1:(n - lag)
    R <- X[idxR, 1]
    S <- X[idxS, 4]
    fn <- function(p, sr){
        v <- .Call("logSRR",
                   log(S),
                   head(p,-1),
                   sr)
        -sum(dnorm(log(R), v, exp(tail(p,1)),log=TRUE))
    }
    p0 <- lapply(srUse, function(sr){
        c0 <- fit$conf
        c0$stockRecruitmentModelCode <- sr
        defpar(fit$data,c0)$rec_pars
    })
    doOne <- function(sr){
        pp <- p0[[match(sr,srUse)]]
        r <- nlminb(c(pp,0),fn,sr=sr, control = list(iter.max = 10000, eval.max = 10000))
        r$AIC <- 2 * length(r$par) + 2 * r$objective
        r
    }
    res <- lapply(srUse, doOne)
    names(res) <- srUse
    res
}
