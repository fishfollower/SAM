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
        alph <- 1/(n+1) * .5
        ## Rule of
        ## Stone, C. (1986). [Generalized Additive Models]: Comment. Statistical Science, 1(3), 312-314
        prb <- tail(head(plogis(seq(logit(1 / (n+1)), logit(n / (n+1)), len = i)),-1),-1)
        kn <- c(quantile(log(ssb),probs = c(0+2/n,prb,1-2/n))) ##c(0+5 * alph,prb,1 - 5 * alph)))
        cnfF0$constRecBreaks <- kn
        dp0 <- defpar(fitRW$data,cnfF0)
        for(nn in names(dp0))
            if(class(dp0[[nn]]) == class(fitRW$pl[[nn]]) &&
               length(dp0[[nn]]) == length(fitRW$pl[[nn]]))
                dp0[[nn]][] <- fitRW$pl[[nn]][]          
        F0 <- try({sam.fit(fitRW$data, cnfF0, dp0, map = map, ...)})
        list(AIC = ifelse(is(F0, "try-error"),NA,AIC(F0)),
             kn = kn,
             fit = F0)
    })
    kk <- which.min(sapply(findknot,function(x)x$AIC))
    if(returnFit)
        return(findknot[[kk]]$fit)
    findknot[[kk]]$kn
}
