## getSplineRecBreaks <- function(dat,
##                                conf,
##                                par = defpar(dat,conf),
##                                srmc = 93,
##                                map = list(),
##                                knotRange = 3:5,
##                                returnFit = TRUE,
##                                ...){
##     if(!match(srmc, c(90,91,92,93), FALSE)){
##         stop("Not a spline recruitment model code")
##     }
##     constRecBreaks <- numeric(0)
##     conf$stockRecruitmentModelCode <- 0
##     if(!is.null(attr(par,"what")) && attr(par,"what") == "Estimate"){
##         par$missing <- NULL
##         attr(par,"what") <- NULL
##     }
##     fitRW <- sam.fit(dat,conf,par, map = map, silent = TRUE)
##     logit <- Vectorize(function(x) log(x/(1-x)))
##     findknot <- lapply(knotRange, function(i){
##         cat("\n\n\n\t",i,"knots\n\n\n")
##         cnfF0 <- fitRW$conf
##         cnfF0$stockRecruitmentModelCode <- srmc
##         ssb <- ssbtable(fitRW)[,1]
##         n <- length(ssb)
##         alph <- 0
##         ## Rule of
##         ## Stone, C. (1986). [Generalized Additive Models]: Comment. Statistical Science, 1(3), 312-314
##         prb <- tail(head(plogis(seq(logit(1 / (n+1)), logit(n / (n+1)), len = i)),-1),-1)
##         F0 <- local({xx <- list(); class(xx) <- "try-error";xx})
##         while(alph < 0.25*n && is(F0, "try-error")){
##             kn <- c(quantile(log(ssb),probs = c(0+alph/n,prb,1-alph/n))) ##c(0+5 * alph,prb,1 - 5 * alph)))
##             cnfF0$constRecBreaks <- kn
##             dp0 <- defpar(fitRW$data,cnfF0)
##             for(nn in names(dp0))
##                 if(class(dp0[[nn]]) == class(fitRW$pl[[nn]]) &&
##                    length(dp0[[nn]]) == length(fitRW$pl[[nn]]))
##                     dp0[[nn]][] <- fitRW$pl[[nn]][]          
##             F0 <- try({sam.fit(fitRW$data, cnfF0, dp0, map = map, ...)})
##             alph <- alph+1
##         }
##         list(AIC = ifelse(is(F0, "try-error"),NA,AIC(F0)),
##              kn = kn,
##              fit = F0)
##     })
##     kk <- which.min(sapply(findknot,function(x)x$AIC))
##     if(returnFit)
##         return(findknot[[kk]]$fit)
##     findknot[[kk]]$kn
## }


getSplineRecBreaks <- function(x, df, Boundary.knots = range(x), method = c("bs","stone"), ...){
    UseMethod("getSplineRecBreaks")
}

getSplineRecBreaks.numeric <- function(x, df, Boundary.knots = range(x), method = c("bs","stone"), ...){
    ## bs way
    if(df < 3)
        stop("df must be at least 3")
    if(df < length(Boundary.knots))
        stop("df can not be lower than the number of Boundary.knots")
    outside <- x < min(Boundary.knots) | x > max(Boundary.knots)
    if(any(outside))
        warning("x has values outside the boundary knots")
    niknots <- df - length(Boundary.knots)
    qi <- seq(from = 0, to = 1, length.out = niknots + 2L)[-c(1L, niknots + 2L)]
    iknots <-  quantile(x[!outside], qi)
    sort(unname(c(iknots,Boundary.knots)))
}

getSplineRecBreaks.sam <- function(x, df, Boundary.knots, method = c("bs","stone"), alpha = 0, ...){
    logssb <- log(ssbtable(x)[,1])
    if(missing(Boundary.knots)){
        n <- length(logssb)        
        Boundary.knots <- quantile(logssb, c(alpha/n,1-alpha/n))
    }
    getSplineRecBreaks(logssb,df,Boundary.knots, method, ...)
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


stockRecruitmentStartingValues <- function(fit, stockRecruitmentModelCode, constRecBreaks, df = 3, ...){
    ##srAll <- c(1,2,3,60,61,63,64,66,67,68,69,90,91,92,93, 201,202,260,261,264,266,26)
    if(missing(stockRecruitmentModelCode)){
        stockRecruitmentModelCode <- c(1,2,3,60,61,63,64,66,67,68,69,90,91,92,93) #srAll
    }## else{
    ##     ii <- match(stockRecruitmentModelCode,srAll)
    ##     if(any(is.na(ii))){
    ##         if(all(is.na(ii)))
    ##             stop("Wrong model codes")
    ##         warning("Unused stockRecruitmentModelCodes")
    ##         ii <- na.omit(ii)
    ##     }
    ##     srUse <- srAll[ii]
    ## }
    X <- summary(fit)
    n <- nrow(X)
    lag <- fit$conf$minAge
    idxR <- (lag + 1):n
    idxS <- 1:(n - lag)
    R <- X[idxR, 1]
    S <- X[idxS, 4]
    year <- as.numeric(rownames(X))
    lastR <- c(NA_real_, head(R,-1))

    getBrks <- function(sr){
        if(sr == 3){
            return(head(seq(min(year),max(year), len = df + 1),-1))
        }else if( sr %in% c(90,91,92,93,290,293,490,493)){ #Add depensatory splines
            return(getSplineRecBreaks(fit, df=df, ...))
        }else{
            return(numeric(0))
        }
    }
    if(missing(constRecBreaks)){
        constRecBreaks <- lapply(stockRecruitmentModelCode, getBrks)
    }
    
    fn <- function(p, sr){
        brk <- constRecBreaks[[match(sr,stockRecruitmentModelCode)]]
        v <- logSRR(log(S),head(p,-1),sr, brk , year, lastR)
        -sum(dnorm(log(R), v, exp(tail(p,1)),log=TRUE))
    }
    p0 <- lapply(stockRecruitmentModelCode, function(sr){
        c0 <- fit$conf
        c0$stockRecruitmentModelCode <- sr
        c0$constRecBreaks <- constRecBreaks[[match(sr,stockRecruitmentModelCode)]]
        defpar(fit$data,c0)$rec_pars
    })
    doOne <- function(sr){
        pp <- p0[[match(sr,stockRecruitmentModelCode)]]
        r <- nlminb(c(pp,0),fn,sr=sr, control = list(iter.max = 10000, eval.max = 10000))
        r$AIC <- 2 * length(r$par) + 2 * r$objective
        r$constRecBreaks <- constRecBreaks[[match(sr,stockRecruitmentModelCode)]]
        r
    }
    res <- lapply(stockRecruitmentModelCode, doOne)
    names(res) <- stockRecruitmentModelCode
    res
}




##' Spline basis for use with formula interface
##'
##' @param x Points to evaluate the basis in
##' @param df Degrees of freedom
##' @param knots Internal knots. If NULL, they are selected from quantiles of x.
##' @param Boundary.knots Boundary knots. Defaults to range of x
##' @param intercept Include an intercept in basis?
##' @return A spline basis
##' @importFrom stats quantile
##' @importFrom utils head tail
##' @author Christoffer Moesgaard Albertsen
##' @export
bc <- function(x, df = 3L, knots = NULL, Boundary.knots = range(x), intercept = FALSE){
    if(is.null(knots)){
        nik <- df-length(Boundary.knots)-intercept
        if(nik < 0){
            warning(gettextf("'df' was too small; have used %d",length(Boundary.knots)-intercept))
            nik <- 0
        }
        if(nik > 0)
            knots <- utils::tail(utils::head(stats::quantile(x,seq(0,1,len=nik+2)),-1),-1)
    }
    Aknots <- sort(c(knots,Boundary.knots))
    v <- .Call("splinebasis_bcR",x,Aknots)
    v[is.na(v) | is.nan(v)] <- 0        # Stability issue quick fix (maybe from dnorm?)
    if(intercept)
        v <- cbind(1,v)
    a <- list(df = df, knots = if (is.null(knots)) numeric(0L) else knots, 
              Boundary.knots = Boundary.knots, intercept = intercept)

    attributes(v) <- c(attributes(v), a)
    class(v) <- c("bc","basis","matrix")
    v
}

##' @importFrom stats makepredictcall
##' @method makepredictcall bc
##' @export
makepredictcall.bc <- function (var, call){
    if (as.character(call)[1L] != "bc") 
        return(call)
    at <- attributes(var)[c("df", "knots", "Boundary.knots", 
        "intercept")]
    xxx <- call[1L:2]
    xxx[names(at)] <- at
    xxx
}



##' Integrated spline basis for use with formula interface
##'
##' @param x Points to evaluate the basis in
##' @param df Degrees of freedom
##' @param knots Internal knots. If NULL, they are selected from quantiles of x.
##' @param Boundary.knots Boundary knots. Defaults to range of x
##' @param intercept Include an intercept in basis?
##' @return A spline basis
##' @importFrom stats quantile
##' @importFrom utils head tail
##' @author Christoffer Moesgaard Albertsen
##' @export
ibc <- function(x, df = 3L, knots = NULL, Boundary.knots = range(x), intercept = FALSE){
    if(is.null(knots)){
        nik <- df-length(Boundary.knots)-intercept
        if(nik < 0){
            warning(gettextf("'df' was too small; have used %d",length(Boundary.knots)-intercept))
            nik <- 0
        }
        if(nik > 0)
            knots <- utils::tail(utils::head(stats::quantile(x,seq(0,1,len=nik+2)),-1),-1)
    }
    Aknots <- sort(c(knots,Boundary.knots))
    v <- .Call("splinebasis_ibcR",x,Aknots)
    if(intercept)
        v <- cbind(1,v)
    a <- list(df = df, knots = if (is.null(knots)) numeric(0L) else knots, 
              Boundary.knots = Boundary.knots, intercept = intercept)

    attributes(v) <- c(attributes(v), a)
    class(v) <- c("ibc","basis","matrix")
    v
}

##' @importFrom stats makepredictcall
##' @method makepredictcall ibc
##' @export
makepredictcall.ibc <- function (var, call){
    if (as.character(call)[1L] != "ibc") 
        return(call)
    at <- attributes(var)[c("df", "knots", "Boundary.knots", 
        "intercept")]
    xxx <- call[1L:2]
    xxx[names(at)] <- at
    xxx
}




##' Double integrated spline basis for use with formula interface
##'
##' @param x Points to evaluate the basis in
##' @param df Degrees of freedom
##' @param knots Internal knots. If NULL, they are selected from quantiles of x.
##' @param Boundary.knots Boundary knots. Defaults to range of x
##' @param intercept Include an intercept in basis?
##' @return A spline basis
##' @importFrom stats quantile
##' @importFrom utils head tail
##' @author Christoffer Moesgaard Albertsen
##' @export
iibc <- function(x, df = 3L, knots = NULL, Boundary.knots = range(x), intercept = FALSE){
    if(is.null(knots)){
        nik <- df-length(Boundary.knots)-intercept
        if(nik < 0){
            warning(gettextf("'df' was too small; have used %d",length(Boundary.knots)-intercept))
            nik <- 0
        }
        if(nik > 0)
            knots <- utils::tail(utils::head(stats::quantile(x,seq(0,1,len=nik+2)),-1),-1)
    }
    Aknots <- sort(c(knots,Boundary.knots))
    v <- .Call("splinebasis_iibcR",x,Aknots)
    if(intercept)
        v <- cbind(1,v)
    a <- list(df = df, knots = if (is.null(knots)) numeric(0L) else knots, 
              Boundary.knots = Boundary.knots, intercept = intercept)

    attributes(v) <- c(attributes(v), a)
    class(v) <- c("iibc","basis","matrix")
    v
}

##' @importFrom stats makepredictcall
##' @method makepredictcall iibc
##' @export
makepredictcall.iibc <- function (var, call){
   if (as.character(call)[1L] != "iibc") 
        return(call)
    at <- attributes(var)[c("df", "knots", "Boundary.knots", 
        "intercept")]
    xxx <- call[1L:2]
    xxx[names(at)] <- at
    xxx
}

