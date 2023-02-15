
.refpointSMethodEnum <- c("Mean"=0,
                          "Median"=1,
                          "Mode"=2,
                          "Quantile"=1)

.SAM_setPB <- function(pb, value, label=""){
    if (!inherits(pb, "txtProgressBar")) 
        stop(gettextf("'pb' is not from class %s", dQuote("txtProgressBar")), 
             domain = NA)
    oldval <- pb$getVal()
    pb$up(value,label)
    invisible(oldval)
}

##' @importFrom utils flush.console
.SAMpb <- function (min = 0, max = 1, initial = 0, char = "=", width = NA, 
    title, label="", format = "| %.3f%%",...) {
    file <- ""
    style <- 3
    .val <- initial
    .killed <- FALSE
    .nb <- 0L
    .pc <- -1L
    nw <- nchar(char, "w")
    if (nw == 0) 
        stop("'char' must have a non-zero width")
    if (max <= min) 
        stop("must have 'max' > 'min'")
    up <- function(value, label="") {
        if (!is.finite(value) || value < min || value > max) 
            return()
        width <- getOption("width") - 10L - nchar(label,"w")  
        if (nw > 1) 
            width <- trunc(width/nw)
        .val <<- value
        nb <- (width * (value - min)/(max - min))
        pc <- (100 * (value - min)/(max - min))
        if (nb == .nb && pc == .pc) 
            return()
        cat(paste0(sprintf("\r %s |",label), strrep(" ", nw * width + 6)), file = file)
        cat(paste(c(sprintf("\r %s |",label), rep.int(char, nb), rep.int(" ", 
            nw * (width - nb)), sprintf(format, pc)), collapse = ""), 
            file = file)
        utils::flush.console()
        .nb <<- nb
        .pc <<- pc
    }
    getVal <- function() .val
    kill <- function() if (!.killed) {
        cat("\n", file = file)
        flush.console()
        .killed <<- TRUE
    }
    up(initial,label)
    structure(list(getVal = getVal, up = up, kill = kill), class = "txtProgressBar")
}




.timeToSeed <- function(){
    tp <- Sys.time()
    tv_sec <- as.numeric(tp)
    tv_nsec <- (as.numeric(format(tp,"%OS6"))-as.numeric(format(tp,"%S")))*10^9
    a <- bitwXor(bitwShiftL(tv_nsec,16), tv_sec)
    b <- bitwShiftL(Sys.getpid(),16)
    bitwXor(a,b)
}

.refpointSMethodParser <- function(x, ...){
    if(length(x) > 1)
        stop("Only one at a time!")
    ## Allow shorthand...
    typePatterns <- list(Mean = "^mean$",
                         Median = "^median$",
                         Mode = "^mode$",
                         Quantile = "^q0\\.[[:digit:]]+$")
    xvalPatterns <- list(Mean = NA,
                         Median = 0.5,
                         Mode = NA,
                         Quantile = c("(q)(0\\.[[:digit:]]+)","\\2"))
    mIndx <- which(sapply(typePatterns, function(p) grepl(p, tolower(x))))
    if(length(mIndx) == 0)
        stop(sprintf("Error in method specification. %s not recognized.",x))
    if(!is.na(xvalPatterns[[mIndx]][1]) & !is.numeric(xvalPatterns[[mIndx]][1])){
        pp <- xvalPatterns[[mIndx]]
        xVal <- as.numeric(gsub(pp[1], pp[2], tolower(x)))
    }else{ # if(!is.na(xvalPatterns[[mIndx]][1]) & is.numeric(xvalPatterns[[mIndx]][1])){
        xVal <- xvalPatterns[[mIndx]][1]
    }
    c(list(methodType = .refpointSMethodEnum[mIndx],
           xVal = xVal),
      list(...))        
}


##' @importFrom stats model.frame model.matrix terms
.refpointSCurveFit <- function(F,C,MT){
    mf <- stats::model.frame(MT$formula, data.frame(F = F))
    X <- stats::model.matrix(MT$formula,mf)    
    if(MT$methodType == 0){       # Mean
        Loss <- function(Obs,Pred) (Obs-Pred)^2
    }else if(MT$methodType == 1){       # Quantiles
        rho <- function(x,q) x * (q - (x<0))
        Loss <- function(Obs,Pred) rho(Obs-Pred, MT$xVal)    
    }else if(MT$methodType == 2){       # Mode
        stop("method not implemented yet")
    }
    fn <- function(par){
        beta <- par[1:(ncol(X))]
        X2 <- X
        pM <- X2 %*% beta
        sum(Loss(C,pM[]))
    }
    opt <- stats::nlminb(numeric(ncol(X)), fn, control = list(iter.max=10000,eval.max=10000))
    attr(opt,"terms") <- stats::terms(mf)
     class(opt) <- "rpscurvefit"
    opt
}

predict.rpscurvefit <- function(x,newF,...){
    X <- model.matrix(attr(x,"terms"),data.frame(F=newF))
    par <- x$par
    beta <- par[1:ncol(X)]
    pM <- X %*% beta
    as.numeric(pM)  
}

.getDoSim <- function(logf1, fit, nYears, aveYears, selYears, pl, logCustomSel = NULL, constraint = "F=%f", deterministicF = TRUE,...){
   if(length(logCustomSel) > 0){
        sel <- exp(logCustomSel)
    }else if(length(selYears) == 0){   
        sel <- NULL
    }else{
        sel <- exp(.logFtoSel(pl$logF, match(selYears, fit$data$years) - 1, fit$conf))
    }
   suppressWarnings(invisible(doSim <- modelforecast(fit,
                                                     rep(sprintf(constraint,exp(logf1)), nYears),
                                                     nosim = 1,
                                                     progress = FALSE,
                                                     returnObj = 2,
                                                     ave.years = aveYears,
                                                     custom_pl = pl,
                                                     cstomSel = sel,
                                                     deterministicF = deterministicF,
                                                     fastFixedF = grepl("[[:space:]]*F[[:space:]]*=[[:space:]]*%f[[:space:]]*",constraint),
                                                     ...)))
   doSim
}

.perRecruitSR <- function(logf, fit, nYears, aveYears, selYears, pl = fit$pl, ct=0, logCustomSel = NULL, nTail = 1,incpb=NULL, doSim = NULL, label="", constraint = "F=%f", deterministicF=TRUE, ...){    
    ## pl$missing <- NULL
    ## attr(pl,"what") <- NULL
    ## f2 <- sam.fit(fit$data,fit$conf,pl, run = FALSE)
    ## f2$opt <- list(par = f2$obj$par, objective=NA)
    ## f2$sdrep <- fit$sdrep[c("estY","covY")]
    ## class(f2) <- "sam"
    if(is.null(doSim))
        doSim <- .getDoSim(logf[1], fit, nYears, aveYears, selYears, pl, logCustomSel, deterministicF)
    if(isTRUE(all.equal(pl,fit$pl))){
        re_pl <- NULL
    }else{
        re_pl <- pl
    }

    ## progress <- TRUE
    ## if(progress){
    ##     pb <- .SAMpb(min = 0, max = length(logf))
    ##     incpb <- function() .SAM_setPB(pb, pb$getVal()+1)
    ## }else{
    ##     incpb <- function(){}
    ## }

    
    do.call("rbind",lapply(logf, function(lf){
        v <- doSim(sprintf(constraint,exp(lf)))
        logRe <- tail(v$logN[1,],nTail)
        logSe <- tail(v$logssb,nTail)
        logSPR <- tail(v$logEmpiricalSPR,nTail)
        if(ct == 0){
            logYe <- tail(v$logCatch,nTail)
            logYPR <- tail(v$logEmpiricalYPR,nTail)
        }else if(ct == 1){
            logYe <- tail(v$logLand,nTail)
            logYPR <- tail(v$logEmpiricalYPR_L,nTail)
        }else{
            logYe <- tail(log(exp(v$logCatch)-(v$logLand)),nTail)
            logYPR <- tail(v$logEmpiricalYPR_D,nTail)
        }
        logLE <- NA_real_ # tail(v$logLifeExpectancy,nTail)
        logYL <- NA_real_ #tail(v$logYLTF,nTail)

        res <- list(logF = rep(lf,length.out = nTail),
                logYPR = logYPR,
                logSPR = logSPR,
                logSe = logSe,
                logRe = logRe,
                logYe = logYe,
                dSR0  = NA_real_,
                logLifeExpectancy = logLE,
                logYearsLost = logYL,
                logDiscYe = NA_real_,
                logDiscYPR = NA_real_)
        if(!is.null(incpb))
            incpb(label)
        as.data.frame(res)
    }))
}

## .perRecruitSR <- function(logf, fit, nYears, aveYears, selYears, pl = fit$pl, ct=0, logCustomSel = numeric(0)){
##     if(length(logCustomSel) > 0){
##         sel <- exp(logCustomSel)
##     }else{
##         sel <- exp(.logFtoSel(pl$logF, selYears, fit$conf))
##     }
##     as.data.frame(.Call(C_perRecruitSR,
##                         logFbar = logf,
##                         tmbdat = fit$obj$env$data,
##                         pl = pl,
##                         sel = sel,
##                         aveYears = aveYears,
##                         nYears = nYears,
##                         CT = ct,
##                         logNinit = as.vector(tail(t(pl$logN),1))
##                         ))
## }

.numDiff <- function(x,y){
    diff(y) / diff(x)    
}

#' @importFrom stats runif predict
.refpointSFitCriteria <- function(rpArgs, pl, MT, fit, nosim, Frange, aveYears, selYears, nYears, catchType, nTail = 1,doSim=NULL,incpb=NULL,label="",constraint="F=%f", deterministicF = TRUE, ...){
    rfv <- function(n,a,b){
        u <- stats::runif(n)
        v1 <- exp(stats::runif(n,log(ifelse(a==0,0.002,a)), log(b)))
        v2 <- stats::runif(n,a,b)
        ifelse(u < 0.25, v1, v2)
    }
    Fvals <- sort(rfv(nosim,Frange[1],Frange[2]))
    PRvals <- .perRecruitSR(log(Fvals),
                            fit=fit,
                            nYears=nYears,
                            aveYears = aveYears,
                            selYears = selYears,
                            pl = pl,
                            ct = catchType,
                            nTail = nTail,
                            incpb = incpb,
                            doSim = doSim,
                            label=sprintf("%s equilibrium simulations",label),
                            constraint = constraint,
                            deterministicF = deterministicF,
                            ...)  
###### Different for different RP's
    getOneRP <- function(rp){
        if(rp$rpType == 1){ ## MSY
            Crit <- exp(PRvals$logYe)
            cutfun <- function(x) x > max(x) * rp$cutoff
            trans <- function(x, report=FALSE, ...){
                v <- exp(x)
                if(report)
                    names(v) <- "MSY"
                v
            }
            fn <- function(x) -predict(CurveFit,trans(x))
            startVals <- function() log(Fseq[which.max(pv)])
        }else if(rp$rpType == 2){ ## MSYRange
            Crit <- exp(PRvals$logYe)
            cutfun <- function(x) x > max(x) * rp$cutoff
            trans <- function(x, report=FALSE, ...){
                dots <- list(...)
                if("keepMSY" %in% names(dots)){
                    keepMSY <- dots$keepMSY
                }else{
                    keepMSY <- FALSE
                }
                xMSY <- exp(utils::head(x,1))
                x2 <- matrix(utils::tail(x,-1),2)
                x3 <- rbind(xMSY + exp(-exp(x2[1,])),
                            xMSY + exp(exp(x2[2,])))
                if(keepMSY)
                    return(c(xMSY, x3))
                if(report){
                    nm <- as.vector(t(outer(paste0(formatC(rp$xVal),"MSYRange"),c("Lower","Upper"),paste)))
                    xOut <- as.vector(x3)
                    names(xOut) <- nm
                    return(xOut)
                }
                x3
            }
            fn <- function(x){
                xx <- trans(x, keepMSY=TRUE)
                p <- predict(CurveFit,as.vector(xx))
                sum((tail(p,-1) - rep(rp$xVal,each=2) * p[1])^2) - p[1]
            }
            startVals <- function(){
                fmsy <- Fseq[which.max(pv)]
                c2 <- sapply(rp$xVal, function(xx) (pv - xx*max(pv))^2)
                f0 <- apply(c2,2,function(cc){
                    fl <- Fseq[Fseq < fmsy][which.min(cc[Fseq < fmsy])]
                    fu <- Fseq[Fseq > fmsy][which.min(cc[Fseq > fmsy])]
                    ##c(log(fmsy-fl),log(fu-fmsy))
                    c(log(-log(log(fl)-log(fmsy))), log(log(log(fu)-log(fmsy))))
                })
                c(log(fmsy),f0)
            }
        }else if(rp$rpType == 3){ ## Max
            Crit <- exp(PRvals$logYPR)
            cutfun <- function(x) rep(TRUE,length(x))
            trans <- function(x, report=FALSE, ...){
                v <- exp(x)
                if(report)
                    names(v) <- "Max"
                v
            }
            fn <- function(x) -predict(CurveFit,trans(x))
            startVals <- function() log(Fseq[which.max(pv)])                        
        }else if(rp$rpType == 4){ ## xdYPR
            stop("Reference point type not implemented yet")
            ## Need derivative!
            Crit <- exp(PRvals$logYPR)
            cutfun <- function(x) rep(TRUE,length(x))
            trans <- function(x, report=FALSE, ...){
                v <- exp(x)
                if(report)
                    names(v) <- paste0(rp$xVal,"dYPR")
                v
            }
            fn <- function(x){
                p0 <- predict(CurveFit,1e-4)
                p <- predict(CurveFit,trans(x))
                sum((p - rp$xVal * p0)^2)
            }
            startVals <- function() sapply(rp$xVal, function(xv) log(Fseq[which.min((pv - xv * pv[1])^2)]))
        }else if(rp$rpType == 5){ ## xSPR
            Crit <- exp(PRvals$logSPR)
            cutfun <- function(x) rep(TRUE,length(x))
            trans <- function(x, report=FALSE, ...){
                v <- exp(x)
                if(report)
                    names(v) <- paste0(rp$xVal,"SPR")
                v
            }
            fn <- function(x){
                p0 <- predict(CurveFit,1e-4)
                p <- predict(CurveFit,trans(x))
                sum((p - rp$xVal * p0)^2)
            }
            startVals <- function() sapply(rp$xVal, function(xv) log(Fseq[which.min((pv - xv * pv[1])^2)]))
        }else if(rp$rpType == 6){ ## xB0
            Crit <- exp(PRvals$logSe)
            cutfun <- function(x) rep(TRUE,length(x))
            trans <- function(x, report=FALSE, ...){
                v <- exp(x)
                if(report)
                    names(v) <- paste0(rp$xVal,"B0")
                v
            }
            fn <- function(x){
                p0 <- predict(CurveFit,1e-4)
                p <- predict(CurveFit,trans(x))
                sum((p - rp$xVal * p0)^2)
            }
            startVals <- function() sapply(rp$xVal, function(xv) log(Fseq[which.min((pv - xv * pv[1])^2)]))
        }else if(rp$rpType == 7){ ## MYPYLdiv
            stop("Reference point type not implemented yet")
            ## Arng <- conf$maxAge - conf$minAge + 1
            ## v <- PRvals$logYe - log(1.0 + exp(PRvals$logYearsLost - log(Arng)))
            ## return(exp(v))
        }else if(rp$rpType == 8){ ## MYPYLprod
            stop("Reference point type not implemented yet")
        }else if(rp$rpType == 9){ ## MDY
            stop("Reference point type not implemented yet")
        }else if(rp$rpType == 10){ ## Crash
            stop("Reference point type not implemented yet")
        }else if(rp$rpType == 11){ ## Ext
            return((exp(PRvals$logSe) - 1)^2)
        }else if(rp$rpType == 12){ ## Lim
            stop("Reference point type not implemented yet")
        }else{
            stop("Reference point type not implemented")
        }
        Frng <- range(Fvals[cutfun(Crit)])
        inRng <- function(x,rng) x > rng[1] & x < rng[2]
        indx <- inRng(Fvals,Frng)
        CurveFit <- .refpointSCurveFit(Fvals[indx], Crit[indx], MT)
        Fseq <- seq(Frng[1],Frng[2],len=200)
        pv <- predict(CurveFit,Fseq)
        opt <- nlminb(startVals(), fn)
        res <- trans(opt$par, report = TRUE)
        attr(res,"curve_fit") <- cbind(F=Fseq,Criteria=pv)
        res
    }
    getDerivedValues <- function(f){
        if(!is.function(MT$derivedSummarizer)){        #Fit
            doOneA <- function(what){
                Crit <- exp(PRvals[[what]])
                if(all(is.na(Crit)))
                   return(NA)
                cutfun <- function(x) x > max(x) * 0.1
                Frng <- range(Fvals[cutfun(Crit)])
                inRng <- function(x,rng) x > rng[1] & x < rng[2]
                indx <- inRng(Fvals,Frng)
                CurveFit <- .refpointSCurveFit(Fvals[indx], Crit[indx], MT)
                predict(CurveFit, f)
            }
            return(sapply(c("logYPR","logSPR","logSe","logRe","logYe","logLifeExpectancy","logYearsLost"), doOneA))
        }else if(is.function(MT$derivedSummarizer)){  #Simulate
            return(sapply(lapply(as.list(.perRecruitSR(rep(log(f),nosim),
                                                         fit=fit,
                                                         nYears=nYears,
                                                         aveYears = aveYears,
                                                         selYears = selYears,
                                                         pl = pl,
                                                         ct = catchType,
                                                       nTail = nTail,
                                                       incpb=incpb,
                                                       doSim = doSim,
                                                       label = sprintf("%s derived values",label),
                                                       constraint = constraint,
                                                       deterministicF = deterministicF,
                                                       ...))[c("logYPR","logSPR","logSe","logRe","logYe","logLifeExpectancy","logYearsLost")],
                          exp),
                          function(x){ if(all(is.na(x))) return(NA); MT$derivedSummarizer}))
        }else{
            stop("Derived type not implemented")
        }
    }
    Forig <- lapply(rpArgs, getOneRP)
    F <- simplify2array(Forig)
    Curves <- lapply(Forig, attr, which = "curve_fit")
    names(Curves) <- names(F)
    D <- sapply(F, getDerivedValues)    
    res <- rbind(logF=unname(F),D)
    colnames(res) <- names(F)
    rownames(res) <- gsub("^log","",rownames(res))
    ## Fseq <- seq(min(Fvals),max(Fvals),len=200)
    ## GraphVals <- rbind(logF=Fseq,sapply(Fseq, getDerivedValues))
    ## colnames(GraphVals) <- Fseq
    ## rownames(GraphVals) <- gsub("^log","",rownames(GraphVals))
    list(Estimates = res,
         ## GraphVals = GraphVals,
         Fvals = Fvals,
         PRvals = PRvals,
         Curves = Curves)
}

## .refpointSGrid <- function(rp, pl, MT, fit){
##     Fvals <- seq(rp$Frange[1], rp$Frange[2],0.02)
##     getC <- function(f){
##         PRvals <- do.call("rbind", replicate(rp$nosim,.perRecruitSR(log(f),
##                                                                     fit=fit,
##                                                                     nYears=rp$nYears,
##                                                                     aveYears = rp$aveYears,
##                                                                     selYears = rp$selYears,
##                                                                     pl = pl,
##                                                                     ct = rp$catchType), simplify=FALSE))
##         exp(PRvals$logYe)
##     }
##     vv <- lapply(Fvals, getC)
  
  
## }

.asympSampleParVec <- function(N,fit, boundary = TRUE, returnList = FALSE){
    C <- t(chol(fit$sdrep$cov.fixed))
    mu <- fit$sdrep$par.fixed

    lower2 <- rep(-Inf, length(mu))
    upper2 <- rep(Inf, length(mu))
    for (nn in names(fit$low)) lower2[names(mu) == nn] = fit$low[[nn]]
    for (nn in names(fit$hig)) upper2[names(mu) == nn] = fit$hig[[nn]]

    doOne <- function(){
        v <- C %*% rnorm(length(mu)) + mu
        if(boundary){
            atLBound <- (fit$opt$par < (lower2 + sqrt(.Machine$double.eps)))
            atUBound <- (upper2 < (fit$opt$par + sqrt(.Machine$double.eps)))
            atBound <- atLBound | atUBound
            v[atBound] <- (atLBound * lower2 + atUBound * upper2)[atBound]
        }
        v[,1]      
    }
    replicate(N,doOne(), !returnList)
}
##' Estimate stochastic reference points
##'
##' The function estimates reference points based on stochastic model forecasts.
##'
##'
##' @details
##' The following reference points are implemented:
##' \describe{
##'    \item{F=x}{F fixed to x, e.g., \code{"F=0.3"} (NOT IMPLEMENTED YET)}
##'    \item{StatusQuo}{F in the last year of the assessment (NOT IMPLEMENTED YET)}
##'    \item{StatusQuo-y}{F in the y years before the last in the assessment, e.g., \code{"StatusQuo-1"} (NOT IMPLEMENTED YET)}
##'    \item{MSY}{F that maximizes yield}
##'    \item{0.xMSY}{Fs that gives 0.x*100\% of MSY, e.g., \code{"0.95MSY"}}
##'    \item{Max}{F that maximizes yield per recruit (NOT IMPLEMENTED YET)}
##'    \item{0.xdYPR}{F such that the derivative of yield per recruit is 0.x times the derivative at F=0, e.g., \code{"0.1dYPR"} (NOT IMPLEMENTED YET)}
##'    \item{0.xSPR}{F such that spawners per recruit is 0.x times spawners per recruit at F=0, e.g., \code{"0.35SPR"} (NOT IMPLEMENTED YET)}
##'    \item{0.xB0}{F such that biomass is 0.x times the biomass at F=0, e.g., \code{"0.2B0"} (NOT IMPLEMENTED YET)}
##' }
##' Reference points can be estimated using these methods:
##' \describe{
##'    \item{Mean}{Use least squares to estimate mean equilibrium values}
##'    \item{Q0.x}{Use quantile regression to estimate the 0.x quantile of equilibrium values}
##'    \item{Median}{Identical to Q0.5}
##'    \item{Mode}{(NOT IMPLEMENTED YET)}
##' }
##' To estimate median equilibrium yield, as required by ICES, the method "Q0.5" should be used.
##' 
##' @examples
##' \dontrun{
##'   stochasticReferencepoints(fit, c("MSY","0.95MSY"))
##' }
##' @param fit a sam fit
##' @param referencepoints a character vector of reference points to estimate (see Details)
##' @param ... Additional arguments passed on
##' @param method estimation method (See Details)
##' @param catchType catch type: catch, landing, discard
##' @param nYears Number of years to forecast
##' @param Frange Range of F values to consider
##' @param nosim Number of simulations for estimation
##' @param aveYears Years to average over for biological input
##' @param constraint Format of forecast constraint. "\%f" is replaced by F values.
##' @param deterministicF If FALSE, modelled logF process noise will be added to target logF in forecasts.
##' @param ... additional parameters that can be passed on
##' @param selYears Years to average over for selectivity
##' @param newton.control List of control parameters for optimization
##' @param seed Seed for simulations
##' @param formula Formula to estimate optimization criteria as a function of F
##' @param nosim_ci Number of simulations for bootstrap confidence intervals
##' @param derivedSummarizer Function to summarize derived per-recruit values
##' @param nTail Number of years from the simulation to include in calculations
##' @return reference point object
##' @export
stochasticReferencepoints <- function(fit,
                                       referencepoints,
                                       ...){
    UseMethod("stochasticReferencepoints")
}




##' @rdname stochasticReferencepoints
##' @method stochasticReferencepoints sam
##' @export
stochasticReferencepoints.sam <- function(fit,
                                          referencepoints,
                                          method = "Q0.5",
                                          catchType = "catch",
                                          nYears = 300,
                                          Frange = c(0,2),
                                          nosim = 1000,
                                          aveYears = c(), #max(fit$data$years)+(-9:0),
                                          selYears = max(fit$data$years),
                                          newton.control = list(),
                                          seed = .timeToSeed(),
                                          formula = ~ibc(F,5),
                                          nosim_ci = 200,
                                          derivedSummarizer = NA,
                                          nTail = 1,
                                          constraint = "F=%f",
                                          deterministicF = TRUE,
                                          ...){


    ## Add some kind of progressbar / messages (and argument to silence)

    oldSeed <- NULL
    if(exists(".Random.seed"))
        oldSeed <- .Random.seed
    on.exit(set.seed(oldSeed))
    set.seed(seed)

    MT <- .refpointSMethodParser(method, formula = formula, derivedSummarizer=derivedSummarizer)

    catchType <- pmatch(catchType,c("catch","landing","discard"))-1
    if(is.na(catchType))
        stop("Invalid catch type")

    aveYearsIn <- aveYears
    ## aveYears <- match(aveYears, fit$data$years) - 1
    ## if(any(is.na(aveYears)))
    ##     stop("aveYears has years without data.")

    selYearsIn <- selYears
    ## selYears <- match(selYears, fit$data$years) - 1
    ## if(any(is.na(selYears)))
    ##     stop("selYears has years without data.")

    if(!all(Frange >= 0) && Frange[1] < Frange[2] && length(Frange) ==2)
        stop("Wrong Frange")
    if(!nosim > 0)
        stop("nosim must be a positive integer")

    ## Get RPs for best fit
    rpArgs <- Reduce(.refpointMerger,
                     lapply(referencepoints, .refpointParser, cutoff = 0.1),
                     list())
    invisible(lapply(rpArgs,.refpointCheckRecruitment,fit=fit))

    doSim <- .getDoSim(logf1=tail(log(fbartable(fit)[,1]),1),
                       fit=fit, nYears = nYears, aveYears = aveYears, selYears = selYears, pl = fit$pl, constraint=constraint,deterministicF=deterministicF,...)
 
    pb <- .SAMpb(min = 0, max = nosim * (nosim_ci + 1 + is.function(derivedSummarizer)*length(rpArgs)))
    incpb <- function(label="") .SAM_setPB(pb, pb$getVal()+1,label)

    v0 <- .refpointSFitCriteria(rpArgs, pl=fit$pl, MT=MT, fit=fit, nosim=nosim, Frange=Frange, aveYears=aveYears, selYears=selYears, nYears=nYears, catchType=catchType, nTail=nTail,incpb=incpb,doSim=doSim,label="Estimation:",constraint=constraint,deterministicF=deterministicF, ...)

    ## Sample to get CIs
    if(nosim_ci > 0){
        plRep <- .asympSampleParVec(nosim_ci,fit, boundary = TRUE, returnList = TRUE)
         vv <- lapply(plRep, function(par){
            oN <- fit$obj
            a <- capture.output(invisible(oN$fn(par)))
            pl <- oN$env$parList(par,oN$env$last.par)
            v <- try({.refpointSFitCriteria(rpArgs,pl=pl, MT=MT, fit=fit, nosim=nosim, Frange=Frange, aveYears=aveYears, selYears=selYears, nYears=nYears, catchType=catchType, nTail=nTail,incpb=incpb,doSim=doSim,label="Confidence intervals:",constraint=constraint,deterministicF=deterministicF, ...)}, silent = TRUE)
 
            v
        })
    ## Get Ye/Re/Se/... (how should they be summarized?)
        ii <- sapply(vv, class) == "try-error"
    }else{
        ii <- logical(0)
    }
    if(sum(!ii) > 0){
        nan2na <- function(x)ifelse(is.nan(x),NA,x)
        resTabs <- lapply(rownames(vv[!ii][[1]]$Estimates), function(nm){
            tab <- cbind(v0$Estimate[nm,],t(apply(do.call("rbind",lapply(vv[!ii], function(x) nan2na(x$Estimates[nm,]))),2,quantile, prob = c(0.025,0.975), na.rm=TRUE)))
            colnames(tab) <- c("Estimate","Low","High")
            rownames(tab) <- colnames(v0$Estimate)
            tab
        })
        names(resTabs) <- rownames(vv[!ii][[1]]$Estimates)
    }else{
        warning("Confidence intervals could not be calculated. Try to increase nosim_ci.")
        resTabs <- lapply(rownames(v0$Estimates), function(nm){
            tab <- cbind(v0$Estimate[nm,], NA, NA)
            colnames(tab) <- c("Estimate","Low","High")
            rownames(tab) <- colnames(v0$Estimate)
            tab
        })
        names(resTabs) <- rownames(v0$Estimate)
    }
    Fseq <- seq(0,2,len=200)
    
    ## Make output tables

      res <- list(tables = list(F = resTabs[["F"]],
                              Yield = resTabs[["Ye"]],
                              YieldPerRecruit = resTabs[["YPR"]],
                              SpawnersPerRecruit = resTabs[["SPR"]],
                              Biomass = resTabs[["Se"]],
                              Recruitment = resTabs[["Re"]],
                              LifeExpectancy = resTabs[["LifeExpectancy"]],
                              LifeYearsLost = resTabs[["YearsLost"]]
                              ),
                graphs = list(F = exp(v0$PRvals$logF),
                              Yield = exp(v0$PRvals$logYe),
                              YieldPerRecruit = exp(v0$PRvals$logYPR),
                              SpawnersPerRecruit = exp(v0$PRvals$logSPR),
                              Biomass = exp(v0$PRvals$logSe),
                              Recruitment = exp(v0$PRvals$logRe),
                              YearsLost = exp(v0$PRvals$logYearsLost),
                              LifeExpectancy = exp(v0$PRvals$logLifeExpectancy)),
                regression = v0$Curves,
                ## opt = NA,
                ## ssdr = sdr,
                fbarlabel = substitute(bar(F)[X - Y], list(X = fit$conf$fbarRange[1], Y = fit$conf$fbarRange[2])),
                stochastic = TRUE
                ## diagonalCorrection = tv
                )
    
    attr(res,"aveYears") <-  aveYearsIn
    attr(res,"selYears") <- selYearsIn
    
    attr(res,"fit") <- fit
    class(res) <- "sam_referencepoints"

    

    res
}
                             
