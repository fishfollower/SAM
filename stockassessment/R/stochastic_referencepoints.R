
.refpointSMethodEnum <- c("Mean"=0,
                          "Median"=1,
                          "Mode"=2,
                          "Quantile"=1)

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
    }else if(!is.na(xvalPatterns[[mIndx]][1]) & is.numeric(xvalPatterns[[mIndx]][1])){
        xVal <- xvalPatters[[mIndx]][1]
    }
    c(list(methodType = .refpointSMethodEnum[mIndx],
           xVal = xVal),
      list(...))        
}



.refpointSCurveFit <- function(F,C,MT){
    mf <- model.frame(MT$formula, data.frame(F = F))
    X <- model.matrix(MT$formula,mf)    
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
    opt <- nlminb(numeric(ncol(X)), fn, control = list(iter.max=10000,eval.max=10000))
    attr(opt,"terms") <- terms(mf)
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



.perRecruitSR <- function(logf, fit, nYears, aveYears, selYears, pl = fit$pl, ct=0, logCustomSel = numeric(0)){
    if(length(logCustomSel) > 0){
        sel <- exp(logCustomSel)
    }else{
        sel <- exp(.logFtoSel(pl$logF, selYears, fit$conf))
    }
    as.data.frame(.Call("perRecruitSR",
                        logFbar = logf,
                        tmbdat = fit$obj$env$data,
                        pl = pl,
                        sel = sel,
                        aveYears = aveYears,
                        nYears = nYears,
                        CT = ct,
                        logNinit = as.vector(tail(t(pl$logN),1))
                        ))
}


.refpointSFitCriteria <- function(rpArgs, pl, MT, fit, nosim, Frange, aveYears, selYears, nYears, catchType){
    rfv <- function(n,a,b){
        u <- runif(n)
        v1 <- exp(runif(n,log(ifelse(a==0,0.002,a)), log(b)))
        v2 <- runif(n,a,b)
        ifelse(u < 0.25, v1, v2)
    }
    Fvals <- sort(rfv(nosim,Frange[1],Frange[2]))
    PRvals <- do.call("rbind",lapply(log(Fvals), .perRecruitSR,
                                     fit=fit,
                                     nYears=nYears,
                                     aveYears = aveYears,
                                     selYears = selYears,
                                     pl = pl,
                                     ct = catchType))  
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
            trans <- function(x, report=FALSE, keepMSY = FALSE, ...){
                xMSY <- exp(head(x,1))
                x2 <- matrix(tail(x,-1),2)
                x3 <- rbind(xMSY - exp(x2[1,]),
                            xMSY + exp(x2[2,]))
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
                xx <- trans(x,keepMSY=TRUE)
                p <- predict(CurveFit,as.vector(xx))
                sum((tail(p,-1) - rep(rp$xVal,each=2) * p[1])^2) - p[1]
            }
            startVals <- function(){
                fmsy <- Fseq[which.max(pv)]
                c2 <- sapply(rp$xVal, function(xx) (pv - xx*max(pv))^2)
                f0 <- apply(c2,2,function(cc){
                    fl <- Fseq[Fseq < fmsy][which.min(cc[Fseq < fmsy])]
                    fu <- Fseq[Fseq > fmsy][which.min(cc[Fseq > fmsy])]
                    c(log(fmsy-fl),log(fu-fmsy))
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
            Crit <- exp(PRvals$logYPR)
            cutfun <- function(x) rep(TRUE,length(x))
            trans <- function(x, report=FALSE, ...){
                v <- exp(x)
                if(report)
                    names(v) <- paste0(rp$xVal,"dYPR")
                v
            }
            ## fn <- function(x) -predict(CurveFit,trans(x))
            ## startVals <- function() log(Fseq[which.max(pv)])            
        }else if(rp$rpType == 5){ ## xSPR

        }else if(rp$rpType == 6){ ## xB0

        }else if(rp$rpType == 7){ ## MYPYLdiv
            Arng <- conf$maxAge - conf$minAge + 1
            v <- PRvals$logYe - log(1.0 + exp(PRvals$logYearsLost - log(Arng)))
            return(exp(v))
        }else if(rp$rpType == 8){ ## MYPYLprod

        }else if(rp$rpType == 9){ ## MDY

        }else if(rp$rpType == 10){ ## Crash

        }else if(rp$rpType == 11){ ## Ext
            return((exp(PRvals$logSe) - 1)^2)
        }else if(rp$rpType == 12){ ## Lim
            stop("Reference point type not implemented yet")
        }else if(rp$rpType == 101){ ## MSYBtrigger
            stop("Reference point type not implemented yet")
        }else if(rp$rpType == 102){ ## Fpa
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
        trans(opt$par, report = TRUE)
    }
    getDerivedValues <- function(f){
        if(!is.function(MT$derivedSummarizer)){        #Fit
            doOneA <- function(what){
                Crit <- exp(PRvals[[what]])
                cutfun <- function(x) x > max(x) * 0.1
                Frng <- range(Fvals[cutfun(Crit)])
                inRng <- function(x,rng) x > rng[1] & x < rng[2]
                indx <- inRng(Fvals,Frng)
                CurveFit <- .refpointSCurveFit(Fvals[indx], Crit[indx], MT)
                predict(CurveFit, f)
            }
            return(sapply(c("logYPR","logSPR","logSe","logRe","logYe","logLifeExpectancy","logYearsLost"), doOneA))
        }else if(is.function(MT$derivedSummarizer)){  #Simulate
            return(sapply(lapply(as.list(do.call("rbind",lapply(rep(log(f),nosim), .perRecruitSR,
                                                         fit=fit,
                                                         nYears=nYears,
                                                         aveYears = aveYears,
                                                         selYears = selYears,
                                                         pl = pl,
                                                         ct = catchType)))[c("logYPR","logSPR","logSe","logRe","logYe","logLifeExpectancy","logYearsLost")],
                          exp),
                          MT$derivedSummarizer))
        }else{
            stop("Derived type not implemented")
        }
    }
    F <- sapply(rpArgs, getOneRP)
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
         PRvals = PRvals)
}

.refpointSGrid <- function(rp, pl, MT, fit){
    Fvals <- seq(rp$Frange[1], rp$Frange[2],0.02)
    getC <- function(f){
        PRvals <- do.call("rbind", replicate(rp$nosim,.perRecruitSR(log(f),
                                                                    fit=fit,
                                                                    nYears=rp$nYears,
                                                                    aveYears = rp$aveYears,
                                                                    selYears = rp$selYears,
                                                                    pl = pl,
                                                                    ct = rp$catchType), simplify=FALSE))
        exp(PRvals$logYe)
    }
    vv <- lapply(Fvals, getC)
  
  
}

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

#' @export
stochasticReferencepoints <- function(fit,
                                       referencepoints,
                                       ...){
    UseMethod("stochasticReferencepoints")
}



#' @export
stochasticReferencepoints.sam <- function(fit,
                                          referencepoints,
                                          method = "Q0.5",
                                          catchType = "catch",
                                          nYears = 300,
                                          Frange = c(0,2),
                                          nosim = 1000,
                                          aveYears = max(fit$data$years)+(-9:0),
                                          selYears = max(fit$data$years),
                                          newton.control = list(),
                                          seed = .timeToSeed(),
                                          formula = ~ibc(F,5),
                                          nosim_ci = 200,
                                          derivedSummarizer = median,
                                          ...){


    ## Add some kind of progressbar / messages (and argument to silence)

    oldSeed <- NULL
    if(exists(".Random.seed"))
        oldSeed <- .Random.seed
    on.exit(set.seed(oldSeed))
    set.seed(seed)

    MT <- .refpointSMethodParser(method, formula = formula)

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

    if(!all(Frange >= 0) && Frange[1] < Frange[2] && length(Frange) ==2)
        stop("Wrong Frange")
    if(!nosim > 0)
        stop("nosim must be a positive integer")

 
    pb <- txtProgressBar(min = 0, max = nosim_ci+2, style = 3)
    incpb <- function() setTxtProgressBar(pb, pb$getVal()+1)
    ## Get RPs for best fit
    rpArgs <- Reduce(.refpointMerger,
                     lapply(referencepoints, .refpointParser, cutoff = 0.1),
                     list())
    a <- capture.output(v0 <- .refpointSFitCriteria(rpArgs, fit$pl, MT, fit, nosim, Frange, aveYears, selYears, nYears, catchType))
   incpb()
    ## Sample to get CIs
    plRep <- .asympSampleParVec(nosim_ci,fit, boundary = TRUE, returnList = TRUE)
    incpb()
    vv <- lapply(plRep, function(par){
        oN <- fit$obj
        a <- capture.output(invisible(oN$fn(par)))
        pl <- oN$env$parList(par,oN$env$last.par)
        a <- capture.output(v <- try({.refpointSFitCriteria(rpArgs,pl, MT, fit, nosim, Frange, aveYears, selYears, nYears, catchType)}, silent = TRUE))
        incpb()
        v
    })
    close(pb)
       
    ## Get Ye/Re/Se/... (how should they be summarized?)
    ii <- sapply(vv, class) == "try-error"
    resTabs <- lapply(rownames(vv[!ii][[1]]$Estimates), function(nm){
        tab <- cbind(v0$Estimate[nm,],t(apply(do.call("rbind",lapply(vv[!ii], function(x) x$Estimates[nm,])),2,quantile, prob = c(0.025,0.975))))
        colnames(tab) <- c("Estimate","Low","High")
        tab
    })
    names(resTabs) <- rownames(vv[!ii][[1]]$Estimates)

    Fseq <- seq(0,2,len=200)
    

      res <- list(tables = list(F = resTabs[["F"]],
                              Yield = resTabs[["Ye"]],
                              YieldPerRecruit = resTabs[["YPR"]],
                              SpawnersPerRecruit = resTabs[["SPR"]],
                              Biomass = resTabs[["Se"]],
                              Recruitment = resTabs[["Re"]],
                              LifeExpectancy = resTabs[["LifeExpectancy"]],
                              LifeYearsLost = resTabs[["YearsLost"]]
                              ),
                graphs = list(F = v0$Fvals,
                              Yield = exp(v0$logYe),
                              YieldPerRecruit = exp(YPRseq),
                              SpawnersPerRecruit = exp(SPRseq),
                              Biomass = exp(Bseq),
                              Recruitment = exp(Rseq),
                              YearsLost = exp(LLseq),
                              LifeExpectancy = exp(LEseq)),
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

    
    ## Make output tables

    list(RP = v0,
         RPsim = vv)
}
                             
