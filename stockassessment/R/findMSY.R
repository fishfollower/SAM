


stockRecruitmentModel <- function(code, rec_pars, covar = NULL){
    final <- function(ssb, f, gr, covar){
        v <- f(ssb)
        if(!is.null(covar)){
            g <- matrix(gr(ssb), 1)
            attr(v,"sd") <- as.vector(sqrt(g %*% covar %*% t(g)))
        }
        return(v)
    }
    if(code == 1){ ## Ricker
        fn <- function(ssb){
            f <- function(thisSSB) exp(rec_pars[1]+log(thisSSB)-exp(rec_pars[2])*thisSSB)
            gr <- function(thisSSB){
                c(exp(rec_pars[1] + log(thisSSB) - exp(rec_pars[2]) * thisSSB),
                  -(exp(rec_pars[1] + log(thisSSB) - exp(rec_pars[2]) * thisSSB) * (exp(rec_pars[2]) * thisSSB))
                  )}
            return(final(ssb, f, gr, covar))
        }
        attr(fn, "SRModel") <- "Ricker"
    }else if(code == 2){ ## Beverton-Holt
        fn <- function(ssb){
            f <- function(thisSSB) exp(rec_pars[1]+log(thisSSB)-log(1.0+exp(rec_pars[2])*thisSSB))
            gr <- function(thisSSB){
                c(exp(rec_pars[1] + log(thisSSB) - log(1 + exp(rec_pars[2]) * thisSSB)),
                  -(exp(rec_pars[1] + log(thisSSB) - log(1 + exp(rec_pars[2]) * thisSSB)) * (exp(rec_pars[2]) * thisSSB/(1 + exp(rec_pars[2]) * thisSSB))))
            }
            return(final(ssb, f, gr, covar))
        }
        attr(fn, "SRModel") <- "BevertonHolt"
    }else if(code == 61){ ## Hockey stick
        fn <- function(ssb){
            f <- function(thisSSB) exp(rec_pars[1] - rec_pars[2] +
                                       log(thisSSB - (0.5 * ((thisSSB - exp(rec_pars[2]))+(0.0)+abs((thisSSB - exp(rec_pars[2]))-(0.0))))))
            gr <- function(thisSSB){
                c(exp(rec_pars[1] - rec_pars[2]) * (thisSSB - 0.5 * (abs(exp(rec_pars[2]) - thisSSB) - exp(rec_pars[2]) + thisSSB)),
                (exp(rec_pars[1] - rec_pars[2]) * (-0.5 * thisSSB * abs(exp(rec_pars[2]) - thisSSB) + 0.5 * (exp(rec_pars[2]) - thisSSB)^2 + 0.5 * exp(rec_pars[2]) * thisSSB - 0.5 * exp(2 * rec_pars[2])))/abs(exp(rec_pars[2]) - thisSSB))
            }
            return(final(ssb, f, gr, covar))
        }
        attr(fn, "SRModel") <- "Hockey"
    }else if(code == 63){ ## Smooth hockey stick
        fn <- function(ssb){
            f <- function(thisSSB) exp(rec_pars[2] +
                                       log(thisSSB + sqrt(exp(2.0 * rec_pars[1]) + (exp(2.0 * rec_pars[3]) / 4.0)) -
                                           sqrt((thisSSB-exp(rec_pars[1]))^2 + (exp(2.0 * rec_pars[3]) / 4.0))))
            gr <- function(thisSSB){ 
                c(exp(rec_pars[2]) * ((exp(rec_pars[1]) * (thisSSB - exp(rec_pars[1])))/sqrt((thisSSB - exp(rec_pars[1]))^2 + 0.25 * exp(2 * rec_pars[3])) + exp(2 * rec_pars[1])/sqrt(exp(2 * rec_pars[1]) + 0.25 * exp(2 * rec_pars[3]))),
                  exp(rec_pars[2]) * (-sqrt((thisSSB - exp(rec_pars[1]))^2 + 0.25 * exp(2 * rec_pars[3])) + sqrt(exp(2 * rec_pars[1]) + 0.25 * exp(2 * rec_pars[3])) + thisSSB),
                  exp(rec_pars[2]) * ((0.25 * exp(2 * rec_pars[3]))/sqrt(exp(2 * rec_pars[1]) + 0.25 * exp(2 * rec_pars[3])) - (0.25 * exp(2 * rec_pars[3]))/sqrt((thisSSB - exp(rec_pars[1]))^2 + 0.25 * exp(2 * rec_pars[3]))))
            }
            return(final(ssb, f, gr, covar))
        }
        attr(fn, "SRModel") <- "SmoothHockey"
    }else if(code == 64){
        fn <- function(ssb){
            f <- function(thisSSB) exp(rec_pars[1] + exp(rec_pars[2]) * log(thisSSB))
            gr <- function(thisSSB){
                c( exp(rec_pars[1]) * thisSSB ^ ( exp(rec_pars[2]) ),
                  exp( rec_pars[1] + rec_pars[2] ) * thisSSB ^ ( exp(rec_pars[2]) ) * log(thisSSB) 
                  )               
            }
            return(final(ssb, f, gr, covar))
        }
        attr(fn, "SRModel") <- "Cushing"
    }else if(code == 65){
        fn <- function(ssb){
            f <- function(thisSSB) exp(rec_pars[1] + log(thisSSB) - log(1.0 + exp(exp(rec_pars[3]) * (log(thisSSB) - rec_pars[2]))))
            gr <- function(thisSSB){
                c( exp(rec_pars[1]) * thisSSB / (1.0 + exp( exp(rec_pars[3]) * ( - rec_pars[2] + log(thisSSB)))),
                (thisSSB * exp(rec_pars[1]+ exp(rec_pars[3]) * (log(thisSSB) - rec_pars[2]) + rec_pars[3]))/(exp(exp(rec_pars[3]) * (log(thisSSB) - rec_pars[2])) + 1)^2,
                -(thisSSB * (log(thisSSB) - rec_pars[2]) * exp(rec_pars[1] + exp(rec_pars[2]) * (log(thisSSB) - rec_pars[2]) + rec_pars[2]))/(exp(exp(rec_pars[3]) * (log(thisSSB) - rec_pars[2])) + 1)^2
                  )
            }
            return(final(ssb, f, gr, covar))
        }
        attr(fn, "SRModel") <- "Shepherd"
    }else{
        stop("Stock recruitment model is not a deterministic function of SSB")
    }
    return(fn)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param fit 
##' @param CI 
##' @param ... 
##' @return 
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

       srfit <- stockRecruitmentModel(fit$conf$stockRecruitmentModelCode,
                                      fit$pl$rec_pars,
                                      covar
                                      )
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
##' @param ... 
##' @return 
##' @author Christoffer Moesgaard Albertsen
##' @export

forecastMSY <- function(fit,
                    nYears = 100,
                    nlminb.control = list(eval.max = 100, iter.max = 100),
                    rec.years = c(),
                    processNoiseF = FALSE,
                    ...){
    UseMethod("findForecastMSY")
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

    ## can we use TMB/AD??
    JacAll <- numDeriv::jacobian(function(x)tail(obj2$gr(x),1),obj2$par)

    dCdTheta <- solve(JacAll[,length(obj2$par),drop=FALSE]) %*% JacAll[,-length(obj2$par),drop=FALSE]

    varLogFScale <- dCdTheta %*% solve(fit$opt$he) %*% t(dCdTheta)

    ## Reuse old fit 
    dG <- rbind(diag(1,length(fit$opt$par)),dCdTheta)
    covAll <- dG %*% solve(fit$opt$he) %*% t(dG)
    covAll[nrow(covAll), ncol(covAll)] <- covAll[nrow(covAll), ncol(covAll)] +  1e-8

    sdr2 <- TMB::sdreport(obj2, obj2$par, solve(covAll))
    
 
    return(list(opt = opt, sdr = sdr2))

}


#' @export
referencepoints <- function(fit,
                            ...){
    UseMethod("referencepoints")
}

#' @export
referencepoints.sam <- function(fit,
                            ...){

    obj0 <- fit$obj
    argsIn <- as.list(obj0$env)[methods::formalArgs(TMB::MakeADFun)[methods::formalArgs(TMB::MakeADFun) != "..."]]
    argsIn$parameters <- fit$pl
    argsIn$random <- unique(names(obj0$env$par[obj0$env$random]))
    ## Add referencepointSet
    argsIn$data$referencepoint <- list(nYears = 100,
                                 aveYears = as.numeric(c(length(argsIn$data$years)-2:1)),
                                 selYears = as.numeric(c(length(argsIn$data$years)-1)),
                                     Fsequence = seq(1e-5,4, len = 200))

    args <- argsIn
    ## Remove referencepoint parameters from map
    map0 <- args$map
    fix <- setdiff(names(args$parameters), args$random)
    args$map <- lapply(args$parameters[fix], function(x)factor(x*NA))

    ## Referencepoints to estimate
    args$map$logScaleFmsy <- NULL
    args$map$logScaleF01 <- NULL
    args$map$logScaleFmax <- NULL
    args$map$logScaleFcrash <- NULL
    args$map$logScaleF35 <- NULL
    if(fit$conf$stockRecruitmentModelCode %in% c(61,63))
        args$map$logScaleFlim <- NULL

    
    obj <- do.call(TMB::MakeADFun, args)
    
    obj$fn(obj$par)
    obj$gr(obj$par)

    opt <- nlminb(obj$par, obj$fn, obj$gr)

    ## Get standard errors
    args <- argsIn
    args$map$logScaleFmsy <- NULL
    args$map$logScaleF01 <- NULL
    args$map$logScaleFmax <- NULL
    args$map$logScaleFcrash <- NULL
    args$map$logScaleF35 <- NULL
    if(fit$conf$stockRecruitmentModelCode %in% c(61,63))
        args$map$logScaleFlim <- NULL

    opl <- obj$env$parList(par = opt$par)
    
    args$parameters$logScaleFmsy <- opl$logScaleFmsy
    args$parameters$logScaleF01 <- opl$logScaleF01
    args$parameters$logScaleFmax <- opl$logScaleFmax
    args$parameters$logScaleFcrash <- opl$logScaleFcrash
    args$parameters$logScaleF35 <- opl$logScaleF35
    if(fit$conf$stockRecruitmentModelCode %in% c(61,63))
        args$map$logScaleFlim <- NULL


    obj2 <- do.call(TMB::MakeADFun, args)

    Sigma <- solve(fit$opt$he)

    ## can we use TMB/AD??
    gridx <- which(names(obj2$par) %in% c("logScaleFmsy",
                                          "logScaleFmax",
                                          "logScaleF01",
                                          "logScaleFcrash",
                                          "logScaleF35",
                                          "logScaleFlim"))
    JacAll <- numDeriv::jacobian(function(x)obj2$gr(x)[gridx],obj2$par)

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
                ssdr = ssdr
                )
                              


    class(res) <- "sam_referencepoints"
    return(res)
}

   
