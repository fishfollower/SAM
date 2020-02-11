
##' Add stock-recruitment curve to srplot
##'
##' @param fit Object to show SR-curve for
##' @param CI Add confidence intervals?
##' @param col Color of fitted line
##' @param cicol Color of confidence intervals
##' @param plot Add the curve to a plot?
##' @param ... not used
##' @seealso srplot
##' @author Christoffer Moesgaard Albertsen
##' @export
addRecruitmentCurve <- function(fit,
                                CI = TRUE,
                                col = rgb(0.6,0,0),
                                cicol = rgb(0.6,0,0,0.3),
                                plot = TRUE,
                                PI = FALSE,
                                picol = rgb(0.6,0,0),
                                pilty = 2,
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
                    plot = TRUE,
                    PI = FALSE,
                    picol = rgb(0.6,0,0),
                    pilty = 2,
                    ...){
       X <- summary(fit)
       R <- X[, 1]
       S <- X[, 4]
       cf <- fit$sdrep$cov.fixed
       covEst <- cf[rownames(cf) %in% c("rec_pars"), colnames(cf) %in% c("rec_pars"), drop = FALSE]
       m <- fit$obj$env$map$rec_pars
       if(is.null(m)){
           covar <- covEst
       }else{
           covar <- covEst[m,m, drop = FALSE]
           covar[is.na(covar)] <- 0
       }

       if(fit$conf$stockRecruitmentModelCode %in% c(0, 3)){
           warning("addRecruitmentCurve is not implemented for this recruitment type.")
       }
       
       srfit <- function(ssb){
           if(fit$conf$stockRecruitmentModelCode %in% c(0, 3)){
               val <- NA
               valsd <- NA
               pisig <- NA
           }else if(fit$conf$stockRecruitmentModelCode %in% c(62)){
               val <- (fit$pl$rec_pars[1])
               g <- matrix(c(1, 0),1)
               valsd <- as.vector(sqrt(g %*% covar %*% t(g)))
               rho <- 2 / ( 1 + exp( -fit$pl$rec_pars[2] ) ) - 1
               pisig <- exp(fit$pl$logSdLogN[fit$conf$keyVarLogN[1]+1]) / sqrt(1 - rho)
           }else if(fit$conf$stockRecruitmentModelCode %in% c(90,91,92)){
               if(fit$conf$stockRecruitmentModelCode == 90){
                   srfun <- function(par){
                       v <- log(ssb) + .Call("ibcdsplineR",
                                             log(ssb),
                                             fit$conf$constRecBreaks,
                                             par)
                       v
                   }
               }else if(fit$conf$stockRecruitmentModelCode == 91){
                   srfun <- function(par){
                       v <- log(ssb) + .Call("ibcsplineR",
                                             log(ssb),
                                             fit$conf$constRecBreaks,
                                             par)
                       v
                   }
               }else{
                   srfun <- function(par){
                       v <- log(ssb) + .Call("bcsplineR",
                                             log(ssb),
                                             fit$conf$constRecBreaks,
                                             par)
                       v
                   }
               }
               val <- srfun(fit$pl$rec_pars)
               g <- matrix(grad(srfun, fit$pl$rec_pars),1)
               valsd <- as.vector(sqrt(g %*% covar %*% t(g)))
               pisig <- exp(fit$pl$logSdLogN[fit$conf$keyVarLogN[1]+1])
           }else{
               v <- .Call("stockRecruitmentModelR",
                          ssb,
                          fit$pl$rec_pars,
                          fit$conf$stockRecruitmentModelCode)
               val <- log(v$Recruits)
               g <- matrix(head((1/v$Recruits) * v$Gradient,-1), 1)
               valsd <- as.vector(sqrt(g %*% covar %*% t(g)))
               pisig <- exp(fit$pl$logSdLogN[fit$conf$keyVarLogN[1]+1])
           }
           res <- exp(val)
           attr(res,"ci_low") <- exp(val - 2 * valsd)
           attr(res,"ci_high") <- exp(val + 2 * valsd)
           attr(res,"pi_low") <- exp(val - 2 * pisig)
           attr(res,"pi_high") <- exp(val + 2 * pisig)
           return(res)
       }

       ssb <- seq(1e-5, max(S), len = 2000)
       tab <- sapply(ssb, function(x) {
           tmp <- srfit(x)
           c(Estimate = as.vector(tmp),
             CIlow = as.vector(attr(tmp,"ci_low")),
             CIhigh = as.vector(attr(tmp,"ci_high")),
             PIlow = as.vector(attr(tmp,"pi_low")),
             PIhigh = as.vector(attr(tmp,"pi_high"))
             )
       })
       if(plot){
           if(CI)
               polygon(c(ssb, rev(ssb)),
                       c(tab["CIlow",],
                         rev(tab["CIhigh",])),
                       col = cicol,
                       border = NA)
           if(fit$conf$stockRecruitmentModelCode %in% c(90,91,92))
               abline(v = exp(fit$conf$constRecBreaks), col = "grey")
           lines(ssb,tab["Estimate",], col = col, lwd = 3)
           if(PI){
               lines(ssb,tab["PIlow",], col = picol, lwd = 3, lty = pilty)
               lines(ssb,tab["PIhigh",], col = picol, lwd = 3, lty = pilty)               
           }
       }
       invisible(srfit)
} 
