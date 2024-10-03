

logSRR <- function(logssb, rec_pars, stockRecruitmentModelCode,
                   constRecBreaks = numeric(0), year = NA_real_, lastR = NA_real_){
    N <- length(logssb)
    year <- rep(year, length.out = N)
    lastR <- rep(lastR, length.out = N)
    .Call(C_logSRR, logssb, rec_pars, stockRecruitmentModelCode, constRecBreaks, year, lastR)
}

stockRecruitmentModelR <- function(logssb, rec_pars, stockRecruitmentModelCode,
                                   constRecBreaks = numeric(0), year = NA_real_, lastR = NA_real_){
    if(length(logssb) != 1 || length(year) != 1 || length(lastR) != 1)
        stop("logssb, year, and lastR must be numeric scalars.")
    v <- .Call(C_stockRecruitmentModelR, logssb, rec_pars, stockRecruitmentModelCode, constRecBreaks, lastR, year)
    list(logRecruits = v$logRecruits,
         Gradient_recpars = utils::head(v$Gradient,-1),
         Gradient_logssb = utils::tail(v$Gradient,1))
}



##' Add stock-recruitment curve to srplot
##'
##' @param fit Object to show SR-curve for
##' @param CI Add confidence intervals?
##' @param col Color of fitted line
##' @param cicol Color of confidence intervals
##' @param plot Add the curve to a plot?
##' @param PI Add prediction intervals?
##' @param picol Color of prediction interval line
##' @param pilty Line type of prediction interval line
##' @param ... not used
##' @seealso srplot
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
##' @param year Show recruitment calculated conditional on this year (for recruitment functions that depend on year)
##' @param lastR Show recruitment calculated conditional on this previous recruitment (for recruitment functions that depend on recruitment the previous year)
##' @importFrom utils head
##' @export
addRecruitmentCurve.sam <- function(fit,
                    CI = TRUE,
                    col = rgb(0.6,0,0),
                    cicol = rgb(0.6,0,0,0.3),
                    plot = TRUE,
                    PI = FALSE,
                    picol = rgb(0.6,0,0),
                    pilty = 2,
                    year = NA_real_,
                    lastR = NA_real_,
                    ...){
       X <- summary(fit)
       R <- X[, 1]
       #S <- X[, 4]
       S <- rbtable(fit)[,1]
       ## cf <- fit$sdrep$cov.fixed
       ## covEst <- cf[rownames(cf) %in% c("rec_pars"), colnames(cf) %in% c("rec_pars"), drop = FALSE]
       covEst <- fit$sdrep$covRecPars
       m <- fit$obj$env$map$rec_pars
       if(is.null(m)){
           covar <- covEst
       }else{
           covar <- covEst[m,m, drop = FALSE]
           covar[is.na(covar)] <- 0
       }

       if(fit$conf$stockRecruitmentModelCode %in% c(-2, -1, 0, 3, 62)){
           warning("addRecruitmentCurve is not intended for time series models.")
       }

       srfit <- function(logssb, year = NA_real_, lastR = NA_real_){
           v <- stockRecruitmentModelR(logssb,
                                       fit$pl$rec_pars,
                                       fit$conf$stockRecruitmentModelCode,
                                       fit$conf$constRecBreaks,
                                       year,
                                       lastR)
           val <- v$logRecruits
           g <- matrix(v$Gradient_recpars, ncol = 1)
           valsd <- as.vector(sqrt(t(g) %*% covar %*% g))
           pisig <- exp(fit$pl$logSdLogN[fit$conf$keyVarLogN[1]+1])
           res <- exp(val)
           attr(res,"ci_low") <- exp(val - 2 * valsd)
           attr(res,"ci_high") <- exp(val + 2 * valsd)
           attr(res,"pi_low") <- exp(val - 2 * pisig)
           attr(res,"pi_high") <- exp(val + 2 * pisig)
           return(res)
       }
       
       getSSBtab <- function(x, year = NA) {
           tmp <- srfit(x, year)
           c(Estimate = as.vector(tmp),
             CIlow = as.vector(attr(tmp,"ci_low")),
             CIhigh = as.vector(attr(tmp,"ci_high")),
             PIlow = as.vector(attr(tmp,"pi_low")),
             PIhigh = as.vector(attr(tmp,"pi_high"))
             )
       }


       if(plot){
           ssbMax <- par("usr")[2]
       }else{
           ssbMax <- max(S)
       }
       ssb <- seq(1e-5,
                  ssbMax, ## max(S),
                  len = 2000)
       tabList <- list(sapply(log(ssb), getSSBtab))
       
       transp <- Vectorize(function(col){
           arg <- as.list(grDevices::col2rgb(col,TRUE)[,1])
           arg$alpha <- 0.1 * 255
           arg$maxColorValue <- 255
           do.call(grDevices::rgb, arg)
       })
       
       if(plot){           
           if(CI)               
               for(i in 1:length(tabList)){
                   polygon(c(ssb, rev(ssb)),
                           c(tabList[[i]]["CIlow",],
                             rev(tabList[[i]]["CIhigh",])),
                           col = ifelse(length(tabList) == 1,cicol,transp(i)),
                           border = NA)
                   }
           if(fit$conf$stockRecruitmentModelCode %in% c(90,91,92,93,290,293))
               abline(v = exp(fit$conf$constRecBreaks), col = "grey")
           for(i in 1:length(tabList))
               lines(ssb,tabList[[i]]["Estimate",], col = ifelse(length(tabList) == 1,col,i), lwd = 3)
           if(length(tabList) > 1){
               for(i in 1:length(tabList)){
                   text(par("usr")[1],tabList[[i]]["Estimate",1],names(tabList)[[i]],pos=4, col = ifelse(length(tabList) == 1,cicol,i))
               }
           }
           if(PI){
               for(i in 1:length(tabList)){               
                   lines(ssb,tabList[[i]]["PIlow",], col = ifelse(length(tabList) == 1,picol,i), lwd = 3, lty = pilty)
                   lines(ssb,tabList[[i]]["PIhigh",], col = ifelse(length(tabList) == 1,picol,i), lwd = 3, lty = pilty)
               }
           }
       }
       invisible(srfit)
} 
