
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
       
       srfit <- function(ssb, year = NA){
           if(fit$conf$stockRecruitmentModelCode %in% c(0)){
               val <- NA
               valsd <- NA
               pisig <- NA
           }else if(fit$conf$stockRecruitmentModelCode %in% c(3)){
               if(is.na(year) && length(year) == 1)
                   stop("A single year must be given")
               indx <- cut(year,c(-Inf,fit$conf$constRecBreaks,Inf))
               val <- fit$pl$rec_pars[indx]
               g <- matrix(0,1,nlevels(indx))
               g[1,indx] <- 1
               valsd <- as.vector(sqrt(g %*% covar %*% t(g)))
               pisig <- exp(fit$pl$logSdLogN[fit$conf$keyVarLogN[1]+1])
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
               g <- matrix(utils::head((1/v$Recruits) * v$Gradient,-1), 1)
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
           ssbMax <- max(c(max(S), par("usr")[2]))
       }else{
           ssbMax <- max(S)
       }
       ssb <- seq(1e-5,
                  ssbMax, ## max(S),
                  len = 2000)
       if(fit$conf$stockRecruitmentModelCode == 3){
           brks <- c(min(fit$data$years),ceiling(fit$conf$constRecBreaks),max(fit$data$years))
           labels <- levels(cut(0,brks,dig.lab=4))
           mid <- head(brks,-1) + diff(brks)/2
           tabList <- lapply(as.list(mid), function(y) sapply(ssb,getSSBtab, year = y))
           names(tabList) <- labels
       }else{
           tabList <- list(sapply(ssb, getSSBtab))
       }

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
           if(fit$conf$stockRecruitmentModelCode %in% c(90,91,92))
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
