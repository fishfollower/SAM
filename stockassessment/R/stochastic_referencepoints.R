
.refpointSMethodEnum <- c("Mean"=0,
                          "Median"=1,
                          "Mode"=2,
                          "Quantile"=1,
                          "LOESS"=3)

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

.refpointSMethodParser <- function(x, monotone = c("n","d","i"), positive = NA, ...){
    if(length(x) > 1)
        stop("Only one at a time!")
    monotone <- match.arg(monotone)
    ## Allow shorthand...
    typePatterns <- list(Mean = "^mean$",
                         Median = "^median$",
                         Mode = "^mode$",
                         Quantile = "^q0\\.[[:digit:]]+$",
                         LOESS = "^loess$")
    xvalPatterns <- list(Mean = NA,
                         Median = 0.5,
                         Mode = NA,
                         Quantile = c("(q)(0\\.[[:digit:]]+)","\\2"),
                         LOESS = 0.2)
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
           xVal = xVal,
           monotone = monotone,
           positive = positive),
      list(...))        
}

removeIntercept <- function(object){
    f <- update(object,~(.)-1)
    if(is(object,"terms")){
        attributes(f) <- attributes(object)
        attr(f,"intercept") <- 0
    }
    return(f)    
}
##' @importFrom stats model.frame model.matrix terms
.refpointSCurveFit <- function(F,C,MT){
    if(MT$useLogF)
        F <- log(F)
    if(MT$methodType == 3){ ## LOESS
        return(loess(C~F,data.frame(C=C,F=F)))
    }
    if(MT$penalize){
        if(!require(RTMB))
            stop("Penalization requires RTMB")
        if(!require(RTMBconvenience))
            stop("Penalization requires RTMBconvenice, available at github.com/calbertsen/RTMBconvenience")
        ## Update formula to give knots
        txt <- as.character(MT$formula)[2]
        splType <- gsub("^([^\\(]+)(\\([[:space:]]*F[[:space:]]*,)([[:space:]]*[[:digit:]]+[[:space:]]*)(\\)$)","\\1",txt)
        Nknots <- as.integer(gsub("^([^\\(]+)(\\([[:space:]]*F[[:space:]]*,)([[:space:]]*[[:digit:]]+[[:space:]]*)(\\)$)","\\3",txt))
        MT$formula <- formula(paste(c("~",sprintf("%s(F,knots=c(%s),Boundary.knots=c(%s,%s),qk=TRUE)",splType,paste(unique(head(tail(seq(min(F),max(F),length.out=Nknots),-1),-1)),collapse=","),min(F),max(F))),collapse = " "))        
        ##MT$formula <- formula(paste(c("~",sprintf("splines::bs(F,knots=c(%s))",paste(seq(min(F),max(F),length.out=Nknots),collapse=","))),collapse = " "))        
    }
    mf <- stats::model.frame(MT$formula, data.frame(F = F,C=C),na.action=stats::na.pass)
    ## mf <- as.data.frame(lapply(mf0, function(x){
    ##     x[!is.finite(x)] <- 0
    ##     x
    ## }))
    ## attributes(mf) <- attributes(mf0)
    X <- stats::model.matrix(MT$formula,mf)
    Xni <- stats::model.matrix(removeIntercept(terms(mf)),mf)
    X[!is.finite(X)] <- 0
    ## Xsd <- stats::model.matrix(~ibc(F,3),data.frame(F=F))
    if(MT$methodType == 0){       # Mean
        Loss <- function(Obs,Pred) ((Obs-Pred))^2 #+ log(Sd)
    }else if(MT$methodType == 1){       # Quantiles
        if(MT$penalize){
            rho <- function(x,q) RTMBconvenience::quantreg_loss(x,q)
        }else{
            rho <- function(x,q) x * (q - (x<0))
        }
        Loss <- function(Obs,Pred) rho((Obs-Pred), MT$xVal) #+ log(Sd)    
    }else if(MT$methodType == 2){       # Mode
        stop("method not implemented yet")
    }
    if(is.na(MT$positive)){
        Cuse <- C
    }else if(MT$positive){
        Cuse <- log(C)
    }else{
        Cuse <- (log(-C))
    }
    p0 <- RTMBconvenience::undim(svd_solve(t(X) %*% X) %*% t(X) %*% Cuse)
    ## Vuse <- log((Cuse-X %*% p0)^2)
    ## p0sd <- RTMBconvenience::undim(svd_solve(t(Xsd) %*% Xsd) %*% t(Xsd) %*% Vuse)    
    if(MT$penalize){
        cat("Fitting penalized GAM-like spline...\n")
        ## Get S matrix (assume first col of X is an intercept)
        S <- local({
            eps <- 1e-4
            grid  <- seq(min(F), max(F), length.out = 1000)
            dx <- diff(grid[1:2])
            X0 <- stats::model.matrix(removeIntercept(terms(mf)),data.frame(F=grid))#[,-1,drop=FALSE]
            X0[!is.finite(X0)] <- 0
            Xp <- stats::model.matrix(removeIntercept(terms(mf)),data.frame(F=grid+eps))#[,-1,drop=FALSE]
            Xp[!is.finite(Xp)] <- 0
            Xm <- stats::model.matrix(removeIntercept(terms(mf)),data.frame(F=grid-eps))#[,-1,drop=FALSE]
            Xm[!is.finite(Xm)] <- 0
            ## Central diffrerences
            X2 <- (Xp - 2 * X0 + Xm) / eps^2
            ## Num int
            t(X2) %*% X2 * dx
        })
        XtX <- t(Xni) %*% Xni
        ## Outer estimation of smoothing penalty using generalized cross validation criteria       
        if(MT$RTMB){
            fnOut <- function(par){
                l <- exp(par$logl)   
                ## mu <- par$mu
                ## Inner estimation of p given mu and penalty
                fn <- MakeTape(function(pIn){
                    l <- exp(pIn[1])
                    p <- pIn[seq_along(p0)+1]
                    ##pSd <- pIn[seq_along(p0sd)+1 + length(p)]
                    if(MT$monotone == "d" || MT$monotone == "pd"){
                        p <- -exp(p[-1])
                    }else if(MT$monotone == "i"){
                        p <- exp(p[-1])
                    }
                    pred <- (RTMBconvenience::undim((RTMB::matrix(X,nrow(X),ncol(X)) %*% p))) ##c(AD(mu),p))))
                    ##predSd <- exp(RTMBconvenience::undim((RTMB::matrix(Xsd,nrow(Xsd),ncol(Xsd)) %*% pSd)))
                    jac <- 0
                    ## if(!is.na(MT$positive)){
                    ##     if(MT$positive){
                    ##         pred <- exp(pred)
                    ##         jac <- log(C)
                    ##     }else{
                    ##         pred <- -exp(pred)
                    ##         jac <- log(-C)
                    ##     }
                    ## }
                    ## if(MT$monotone == "n"){
                    pen <- (RTMBconvenience::undim(t(p[-1]) %*% RTMB::matrix(S,nrow(S),ncol(S)) %*% p[-1]))
                    ## }else{
                    ##     pen <- (undim(t(p[-(1]) %*% RTMB::matrix(S[-1,-1],nrow(S)-1,ncol(S)-1) %*% p[-1]))
                    ## }
                    sum(Loss((Cuse),pred)) + l * pen
                },c(0,p0*0))
                fnnewt <- fn$newton(2:(length(p0)+1), trace=0,on_failure_return_nan=FALSE,sparse=TRUE, simplify=TRUE,maxit=10000,on_failure_give_warning=FALSE)
                ## numeric(ncol(X)+1)
                beta0 <- c(fnnewt(par$logl)) ## Does not include logl
                beta <- beta0[seq_along(p0)]
                REPORT(beta) ## Needs to be reported before transformations!
                ##betaSd <- beta0[seq_along(p0sd) + length(p0)]
                nll <- fn(c(par$logl,beta0))
                REPORT(nll)
                jac <- 0
                if(MT$monotone == "d" || MT$monotone == "pd"){
                    beta[-1] <- -exp(beta[-1])
                    jac <- log(C)
                }else if(MT$monotone == "i"){
                    beta[-1] <- exp(beta[-1])
                    jac <- log(-C)
                }
                pred <- RTMBconvenience::undim((RTMB::matrix(X,nrow(X),ncol(X)) %*% beta))
                ##predSd <- exp(RTMBconvenience::undim((RTMB::matrix(Xsd,nrow(Xsd),ncol(Xsd)) %*% betaSd)))
                ## if(!is.na(MT$positive)){
                ##         if(MT$positive){
                ##             pred <- exp(pred)
                ##         }else{
                ##             pred <- -exp(pred)
                ##         }
                ## }
                if(MT$monotone == "n"){               
                    A   <- solve(XtX + l * S, XtX)
                }else{
                    A <- solve(XtX[-1,-1] + l * S[-1,-1], XtX[-1,-1])
                }
                edf <- sum(diag(A))   # = trace(H)
                F <- sum(diag(A %*% A)) / length(C)
                ## REPORT(betaSd)
                ##spl <- cbind(seq(0,1,len=1000),RTMB::matrix(basis_fun(seq(0,1,len=1000),knots),nrow=1000) %*% beta + par$mu)
                ##REPORT(spl)
                REPORT(edf)            
                ## length(C) * sum((C-pred)^2) / (length(C) - edf)^2
                G <- length(C) * sum(Loss(Cuse,pred)) / (length(C) - MT$rho * edf)^2
                MT$gamma * G + (1 - MT$gamma) * F
            }
            ## if(is.na(MT$positive)){
            ##     mu0 <- mean(C)
            ## }else if(MT$positive){
            ##     mu0 <- log(mean(range((C))))
            ## }else{
            ##     mu0 <- log(mean(range((-C)))) 
            ## }
            par0 <- list(logl=0)#,mu=mu0 * 0)
            obj <- RTMB::MakeADFun(fnOut,par0)
            ll <- (-20:20)
            vv <- sapply(ll, obj$fn)          
            ## Grid search
            ## Check gradient
            gg <- obj$gr(ll[which.min(vv)])
            if(!is.finite(gg)){
                opt0 <- nlminb(ll[which.min(vv)],obj$fn)
            }else{
                opt0 <- nlminb(ll[which.min(vv)],obj$fn,obj$gr,control=list(eval.max=10000,iter.max=10000,trace=0))
                opt0 <- nlminb(opt0$par,obj$fn,obj$gr,obj$he,control=list(eval.max=10000,iter.max=10000,trace=0))
            }
            ## opt0 <- nlminb(obj$par,obj$fn,obj$gr, obj$he)
            rp <- obj$report(opt0$par)
            opt <- opt0
            opt$par <- rp$beta
            ## opt$parSd <- rp$betaSd
            ##opt$spl <- rp$spl
            opt$edf <- rp$edf
            opt$objective <- opt0$objective
            opt$p0 <- p0
            attr(opt,"orig") <- opt0
            opt$gridSearch <- cbind(F=ll,Cuse=vv)
        }else{
            fnOut <- function(par){
                l <- exp(par[1])   
                ## mu <- par$mu
                ## Inner estimation of p given mu and penalty
                fn <- function(pIn){
                    p <- pIn[seq_along(p0)]
                    ## pSd <- pIn[seq_along(p0sd)+length(p0)]
                    ## l <- exp(pIn[1])
                    if(MT$monotone == "d" || MT$monotone == "pd"){
                        p <- -exp(p[-1])
                    }else if(MT$monotone == "i"){
                        p <- exp(p[-1])
                    }
                    pred <- RTMBconvenience::undim(X %*% p) ##c(AD(mu),p))))
                    ## predSd <- exp(RTMBconvenience::undim(Xsd %*% pSd))
                    ## jac <- 0
                    ## if(!is.na(MT$positive)){
                    ##     if(MT$positive){
                    ##         pred <- exp(pred)
                    ##         jac <- log(C)
                    ##     }else{
                    ##         pred <- -exp(pred)
                    ##         jac <- log(-C)
                    ##     }
                    ## }
                    ## if(MT$monotone == "n"){
                    pen <- (RTMBconvenience::undim(t(p[-1]) %*% RTMB::matrix(S,nrow(S),ncol(S)) %*% p[-1]))
                    ## }else{
                    ##     pen <- (undim(t(p[-(1]) %*% RTMB::matrix(S[-1,-1],nrow(S)-1,ncol(S)-1) %*% p[-1]))
                    ## }
                    sum(Loss((Cuse),pred)) + l * pen
                }
                optIn <- nlminb(c(p0),fn,control=list(eval.max=10000,iter.max=10000,trace=1))
                beta0 <- optIn$par
                beta <- beta0[seq_along(p0)]
                betaNonTransformed <- beta
                ## betaSd <- beta0[seq_along(p0sd)+length(p0)]
                nll <- optIn$objective                
                if(MT$monotone == "d" || MT$monotone == "pd"){
                    beta[-1] <- -exp(beta[-1])
                }else if(MT$monotone == "i"){
                    beta[-1] <- exp(beta[-1])
                }
                pred <- as.vector(X %*% beta)
                ## predSd <- exp(RTMBconvenience::undim(Xsd %*% betaSd))
                ## if(!is.na(MT$positive)){
                ##         if(MT$positive){
                ##             pred <- exp(pred)
                ##         }else{
                ##             pred <- -exp(pred)
                ##         }
                ## }
                if(MT$monotone == "n"){               
                    A   <- solve(XtX + l * S, XtX)
                }else{
                    A <- solve(XtX[-1,-1] + l * S[-1,-1], XtX[-1,-1])
                }
                edf <- sum(diag(A))   # = trace(H)
                F <- sum(diag(A %*% A)) / length(C)                
                ## length(C) * sum((C-pred)^2) / (length(C) - edf)^2
                ## r <- length(C) * sum(Loss(Cuse,pred,predSd)) / (length(C) - edf)^2
                G <- length(C) * sum(Loss(Cuse,pred)) / (length(C) - MT$rho * edf)^2
                r <- (MT$gamma * G + (1 - MT$gamma) * F)
                attr(r,"edf") <- edf
                attr(r,"beta") <- betaNonTransformed
                ## attr(r,"betaSd") <- betaSd
                attr(r,"nll") <- nll
                r
            }
            ## Grid search
            ll <- (-5:30)
            vv <- sapply(ll, fnOut)
            opt0 <- nlminb(ll[which.min(vv)],fnOut,control=list(eval.max=10000,iter.max=10000,trace=1))
            ##opt0 <- nlminb(0, fnOut)
            opt <- opt0
            vx <- fnOut(opt$par)
            opt$par <- attr(vx,"beta")
            ## opt$parSd<- attr(vx,"betaSd")
            ##opt$spl <- rp$spl
            opt$edf <- attr(vx,"edf")
            opt$objective <- opt0$objective
            opt$p0 <- p0
            attr(opt,"orig") <- opt0
        }
    }else{
        cat("Fitting spline...\n")
        fn <- function(par){
            beta <- par[1:(ncol(X))]
            ## betaSd <- par[1:ncol(Xsd) + (ncol(X))]
            if(MT$monotone == "d" || MT$monotone == "pd"){
                beta[-1] <- -exp(beta[-1])
            }else if(MT$monotone == "i"){
                beta[-1] <- exp(beta[-1])
            }
            X2 <- X
            pM <- X2 %*% beta
            ## predSd <- exp((Xsd %*% betaSd))
            ## if(!is.na(MT$positive)){
            ##     if(MT$positive){
            ##         pM <- exp(pM)
            ##     }else{
            ##         pM <- -exp(pM)
            ##     }
            ## }
            sum(Loss(Cuse,pM[]))
        }
        opt <- stats::nlminb(numeric(ncol(X)), fn, control = list(iter.max=10000,eval.max=10000,trace=0))
        opt$parAll <- opt$par
        opt$par <- opt$parAll[1:ncol(X)]
        ##opt$parSd <- opt$parAll[1:ncol(Xsd) + ncol(X)]
    }
    cat("Done fitting spline...\n")
    attr(opt,"terms") <- stats::terms(mf)
    attr(opt,"MT") <- MT
    class(opt) <- "rpscurvefit"
    opt
}

##' @importFrom stats predict
##' @method predict rpscurvefit
##' @export
predict.rpscurvefit <- function(x,newF,...){
    MT <- attr(x,"MT")
    if(MT$useLogF)
        newF <- log(newF)
    X <- model.matrix(attr(x,"terms"),data.frame(F=newF))
    par <- x$par
    beta <- par[1:ncol(X)]
    if(MT$monotone == "d"){
        beta[-1] <- -exp(beta[-1])
    }else if(MT$monotone == "i"){
        beta[-1] <- exp(beta[-1])
    }
    pM <- X %*% beta
    if(!is.na(MT$positive)){
        if(MT$positive){
            pM <- exp(pM)
        }else{
            pM <- -exp(pM)
        }
    }
    as.numeric(pM)  
}

predict.rpsnocurve <- function(x,newF,...){
    rep(NA,length(newF))
}

.getDoSim <- function(logf1, fit, nYears, aveYears, selYears, pl, logCustomSel = NULL, constraint = "F=%f", deterministicF = TRUE,...){
   if(length(logCustomSel) > 0){
        sel <- exp(logCustomSel)
    }else if(length(selYears) == 0){   
        sel <- NULL
    }else{
        sel <- exp(.logFtoSel(pl$logF, match(selYears, fit$data$years) - 1, fit$conf))
    }
   cat("Using fastFixedF:", grepl("[[:space:]]*F[[:space:]]*=[[:space:]]*%f[[:space:]]*",constraint),"\n")
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

.perRecruitSR_Calc <- function(logf, fit, nYears, aveYears, selYears, pl = fit$pl, ct = 0, logCustomSel = numeric(0), logNinit = fit$pl$logN[,ncol(fit$pl$logN)], DT=0){
    if(length(logCustomSel) > 0){
        sel <- exp(logCustomSel)
    }else{
        sel <- exp(.logFtoSel(pl$logF, selYears, fit$conf))
    }
    .Call(C_perRecruitSR_Calc,
          logFbar = logf,
          tmbdat = fit$obj$env$data,
          pl = pl,
          sel = sel,
          aveYears = aveYears,
          nYears = ifelse(nYears==0,150,nYears),
          CT = ct,
          logNinit = logNinit,
          logNFY = log(ntable(fit)[1,]),
          DT = DT
          )
}

.perRecruitSR_Calc_MakePtr <- function(logf, fit, nYears, aveYears, selYears, pl = fit$pl, ct = 0, logCustomSel = numeric(0), logNinit = fit$pl$logN[,ncol(fit$pl$logN)], DT=0){
    if(length(logCustomSel) > 0){
        sel <- exp(logCustomSel)
    }else{
        sel <- exp(.logFtoSel(pl$logF, selYears, fit$conf))
    }
    .Call(C_MakePtr_perRecruitSR_Calc,
          logFbar = logf,
          tmbdat = fit$obj$env$data,
          pl = pl,
          sel = sel,
          aveYears = aveYears,
          nYears = ifelse(nYears==0,150,nYears),
          CT = ct,
          logNinit = logNinit,
          logNFY = log(ntable(fit)[1,]),
          DT = DT
          )
}

.perRecruitSR_Calc_EvalPtr <- function(ptr, logf, outType = 0, Ntail = 0, q = NA_real_){
     .Call(C_EvalPtr_perRecruitSR_Calc,
          r_ptr = ptr,
          logFbar = logf,
          outType = outType,
          Ntail = Ntail,
          q = q          
          )
}




.perRecruitSR <- function(logf, fit, nYears, aveYears, selYears, pl = fit$pl, ct=0, logCustomSel = NULL, nTail = 1,incpb=NULL, doSim = NULL, label="", constraint = "F=%f", deterministicF=TRUE, tailSummarizer=mean,ncores=1,mc.type, ...){    
    ## pl$missing <- NULL
    ## attr(pl,"what") <- NULL
    ## f2 <- sam.fit(fit$data,fit$conf,pl, run = FALSE)
    ## f2$opt <- list(par = f2$obj$par, objective=NA)
    ## f2$sdrep <- fit$sdrep[c("estY","covY")]
    ## class(f2) <- "sam"
    if(is.null(doSim))
        doSim <- .getDoSim(max(logf)+0.2, fit, nYears, aveYears, selYears, pl, logCustomSel, deterministicF, ...)
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

    
    do.call("rbind",.SAM_lapply(logf, function(lf){
        v <- doSim(sprintf(constraint,exp(lf)))
        logRe <- tailSummarizer(tail(v$logN[1,],nTail))
        logSe <- tailSummarizer(tail(v$logssb,nTail))
        logSPR <- tailSummarizer(tail(v$logEmpiricalSPR,nTail))
        if(ct == 0){
            logYe <- tailSummarizer(tail(v$logCatch,nTail))
            logYPR <- tailSummarizer(tail(v$logEmpiricalYPR,nTail))
        }else if(ct == 1){
            logYe <- tailSummarizer(tail(v$logLand,nTail))
            logYPR <- tailSummarizer(tail(v$logEmpiricalYPR_L,nTail))
        }else{
            logYe <- tailSummarizer(tail(log(exp(v$logCatch)-(v$logLand)),nTail))
            logYPR <- tailSummarizer(tail(v$logEmpiricalYPR_D,nTail))
        }
        logLE <- NA_real_ # tail(v$logLifeExpectancy,nTail)
        logYL <- NA_real_ #tail(v$logYLTF,nTail)

        res <- list(logF = lf,#mean(rep(lf,length.out = nTail)),
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
    }, ncores = ncores, type = mc.type))
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


## TODO:
## [ ] Restructure to do curves for all values in the beginning and only once, save them in regression, and use them to get RPs and derived values
#' @importFrom stats runif predict
.refpointSFitCriteria <- function(rpArgs, pl, MT, fit, nosim, Frange, aveYears, selYears, nYears, catchType, nTail = 1,doSim=NULL,incpb=NULL,label="",constraint="F=%f", deterministicF = TRUE, randomF = TRUE, knots = 5, tailSummarizer = mean, gridSummarizer = median, ncores=1,mc.type="mclapply",defaultSpline="ibc", ...){
    rfv <- function(n,a,b){
        u <- stats::runif(n)
        v1 <- exp(stats::runif(n,log(ifelse(a==0,0.002,a)), log(b)))
        v2 <- stats::runif(n,a,b)
        ifelse(u < 0.25, v1, v2)
    }
    if(randomF){
        ##Fvals <- sort(c(rep(1e-6,100), rfv(nosim,Frange[1],Frange[2])))
        Fvals <- sort(rfv(nosim,Frange[1],Frange[2]))
    }else{
        ##Fvals <- pmax(1e-6,rep(Frange,each = nosim))
        Fvals <- pmax(1e-6,rep(Frange,each = nosim))
    }
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
                            tailSummarizer = tailSummarizer,
                            ncores=ncores,
                            mc.type=mc.type,
                            ...)
    ## Fit curves
    if(!randomF){
        PRvals <- do.call(rbind,lapply(split(PRvals,PRvals$logF),function(x) as.data.frame(lapply(as.list(x),gridSummarizer))))
        ## PRcurves <- 
    }else{
        cat("\nFitting curves...\n")
        doOneCurve <- function(what){
            cat(what,"\n")
            Crit <- exp(PRvals[[what]])
            if(all(is.na(Crit))){
                r <- NA
                class(r) <- "rpsnocurve"
                return(r)
            }
            cutfun <- function(x) x > max(x) * 0#rp$cutoff
            Frng <- range(Fvals[cutfun(Crit)])
            inRng <- function(x,rng) x > rng[1] & x < rng[2]
            ##indx <- inRng(Fvals,Frng)
            indx <- seq_along(Fvals)
            .getCurve <- function(k, Fv, Cv, MT){
                MT$positive <- TRUE
                MT$monotone <- "n"
                MT$formula <- as.formula(bquote(~.(defaultSpline)(F,.(k))))
                if(!MT$penalize){
                    if(what %in% c("logSPR","logSe","logLifeExpectancy")){
                        MT$formula <- as.formula(bquote(~iibc(F,.(k))))
                        MT$monotone <- "d"
                    }else if(what %in% c("logYearsLost")){
                        MT$formula <- as.formula(bquote(~iibc(F,.(k))))
                        MT$monotone <- "i"
                    }
                }else{
                    MT$formula <- as.formula(bquote(~.(defaultSpline)(F,.(k))))                   
                    ## MT$formula <- as.formula(bquote(~splines::bs(F,.(k))))
                }
                CurveFit <- .refpointSCurveFit(Fv, Cv, MT)
                if(MT$methodType == 1){
                    logLik <- length(Fv)*(log(MT$xVal*(1-MT$xVal)) - 1 - log(CurveFit$objective/length(Fv)))
                    npar <- length(CurveFit$par)
                }else if(MT$methodType == 0){
                    logLik <- -CurveFit$objective
                    npar <- length(CurveFit$par)
                }else{
                    logLik <- NA
                    npar <- NA
                }
                attr(CurveFit,"AIC") <- -2 * logLik  + 2 * npar
                attr(CurveFit,"knots") <- k
                CurveFit
            }
            if(MT$methodType == 3){
                CurveFit <- .refpointSCurveFit(Fvals[indx], Crit[indx], MT)
            }else{
                if(is.na(knots)){            
                    candidateCurves <- lapply(3:20,function(k) .getCurve(k, Fv=Fvals[indx], Cv=Crit[indx], MT=MT))
                    canAIC <- sapply(candidateCurves, attr, which = "AIC")
                    cat(canAIC)
                    CurveFit <- candidateCurves[[which.min(canAIC)]]
                }else{
                    CurveFit <- .getCurve(knots, Fv=Fvals[indx], Cv=Crit[indx], MT=MT)
                }
            }
            CurveFit
        }
        PRcurves <- lapply(c("logYPR","logSPR","logSe","logRe","logYe","logLifeExpectancy","logYearsLost"), doOneCurve)
        names(PRcurves) <- c("YPR","SPR","Se","Re","Ye","LifeExpectancy","YearsLost")
        Fseq <- pmax(1e-5,seq(min(Frange),max(Frange),len=1000))
        cat("Calculate regressions...\n")
        regressions <- c(list(F=Fseq),lapply(PRcurves,function(x) predict(x,Fseq)))
    }
    
###### Different for different RP's
    getOneRP <- function(rp){
        Fvals <- exp(PRvals$logF)
        if(rp$rpType == 1){ ## MSY
            ## Crit <- exp(PRvals$logYe)
            ## cutfun <- function(x) x > max(x) * rp$cutoff
            trans <- function(x, report=FALSE, ...){
                v <- exp(x)
                if(report)
                    names(v) <- "MSY"
                v
            }
            fn <- function(x){
                Ye <- predict(PRcurves$Ye, trans(x))
                -log(Ye)
            }
            startVals <- function(Fseq){
                pv <- regressions$Ye
                log(Fseq[which.max(pv)])
            }
        }else if(rp$rpType == 2){ ## MSYRange
            ## Crit <- exp(PRvals$logYe)
            ## cutfun <- function(x) x > max(x) * rp$cutoff
            trans <- function(x, report=FALSE, ...){
                dots <- list(...)
                if("keepMSY" %in% names(dots)){
                    keepMSY <- dots$keepMSY
                }else{
                    keepMSY <- FALSE
                }
                xMSY <- exp(utils::head(x,1))
                x2 <- matrix(utils::tail(x,-1),2)
                x3 <- rbind(exp(log(xMSY) - exp(-x2[1,])),
                            exp(log(xMSY) + exp(x2[2,])))
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
                logYe <- log(predict(PRcurve$Ye,as.vector(xx)))
                sum((tail(p,-1) - rep(rp$xVal,each=2) * p[1])^2) - p[1]
            }
            startVals <- function(Fseq){
                pv <- regressions$Ye
                fmsy <- Fseq[which.max(pv)]
                c2 <- sapply(rp$xVal, function(xx) (pv - xx*max(pv))^2)
                f0 <- apply(c2,2,function(cc){
                    fl <- Fseq[Fseq < fmsy][which.min(cc[Fseq < fmsy])]
                    fu <- Fseq[Fseq > fmsy][which.min(cc[Fseq > fmsy])]
                    ##c(log(fmsy-fl),log(fu-fmsy))
                    c(-log(log(fmsy)-log(fl)), log(log(fu)-log(fmsy)))
                })
                c(log(fmsy),f0)
            }
        }else if(rp$rpType == 3){ ## Max
            stop("Not ready")
            Crit <- exp(PRvals$logYPR)
            cutfun <- function(x) rep(TRUE,length(x))
            trans <- function(x, report=FALSE, ...){
                v <- exp(x)
                if(report)
                    names(v) <- "Max"
                v
            }
            fn <- function(x) -predict(CurveFit,trans(x))
            startVals <- function(Fseq,pv) log(Fseq[which.max(pv)])                        
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
            startVals <- function(Fseq,pv) sapply(rp$xVal, function(xv) log(Fseq[which.min((pv - xv * pv[1])^2)]))
        }else if(rp$rpType == 5){ ## xSPR
            stop("Not ready")
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
            startVals <- function(Fseq,pv) sapply(rp$xVal, function(xv) log(Fseq[which.min((pv - xv * pv[1])^2)]))
        }else if(rp$rpType == 6){ ## xB0
            stop("Not ready")
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
            startVals <- function(Fseq,pv) sapply(rp$xVal, function(xv) log(Fseq[which.min((pv - xv * pv[1])^2)]))
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
            return(((PRvals$logSe) - 0)^2)
        }else if(rp$rpType == 12){ ## Lim
            stop("Reference point type not implemented yet")
        }else{
            stop("Reference point type not implemented")
        }
        opt <- nlminb(startVals(regressions$F), fn,control=list(iter.max=10000,eval.max=10000))
        res <- trans(opt$par, report = TRUE)
        ## attr(res,"curve_fit_list") <- candidateCurves
        ## attr(res,"curve_fit_aiclist") <- canAIC
        ## attr(res,"curve_fit_opt") <- CurveFit
        ## attr(res,"curve_fit") <- cbind(F=Fseq,Criteria=pv)
        res
    }
    getDerivedValues <- function(f){
        if(!is.function(MT$derivedSummarizer)){        #Fit
            v <- sapply(PRcurves[c("YPR","SPR","Se","Re","Ye","LifeExpectancy","YearsLost")], function(x) predict(x, f))
            names(v) <- paste0("log",names(v))
            return(v)
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
    F <- Reduce("c",(Forig))
    Curves <- lapply(Forig, attr, which = "curve_fit")
    CurveFits <- lapply(Forig, attr, which = "curve_fit_opt")
    CurveFitsList <- lapply(Forig, attr, which = "curve_fit_list")
    CurveFitsAIC <- lapply(Forig, attr, which = "curve_fit_aiclist")
    ##names(Curves) <- names(rpArgs)
    D <- sapply(F, getDerivedValues)    
    res <- rbind(logF=unname(F),D)
    #colnames(res) <- Reduce("c",lapply(Forig,names))
    rownames(res) <- gsub("^log","",rownames(res))
    ## Fseq <- seq(min(Fvals),max(Fvals),len=200)
    ## GraphVals <- rbind(logF=Fseq,sapply(Fseq, getDerivedValues))
    ## colnames(GraphVals) <- Fseq
    ## rownames(GraphVals) <- gsub("^log","",rownames(GraphVals))
    list(Estimates = res,
         ## GraphVals = GraphVals,
         Fvals = Fvals,
         PRvals = PRvals,
         Curves = PRcurves,
         regressions = regressions
         ## CurveFits = CurveFits,
         ## CurveFitsList = CurveFitsList,
         ## CurveFitsAIC = CurveFitsAIC
         )
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
##' Note that this function is highly experimental.
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
##' @param knots Number of knots to use
##' @param nosim_ci Number of simulations for bootstrap confidence intervals
##' @param derivedSummarizer Function to summarize derived per-recruit values
##' @param nTail Number of years from the simulation to include in calculations
##' @param Fsequence F sequence to explore
##' @param run run it?
##' @param DT ...
##' @param equilibriumMethod method to use
##' @param ncores Number of cores
##' @return reference point object
##' @export
stochasticReferencepoints <- function(fit,
                                       referencepoints,
                                       ...){
    UseMethod("stochasticReferencepoints")
}




##' @rdname stochasticReferencepoints
##' @method stochasticReferencepoints sam
##' @importFrom stats lm
##' @export
stochasticReferencepoints.sam <- function(fit,
                                          referencepoints,
                                          method = "Median",
                                          catchType = "catch",
                                          nYears = 100,
                                          Frange = c(0,2),
                                          nosim = 200,
                                          aveYears = max(fit$data$years)+(-4:0),
                                          selYears = max(fit$data$years)+(-4:0),
                                          newton.control = list(),
                                          seed = .timeToSeed(),
                                          knots = NA,
                                          penalizeSpline = TRUE,
                                          defaultSpline = "ibc",
                                          logFspline = FALSE,
                                          splineGamma = ifelse(nosim < 100,0.2,0.3),
                                          splineRho = ifelse(nosim < 100,1.3,2),
                                          useRTMB = TRUE,
                                          nosim_ci = 0,
                                          derivedSummarizer = NA,
                                          nTail = 1,
                                          constraint = "F=%f",
                                          deterministicF = FALSE,
                                          processNoiseF = FALSE,
                                          Fsequence = seq(min(Frange),max(Frange), len = 50),
                                          run = TRUE,
                                          DT = 0,
                                          equilibriumMethod = c("ES","EC","AD"),
                                          ncores = 1,
                                          tailSummarizer = mean,
                                          gridSummarizer = median,
                                          mc.type = "mclapply",
                                          ...){

    equilibriumMethod <- match.arg(equilibriumMethod)
    if(equilibriumMethod == "EC"){

        if(!all(diff(Fsequence) > 0) || !all(Fsequence >= 0))
            stop("Values of Fsequence must be positive and increasing.")
        if(!isTRUE(all.equal(Fsequence[1],0, check.attributes = FALSE, use.names = FALSE)))
            warning("The first value of Fsequence should be 0.")

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

        typePatterns <- list(Median = "^median$",
                             Mean = "^mean$",
                             Mode = "^mode$",
                             Quantile = "^q0\\.[[:digit:]]+$")
        mIndx <- which(sapply(typePatterns, function(p) grepl(p, tolower(method))))
        if(length(mIndx) == 0)
            stop(sprintf("Error in method specification. %s not recognized.",typePatterns))
        mIndx <- mIndx[1]
        if(mIndx == 4){            
            q <- as.numeric(gsub("(q)(0\\.[[:digit:]]+)","\\2",tolower(method)))
        }else{
            q <- NA_real_
        }
        
    ## Parse input reference points
     rpArgs <- Reduce(.refpointMerger,
                         lapply(referencepoints, .refpointParser, nYears = nYears, aveYears = aveYears, selYears = selYears, logCustomSel = numeric(0), catchType = catchType - 1,logN0=fit$pl$logN[,ncol(fit$pl$logN)],stochasticType=mIndx,q=q, DT=DT, Ntail = nTail),
                         list())
    ## Add starting values    
    rpArgs <- lapply(rpArgs, .refpointStartingValue, fit = fit, Fsequence = Fsequence)
    ## Add Fsequence for plotting
    rp0 <- list(rpType = -1,
                xVal = log(Fsequence),
                nYears = nYears,
                aveYears = aveYears,
                selYears = selYears,
                logCustomSel = numeric(0),
                catchType = catchType - 1,
                logF0 = log(Fsequence),
                logN0=numeric(0),
                stochasticType=mIndx,
                q=q,
                DT = DT,
                Ntail = nTail
                )
        rpArgs <- c(list(rp0),rpArgs)
       pb <- .SAMpb(min = 0, max = nosim_ci + 1, label="Point estimate")
        incpb <- function(label="") .SAM_setPB(pb, pb$getVal()+1,label)
        ssdr <- .refpointOptimizer(fit, rpArgs, nosim_ci,ncores, incpb)

        res <- .refpointOutput(ssdr,rpArgs, fit, biasCorrect, aveYearsIn, selYearsIn,

                               c(), #Fsequence,
                               referencepoints,
                               TRUE)
         attr(res,"equilibriumMethod") <- equilibriumMethod
        return(res)
        
    }else if(equilibriumMethod == "AD"){
##########################################################################################
################################## Approximation #########################################
##########################################################################################
        cat("Starting...\n")
        typePatterns <- list(Median = "^median$",
                             Mean = "^mean$",
                             Mode = "^mode$",
                             Quantile = "^q0\\.[[:digit:]]+$")
        mIndx <- which(sapply(typePatterns, function(p) grepl(p, tolower(method))))
        if(length(mIndx) == 0)
            stop(sprintf("Error in method specification. %s not recognized.",x))
        mIndx <- mIndx[1]
        if(mIndx == 4){            
            q <- as.numeric(gsub("(q)(0\\.[[:digit:]]+)","\\2",tolower(method)))
        }else{
            q <- NA_real_
        }
        catchType <- pmatch(catchType,c("catch","landing","discard"))
        if(is.na(catchType))
            stop("Invalid catch type")

        aveYearsIn <- aveYears
        if(length(aveYears) == 0)
            aveYears <- max(fit$data$years)+(-9:0)
        aveYears <- match(aveYears, fit$data$years) - 1
        if(any(is.na(aveYears)))
            stop("aveYears has years without data.")

        selYearsIn <- selYears
        selYears <- match(selYears, fit$data$years) - 1
        if(any(is.na(selYears)))
            stop("selYears has years without data.")

        cat("Parsing ref points...\n")
        rpArgs <- Reduce(.refpointMerger,
                         lapply(referencepoints, .refpointParser, nYears = nYears, aveYears = aveYears, selYears = selYears, logCustomSel = numeric(0), catchType = catchType - 1,logN0=fit$pl$logN[,ncol(fit$pl$logN)],stochasticType=mIndx,q=q, DT=DT),
                         list())
        cat("Adding starting values...\n")
        ## Add starting values    
        rpArgs <- lapply(rpArgs, .refpointStartingValue, fit = fit, Fsequence = Fsequence,stochasticType=0, q=q)
        ## Add Fsequence for plotting
        ## cat("Adding F sequence...\n")
        ## rp0 <- list(rpType = -1,
        ##             xVal = log(Fsequence[1]),
        ##             nYears = nYears,
        ##             aveYears = aveYears,
        ##             selYears = selYears,
        ##             logCustomSel = numeric(0),
        ##             catchType = catchType - 1,
        ##             logF0 = log(Fsequence[1]),
        ##             logN0=fit$pl$logN[,ncol(fit$pl$logN)],
        ##             stochasticType=mIndx,
        ##             q=q)
        cat("Prepare for TMB object...\n")
        ## Make list for TMB
        obj0 <- fit$obj
        argsIn <- as.list(obj0$env)[setdiff(methods::formalArgs(TMB::MakeADFun),"...")]
        argsIn$silent <- obj0$env$silent
        argsIn$parameters <- fit$pl
        argsIn$random <- unique(names(obj0$env$par[obj0$env$random]))
        argsIn$data$reportingLevel <- -1

        argsIn$data$referencepoints <- rpArgs ##c(list(rp0), rpArgs)
        attr(argsIn$data$referencepoints,"newton_config") <- newton.control
        args <- argsIn

        if(!run) return(args)
        cat("Make TMB Object...\n")
        objSDR <- do.call(TMB::MakeADFun, args)
        cat("Run $fn...\n")        
        objSDR$fn(fit$opt$par)
        cat("sdreport...\n")
        biasCorrect <- FALSE
        sdr <- TMB::sdreport(objSDR, objSDR$par, fit$opt$he,
                             bias.correct= biasCorrect,
                             skip.delta.method = biasCorrect,
                             bias.correct.control = list(sd = TRUE,
                                                         split = objSDR$env$ADreportIndex()[grepl("referencepoint_[[:digit:]]+_.+",names(objSDR$env$ADreportIndex()))]
                                                         ))
        ssdr <- summary(sdr)
        cat("Make tables...\n")        

        ## Make tables
        res <- .refpointOutput(ssdr,rpArgs, fit, biasCorrect, aveYearsIn, selYearsIn, numeric(0), referencepoints,TRUE)
        attr(res,"equilibriumMethod") <- equilibriumMethod
   
        cat("Return...\n")                
        return(res)
        
    }else if(equilibriumMethod == "ES"){
##########################################################################################
################################## Simulation ############################################
##########################################################################################
        
        ## Add some kind of progressbar / messages (and argument to silence)

        oldSeed <- NULL
        if(exists(".Random.seed"))
            oldSeed <- .Random.seed
        on.exit(set.seed(oldSeed))
        set.seed(seed)

        MT <- .refpointSMethodParser(method, formula = NA, derivedSummarizer=derivedSummarizer, positive = TRUE, penalize = penalizeSpline, useLogF = logFspline, RTMB = useRTMB, gamma = splineGamma, rho = splineRho)

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

        if(!all(Frange >= 0) && ((Frange[1] < Frange[2] && length(Frange) ==2) || !randomF))
            stop("Wrong Frange")
        if(!nosim > 0)
            stop("nosim must be a positive integer")

        ## Get RPs for best fit
        rpArgs <- Reduce(.refpointMerger,
                         lapply(referencepoints, .refpointParser, cutoff = 0),
                         list())
        invisible(lapply(rpArgs,.refpointCheckRecruitment,fit=fit))

        if (!.checkFullDerived(fit) && any(sapply(rpArgs, function(x) x$rpType %in% c(3,4,5))))
            stop("The reference points specified needs a fit with all derived values. Fit with `fullDerived=TRUE` or update with `getAllDerivedValues`.")

        doSim <- .getDoSim(logf1= max(Frange)+0.2,##tail(log(fbartable(fit)[,1]),1),
                           fit=fit, nYears = nYears, aveYears = aveYears, selYears = selYears, pl = fit$pl, constraint=constraint,deterministicF=deterministicF, processNoiseF=processNoiseF,...)
        
        pb <- .SAMpb(min = 0, max = nosim * (nosim_ci + 1 + is.function(derivedSummarizer)*length(rpArgs)))
        incpb <- function(label="") .SAM_setPB(pb, pb$getVal()+1,label)

        v0 <- .refpointSFitCriteria(rpArgs, pl=fit$pl, MT=MT, fit=fit, nosim=nosim, Frange=Frange, aveYears=aveYears, selYears=selYears, nYears=nYears, catchType=catchType, nTail=nTail,incpb=incpb,doSim=doSim,label="Estimation:",constraint=constraint,deterministicF=deterministicF, processNoiseF=processNoiseF, knots=knots, tailSummarizer = tailSummarizer, gridSummarizer = gridSummarizer, ncores=ncores,mc.type=mc.type, defaultSpline=defaultSpline, ...)

        ## Sample to get CIs
        if(nosim_ci > 0){
            plRep <- .asympSampleParVec(nosim_ci,fit, boundary = TRUE, returnList = TRUE)
            vv <- lapply(plRep, function(par){
                oN <- fit$obj
                a <- capture.output(invisible(oN$fn(par)))
                pl <- oN$env$parList(par,oN$env$last.par)
                v <- try({.refpointSFitCriteria(rpArgs,pl=pl, MT=MT, fit=fit, nosim=nosim, Frange=Frange, aveYears=aveYears, selYears=selYears, nYears=nYears, catchType=catchType, nTail=nTail,incpb=incpb,doSim=doSim,label="Confidence intervals:",constraint=constraint,deterministicF=deterministicF,knots=knots, tailSummarizer = tailSummarizer, ncores=ncores,mc.type=mc.type, ...)}, silent = TRUE)
                
                v
            })
            ## Get Ye/Re/Se/... (how should they be summarized?)
            ii <- sapply(vv, class) == "try-error"
        }else{
            ii <- logical(0)
        }
        if(sum(!ii) > 0){
            nan2na <- function(x)ifelse(is.nan(unlist(x)),NA,unlist(x))
            resTabs <- lapply(rownames(vv[!ii][[1]]$Estimates), function(nm){
                ci_deltaMethod <- FALSE
                if(ci_deltaMethod){
                    ## DOI:10.1109/WSC.2006.323107
                    Sigma <- fit$sdrep$cov.fixed
                    doOneCI <- function(jj){
                        Theta <- do.call("rbind",lapply(vv[!ii], function(x) nan2na(x$Estimates[nm,jj])))
                        dPar <- as.data.frame(do.call("rbind",plRep[!ii]))
                        dParM <- sapply(dPar,mean)
                        colnames(dPar) <- paste0(seq_len(ncol(dPar)),"_",colnames(dPar))
                        d0 <- cbind(data.frame(Theta = log(Theta)), dPar)
                        if(nrow(dPar) > ncol(dPar)){
                            gr0 <- tail(coef(stats::lm(Theta~.,data=d0)),-1)
                        }else{
                            gr0 <- sapply(seq_len(ncol(dPar)), function(kk){ sum((dPar[,kk]-dParM[kk])*(Theta-mean(Theta))) / sum((dPar[,kk]-dParM[kk])^2) })
                        }
                        sdErr <- as.numeric(t(gr0) %*% Sigma %*% (gr0))
                        exp(log(v0$Estimates[nm,jj])+ 2 * c(-1,1) * sdErr)
                    }
                    CIest <- t(sapply(seq_len(ncol(v0$Estimates)), doOneCI))
                    rownames(CIest) <- colnames(v0$Estimates)
                }else{
                    CIest <- t(apply(do.call("rbind",lapply(vv[!ii], function(x) nan2na(x$Estimates[nm,]))),2,quantile, prob = c(0.025,0.975), na.rm=TRUE))
                }
                tab <- cbind(v0$Estimate[nm,],CIest)
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
                    regression = v0$regressions,
                    curve_opt = v0$Curves ,
                    ## curve_list = v0$CurveFitsList,
                    ## curve_aics = v0$CurveFitsAIC,
                    ## opt = NA,
                    ## ssdr = sdr,
                    fbarlabel = substitute(bar(F)[X - Y], list(X = fit$conf$fbarRange[1], Y = fit$conf$fbarRange[2])),
                    stochastic=TRUE
                    ## diagonalCorrection = tv
                    )

             attr(res,"stochasticType") <- sapply(rpArgs,function(x) x$stochasticType)
        attr(res,"stochasticQ") <- sapply(rpArgs,function(x) x$q)
        attr(res,"nYears") <- sapply(rpArgs,function(x) x$nYears)
        attr(res,"equilibriumMethod") <- equilibriumMethod
   
        attr(res,"aveYears") <-  aveYearsIn
        attr(res,"selYears") <- selYearsIn
        
        attr(res,"fit") <- fit
        class(res) <- "sam_referencepoints"

        

        return(res)
    }
}
    
