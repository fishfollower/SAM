spsa <- function(p, fn, twosided=TRUE, maxit = 20, tol = 1e-7, ...){
    tol <- rep(tol, length.out = 6)
    matsqrt <- function(X){
        d <- svd(X)
        d$v %*% diag(sqrt(d$d),length(d$d),length(d$d)) %*% t(d$v)
    }
    np <- length(p)
    k <- 1
    ak <- 1; ghat <- 1
    Hbar <- diag(np)
    XAc <- Inf
    FAc <- Inf
    XRc <- Inf
    FRc <- Inf
    Gc <- Inf
    ylast <- Inf
    a <- 1
    alpha <- 0.602
    gamma <- 0.101
    A <- 0.05 * 200
    ###vvv <- c(fn(p,...),fn(p,...),fn(p,...),fn(p,...))
    vvv <- unlist(lapply(replicate(5,p,FALSE),fn,...))
    cat(mean(vvv))
    c0 <- 2 * sd(as.vector(vvv))
    c1 <- 2 * c0
    yplus <- yminus <- mean(vvv)
    while(k <= maxit && XAc > tol[1] && XRc > tol[2] && FAc > tol[3] && FRc > tol[4]  && Gc > tol[5]){
        ak <- a / (k+A) ^ alpha
        ## vvv <- unlist(lapply(replicate(5,p,FALSE),fn,...))
        ## c0 <- 3 * sd(as.vector(vvv))
        ## c1 <- 2 * c0    
        ck <- c0 / k ^ gamma
        ck2 <- c1 / k ^ gamma
        Delta_k <- 2 * rbinom(np,1,0.5) - 1
        Delta_k2 <- 2 * rbinom(np,1,0.5) - 1
        pplus <- p + ck * Delta_k
        pminus <- p - ck * Delta_k
        pplusplus <- pplus + ck2 * Delta_k2
        pminusplus <- pminus + ck2 * Delta_k2
        lastF <- 0.5 * (yplus + yminus)
        yplus <- fn(pplus,...)
        yminus <- fn(pminus,...)
        yplusplus <- fn(pplusplus,...)
        yminusplus <- fn(pminusplus,...)
        ghat <- (yplus-yminus) / (2.0 * ck * Delta_k)
        Gplus <- (yplusplus - yplus) / (ck2 * Delta_k2)
        Gminus <- (yminusplus - yminus) / (ck2 * Delta_k2)
        dG <- Gplus - Gminus
        Hx <- dG %o% (1 / (2 * ck * Delta_k))
        Hhat <- 0.5 * (Hx + t(Hx))
        ## p <- p - ak * ghat
        Hbar <- k/(k+1) * Hbar + 1/(k+1) * Hhat
        Huse <- matsqrt( Hbar %*% Hbar) + diag(1/k,np,np)
        lastP <- p
        step <- as.vector(ak * solve(Huse) %*% ghat)        
        if(norm(step,"2") > 0.1 * norm(p,"2")){
            step <- step / norm(step,"2") * 0.1 * norm(p,"2")
        }
        ## ynext <- replicate(5, fn(p + step))
        ## if(mean(ynext) < ylast - sd(ynext)){
            p <- p - step
            XAc <- norm(p - lastP,"2")
            XRc <- XAc / norm(lastP,"2")
            FAc <- norm(0.5 * (yplus + yminus) - lastF,"2")
            FRc <- FAc / norm(lastF,"2")
            Gc <- max(abs(ghat))
        ## }
        cat(k,mean(c(yplus,yminus)),"\t",paste(p,collapse=" "),"\n")
        k <- k+1
    }
    list(par = p, objective = fn(p,...), iter=k, H = Huse,
         convergence = c(XAc > tol[1],
                         XRc > tol[2],
                         FAc > tol[3],
                         FRc > tol[4],
                         Gc > tol[3]))
}

stochastic_referencepoints <- function(fit, FUN, ...){
    UseMethod("stochastic_referencepoints")
}



stochastic_referencepoints.sam <- function(fit,
                                           FUN = median,
                                           nYears,
                                           Fsequence = seq(0,4,len = 200),
                                           aveYears = max(fit$data$years)+(-9:0),
                                           selYears = max(fit$data$years),
                                           SPRpercent = c(0.35),
                                           dYPRpercent = c(0.1),
                                           B0percent = c(0.2),
                                           catchType = "catch",
                                           MSYreduction = c(0.05),
                                           seed,
                                           ...){

    if(missing(seed)){
        set.seed(NULL)
        seed <- .Random.seed
    }
    
    getStochPR <- function(logf, n, nYears, sel, pl){
        set.seed(seed)
        v <- do.call("rbind",replicate(n,.Call("stochPerRecruitR",
                                          logFbar = logf,
                                          tmbdat = fit$obj$env$data,
                                          pl = pl,
                                          sel = sel,
                                          aveYears = aveYearsI,
                                          nYears = nYears,
                                          CT = 0,
                                          logNinit = as.vector(tail(t(pl$logN),1))
                                          ), simplify = FALSE)
                     )
        storage.mode(v) <- "double"
        v
    }
    getSummary <- function(v, what){
        v2 <- v[,what]
        c(Estimate = FUN(v2),
          Low = quantile(v2,0.025),
          High = quantile(v2,0.975))
    }

    ## skeleton <- list(logFmsy = 0,
    ##                  logFmsyRange = matrix(0,2,length(MSYreduction)),
    ##                  logFmypyl = 0,
    ##                  logFmdy = 0,
    ##                  logFmax = 0,
    ##                  logFxdYPR = numeric(length(dYPRpercent)),
    ##                  logFcrash = 0,
    ##                  logFext = 0,
    ##                  logFxPercent = numeric(length(SPRpercent)),
    ##                  logFxB0 = numeric(length(B0percent)),
    ##                  logFlim = 0)
    ## ## Assume everything is there
    ## ##map <- list()
    ## ## Functions to optimize
    ## yield <- function(f, n){
    ##     v <- getStochPR(log(f), n, nYears)[,"logYe"]
    ##     FUN(exp(v))
    ## }
    ## SPR <- function(f, n){
    ##     v <- getStochPR(log(f), n, nYears)[,"logSPR"]
    ##     FUN(exp(v))
    ## }
    ## Se <- function(f, n){
    ##     v <- getStochPR(log(f), n, nYears)[,"logSe"]
    ##     FUN(exp(v))
    ## }
    ## YPR <- function(f, n){
    ##     v <- getStochPR(log(f), n, nYears)[,"logYPR"]
    ##     FUN(exp(v))
    ## }
    ## MSYRange <- c(0.95)
    ## fn <- function(p0, n){
    ##     par <- relist(p0, skeleton)
    ##     ## Handle map
    ##     nll <- skeleton
    ##     if(!is.na(par$logFmsy)){
    ##         tmpYmsy <- yield(exp(par$logFmsy),n)
    ##         nll$logFmsy <- -log(tmpYmsy)
    ##         for(i in 1:ncol(par$logFmsyRange)){
    ##             if(!is.na(par$logFmsyRange[1,i])){
    ##                 logF <- par$logFmsy - exp(-par$logFmsyRange[1,i])
    ##                 tmp <- log(yield(exp(logF), n)) - (log(MSYRange[i]) + log(tmpYmsy))
    ##                 nll$logFmsyRange[1,i] <- tmp^2
    ##             }
    ##             if(!is.na(par$logFmsyRange[2,i])){
    ##                 logF <- par$logFmsy + exp(par$logFmsyRange[2,i])
    ##                 tmp <- log(yield(exp(logF), n)) - (log(MSYRange[i]) + log(tmpYmsy))
    ##                 nll$logFmsyRange[2,i] <- tmp^2
    ##             }
    ##         }
    ##     }
    ##     sum(unlist(nll), na.rm = TRUE)
    ## }
    
    ## opt <- spsa(log(opt0$par),fn, N=50, maxit = 20000))
    ## ## Fixed RP's
    ## F <- c(sq = unname(tail(fbartable(fit)[,1],1)),
    ##        zero = 0)

    aveYearsI <- as.numeric(match(aveYears, fit$data$years))
    selYearsI <- match(selYears, fit$data$years)
    toSeq <- function(x) seq(x[1], x[2], by = 1)
    logF2AvgSel <- function(logF){
        fbar <- colMeans(exp(logF)[fit$conf$keyLogFsta[1,toSeq(fit$conf$fbarRange)]+1,,drop=FALSE])
        rowMeans(exp(logF) / matrix(fbar, nrow(logF), ncol(logF), byrow=TRUE))
    }
    joint_fn <- function(par, rp_logf){
        ## ONLY MSY FOR NOW
        a <- capture.output(fit$obj$fn(par))
        pl <- fit$obj$env$parList(par = fit$obj$env$last.par)
        logF <- pl$logF[,selYearsI,drop=FALSE]
        sel <- logF2AvgSel(logF)
        v <- getStochPR(rp_logf,100, 0, sel, pl)
        -log(median(exp(v[,"logYe"])))
    }

    opt <- nlminb(log(0.5), function(x)joint_fn(fit$opt$par,x), control = list(trace=1))
    
    ## Fsequence
    allVal <- lapply(Fsequence, function(f) getStochPR(log(f), 100, 1000,logF2AvgSel(fit$pl$logF[,selYearsI,drop=FALSE]),fit$pl))
    vFun1000 <- do.call("rbind",lapply(allVal, function(x) apply(exp(x),2,FUN)))
    cfL <- lapply(colnames(vFun)[-1], function(nm) covafillr::covafill(vFun[,1], vFun[,nm], p = 7L))
    names(cfL) <- colnames(vFun)[-1]

    ## Fmsy
    opt_fmsy<- nlminb((Fsequence[which.max(vv)]),
               objective = function(x) -cf$predict(x)[1,1],
               gradient = function(x) -cf$predict(x)[1,2],
               hessian = function(x) matrix(-cf$predict(x)[1,3],1,1))

    ## Uncertainty

}
                             
