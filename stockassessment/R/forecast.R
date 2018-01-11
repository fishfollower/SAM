##' rmvnorm helper function to draw multivariate normal samples
##' @param n the number of samples.
##' @param mu the mean vector.
##' @param Sigma a positive-definite symmetric matrix specifying the covariance matrix.
##' @details Generates samples via the Cholesky decomposition, which is less platform dependent than eigenvalue decomposition.
##' @return If n = 1 a vector of the same length as mu, otherwise an n by length(mu) matrix with one sample in each row.
##' @export
rmvnorm <- function(n = 1, mu, Sigma){
  p <- length(mu)
  if(!all(dim(Sigma) == c(p, p))){
    stop("incompatible arguments")
  }
  idx <- diag(Sigma) > .Machine$double.xmin 
  L <- matrix(0,p,p)
  if(any(idx)){
    L[idx,idx] <- chol(Sigma[idx,idx])
  }
  X <- matrix(rnorm(p * n), n)
  X <- drop(mu) + t(X%*%L)
  if(n == 1){
    drop(X)
  }else{
    t(X)
  }
}

##' forecast function to do shortterm
##' @param fit an assessment object of type sam, as returned from the function sam.fit
##' @param fscale a vector of f-scales. See details.  
##' @param catchval a vector of target catches. See details.   
##' @param fval a vector of target f values. See details.  
##' @param nosim number of simulations default is 1000
##' @param year.base starting year default last year in assessment. Currently it is only supported to use last assessment year or the year before  
##' @param ave.years vector of years to average for weights, maturity, M and such  
##' @param rec.years vector of years to use to resample recruitment from
##' @param label optional label to appear in short table
##' @param overwriteSelYears if a vector of years is specified, then the average selectivity of those years is used (not recommended)
##' @param deterministic option to turn all process noise off (not recommended, as it will likely cause bias)
##' @details There are three ways to specify a scenario. If e.g. four F values are specified (e.g. fval=c(.1,.2,.3,4)), then the first value is used in the last assessment year (base.year), and the three following in the three following years. Alternatively F's can be specified by a scale, or a target catch. Only one option can be used per year. So for instance to set a catch in the first year and an F-scale in the following one would write catchval=c(10000,NA,NA,NA), fscale=c(NA,1,1,1). The length of the vector specifies how many years forward the scenarios run. 
##' @return an object of type samforecast
##' @importFrom stats median uniroot quantile
##' @importFrom Matrix bdiag
##' @export
forecast <- function(fit, fscale=NULL, catchval=NULL, fval=NULL, nosim=1000, year.base=max(fit$data$years), ave.years=max(fit$data$years)+(-4:0), rec.years=max(fit$data$years)+(-9:0), label=NULL, overwriteSelYears=NULL, deterministic=FALSE){

  idxN <- 1:nrow(fit$rep$nvar)
    
  idxF <- 1:nrow(fit$rep$fvar)+nrow(fit$rep$nvar)
    
  resample <- function(x, ...){
    if(deterministic){
      ret <- mean(x)
    }else{
      ret <- x[sample.int(length(x), ...)]
    }
    return(ret)
  }
    
  if(missing(fscale)&missing(fval)&missing(catchval))stop("No scenario is specified")    
  if(missing(fscale)&!missing(fval))fscale<-rep(NA,length(fval))
  if(missing(fscale)&!missing(catchval))fscale<-rep(NA,length(catchval))
  if(missing(fval)&!missing(fscale))fval<-rep(NA,length(fscale))
  if(missing(fval)&!missing(catchval))fval<-rep(NA,length(catchval))
  if(missing(catchval)&!missing(fval))catchval<-rep(NA,length(fval))
  if(missing(catchval)&!missing(fscale))catchval<-rep(NA,length(fscale))
        
  if(!all(rowSums(!is.na(cbind(fscale, catchval, fval)))==1)){
    stop("For each forecast year exactly one of fscale, catchval or fval must be specified (all others must be set to NA)")
  }

  if(!is.null(overwriteSelYears)){
    fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)  
    Ftab <- faytable(fit)
    fixedsel <- colMeans(Ftab[as.integer(rownames(Ftab))%in%overwriteSelYears,,drop=FALSE])
    fixedsel <- fixedsel/mean(fixedsel[fromto[1]:fromto[2]])
  }
    
  getF <- function(x, fleet=which(fit$data$fleetTypes==0)){
    getfleet <- function(f){
      idx <- fit$conf$keyLogFsta[f,]+2    
      ret <- c(NA,exp(x[idxF]))[idx]
      ret[is.na(ret)] <- 0
      ret
    }
    FF <- lapply(fleet,getfleet)
    ret <- Reduce("+",FF) 
    if(!is.null(overwriteSelYears)){
      fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)    
      thisfbar<-mean(ret[fromto[1]:fromto[2]])
      ret<-fixedsel*thisfbar
      FF[]<-NA
    }
    attr(ret,"byFleet") <- do.call(cbind,FF)
    ret
  }

  fbar<-function(x){
    fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)  
    mean(getF(x)[fromto[1]:fromto[2]])
  }    
    
  getN <- function(x){
    idx <- fit$conf$keyLogFsta[1,]+1
    nsize <- length(idx)
    ret <- exp(x[1:nsize])
    ret
  }

  getProcessVar <- function(fit, set.rec.var.zero=TRUE){
    cov <- bdiag(fit$rep[c("nvar","fvar")])
    if(set.rec.var.zero){
      cov[1,1] <- 0
    }
    as.matrix(cov)
  }
    
  step <- function(x, nm, recpool, scale, inyear=FALSE){
    F <- getF(x)
    N <- getN(x)
    if(!inyear){
      Z <- F+nm
      n <- length(N)
      N <- c(resample(recpool,1),N[-n]*exp(-Z[-n])+c(rep(0,n-2),N[n]*exp(-Z[n])))
    }
    xx <- rep(NA,length=length(x))
    xx[idxN] <- log(N)
    xx[idxF] <- x[idxF]+log(scale)
    return(xx)
  }
    
  scaleF <- function(x, scale){
    x[idxF] <- x[idxF]+log(scale)
    return(x)
  }
    
  sel<-function(x){
    getF(x)/fbar(x)
  }

  catch <- function(x, nm, cw){
    F <- getF(x)
    FF <- attr(F,"byFleet")
    Z <- F+nm
    N <- getN(x)
    C <- FF/Z*(1-exp(-Z))*N
    TCW <- sum(C*cw)
    attr(TCW, "byFleet") <- colSums(C*cw)
    return(TCW)
  }

  ssb <- function(x, nm, sw, mo, pm, pf){
    F <- getF(x)
    FF <- attr(F,"byFleet")
    ZZ <- pm*nm+rowSums(pf*FF)
    N <- getN(x)*exp(-ZZ)
    return(sum(N*mo*sw))
  }

  rectab<-rectable(fit)
  recpool<-rectab[rownames(rectab)%in%rec.years,1]

  # Get final state
  if(year.base==max(fit$data$years)){
    est <- fit$sdrep$estY
    cov <- fit$sdrep$covY
  }
  if(year.base==(max(fit$data$years)-1)){  
    est <- fit$sdrep$estYm1
    cov <- fit$sdrep$covYm1
  }
  if(year.base<(max(fit$data$years)-1)){
    stop("State not saved, so cannot proceed from this year")
  }
  if(deterministic){
    cov<-cov*0
  }
  sim <- rmvnorm(nosim, mu=est, Sigma=cov)
    
  doAve <- function(x){
    if(length(dim(x))==2){
      ret <- colMeans(x[rownames(x)%in%ave.years,,drop=FALSE])
    }
    if(length(dim(x))==3){
      ret <- apply(x[rownames(x)%in%ave.years,,,drop=FALSE],c(2,3),mean)
    }
    ret
  }

  ave.sw <- doAve(fit$data$stockMeanWeight)
  ave.cw <- doAve(fit$data$catchMeanWeight)
  ave.mo <- doAve(fit$data$propMat)
  ave.nm <- doAve(fit$data$natMor)
  ave.lf <- doAve(fit$data$landFrac)
  ave.lw <- doAve(fit$data$landMeanWeight)
  ave.pm <- doAve(fit$data$propM)
  ave.pf <- doAve(fit$data$propF)
  getThisOrAve <- function(x,y, ave){
    if(y %in% rownames(x)){
      if(length(dim(x))==2){  
        ret <- x[which(rownames(x)==y),]
      }else{
        ret <- x[which(rownames(x)==y),,]
      }
    }else{
      ret <- ave
    }
    ret
  }
  procVar<-getProcessVar(fit)  
  simlist<-list()
  for(i in 0:(length(fscale)-1)){
    y<-year.base+i
    
    sw<-getThisOrAve(fit$data$stockMeanWeight, y, ave.sw)
    cw<-getThisOrAve(fit$data$catchMeanWeight, y, ave.cw)
    mo<-getThisOrAve(fit$data$propMat, y, ave.mo)
    nm<-getThisOrAve(fit$data$natMor, y, ave.nm)
    lf<-getThisOrAve(fit$data$landFrac, y, ave.lf)    
    lw<-getThisOrAve(fit$data$landMeanWeight, y, ave.lw)
    pm<-getThisOrAve(fit$data$propM, y, ave.pm)
    pf<-getThisOrAve(fit$data$propF, y, ave.pf)

    sim <- t(apply(sim, 1, function(s)step(s, nm=nm, recpool=recpool, scale=1, inyear=(i==0))))
    if(i!=0){
      if(deterministic)procVar<-procVar*0  
      sim <- sim + rmvnorm(nosim, mu=rep(0,nrow(procVar)), Sigma=procVar)
    }
    
    if(!is.na(fscale[i+1])){
      sim<-t(apply(sim, 1, scaleF, scale=fscale[i+1]))    
    }

    if(!is.na(fval[i+1])){
      curfbar<-median(apply(sim, 1, fbar))
      adj<-fval[i+1]/curfbar
      sim<-t(apply(sim, 1, scaleF, scale=adj))    
    }

    if(!is.na(catchval[i+1])){
      simtmp<-NA
      fun<-function(s){
        simtmp<<-t(apply(sim, 1, scaleF, scale=s))      
        simcat<-apply(simtmp, 1, catch, nm=nm, cw=cw)
        return(catchval[i+1]-median(simcat))
      }
      ff <- uniroot(fun, c(0,100))$root
      sim <- simtmp
    }
    
    fbarsim <- apply(sim, 1, fbar)
    catchsim <- apply(sim, 1, catch, nm=nm, cw=cw)
    ssbsim <- apply(sim, 1, ssb, nm=nm, sw=sw, mo=mo, pm=pm, pf=pf)
    recsim <- exp(sim[,1])
    catchbysim <- apply(sim, 1, function(x)attr(catch(x, nm=nm, cw=cw), "byFleet"))
    simlist[[i+1]] <- list(sim=sim, fbar=fbarsim, catch=catchsim, ssb=ssbsim, rec=recsim, year=y, catchby=catchbysim)
  }
  attr(simlist, "fit")<-fit

  collect <- function(x){
    quan <- quantile(x, c(.50,.025,.975))
    c(median=quan[1], low=quan[2], high=quan[3])
  }
  fbar <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$fbar))),3)
  rec <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$rec))))
  ssb <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$ssb))))
  catch <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$catch))))
  if(sum(fit$data$fleetTypes==0)==1){
    catchby <- catch
  }else{
    catchby <- round(do.call(rbind, lapply(simlist,function(xx)as.vector(apply(xx$catchby,1,collect)))))
  }
  tab <- cbind(fbar, rec, ssb, catch)
  rownames(tab) <- unlist(lapply(simlist, function(xx)xx$year))
  rownames(catchby) <- rownames(tab)
  nam <- c("median","low","high")
  colnames(tab) <- paste0(rep(c("fbar:","rec:","ssb:","catch:"), each=length(nam)), nam)
  colnames(catchby) <- paste0(rep(paste0("F",1:sum(fit$data$fleetTypes==0),":"), each=length(nam)), nam)
  attr(simlist, "tab") <- tab
  attr(simlist, "catchby") <- catchby
  shorttab <- t(tab[,grep("median",colnames(tab))])
  rownames(shorttab) <- sub(":median","",paste0(label,if(!is.null(label))":",rownames(shorttab)))
  attr(simlist, "shorttab") <- shorttab
  attr(simlist, "label") <- label  
  class(simlist) <- "samforecast"
  simlist
}
