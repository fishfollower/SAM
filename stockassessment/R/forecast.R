##' forecast function to do shortterm
##' @param fit an assessment object of type sam, as returned from the function sam.fit
##' @param fscale a vector of f-scales. See details.  
##' @param catchval a vector of target catches. See details.   
##' @param fval a vector of target f values. See details.  
##' @param nosim number of simulations default is 1000
##' @param year.base starting year default last year in assessment. Currently it is only supported to use last assessment year or the year before  
##' @param ave.years vector of years to average for weights, maturity, M and such  
##' @param rec.years vector of years to use to resample recruitment from 
##' @details There are three ways to specify a scenario. If e.g. four F values are specified (e.g. fval=c(.1,.2,.3,4)), then the first value is used in the last assessment year (base.year), and the three following in the three following years. Alternatively F's can be specified by a scale, or a target catch. Only one option can be used per year. So for instance to set a catch in the first year and an F-scale in the following one would write catchval=c(10000,NA,NA,NA), fscale=c(NA,1,1,1). The length of the vector specifies how many years forward the scenarios run. 
##' @return an object of type samforecast
##' @importFrom stats median uniroot quantile
##' @importFrom MASS mvrnorm
##' @export
forecast <- function(fit, fscale=NULL, catchval=NULL, fval=NULL, nosim=1000, year.base=max(fit$data$years), ave.years=max(fit$data$years)+(-4:0), rec.years=max(fit$data$years)+(-9:0)){

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
    
  getF <- function(x){
    idx <- fit$conf$keyLogFsta[1,]+1
    nsize <- length(idx)
    ret <- exp(x[nsize+idx])
    ret[idx==0] <- 0
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

  getState <- function(N,F){
    k <- fit$conf$keyLogFsta[1,]
    F <- F[k>=0]
    k <- unique(k[k>=0]+1)
    x <- log(c(N,F[k]))
    x
  }    

  getProcessVar <- function(fit){
    cof <- coef(fit)
    sdN <- exp(cof[names(cof)=="logSdLogN"][fit$conf$keyVarLogN+1])
    sdN[1]<-0
    nN <- length(sdN)
    sdF <- exp(cof[names(cof)=="logSdLogFsta"][fit$conf$keyVarF+1])
    k<-fit$conf$keyLogFsta[1,]
    sdF <- sdF[k>=0]
    k <- unique(k[k >= 0] + 1)
    sdF <-sdF[k]
    nF <- length(sdF)
    if(fit$conf$corFlag==0){
      corr <- diag(nF)    
    }
    if(fit$conf$corFlag==1){
      y <- cof[names(cof)=="itrans_rho"]
      rho <- 2/(1+exp(-2*y))-1
      corr <- matrix(rho, nrow=nF, ncol=nF)
      diag(corr) <- 1
    }
    if(fit$conf$corFlag==2){
      y <- cof[names(cof)=="itrans_rho"]
      rho <- 2/(1+exp(-2*y))-1
      corr <- diag(sdF)
      corr <- rho^abs(row(corr)-col(corr))
    }
    cov <- matrix(0,nrow=nN+nF,ncol=nN+nF)
    cov[1:nN,1:nN] <- diag(sdN^2)
    cov[nN+1:nF,nN+1:nF] <- (sdF%*%t(sdF))*corr
    cov
  }
    
  step <- function(x, nm, recpool, scale, inyear=FALSE){
    F <- getF(x)
    N <- getN(x)
    if(!inyear){
      Z <- F+nm
      n <- length(N)
      N <- c(sample(recpool,1),N[-n]*exp(-Z[-n])+c(rep(0,n-2),N[n]*exp(-Z[n])))
    }
    F <- F*scale
    xx <- getState(N,F)
    return(xx)
  }
    
  scaleF <- function(x, scale){
    F <- getF(x)*scale
    N <- getN(x)
    xx <- getState(N,F)
    return(xx)
  }
    
  sel<-function(x){
    getF(x)/fbar(x)
  }

  catch <- function(x, nm, cw){
    F <- getF(x)
    Z <- F+nm
    N <- getN(x)
    C <- F/Z*(1-exp(-Z))*N
    return(sum(cw*C))
  }

  ssb <- function(x, nm, sw, mo, pm, pf){
    F <- getF(x)
    N <- getN(x)*exp(-pm*nm-pf*F)
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
  sim<-MASS::mvrnorm(nosim, mu=est, Sigma=cov)
    
  if(!all.equal(est,getState(getN(est),getF(est))))stop("Sorry somthing is wrong here (check code for getN, getF, and getState)")  

  doAve<-function(x,y)colMeans(x[rownames(x)%in%ave.years,,drop=FALSE]) 
  ave.sw<-doAve(fit$data$stockMeanWeight)
  ave.cw<-doAve(fit$data$catchMeanWeight)
  ave.mo<-doAve(fit$data$propMat)
  ave.nm<-doAve(fit$data$natMor)
  ave.lf<-doAve(fit$data$landFrac)
  ave.lw<-doAve(fit$data$landMeanWeight)
  ave.pm<-doAve(fit$data$propM)
  ave.pf<-doAve(fit$data$propF)
  getThisOrAve<-function(x,y, ave){
    if(y %in% rownames(x)){
      ret<-x[which(rownames(x)==y),]
    }else{
      ret<-ave
    }
    return(ret)
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
      sim <- sim + MASS::mvrnorm(nosim, mu=rep(0,nrow(procVar)), Sigma=procVar)
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
    simlist[[i+1]] <- list(sim=sim, fbar=fbarsim, catch=catchsim, ssb=ssbsim, rec=recsim, year=y)
  }
    
  attr(simlist, "fit")<-fit

  collect <- function(x){
    quan <- quantile(x, c(.50,.025,.975))
    c(median=quan[1], low=quan[2], hig=quan[3])
  }
  fbar <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$fbar))),3)
  rec <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$rec))))
  ssb <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$ssb))))
  catch <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$catch))))
  tab <- cbind(fbar, rec,ssb,catch)
  rownames(tab) <- unlist(lapply(simlist, function(xx)xx$year))
  nam <- c("median","low","hig")
  colnames(tab) <- paste0(rep(c("fbar:","rec:","ssb:","catch:"), each=length(nam)), nam)
  attr(simlist, "tab")<-tab
    
  class(simlist) <- "samforecast"
  simlist
}
