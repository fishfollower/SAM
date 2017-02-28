##' Table helper 
##' @param fit returned object from sam.fit
##' @param quoted name of what to extract
##' @param x rownames of table
##' @param trans function to be applied
##' @details ...
.tableit <-function (fit, what, x=fit$data$years, trans=function(x)x){
   idx<-names(fit$sdrep$value)==what
   y<-fit$sdrep$value[idx]
   ci<-y+fit$sdrep$sd[idx]%o%c(-2,2)
   ret<-trans(cbind(y,ci))
   rownames(ret)<-x
   colnames(ret)<-c("Estimate","Low","Hig")
   return(ret)
}

##' SSB table 
##' @param  fit ...
##' @details ...
##' @export
ssbtable<-function(fit){
   ret<-.tableit(fit, "logssb", trans=exp) 
   return(ret)
}

##' TSB table 
##' @param  fit ... 
##' @details ...
##' @export
tsbtable<-function(fit){
   ret<-.tableit(fit, "logtsb", trans=exp) 
   return(ret)
}

##' Fbar table 
##' @param  fit ... 
##' @details ...
##' @export
fbartable<-function(fit){
   ret<-.tableit(fit, "logfbar", trans=exp) 
   return(ret)
}

##' Recruit table 
##' @param  fit ... 
##' @details ...
##' @export
rectable<-function(fit){
   ret<-.tableit(fit, "logR", trans=exp) 
   return(ret)
}

##' Catch table 
##' @param  fit ... 
##' @details ...
##' @export
catchtable<-function(fit){
   xx <- as.integer(rownames(fit$data$catchMeanWeight))
   ret <- .tableit(fit, x=xx, "logCatch", trans=exp) 
   return(ret)
}

##' Yield per recruit calculation
##' @param fit the object returned from sam.fit
##' @param Flimit Upper limit for Fbar
##' @param Fdelta increments on the Fbar axis 
##' @param aveYears Number of years back to use when calculating averages (selection, weights, ...)
##' @param ageLimit Oldest age used (should be high)
##' @export
ypr<-function(fit, Flimit=2, Fdelta=0.01, aveYears=15, ageLimit=100){
  idxF <- fit$conf$keyLogFsta[1,]+1
  barAges <- do.call(":",as.list(fit$conf$fbarRange))+(1-fit$conf$minAge) 
  last.year.used=max(fit$data$years)
  idxno<-which(fit$data$years==last.year.used)
  #dim<-fit.current$stateDim
  #idxN<-1:ncol(stock.mean.weight) 
  F <- exp(fit$pl$logF[idxF,])
  
  sel<-function(){
    Sa<-rep(0,nrow(F))
    K<-0
    for(i in 0:(aveYears-1)){
      thisF<-F[,idxno-i]
      Sa<-Sa+thisF
      K<-K+fbartable(fit)[idxno-i]
    }
    return(Sa/K)
  }

  extend<-function(x,len=100){
    ret<-numeric(len)
    ret[1:length(x)]<-x
    ret[-c(1:length(x))]<-x[length(x)]
    ret
  }

  ave.sl<-sel()
  ave.sw<-colMeans(fit$data$stockMeanWeight[(idxno-aveYears+1):idxno,,drop=FALSE])
  ave.cw<-colMeans(fit$data$catchMeanWeight[(idxno-aveYears+1):(idxno-1),,drop=FALSE])
  ave.pm<-colMeans(fit$data$propMat[(idxno-aveYears+1):idxno,,drop=FALSE])
  ave.nm<-colMeans(fit$data$natMor[(idxno-aveYears+1):idxno,,drop=FALSE])
  ave.lf<-colMeans(fit$data$landFrac[(idxno-aveYears+1):(idxno-1),,drop=FALSE])
  ave.cw.land<-colMeans(fit$data$landMeanWeight[(idxno-aveYears+1):(idxno-1),,drop=FALSE])

  N<-numeric(ageLimit)
  N[1]<-1.0
  M<-extend(ave.nm)
  sw<-extend(ave.sw)
  cw<-extend(ave.cw.land)
  pm<-extend(ave.pm)
  lf<-extend(ave.lf)
 
  deltafirst <- 0.00001
  delta <- Fdelta
  scales<-c(0, deltafirst, seq(0.01, Flimit, by=delta))
  yields<-numeric(length(scales))
  ssbs<-numeric(length(scales))
  for(i in 1:length(scales)){
    scale<-scales[i]
    F<-extend(ave.sl*scale)
    Z<-M+F
    for(a in 2:length(N)){
      N[a]<-N[a-1]*exp(-Z[a-1])  
    }
    C<-F/Z*(1-exp(-Z))*N*lf  
    Y<-sum(C*cw)
    yields[i]<-Y
    ssbs[i]<-sum(N*pm*sw)
  }

  fmaxidx<-which.max(yields)
  fmax<-scales[fmaxidx]

  deltaY<-diff(yields)
  f01idx<-which.min((deltaY/delta-0.1*deltaY[1]/deltafirst)^2)+1
  f01<-scales[f01idx]

  f35spridx<-which.min((ssbs-0.35*ssbs[1])^2)+1
  f35<-scales[f35spridx]
  
  fbarlab <- substitute(bar(F)[X - Y], list(X = fit$conf$fbarRange[1], Y = fit$conf$fbarRange[2]))
  ret<-list(fbar=scales, ssb=ssbs, yield=yields, fbarlab=fbarlab, f35=f35, f01=f01, fmax=fmax, 
            f35Idx=f35spridx, f01Idx=f01idx, fmaxIdx=fmaxidx)
  class(ret)<-"samypr"
  return(ret)
}
