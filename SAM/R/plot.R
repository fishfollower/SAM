##' Plot helper 
##' @param  fit ... 
##' @details ...
.plotit <-function (what, x=fit$data$years, ylab=what, xlab="Years", ex=numeric(0), trans=function(x)x ,...){
   idx<-names(fit$sdrep$value)==what
   y<-fit$sdrep$value[idx]
   ci<-y+fit$sdrep$sd[idx]%o%c(-2,2)
   plot(x,trans(y), xlab=xlab, ylab=ylab, type="l", lwd=3, ylim=range(c(trans(ci),0,ex)), las=1,...)
   polygon(c(x,rev(x)), y = c(trans(ci[,1]),rev(trans(ci[,2]))), border = gray(.5,alpha=.5), col = gray(.5,alpha=.5))
   grid(col="black")
 }

##' SAM Fbar plot 
##' @param  fit ...
##' @param  partial ...
##' @details ...
##' @export
fbarplot<-function(fit,partial=TRUE,...){
  fbarRange<-fit$conf$fbarRange
  fbarlab=substitute(bar(F)[X-Y],list(X=fbarRange[1],Y=fbarRange[2]))
  fmat<-fit$pl$logF[fit$conf$keyLogFsta[1,]+1,]
  idx<-which(fit$conf$minAge:fit$conf$maxAge %in% fbarRange[1]:fbarRange[2])
  exx <- if(partial){exp(fmat[idx,])}else{numeric(0)}
  .plotit("logfbar", ylab=fbarlab, trans=exp, ex=exx,...)
  if(partial){
    matplot(fit$data$years, t(exp(fmat[idx,])), add=TRUE, type="b", col="lightblue", pch=as.character(fbarRange[1]:fbarRange[2]))
  }  
}


##' SAM SSB plot 
##' @param  fit ...
##' @param  partial ...
##' @details ...
##' @export
ssbplot<-function(fit,...){
  .plotit("logssb", ylab="SSB", trans=exp,...)
}

##' SAM Recruits plot 
##' @param  fit ...
##' @param  partial ...
##' @details ...
##' @export
recplot<-function(fit,...){
  .plotit("logR", ylab="Recruits", trans=exp,...)
}
