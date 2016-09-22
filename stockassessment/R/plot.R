##' Plot helper
##' @param  fit the fitted object from sam.fit
##' @param  what quoted name of object to extract
##' @param  x x-alues
##' @param  ymab label on y-axis
##' @param  xlab label on x-axis
##' @param  ex extra y's to make room for
##' @param  trans function to transform values by
##' @param  ... extra arguments transferred to plot
##' @importFrom graphics plot polygon grid
##' @importFrom grDevices gray
##' @details ...
.plotit <-function (fit, what, x=fit$data$years, ylab=what, xlab="Years", ex=numeric(0), trans=function(x)x ,...){
   idx<-names(fit$sdrep$value)==what
   y<-fit$sdrep$value[idx]
   ci<-y+fit$sdrep$sd[idx]%o%c(-2,2)
   plot(x,trans(y), xlab=xlab, ylab=ylab, type="l", lwd=3, ylim=range(c(trans(ci),0,ex)), las=1,...)
   polygon(c(x,rev(x)), y = c(trans(ci[,1]),rev(trans(ci[,2]))), border = gray(.5,alpha=.5), col = gray(.5,alpha=.5))
   grid(col="black")
 }

##' SAM Fbar plot 
##' @param fit the object returned from sam.fit 
##' @param partial true if included partial F's are to be plotted 
##' @param ... extra arguments transferred to plot
##' @importFrom graphics matplot
##' @details ...
##' @export
fbarplot<-function(fit,partial=TRUE,...){
  fbarRange<-fit$conf$fbarRange
  fbarlab=substitute(bar(F)[X-Y],list(X=fbarRange[1],Y=fbarRange[2]))
  fmat<-fit$pl$logF[fit$conf$keyLogFsta[1,]+1,]
  idx<-which(fit$conf$minAge:fit$conf$maxAge %in% fbarRange[1]:fbarRange[2])
  exx <- if(partial){exp(fmat[idx,])}else{numeric(0)}
  .plotit(fit, "logfbar", ylab=fbarlab, trans=exp, ex=exx,...)
  if(partial){
    matplot(fit$data$years, t(exp(fmat[idx,])), add=TRUE, type="b", col="lightblue", pch=as.character(fbarRange[1]:fbarRange[2]))
  }  
}


##' SAM SSB plot 
##' @param fit the object returned from sam.fit
##' @param ... extra arguments transferred to plot
##' @details ...
##' @export
ssbplot<-function(fit, ...){
  .plotit(fit, "logssb", ylab="SSB", trans=exp,...)
}

##' SAM TSB plot 
##' @param fit the object returned from sam.fit
##' @param ... extra arguments transferred to plot
##' @details ...
##' @export
tsbplot<-function(fit, ...){
  .plotit(fit, "logtsb", ylab="TSB", trans=exp,...)
}

##' SAM Recruits plot 
##' @param fit the object returned from sam.fit
##' @param ... extra arguments transferred to plot
##' @details ...
##' @export
recplot<-function(fit,...){
   lab<-paste("Recruits (age ", fit$conf$minAge, ")", sep="")
  .plotit(fit, "logR", ylab=lab, trans=exp,...)
}
