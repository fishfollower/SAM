##' Plot SAM model 
##' @param  ... ... 
##' @details ...
##' @export
plot.sam<-function(x,...){
 .plotit <-function (what, x=fit$data$years, ylab=what, xlab="Years", trans=function(x)x ,...){
    idx<-names(fit$sdrep$value)==what
    y<-fit$sdrep$value[idx]
    ci<-y+fit$sdrep$sd[idx]%o%c(-2,2)
    plot(x,trans(y), xlab=xlab, ylab=ylab, type="l", lwd=3, ylim=range(c(trans(ci),0)), las=1,...)
    polygon(c(x,rev(x)), y = c(trans(ci[,1]),rev(trans(ci[,2]))), border = gray(.5,alpha=.5), col = gray(.5,alpha=.5))
    grid(col="black")
  }
  par(mfrow=c(3,1))
  .plotit("logssb", ylab="SSB", trans=exp)
  .plotit("logfbar", ylab="Fbar", trans=exp)
  .plotit("logR", ylab="R", trans=exp)
}
