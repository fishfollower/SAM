##' Plot SAM model 
##' @param  x ...
##' @param  ... extra arguments
##' @details ...
##' @export
plot.sam<-function(x, ...){
  par(mfrow=c(3,1))
  ssbplot(x,...)
  fbarplot(x,...)
  recplot(x,...)
  par(mfrow=c(1,1))
}

##' Print SAM fit 
##' @param  x ...
##' @param  ... extra arguments
##' @details ...
##' @export
print.sam<-function(x, ...){
 cat("SAM model: log likelihood is", logLik(x,...),"Convergence", ifelse(0==x$opt$convergence, "OK\n", "failed\n"))
}

##' Log likelihood SAM fit 
##' @param  object ...
##' @param  ... extra arguments
##' @details ...
##' @export
logLik.sam<-function(object, ...){
 ret<- -object$opt$objective
 attr(ret,"df")<-length(object$opt$par)
 class(ret)<-"logLik"
 ret
}
