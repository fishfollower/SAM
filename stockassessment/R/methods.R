##' Plot sam object 
##' @method plot sam
##' @param  x ...
##' @param  ... extra arguments
##' @importFrom graphics par
##' @details ...
##' @export
plot.sam<-function(x, ...){
  par(mfrow=c(3,1))
  ssbplot(x,...)
  fbarplot(x,...)
  recplot(x,...)
  par(mfrow=c(1,1))
}

##' Print sam object 
##' @method print sam 
##' @param  x ...
##' @param  ... extra arguments
##' @details ...
##' @export
print.sam<-function(x, ...){
 cat("SAM model: log likelihood is", logLik.sam(x,...),"Convergence", ifelse(0==x$opt$convergence, "OK\n", "failed\n"))
}

##' Log likelihood of sam object 
##' @method logLik sam 
##' @param  object sam fitted object (result from sam.fit)
##' @param  ... extra arguments
##' @details ...
##' @export
logLik.sam<-function(object, ...){
 ret<- -object$opt$objective
 attr(ret,"df")<-length(object$opt$par)
 class(ret)<-"logLik"
 ret
}

##' Extract number of observations from sam object 
##' @method nobs sam 
##' @param object sam fitted object (result from sam.fit)
##' @param ... extra arguments
##' @importFrom stats nobs
##' @details ...
##' @export
nobs.sam<-function(object, ...){
 as.integer(object$data$nobs)
}
