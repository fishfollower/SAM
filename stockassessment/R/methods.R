##' Plot SAM model 
##' @param  fit ... 
##' @details ...
##' @export
plot.sam<-function(fit,...){
  par(mfrow=c(3,1))
  ssbplot(fit,...)
  fbarplot(fit,...)
  recplot(fit,...)
  par(mfrow=c(1,1))
}

##' Print SAM fit 
##' @param  fit ... 
##' @details ...
##' @export
print.sam<-function(fit,...){
 cat("SAM model: log likelihood is", logLik(fit,...),"Convergence", ifelse(0==fit$opt$convergence, "OK\n", "failed\n"))
}

##' Log likelihood SAM fit 
##' @param  fit ... 
##' @details ...
##' @export
logLik.sam<-function(fit,...){
 ret<- -fit$opt$objective
 attr(ret,"df")<-length(fit$opt$par)
 class(ret)<-"logLik"
 ret
}
