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
