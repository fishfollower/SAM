##' Plot SAM model 
##' @param  fit ... 
##' @details ...
##' @export
plot.sam<-function(fit,...){
  par(mfrow=c(3,1))
  .plotit("logssb", ylab="SSB", trans=exp)
  .plotit("logfbar", ylab="Fbar", trans=exp)
  .plotit("logR", ylab="R", trans=exp)
}
