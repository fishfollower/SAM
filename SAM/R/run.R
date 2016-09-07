##' Fit SAM model 
##' @param data ...
##' @param parameters ... 
##' @details ...
##' @useDynLib SAM
##' @export
##' @examples
##' data(nscodData)
##' data(nscodParameters)
##' fit <- sam.fit(nscodData, nscodParameters)


sam.fit <- function(data, parameters){
  obj <- MakeADFun(data,parameters,random=c("logN", "logF"), DLL="SAM")
  opt<-nlminb(obj$par,obj$fn,obj$gr,control=list(trace=1, eval.max=1200, iter.max=900))
  sdrep<-sdreport(obj)
  pl <- as.list(sdrep,"Est")
  plsd <- as.list(sdrep,"Std")
  return(list(sdrep=sdrep, pl=pl, plsd=plsd))
}


