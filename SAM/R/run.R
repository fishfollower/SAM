##' Fit SAM model 
##' @param data ...
##' @param conf ...
##' @param parameters ...
##' @importFrom TMB MakeADFun sdreport
##' @details ...
##' @useDynLib SAM
##' @export
##' @examples
##' data(nscodData)
##' data(nscodConf)
##' data(nscodParameters)
##' fit <- sam.fit(nscodData, nscodConf, nscodParameters)
sam.fit <- function(data, conf, parameters){
  tmball <- c(data, conf)
  obj <- MakeADFun(tmball,parameters,random=c("logN", "logF"), DLL="SAM")
  opt<-nlminb(obj$par,obj$fn,obj$gr,control=list(trace=1, eval.max=1200, iter.max=900))
  sdrep<-sdreport(obj)
  pl <- as.list(sdrep,"Est")
  plsd <- as.list(sdrep,"Std")
  ret <- list(sdrep=sdrep, pl=pl, plsd=plsd, data=data, conf=conf, opt=opt)
  class(ret)<-"sam"
  return(ret)
}


