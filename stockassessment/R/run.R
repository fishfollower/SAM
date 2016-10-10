##' Fit SAM model 
##' @param data ...
##' @param conf ...
##' @param parameters ...
##' @param newtonsteps ...
##' @param ... extra arguments to MakeADFun
##' @importFrom TMB MakeADFun sdreport
##' @importFrom stats nlminb optimHess
##' @details ...
##' @useDynLib stockassessment
##' @export
##' @examples
##' data(nscodData)
##' data(nscodConf)
##' data(nscodParameters)
##' fit <- sam.fit(nscodData, nscodConf, nscodParameters)
sam.fit <- function(data, conf, parameters, newtonsteps=3,...){
  tmball <- c(data, conf)
  obj <- MakeADFun(tmball, parameters, random=c("logN", "logF"), DLL="stockassessment", ...)
  opt <- nlminb(obj$par, obj$fn,obj$gr ,control=list(trace=1, eval.max=2000))
  for(i in seq_len(newtonsteps)) { # Take a few extra newton steps 
    g <- as.numeric( obj$gr(opt$par) )
    h <- optimHess(opt$par, obj$fn, obj$gr)
    opt$par <- opt$par - solve(h, g)
    opt$objective <- obj$fn(opt$par)
  }
  sdrep <- sdreport(obj,opt$par)
  pl <- as.list(sdrep,"Est")
  plsd <- as.list(sdrep,"Std")
  ret <- list(sdrep=sdrep, pl=pl, plsd=plsd, data=data, conf=conf, opt=opt)
  class(ret)<-"sam"
  return(ret)
}
