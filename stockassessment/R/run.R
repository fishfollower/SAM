##' Fit SAM model 
##' @param data data for the sam model as returned from the setup.sam.data function
##' @param conf model configuration for the model in a format similar to what is returned from the defcon function
##' @param parameters initial values for the model in a format similar to what is returned from the defpar function
##' @param newtonsteps optional extra true extra newton steps
##' @param rm.unidentified option to eliminate unidentified model parameters based on gradient in initial value (somewhat experimental)
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
sam.fit <- function(data, conf, parameters, newtonsteps=3, rm.unidentified=FALSE,...){
  tmball <- c(data, conf)
  ran <- c("logN", "logF")
  obj <- MakeADFun(tmball, parameters, random=ran, DLL="stockassessment", ...)
  if(rm.unidentified){
    skel <- parameters[!names(parameters)%in%ran]
    gr <- obj$gr()
    safemap <- relist(gr,skel)
    safemap <- lapply(safemap, function(x)factor(ifelse(abs(x)>1.0e-15,1:length(x),NA)))
    obj <- MakeADFun(tmball, parameters, random=ran, map=safemap, DLL="stockassessment", ...)
  }
  opt <- nlminb(obj$par, obj$fn,obj$gr ,control=list(trace=1, eval.max=2000, iter.max=1000))
  for(i in seq_len(newtonsteps)) { # Take a few extra newton steps 
    g <- as.numeric( obj$gr(opt$par) )
    h <- optimHess(opt$par, obj$fn, obj$gr)
    opt$par <- opt$par - solve(h, g)
    opt$objective <- obj$fn(opt$par)
  }
  sdrep <- sdreport(obj,opt$par)
  pl <- as.list(sdrep,"Est")
  plsd <- as.list(sdrep,"Std")
  ret <- list(sdrep=sdrep, pl=pl, plsd=plsd, data=data, conf=conf, opt=opt, obj=obj)
  class(ret)<-"sam"
  return(ret)
}
