##' Fit SAM model 
##' @param data data for the sam model as returned from the setup.sam.data function
##' @param conf model configuration which can be set up using the \code{\link{defcon}} function and then modified.
##' @param parameters initial values which can be set up using the \code{\link{defpar}} function and then modified.
##' @param newtonsteps optional extra true extra newton steps
##' @param rm.unidentified option to eliminate unidentified model parameters based on gradient in initial value (somewhat experimental)
##' @param run if FALSE return AD object without running the optimization
##' @param lower named list with lower bounds for optimization (only met before extra newton steps)
##' @param upper named list with upper bounds for optimization (only met before extra newton steps)
##' @param sim.condRE logical with default \code{TRUE}. Simulated observations will be conditional on estimated values of F and N, rather than also simulating F and N forward from their initial values.
##' @param ... extra arguments to MakeADFun
##' @return an object of class \code{sam}
##' @details ...
##' @importFrom TMB MakeADFun sdreport
##' @importFrom stats nlminb optimHess
##' @importFrom utils relist
##' @useDynLib stockassessment
##' @export
##' @examples
##' data(nscodData)
##' data(nscodConf)
##' data(nscodParameters)
##' fit <- sam.fit(nscodData, nscodConf, nscodParameters)
sam.fit <- function(data, conf, parameters, newtonsteps=3, rm.unidentified=FALSE,run=TRUE, lower=getLowerBounds(parameters), upper=getUpperBounds(parameters), sim.condRE=TRUE, ...){
  tmball <- c(data, conf, simFlag=as.numeric(sim.condRE))    
  nmissing <- sum(is.na(data$logobs))
  parameters$missing <- numeric(nmissing)
  ran <- c("logN", "logF", "missing")
  obj <- MakeADFun(tmball, parameters, random=ran, DLL="stockassessment", ...)
  if(rm.unidentified){
    skel <- parameters[!names(parameters)%in%ran]
    gr <- obj$gr()
    safemap <- relist(gr,skel)
    safemap <- lapply(safemap, function(x)factor(ifelse(abs(x)>1.0e-15,1:length(x),NA)))
    obj <- MakeADFun(tmball, parameters, random=ran, map=safemap, DLL="stockassessment", ...)
  }
  
  lower2<-rep(-Inf,length(obj$par))
  upper2<-rep(Inf,length(obj$par))
  for(nn in names(lower)) lower2[names(obj$par)==nn]=lower[[nn]]
  for(nn in names(upper)) upper2[names(obj$par)==nn]=upper[[nn]]

  if(!run) return( list(sdrep=NA, pl=parameters, plsd=NA, data=data, conf=conf, opt=NA, obj=obj) )
  
  opt <- nlminb(obj$par, obj$fn,obj$gr ,control=list(trace=1, eval.max=2000, iter.max=1000),lower=lower2,upper=upper2)
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

##' Bounds
##' @param parameters initial values for the model in a format similar to what is returned from the defpar function
##' @return a named list
getLowerBounds<-function(parameters){
    list(sigmaObsParUS=rep(-10,length(parameters$sigmaObsParUS)))
}

##' Bounds
##' @param parameters initial values for the model in a format similar to what is returned from the defpar function
##' @return a named list
getUpperBounds<-function(parameters){
    list(sigmaObsParUS=rep(10,length(parameters$sigmaObsParUS)))
}
    
