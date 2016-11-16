##' Plot sam object 
##' @method plot sam
##' @param  x ...
##' @param  ... extra arguments
##' @importFrom graphics par
##' @details ...
##' @export
plot.sam<-function(x, ...){
  op<-par(mfrow=c(3,1))
  ssbplot(x,...)
  fbarplot(x,...)
  recplot(x,...)
  par(op)
}

##' Plot sam object 
##' @method plot samset
##' @param  x ...
##' @param  ... extra arguments
##' @importFrom graphics par
##' @details ...
##' @export
plot.samset<-function(x, ...){
  op<-par(mfrow=c(3,1))
  ssbplot(x,...)
  fbarplot(x,...)
  recplot(x,...)
  par(op)
}

##' Plot sam residuals 
##' @method plot samres 
##' @param  x ...
##' @param  ... extra arguments
##' @details ...
##' @export
##' @examples
##' \dontrun{
##' data(nscodData)
##' data(nscodConf)
##' data(nscodParameters)
##' fit <- sam.fit(nscodData, nscodConf, nscodParameters)
##' par(ask=FALSE)
##' plot(residuals(fit))
##' }
plot.samres<-function(x, ...){
  plotby(x$year, x$age, x$residual, by=x$fleet, xlab="Year", ylab="Age", ...)
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

##' Print samres object 
##' @method print samres 
##' @param  x ...
##' @param  ... extra arguments
##' @details ...
##' @export
print.samres<-function(x, ...){
 class(x)<-NULL
 print(as.data.frame(x))
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

##' Extract fixed coefficients of sam object 
##' @method coef sam 
##' @param  object sam fitted object (result from sam.fit)
##' @param  ... extra arguments
##' @details ...
##' @importFrom stats coef
##' @export
coef.sam <- function(object, ...){
 ret <- object$sdrep$par.fixed
 attr(ret,"cov") <- object$sdrep$cov.fixed
 attr(ret,"sd") <- sqrt(diag(object$sdrep$cov.fixed))
 class(ret)<-"samcoef"
 ret
}

##' Print samcoef object 
##' @method print samcoef 
##' @param  x ...
##' @param  ... extra arguments
##' @details ...
##' @export
print.samcoef<-function(x, ...){
 y<-as.vector(x)
 names(y)<-names(x)
 print(y)
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

##' Extract residuals from sam object 
##' @method residuals sam 
##' @param object sam fitted object (result from sam.fit)
##' @param discrete logical if model contain discrete observations  
##' @param ... extra arguments for TMB's oneStepPredict
##' @importFrom stats residuals
##' @importFrom TMB oneStepPredict
##' @details ...
##' @export
residuals.sam<-function(object, discrete=FALSE, ...){
 res <- oneStepPredict(object$obj, observation.name="logobs", data.term.indicator="keep", discrete=discrete,...)
 ret <- cbind(object$data$aux, res)
 class(ret)<-"samres"
 ret
}

##' Summary of sam object 
##' @method summary sam 
##' @param object sam fitted object (result from sam.fit)
##' @param ... extra arguments 
##' @details ...
##' @export
summary.sam<-function(object, ...){
 ret <- cbind(round(rectable(object)), round(ssbtable(object)), round(fbartable(object),2))
 colnames(ret)[1] <- paste("R(age ", object$conf$minAge, ")", sep="") 
 colnames(ret)[4] <- "SSB"
 colnames(ret)[7] <- paste("Fbar(",object$conf$fbarRange[1], "-", object$conf$fbarRange[2], ")", sep="")
 ret
}
##' Simulate from a sam object 
##' @method simulate sam 
##' @param object sam fitted object (result from sam.fit)
##' @param nsim number of response lists to simulate. Defaults to 1.
##' @param seed random number seed
##' @param ... extra arguments 
##' @details ...
##' @return returns a list of lists. The outer list has length \code{nsim}. Each inner list contains simulated values of \code{logF}, \code{logN}, and \code{obs} with dimensions equal to those parameters.
##' @export
simulate.sam<-function(object, nsim=1, seed=NULL, ...){
 ret <- list()
 if(!is.null(seed)) set.seed(seed)
 for (i in 1:nsim)
 {
   ret[[i]] <- object$obj$simulate()
 }
 ret
}
