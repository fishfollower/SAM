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

##' Extract residuals from sam object 
##' @method residuals sam 
##' @param object sam fitted object (result from sam.fit)
##' @param ... extra arguments for TMB's oneStepPredict
##' @importFrom stats residuals
##' @importFrom TMB oneStepPredict
##' @details ...
##' @export
residuals.sam<-function(object, ...){
 res <- oneStepPredict(object$obj, observation.name="logobs", data.term.indicator="keep",...)
 ret <- cbind(object$data$obs, res)
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


