##' Fit SAM model 
##' @param data data for the sam model as returned from the setup.sam.data function
##' @param conf model configuration which can be set up using the \code{\link{defcon}} function and then modified.
##' @param parameters initial values which can be set up using the \code{\link{defpar}} function and then modified.
##' @param newtonsteps optional extra true newton steps
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
##' @importFrom utils relist packageDescription
##' @useDynLib stockassessment
##' @export
##' @examples
##' data(nscodData)
##' data(nscodConf)
##' data(nscodParameters)
##' fit <- sam.fit(nscodData, nscodConf, nscodParameters)
sam.fit <- function(data, conf, parameters, newtonsteps=3, rm.unidentified=FALSE, run=TRUE, lower=getLowerBounds(parameters), upper=getUpperBounds(parameters), sim.condRE=TRUE, ...){
  definit <- defpar(data, conf)
  if(!identical(parameters,relist(unlist(parameters), skeleton=definit))){
    warning("Initial values are not consistent, so running with default init values from defpar()")
    parameters<-definit
  }
  data<-clean.void.catches(data,conf)
  tmball <- c(data, conf, simFlag=as.numeric(sim.condRE))
  if(is.null(tmball$resFlag)){tmball$resFlag <- 0}  
  nmissing <- sum(is.na(data$logobs))
  parameters$missing <- numeric(nmissing)
  extraArgs <- list(...)
  map <- list()
  if(length(extraArgs)>0){
    if("map" %in% names(extraArgs)){
      map <- extraArgs[["map"]]
    }
  }
  if(dim(parameters$logP)[1]>1){
    ran <- c("logN", "logF","logP", "missing")
    if(length(map)>0){
      obj <- MakeADFun(tmball, parameters, random=ran, DLL="stockassessment",map=map,...)
    } else {
      obj <- MakeADFun(tmball, parameters, random=ran, DLL="stockassessment",...)
    }
  } else {
    ran <- c("logN", "logF", "missing")
    if(length(extraArgs)>0){
      if("map" %in% names(extraArgs)){
        map <- c(map,list(logP=as.factor(0)))
      }
    } else {
      map <- list(logP=as.factor(0))
    }
    obj <- MakeADFun(tmball, parameters, random=ran, DLL="stockassessment", map=map,...)
  }


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
    mapped <- which(names(opt$par) %in% names(map))
    g <- as.numeric( obj$gr(opt$par) )
    h <- optimHess(opt$par, obj$fn, obj$gr)
    if(length(mapped)>0){
      opt$par[-mapped] <- opt$par[-mapped] - solve(h[-mapped,-mapped], g[-mapped])
    } else {
      opt$par <- opt$par - solve(h, g)
    }
    opt$objective <- obj$fn(opt$par)
  }
  rep <- obj$report()
  sdrep <- sdreport(obj,opt$par)

  # Last two states
  idx <- c(which(names(sdrep$value)=="lastLogN"),which(names(sdrep$value)=="lastLogF"))
  sdrep$estY <- sdrep$value[idx]
  sdrep$covY <- sdrep$cov[idx,idx]

  idx <- c(which(names(sdrep$value)=="beforeLastLogN"),which(names(sdrep$value)=="beforeLastLogF"))
  sdrep$estYm1 <- sdrep$value[idx]
  sdrep$covYm1 <- sdrep$cov[idx,idx]

  pl <- as.list(sdrep,"Est")
  plsd <- as.list(sdrep,"Std")

  sdrep$cov<-NULL # save memory

  ret <- list(sdrep=sdrep, pl=pl, plsd=plsd, data=data, conf=conf, opt=opt, obj=obj, rep=rep, low=lower2, hig=upper2)
  attr(ret, "RemoteSha") <- substr(packageDescription("stockassessment")$RemoteSha, 1, 12)
  attr(ret, "Version") <- packageDescription("stockassessment")$Version
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
    
##' remove void catches
##' @param dat data for the sam model as returned from the setup.sam.data function
##' @param conf model configuration which can be set up using the \code{\link{defcon}} function and then modified
##' @return an updated dataset without the catches where F is fixed to zero
clean.void.catches<-function(dat, conf){
  cfidx <- which(dat$fleetTypes==0)
  aidx <- unique(dat$aux[dat$aux[,2]%in%cfidx,3]-conf$minAge+1)
  faidx <- as.matrix(expand.grid(cfidx, aidx))
  faidx <- faidx[which(conf$keyLogFsta[faidx]== -1),,drop=FALSE]
  rmidx <- paste0(dat$aux[,2],"x",dat$aux[,3]-conf$minAge+1) %in%  paste0(faidx[,1],"x",faidx[,2])
  dat$aux <- dat$aux[!rmidx,]
  dat$logobs <- dat$logobs[!rmidx]
  dat$weight <- dat$weight[!rmidx]
  dat$nobs<-sum(!rmidx)
  dat$minAgePerFleet<-as.integer(tapply(dat$aux[,"age"], INDEX=dat$aux[,"fleet"], FUN=min))
  dat$maxAgePerFleet<-as.integer(tapply(dat$aux[,"age"], INDEX=dat$aux[,"fleet"], FUN=max))
  newyear<-min(as.numeric(dat$aux[,"year"])):max(as.numeric(dat$aux[,"year"]))
  newfleet<-min(as.numeric(dat$aux[,"fleet"])):max(as.numeric(dat$aux[,"fleet"]))
  mmfun<-function(f,y, ff){idx<-which(dat$aux[,"year"]==y & dat$aux[,"fleet"]==f); ifelse(length(idx)==0, NA, ff(idx)-1)}
  dat$idx1<-outer(newfleet, newyear, Vectorize(mmfun,c("f","y")), ff=min)
  dat$idx2<-outer(newfleet, newyear, Vectorize(mmfun,c("f","y")), ff=max)
  dat
}


##' Jitter runs 
##' @param fit a fitted model object as returned from sam.fit
##' @param nojit a list of vectors. Each element in the list specifies a run where the fleets mentioned are omitted
##' @param par initial values to jitter around. The defaule ones are returned from the defpar function
##' @param sd the standard deviation used to jitter the initial values (most parameters are on a log scale, so similar to cv) 
##' @param ncores the number of cores to attemp to use
##' @return A "samset" object, which is basically a list of sam fits
##' @details ...
##' @importFrom parallel detectCores makeCluster clusterEvalQ parLapply stopCluster
##' @importFrom stats rnorm
##' @export
jit <- function(fit, nojit=10, par=defpar(fit$data, fit$conf), sd=.25, ncores=detectCores()){
  parv <- unlist(par)
  pars <- lapply(1:nojit, function(i)relist(parv+rnorm(length(parv),sd=sd), par))
  cl <- makeCluster(ncores) #set up nodes
  clusterEvalQ(cl, {library(stockassessment)}) #load the package to each node
  fits <- parLapply(cl, pars, function(p)sam.fit(fit$data, fit$conf, p, silent = TRUE))
  stopCluster(cl) #shut it down
  attr(fits,"fit") <- fit
  attr(fits,"jitflag") <- 1
  class(fits) <- c("samset")
  fits
}
