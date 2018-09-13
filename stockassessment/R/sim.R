##' Simulate data from fitted model and re-estimate from each run  
##' @param fit a fitted model returned from sam.fit 
##' @param nsim number of simulations 
##' @param ncores number of cores to be used
##' @importFrom parallel detectCores makeCluster parLapply stopCluster clusterEvalQ
##' @export
simstudy <- function(fit, nsim, ncores = detectCores()){
  simdata <- simulate(fit, nsim=nsim,  full.data=TRUE)
  cl <- makeCluster(ncores) #set up nodes
  clusterEvalQ(cl, {library(stockassessment)}) #load the package to each node
  runs <- parLapply(cl, simdata, function(x)sam.fit(x, fit$conf, fit$obj$env$par))
  stopCluster(cl) #shut it down
  attr(runs, "fit") <- fit
  class(runs)<-"samset"
  runs
}
