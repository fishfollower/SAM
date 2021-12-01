##' Harvest control rule forecast
##'
##' @param fit A SAM fit
##' @param ... other variables used by the methods
##' @seealso modelforecast
##' @export
hcr <- function(fit, ...){
    UseMethod("hcr")
}



##' Harvest control rule forecast
##'
##' The formula below is used to determine a new F based on the previous SSB.
##' \deqn{F = \left\{
##'               \begin{array}{ll}
##'                   F_{low} & SSB < B_{low} \\
##'                   min\left(Ftarget, \max\left( F_{lim}, (SSB - B_{lim}) \cdot (F_{target} - F_{lim}) / (B_{trigger}-B_{lim}) \right)\right) & SSB \ge B_{low}
##'               \end{array}\right.
##' }
##' 
##' @param fit A SAM fit
##' @param nYears Number of years to forecast
##' @param Ftarget Target F for high SSB
##' @param Btrigger SSB that triggers the control rule
##' @param Flim F used for SSB = Blim
##' @param Blim Between Blim and Btrigger, F values are selected based on linear interpolation from Flim to Ftarget
##' @param Flow F used for SSB < Blow
##' @param Blow SSB for which Flow is used below
##' @param nosim Number of simulations to do. If NULL a model forecast based on the Laplace approximation is used
##' @param ave.years vector of years to average for weights, maturity, M and such  
##' @param rec.years vector of years to use to resample recruitment from. If an empty vector is given, recruitment is based on the fitted model.
##' @param preForecast list of forecast parameters (i.e., fval, fscale, catchval, landval, or nextssb) to use before the HCR
##' @param currentSSB if TRUE, SSB at the begining of the control rule year is used. If FALSE, SSB at the begining of the previous year is used.
##' @param ... additional arguments passed to \link{modelforecast}
##' @return hcr object 
##' @author Christoffer Moesgaard Albertsen
##' @rdname hcr
##' @method hcr sam
##' @export
hcr.sam <- function(fit,
                nYears = 20,
                Ftarget,
                Btrigger,
                Flim = 0,
                Blim = 0,
                Flow = 0,
                Blow = 0,
                nosim = 10000,
                ave.years = max(fit$data$years)+(-4:0),
                rec.years = numeric(0), #max(fit$data$years)+(-9:0),
                preForecast = list(),
                currentSSB = FALSE,
                ...
                ){

    xargs <- list(...)
    if(length(preForecast) > 0){
        preYears <- unique(sapply(preForecast, length))
        if(length(preYears) > 1)
            stop("Forecast variables in preForecast must have the same length.")
        preForecast <- lapply(preForecast, function(x) c(x, rep(NA,nYears)))
    }else{
        preYears <- 0
    }

    args <- c(list(fit = fit,
                 hcr = rep(c(NA,1),c(preYears,nYears)),
                 hcrConf = c(Ftarget[1], Flim[1], Flow[1], Blim[1], Blow[1], Btrigger[1]),
                 hcrCurrentSSB = as.integer(currentSSB),
                 nosim = nosim,
                 ave.years = ave.years,
                 rec.years = rec.years),
              preForecast,
              list(...)
              )
    f <- do.call(modelforecast, args)
    r <- list(Ftarget = Ftarget,
              Flim = Flim,
              Flow = Flow,
              Blim = Blim,
              Blow = Blow,
              Btrigger = Btrigger,
              hcr = function(ssb) .Call("hcrR", ssb, c(Ftarget[1], Flim[1], Flow[1], Blim[1], Blow[1], Btrigger[1]), PACKAGE = "stockassessment"),
              forecast = f)
    attr(r,"fit") <- fit
    class(r) <- "hcr"
    r
}
