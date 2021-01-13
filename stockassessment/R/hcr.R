##' Harvest control rule forecast
##'
##' @details
##' The formula below is used to determine a new F based on the previous SSB.
##' \deqn{F = \left\{
##'               \begin{array}{ll}
##'                   F_{low} & SSB < B_{low} \\
##'                   min\left(Ftarget, \max\left( F_{lim}, (SSB - B_{lim}) \cdot (F_{target} - F_{lim}) / (B_{trigger}-B_{lim}) \right)\right) & SSB \qe B_{low}
##'               \end{array}\right.
##' }
##' 
##' @param fit A SAM fit
##' @param nYears Number of years to forecast
##' @param Ftarget Target F for high SSB
##' @param Btrigger SSB that triggers the control rule
##' @param Flim 
##' @param Blim 
##' @param Flow F used for SSB < Blow
##' @param Blow SSB for which Flow is used below
##' @param nosim 
##' @param ave.years 
##' @param rec.years 
##' @return hcr object 
##' @author Christoffer Moesgaard Albertsen
##' @export
hcr <- function(fit,
                nYears = 20,
                Ftarget,
                Btrigger,
                Flim = 0,
                Blim = 0,
                Flow = 0,
                Blow = 0,
                nosim = NULL,
                ave.years = max(fit$data$years)+(-4:0),
                rec.years = c(), #max(fit$data$years)+(-9:0),
                nCatchAverageYears = 1,
                customSel = NULL,
                overwriteSelYears = NULL,
                lagR = FALSE,
                splitLD = FALSE,
                addTSB = FALSE,
                biasCorrect = FALSE,
                returnAllYears = FALSE,                    
                ...
                ){

    f <- modelforecast(fit,
                       hcr = rep(1,nYears),
                       hcrConf = c(Ftarget[1], Flim[1], Flow[1], Blim[1], Blow[1], Btrigger[1]),
                       nosim = nosim,
                       ave.years = ave.years,
                       rec.years = rec.years,
                       nCatchAverageYears = nCatchAverageYears,
                       customSel = customSel,
                       overwriteSelYears = overwriteSelYears,
                       lagR = lagR,
                       splitLD = splitLD,
                       addTSB = addTSB,
                       biasCorrect = biasCorrect,
                       returnAllYears = returnAllYears,
                       ...
                       )
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

