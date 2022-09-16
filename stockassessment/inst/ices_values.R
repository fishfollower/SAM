
##' Function to estimate ICES values
##'
##' The function calculates MSYBtrigger and Bpa from estimated reference points through a simulation forecast. If a sam object is used, other reference points (e.g., Fmsy) are estimated first.
##' 
##' @param x object to calculate for
##' @param nosim Number of simulations
##' @param nyears Number of years to forecast before equilibrium
##' @param ntail Number of years to use for calculations
##' @param ... Other parameters passed to \link{referencepoints}
##' @return ICES values
##' @author Christoffer Moesgaard Albertsen
##' @export
ICESvalues <- function(x,
                        nosim,
                       nyears,
                       ntail,
                        ...){
    UseMethod("ICESvalues")
}


##' @rdname ICESvalues
##' @method ICESvalues sam
##' @export
ICESvalues.sam <- function(x,
                            nosim,
                            nyears,
                            ...){
    rp <- referencepoints(x,...)
    ICESvalues(rp, nosim, nyears)
}

stat_mode <- function(x, ...){
    d <- density(x, n = 100000, ...)
    d$x[which.max(d$y)]
}


getFp05 <- function(fit, rp, Btrigger, nosim, nyears, ntail){
    
    getP <- function(Fval){
        h1 <- hcr(fit, Ftarget = Fval, Btrigger = as.numeric(Btrigger), nosim = nosim, nYears = nyears)
        mean(unlist(lapply(tail(h1$forecast, ntail),function(x)x$ssb)) < rp$tables$Biomass["lim","Estimate"], na.rm = TRUE)
    }
    getObj <- function(F){
        getP(F) - 0.05
    }

    F1 <- rp$tables$F["MSY","Estimate"]
    p1 <- getP(F1)
    while(p1 > 0.05){
        F2 <- F1
        F1 <- 0.5 * F1
        p1 <- getP(F1)
    }
    F2 <- rp$tables$F["lim","Estimate"]
    p2 <- getP(F2)
    while(p2 < 0.05){
        F1 <- F2
        F2 <- 1.5 * F2
        p2 <- getP(F2)
    }    
    opt <- uniroot(getObj, c(F1,F2))
    h <- hcr(fit, Ftarget = opt$root, Btrigger = as.numeric(Btrigger), nosim = nosim, nYears = nyears)
    list(Fp05 = opt$root,
         HCR = h,
         opt <- opt)
}


##' @rdname ICESvalues
##' @method ICESvalues sam_referencepoints
##' @export
ICESvalues.sam_referencepoints <- function(x,
                                           nosim = 1000,
                                           nyears = 100,
                                           ntail = 10,
                                           quantile_buffer_Fmsy = 0.5,
                                           calculate_Fp05 = TRUE,
                                            ...){
    fit <- attr(x,"fit")
    Fmsy <- summary(x$sdr)["referencepoint.logFmsy",1:2]
    FmsyUse <- Fmsy[1]
    if(is.finite(Fmsy[2]) && quantile_buffer_Fmsy > 0 && quantile_buffer_Fmsy < 1)
        FmsyUse <- qnorm(quantile_buffer_Fmsy, Fmsy[1], Fmsy[2])
    cat("\nCalculating MSYBtrigger\n")
    forecast_MSYBtrigger = modelforecast(fit,
                                         fval = exp(rep(FmsyUse, nyears)),
                                         nosim = nosim,
                                         processNoiseF=FALSE,
                                         deterministicF=TRUE,
                                         resampleFirst=FALSE,                                     
                                         rec.years=c(),
                                         ave.years=attr(x,"aveYears"),
                                         overwriteSelYears = attr(x,"selYears"),
                                         biasCorrect=FALSE)
    Blim <- x$tables$Biomass["lim","Estimate"]
    SigmaB <- diff(log(tail(ssbtable(fit),1))[c(1,3)])/2
    Bpa <- Blim * exp(1.645 * SigmaB)
    Flim <- x$tables$F["lim","Estimate"]
    SigmaF <- diff(log(tail(fbartable(fit),1))[c(1,3)])/2
    Fpa2017 <- Flim * exp(1.645 * SigmaF)
    SSB_5pct <- unname(quantile(unlist(lapply(tail(forecast_MSYBtrigger,ntail),
                                              function(x) x$ssb)), 0.05, na.rm = TRUE))
    SSB_mean <- unname(mean(unlist(lapply(tail(forecast_MSYBtrigger,ntail),
                                                 function(x) x$ssb)), na.rm = TRUE))
    SSB_median <- unname(quantile(unlist(lapply(tail(forecast_MSYBtrigger,ntail),
                                                 function(x) x$ssb)), 0.5, na.rm = TRUE))
    SSB_mode <- unname(stat_mode(unlist(lapply(tail(forecast_MSYBtrigger,ntail),
                                               function(x) x$ssb)), na.rm = TRUE))    
    MSYBt <- pmax(Bpa, SSB_5pct, na.rm = TRUE)
    attr(MSYBt,"SigmaB") <- SigmaB
    attr(MSYBt,"SSB_5pct") <- SSB_5pct
    attr(MSYBt,"Bpa") <- Bpa
    class(MSYBt) <- "sam_icesval_MSYBtrigger"
    
    r <- list(MSYBtrigger = MSYBt,
              SSB_summary = c(Mean = SSB_mean,
                              Mode = SSB_mode,
                              Median = SSB_median),
              Bpa = Bpa,
              Fpa2017 = Fpa2017)

    class(r) <- "sam_icesval"    
    attr(r,"fit") <- fit
    attr(r,"referencepoints") <- x
    attr(r,"ssb_samples") <- unlist(lapply(tail(forecast_MSYBtrigger,ntail),
                                           function(x) x$ssb))
    attr(r,"catch_samples") <- unlist(lapply(tail(forecast_MSYBtrigger,ntail),
                                             function(x) x$catch))
    if(calculate_Fp05){
        cat("\nCalculating Fp.05\n")
        Fp050 <- getFp05(fit, x, as.numeric(MSYBt), nosim, nyears, ntail)
        r$Fp05 <- Fp050$F
        attr(r,"HCR_Fp05") <- Fp050$F$HCR
        r$Fpa2021 <- Fp050$F        
    }else{
        r$Fp05 <- NA
        attr(r,"HCR_Fp05") <- NULL
        r$Fpa2021 <- NA        
    }
    r
}

#' @export
print.sam_icesval <- function(x, ...){
    cat("\nEstimated ICES values\n\n")
    v1 <- c("MSYBtrigger:","Bpa:")
    v2 <- sapply(c(x$MSYBtrigger,x$Bpa),function(x)formatC(x,format="f",...))
    nc <- pmax(40, max(nchar(v1))+max(nchar(v2))+1)
    for(i in seq_along(v1))
        cat(v1[i], sprintf(sprintf("%%%ds",30-nchar(v1[i])),v2[i]),"\n")
}
