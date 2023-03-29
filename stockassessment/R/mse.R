
## Function to add simulated years to fit
addSimulatedYears <- function(fit, constraints,resampleFirst=FALSE, ...){
    doSim <- modelforecast(fit, constraints, nosim=1, returnObj=2,addDataYears=TRUE,resampleFirst=resampleFirst, useModelLastN = FALSE, ...)    
    v <- doSim()
    obj <- environment(doSim)$obj
    names(v) <- gsub("dat\\.","",names(v))
    dat <- fit$data
    dat <- dat[!duplicated(names(dat))]
    nms <- intersect(names(dat), names(v))
    dat[nms] <- v[nms]
    nms1 <- c("aux","auxData","idx1","idx2","idxCor","weight")
    dat[nms1] <- obj$env$data[nms1]
    dat$years <- min(as.numeric(dat$aux[,"year"])):max(as.numeric(dat$aux[,"year"]))
    ## Fix dimensionnames
    dmnm <- list(dat$years, fit$conf$minAge:fit$conf$maxAge, NULL)
    dimnames(dat$propMat) <- dimnames(dat$stockMeanWeight) <- dimnames(dat$natMor) <- dimnames(dat$propM) <- dmnm[1:2]
    dimnames(dat$catchMeanWeight) <- dimnames(dat$landFrac) <- dimnames(dat$disMeanWeight) <- dimnames(dat$landMeanWeight) <- dimnames(dat$propF) <- dmnm
    dat$noYears <- length(dat$years)
    cnf <- fit$conf
    pl <- fit$pl
    pl$missing <- NULL
    attr(pl,"what") <- NULL
    nms2 <- intersect(names(pl), names(v))
    pl[nms2] <- v[nms2]    
    map <- obj$env$map
    newFit <- suppressWarnings(sam.fit(dat,cnf,pl,#map=map,
                                       run=FALSE, check.parameters=FALSE))
    newFit$opt <- list(par = newFit$obj$par,
                       objective = NA,
                       convergence = 0)
    plMap <- newFit$pl
    map <- newFit$obj$env$map
    with.map <- intersect(names(plMap), names(map))
    applyMap <- function(par.name) {
        tapply(plMap[[par.name]], map[[par.name]], mean)
    }
    plMap[with.map] <- sapply(with.map, applyMap, simplify = FALSE)         
    p <- unlist(plMap)
    names(p) <- names(newFit$obj$env$last.par)          
    newFit$obj$env$last.par <- newFit$obj$env$last.par.best <- p
    ## SDREPORT
    obj2 <- TMB::MakeADFun(newFit$obj$env$data, newFit$obj$env$parameters, type = "ADFun", 
                      ADreport = TRUE, DLL = newFit$obj$env$DLL, silent = newFit$obj$env$silent)
    
    newFit$rep <- newFit$obj$report(p)
    sdv <- obj2$fn(p)
    sdrep <- list(value = sdv,
                         sd = rep(0,length(sdv)))
    idx <- c(which(names(sdrep$value) == "lastLogN"), which(names(sdrep$value) == 
                                                            "lastLogF"))
    sdrep$estY <- sdrep$value[idx]
    sdrep$covY <- matrix(0,length(idx),length(idx))
    idx <- c(which(names(sdrep$value) == "beforeLastLogN"), which(names(sdrep$value) == 
                                                                  "beforeLastLogF"))
    sdrep$estYm1 <- sdrep$value[idx]
    sdrep$covYm1 <- matrix(0,length(idx),length(idx))
    newFit$sdrep <- sdrep
    class(newFit) <- "sam"
    newFit
}




updateAssessment <- function(OM, EM, knotRange, AdviceLag, intermediateFleets){
    datNew <- OM$data
    confNew <- EM$conf
    if(AdviceLag > 0){
        if(length(intermediateFleets) == 0){ ## Reduce fully
            datNew <- reduce(datNew,
                             year = tail(datNew$years, AdviceLag),
                             conf = confNew)
            confNew <- attr(datNew,"conf")
        }else{ ## Keep some fleets in first year
            ## First remove intermediate years 2+
            if(AdviceLag > 1){
                datNew <- reduce(datNew,
                                 year = tail(datNew$years, AdviceLag-1),
                                 conf = confNew)
                confNew <- attr(datNew,"conf")
            }
            ## Remove fleets not in intermediateFleets for first intermediate year
            if(is.character(intermediateFleets)){
                intermediateFleets <- match(intermediateFleets, attr(datNew,"fleetNames"))
                if(any(is.na(intermediateFleets)))
                    stop("intermediateFleets names does not match model fleet names")
            }
            datNew <- reduce(datNew,
                             year = tail(datNew$years, 1),
                             fleet = setdiff(seq_along(datNew$fleetTypes), intermediateFleets),
                             conf = confNew)
            confNew <- attr(datNew,"conf")            
        }
    }
    ## Prepare parameters
    dp <- defpar(datNew, confNew)
    pl <- EM$pl
    for (nn in intersect(names(dp), names(pl)))
        if (length(dp[[nn]]) == length(pl[[nn]])){
            dp[[nn]][] <- pl[[nn]][]
        }## else if(nn %in% c("logN","logF")){
        ##     nnIndx <- pmin(ncol(dp[[nn]]),ncol(pl[[nn]]))
        ##     dp[[nn]][,seq_len(nnIndx)] <- pl[[nn]][,seq_len(nnIndx),drop=FALSE]
        ##     if(ncol(pl[[nn]]) < ncol(dp[[nn]]))
        ##         dp[[nn]][,-seq_len(nnIndx)] <- do.call("cbind", replicate(ncol(dp[[nn]]) - ncol(pl[[nn]]),pl[[nn]][,ncol(pl[[nn]]),drop=FALSE],simplify=FALSE))
        ## }

    if(length(knotRange) > 0 &&
       (EM$conf$stockRecruitmentModelCode == 90 ||
        EM$conf$stockRecruitmentModelCode == 91 ||
        EM$conf$stockRecruitmentModelCode == 92)){
        capture.output(EM_New <- stockassessment:::getSplineRecBreaks(datNew, confNew, dp, srmc=EM$conf$stockRecruitmentModelCode, knotRange = knotRange, returnFit=TRUE, sim.condRE=c(TRUE,FALSE), silent = TRUE, newtonsteps=0))
    }else{
        capture.output(EM_New <- sam.fit(datNew, confNew, dp, silent = TRUE, newtonsteps=0))
    }
    return(EM_New)
}


##' Management strategy evaluation using SAM models
##'
##' @param OM sam.fit that will work as operating model
##' @param EM sam.fit that will work as estimation model
##' @param nYears Number of years to run simulation
##' @param forecastSettings Settings to do forecast that determines advice
##' @param AdviceYears Number of years advice given at a time. How advice is given is determined by forecastSettings
##' @param AdviceLag Lag between assessment and advice 
##' @param initialAdvice Advice in the first AdviceLag years
##' @param implementationError Function to add implementation error (i.e, transform advice to target catch)
##' @param knotRange Range of spline knot values to try
##' @param intermediateFleets Fleets that are available in the (first) intermediate year
##' @param ... arguments passed on to addSimulatedYears
##' @return a list with MSE result
MSE <- function(OM,
                EM,
                nYears,
                forecastSettings,
                AdviceYears = 1,
                AdviceLag = 0,
                initialAdvice = NA,
                implementationError = function(x) x,
                knotRange = 3,
                intermediateFleets = numeric(0),
                ...){

    if(max(OM$data$aux[OM$data$aux[,"fleet"] %in% which(OM$data$fleetTypes==0),"year"]) != max(OM$data$years))
        stop("Operating model must have catches in the final year")

    AdviceLag <- pmax(0, AdviceLag)
    OM_pl <- OM$pl
    trueSel <- as.numeric(exp(t(OM_pl$logF[OM$conf$keyLogFsta[1,]+1,ncol(OM_pl$logF)])) / tail(fbartable(OM)[,1],1))
    EM_pl <- EM$pl
    
    nYOld <- nYears
    nYears <- max(seq(1,nYears + AdviceYears-1, by = AdviceYears)) + (AdviceYears - 1)
    if(nYOld < nYears)
        message(sprintf("nYears changed to %d to fit the AdviceYears increments.",nYears))
    ssb <- matrix(NA,nYears+AdviceLag,5)
    fbar <- matrix(NA,nYears+AdviceLag,5)
    rec <- matrix(NA,nYears+AdviceLag,5)
    catch <- matrix(NA,nYears+AdviceLag,5)

    rownames(ssb) <- rownames(fbar) <- rownames(rec) <- rownames(catch) <- seq(max(OM$data$years) + 1,len = nYears + AdviceLag)
    colnames(ssb) <- colnames(fbar) <- colnames(rec) <- colnames(catch) <- c("Advice","True","Estimate","Low","High")

    OM_update <- OM
    OM_data <- OM$data
    EM_update <- EM
    msg <- "OK"

    catch[seq_len(AdviceLag),"Advice"] <- rep(initialAdvice,length.out = AdviceLag)

    AdviceToCatchConstraint <- Vectorize(function(x){
        if(is.na(x)) return(NA)
        sprintf("C=%f",implementationError(x))
    })

    if(length(forecastSettings$constraints) != AdviceLag + AdviceYears + 1){
        warning("Length of forecastSettings$constraints should equal AdviceLag + AdviceYears + 1. Modifying to match.")
        if(is.null(forecastSettings$constraints)){
            forecastSettings$constraints <- rep(NA, length.out = AdviceLag + AdviceYears + 1)
        }else{
            forecastSettings$constraints <- rep(forecastSettings$constraints, length.out = AdviceLag + AdviceYears + 1)            
        }        
    }

    if(AdviceLag > 0){
        capture.output(OM_update <- try({addSimulatedYears(OM_update, 
                                                          constraints = AdviceToCatchConstraint(rep(initialAdvice,length.out = AdviceLag)),
                                                          ...)}))
        iy <- head(rownames(ssb),AdviceLag)
        ssb[iy,"True"] <- ssbtable(OM_update)[iy,1]
        fbar[iy,"True"] <- fbartable(OM_update)[iy,1]
        rec[iy,"True"] <- rectable(OM_update)[iy,1]
        catch[iy,"True"] <- catchtable(OM_update)[iy,1]
    }
    for(i in seq(1,nYears-(AdviceYears-1), by = AdviceYears)){ # Index over assessment year
        yr <- rownames(ssb)[seq(i,len=AdviceYears)]
        yr_tac <- rownames(ssb)[seq(i+AdviceLag,len=AdviceYears)]
        
        cat("\n\n\nSimulation year",i,"\n")
        cat("\tAdvice year",yr_tac,"\n")
        cat("\tLast true data year",tail(rownames(fbartable(OM_update)),1),"\n")
        cat("\tLast observed data year",tail(rownames(fbartable(EM_update)),1),"\n")
        cat("\tAssessment year",yr[1],"\n")
        cat("\tCurrent SSB",ssb[yr[1],"True"],"\n")

        fcThisYear <- forecastSettings
        fcThisYear$constraints <- gsub("%ADVICE%",catch[yr,"Advice"],fcThisYear$constraints)
        fcThisYear$constraints[grepl("C=NA",fcThisYear$constraints)] <- NA
        adviceForecast <- try({do.call(modelforecast, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
        if(class(adviceForecast) == "try-error"){
            msg <- "Advice forecast error"
            break;
        }
        afFTab <- attr(adviceForecast,"tab")        
        tabLab <- attr(adviceForecast, "estimateLabel")
        cat("\tAdvice",paste(yr_tac,afFTab[yr_tac,sprintf("catch:%s",tabLab)],sep=": ",collapse="; "),"\n\n\n")
        
        catch[yr_tac,"Advice"] <- afFTab[yr_tac,sprintf("catch:%s",tabLab)]
        fbar[yr_tac,"Advice"] <- afFTab[yr_tac,sprintf("fbar:%s",tabLab)]
        ssb[yr_tac,"Advice"] <- afFTab[yr_tac,sprintf("ssb:%s",tabLab)]
        rec[yr_tac,"Advice"] <- afFTab[yr_tac,sprintf("rec:%s",tabLab)]
        
        ## Simulate next year (implement different rules for AdviceYears > 0)
        capture.output(OM_update <- try({addSimulatedYears(OM_update, 
                                                          constraints = AdviceToCatchConstraint(rep(catch[yr_tac,"Advice"],length.out = AdviceYears)),
                                                          ...)}))
        if(class(OM_update) == "try-error"){
            msg <- "Adding simulated year error"
            break;
        }
        ssb[yr_tac,"True"] <- ssbtable(OM_update)[yr_tac,1]
        fbar[yr_tac,"True"] <- fbartable(OM_update)[yr_tac,1]
        rec[yr_tac,"True"] <- rectable(OM_update)[yr_tac,1]
        catch[yr_tac,"True"] <- catchtable(OM_update)[yr_tac,1]
        
        ## New assessment
        EM_update <- try({updateAssessment(OM_update, EM_update, knotRange, AdviceLag, intermediateFleets)})
        if(class(EM_update) != "sam"){
            msg <- "Assessment error"
            break;
        }
        ssb[yr,c("Estimate","Low","High")] <- ssbtable(EM_update)[yr,,drop=FALSE]
        fbar[yr,c("Estimate","Low","High")] <- fbartable(EM_update)[yr,,drop=FALSE]
        rec[yr,c("Estimate","Low","High")] <- rectable(EM_update)[yr,,drop=FALSE]
        catch[yr,c("Estimate","Low","High")] <- catchtable(EM_update)[yr,,drop=FALSE]        
    }
    ## Last year assessments
    if(msg == "OK")
        for(i0 in seq_len(AdviceLag)){
            yr <- rownames(ssb)[nYears + i0]       
            cat("Assessment for last intermediate year",i0,"/",AdviceLag,"\n")
            cat("\tLast true data year",tail(rownames(fbartable(OM_update)),1),"\n")
            cat("\tLast observed data year",tail(rownames(fbartable(OM_update)),1),"\n")
            cat("\tAssessment year",yr[1],"\n")
            EM_update <- try({updateAssessment(OM_update, EM_update, knotRange, AdviceLag - i0,intermediateFleets)})
            if(class(EM_update) != "sam"){
                msg <- "Assessment error"
                break;
            }
            ssb[yr,c("Estimate","Low","High")] <- ssbtable(EM_update)[yr,,drop=FALSE]
            fbar[yr,c("Estimate","Low","High")] <- fbartable(EM_update)[yr,,drop=FALSE]
            rec[yr,c("Estimate","Low","High")] <- rectable(EM_update)[yr,,drop=FALSE]
            catch[yr,c("Estimate","Low","High")] <- catchtable(EM_update)[yr,,drop=FALSE]        
        }
    

    ## Return
    list(OM = OM_update,
         EM = EM_update,
         ssb = ssb,
         fbar = fbar,
         rec = rec,
         catch = catch,
         msg = msg)

}
