
updateAssessment <- function(fitT, fit0, knotRange, intermediateYears){
    datNew <- fitT$data
    confNew <- fit0$conf
    if(intermediateYears > 0){
        datNew <- reduce(datNew,
                         year = tail(datNew$years, intermediateYears),
                         conf = confNew)
        confNew <- attr(datNew,"conf")
    }
    ## Prepare parameters
    dp <- defpar(datNew, confNew)
    pl <- fit0$pl
    for (nn in intersect(names(dp), names(pl)))
        if (length(dp[[nn]]) == length(pl[[nn]])){
            dp[[nn]] <- pl[[nn]]
        }else if(nn %in% c("logN","logF")){
            nnIndx <- pmin(ncol(dp[[nn]]),ncol(pl[[nn]]))
            dp[[nn]][,seq_len(nnIndx)] <- pl[[nn]][,seq_len(nnIndx),drop=FALSE]
            if(ncol(pl[[nn]]) < ncol(dp[[nn]]))
                dp[[nn]][,-seq_len(nnIndx)] <- do.call("cbind", replicate(ncol(dp[[nn]]) - ncol(pl[[nn]]),pl[[nn]][,ncol(pl[[nn]]),drop=FALSE],simplify=FALSE))
        }

    if(length(knotRange) > 0 &&
       (fit0$conf$stockRecruitmentModelCode == 90 ||
        fit0$conf$stockRecruitmentModelCode == 91 ||
        fit0$conf$stockRecruitmentModelCode == 92)){
        capture.output(f0upN <- try({stockassessment:::getSplineRecBreaks(datNew, confNew, dp, srmc=fit0$conf$stockRecruitmentModelCode, knotRange = knotRange, returnFit=TRUE, sim.condRE=c(TRUE,FALSE), silent = TRUE, newtonsteps=0)}))
    }else{
        capture.output(f0upN <- try({sam.fit(datNew, confNew, dp, silent = TRUE, newtonsteps=0)}))
    }
    return(f0upN)
}



MSE <- function(OM,
                EM,
                nYears,
                forecastSettings,
                AdviceYears = 1,
                AdviceLag = 0,
                initialAdvice = NA,
                implementationError = function(x) x,
                knotRange = 3,
                ...){

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
        sprintf("C = %f",implementationError(x))
    })

    if(length(forecastSettings$constraints) != AdviceLag + AdviceYears){
        warning("Length of forecastSettings$constraints should equal AdviceLag + AdviceYears. Modifying to match.")
        if(is.null(forecastSettings$constraints)){
            forecastSettings$constraints <- rep(NA, length.out = AdviceLag + AdviceYears)
        }else{
            forecastSettings$constraints <- rep(forecastSettings$constraints, length.out = AdviceLag + AdviceYears)
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

        adviceForecast <- do.call(modelforecast, c(list(fit = EM_update, progress=FALSE), forecastSettings))
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
        
        ## Simulate next year
        capture.output(OM_update <- try({addSimulatedYears(OM_update, 
                                                          constraints = AdviceToCatchConstraint(rep(initialAdvice,length.out = AdviceLag)),
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
        EM_update <- updateAssessment(OM_update, EM_update, knotRange, AdviceLag)
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
            EM_update <- updateAssessment(OM_update, EM_update, knotRange, AdviceLag - i0)
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
