##' Model based forecast function
##' @rdname modelforecast 
##' @param fit an assessment object of type sam, as returned from the function sam.fit
##' @param ... other variables used by the methods
##' @export
modelforecast <- function(fit, ...){
    UseMethod("modelforecast")
}

##' Model based forecast function
##' @param fscale a vector of f-scales. See details.  
##' @param catchval a vector of target catches. See details.
##' @param fval a vector of target f values. See details.
##' @param nextssb a vector target SSB values the following year. See details
##' @param landval a vector of target catches. See details.   
##' @param findMSY Should not be used. See \link{forecastMSY}.
##' @param hcr Should not be used. See \link{hcr}.
##' @param nosim number of simulations. Not used.
##' @param year.base starting year default last year in assessment. Currently it is only supported to use last assessment year or the year before  
##' @param ave.years vector of years to average for weights, maturity, M and such  
##' @param rec.years vector of years to use to resample recruitment from. If the vector is empty, the stock recruitment model is used.
##' @param label optional label to appear in short table
##' @param overwriteSelYears if a vector of years is specified, then the average selectivity of those years is used (not recommended)
##' @param deterministicF option to set F variance to (almost) zero (not recommended)
##' @param processNoiseF option to turn off process noise in F
##' @param resampleFirst Resample base year when nosim > 0?
##' @param customSel supply a custom selection vector that will then be used as fixed selection in all years after the final assessment year (not recommended)
##' @param lagR if the second youngest age should be reported as recruits
##' @param splitLD if TRUE the result is split in landing and discards
##' @param addTSB if TRUE the total stock biomass (TSB) is added
##' @param biasCorrect Do bias correction of reported variables. Can be turned off to reduce running time (not recommended).
##' @param returnAllYears If TRUE, all years are bias corrected. Otherwise, only forecast years are corrected.
##' @param useUniroot Use uniroot to find catchval, landval, and nextssb?
##' @param nCatchAverageYears Should not be used. See \link{forecastMSY}.
##' @param returnObj Only return TMB object?
##' @param hcrConf Should not be used. See \link{hcr}.
##' @param hcrCurrentSSB Should not be used. See \link{hcr}.
##' @param progress Show progress bar for simulations?
##' @param estimate the summary function used (typically mean or median) for simulations 
##' @details There are four ways to specify a scenario. If e.g. four F values are specified (e.g. fval=c(.1,.2,.3,4)), then the first value is used in the year after the last assessment year (base.year + 1), and the three following in the three following years. Alternatively F's can be specified by a scale, or a target catch. Only one option can be used per year. So for instance to set a catch in the first year and an F-scale in the following one would write catchval=c(10000,NA,NA,NA), fscale=c(NA,1,1,1). If only NA's are specified in a year, the F model is used for forecasting. The length of the vector specifies how many years forward the scenarios run. 
##' @return an object of type samforecast
##' @seealso forecast
##' @importFrom methods formalArgs
##' @importFrom stats var
##' @importFrom utils tail
##' @rdname modelforecast
##' @method modelforecast sam
##' @export
modelforecast.sam <- function(fit,
                     fscale = NULL,
                     catchval = NULL,
                     fval = NULL,
                     nextssb = NULL,
                     landval = NULL,
                     findMSY = NULL,
                     hcr = NULL,
                     nosim = NULL,
                     year.base = max(fit$data$years),
                     ave.years = c(),
                     rec.years = c(), #max(fit$data$years)+(-9:0),
                     label = NULL,
                     overwriteSelYears = NULL,
                     deterministicF = FALSE,
                     processNoiseF = FALSE,
                     resampleFirst = !is.null(nosim) && nosim > 0,
                     customSel = NULL,
                     lagR = FALSE,
                     splitLD = FALSE,
                     addTSB = FALSE,
                     biasCorrect = FALSE,
                     returnAllYears = FALSE,
                     useUniroot = FALSE,
                     nCatchAverageYears = 1,
                     returnObj = FALSE,
                     hcrConf = numeric(0),
                     hcrCurrentSSB = 0,
                     progress = TRUE,
                     estimate = median,
                     ...
                     ){

    if(progress && !returnObj && !is.null(nosim) && nosim > 0){
        pb <- txtProgressBar(min = 0, max = nosim+3, style = 3)
        incpb <- function() setTxtProgressBar(pb, pb$getVal()+1)
    }else{
        incpb <- function(){ return(invisible(NULL)) }
    }


    ## Handle year.base < max(fit$data$years)
    if(year.base > max(fit$data$years)){
        stop("")
    }else if(year.base < max(fit$data$years)){
        warning("year.base is ignored for now")
    }

    if(length(ave.years) == 0){
        useModelBio <- TRUE
        ave.years = max(fit$data$years)+(-4:0)
    }else{
        useModelBio <- FALSE
    }
    
    ## Checks
    if(deterministicF && length(fscale) > 0 && any(!is.na(fscale)) && is.null(customSel))
        warning("Forecasted F values may depend on the last estimated F vector and can therefore have uncertainty. Provide a custom selectivity to circumvent this.")
  
    
    ## Get number of forecast years
    lengthVec <- c(length(fscale),
                   length(catchval),
                   length(fval),length(nextssb),length(landval),
                   length(hcr))
    if(any(lengthVec > 0 & lengthVec < max(lengthVec)))
        stop("All target vectors must have the same length")    
    nYears <- max(length(fscale),length(catchval),length(fval),length(nextssb),length(landval),length(findMSY), length(hcr))

    ## Convert input to an F model code and a target value
    if(is.null(fscale))
        fscale <- rep(NA_real_, nYears)
    if(is.null(catchval))
        catchval <- rep(NA_real_, nYears)
    if(is.null(fval))
        fval <- rep(NA_real_, nYears)
    if(is.null(nextssb))
        nextssb <- rep(NA_real_, nYears)
    if(is.null(landval))
        landval <- rep(NA_real_, nYears)
    if(is.null(findMSY))
        findMSY <- rep(NA_real_, nYears)
    if(is.null(hcr))
        hcr <- rep(NA_real_, nYears)

    if(any(!is.na(nextssb)))
        warning("The nextssb target is not currently supported.")
    if(any(!is.na(nextssb)) && (any(fit$data$propM != 0) || any(fit$data$propM != 0)))
        warning("The nextssb target is at the beginning of the year.")
    
    tab <- rbind(fscale,fval,catchval,nextssb,landval, findMSY, hcr)
    FModel <- apply(tab,2, function(x){
        y <- which(!is.na(x))
        switch(as.character(length(y)),
               "0"=0,
               "1"=y,
               stop("At most one target can be specified per year."))
    })
    target <- apply(tab,2, function(x){
        y <- which(!is.na(x))
        switch(as.character(length(y)),
               "0"=NA_real_,
               "1"=x[y],
               stop("At most one target can be specified per year."))
    })

    ## Use custom selectivity?
    if(is.null(customSel)){
        customSel <- numeric(0)
        if(!is.null(overwriteSelYears)){
            fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)  
            Ftab <- faytable(fit)
            customSel <- colMeans(Ftab[as.integer(rownames(Ftab)) %in% overwriteSelYears, , drop=FALSE])
            customSel <- customSel/mean(customSel[fromto[1]:fromto[2]])
        }
    }else{
        if(!is.null(overwriteSelYears))
            warning("overwriteSelYears is ignored when customSel is given.")
        fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)  
        customSel <- customSel/mean(customSel[fromto[1]:fromto[2]])
    }

    ## Get recruitment model
    if(length(rec.years) == 0){
        recModel <- rep(0,nYears)
        logRecruitmentMedian <- NA_real_
        logRecruitmentVar <- NA_real_
    }else{
        rectab <- rectable(fit)
        recpool <- rectab[rownames(rectab)%in%rec.years,1]
        recModel <- rep(1,nYears)
        logRecruitmentMedian <- median(log(recpool))
        logRecruitmentVar <- stats::var(log(recpool))
    }
    
    ## Get F process time scale model
    ## By default, scale as random walk
    fsdTimeScaleModel <- rep(0,nYears)
    if(deterministicF){ ## 'Zero' variance of F process
        fsdTimeScaleModel <- rep(2,nYears)
    }else if(!processNoiseF){ ## Constant variance of F process
        fsdTimeScaleModel <- rep(1,nYears)
    }
    ## When F model is used, fsdTimeScaleModel should be 1
    fsdTimeScaleModel[FModel == 0] <- 1

    
    ## Convert average years to indices
    ave.yearsIn <- ave.years
    ave.years <- match(ave.years, fit$data$years) - 1
    if(any(is.na(ave.years)))
        stop("ave.years has years without data.")
    
    ## Prepare forecast
    obj0 <- fit$obj
    obj0$fn(fit$opt$par)
    obj0$gr(fit$opt$par)
   
    args <- as.list(obj0$env)[methods::formalArgs(TMB::MakeADFun)[methods::formalArgs(TMB::MakeADFun) != "..."]]
    pl <- fit$pl
    pl$logF <- cbind(pl$logF,matrix(pl$logF[,ncol(pl$logF)],nrow(pl$logF),nYears))
    pl$logN <- cbind(pl$logN,matrix(pl$logN[,ncol(pl$logN)],nrow(pl$logN),nYears))
    if(useModelBio){
        extendBio <- function(x) rbind(x,matrix(x[nrow(x)],nYears,ncol(x)))
        if(nrow(pl$logitMO) > 0){            
            pl$logitMO <- extendBio(pl$logitMO)
        }else{
            warning(sprintf("No MO random effects. Using data average over %s.",paste(ave.yearsIn,collapse=", ")))
        }
        if(nrow(pl$logNM) > 0){            
            pl$logNM <- extendBio(pl$logNM)
        }else{
            warning(sprintf("No NM random effects. Using data average over %s.",paste(ave.yearsIn,collapse=", ")))
        }
        if(nrow(pl$logSW) > 0){            
            pl$logSW <- extendBio(pl$logSW)
        }else{
            warning(sprintf("No SW random effects. Using data average over %s.",paste(ave.yearsIn,collapse=", ")))
        }
        if(nrow(pl$logCW) > 0){            
            pl$logCW <- extendBio(pl$logCW)
        }else{
            warning(sprintf("No CW random effects. Using data average over %s.",paste(ave.yearsIn,collapse=", ")))
        }
    }
    args$parameters <- pl
    args$random <- unique(names(obj0$env$par[obj0$env$random]))
    args$data$forecast <- list(nYears = as.numeric(nYears),
                          nCatchAverageYears = as.numeric(nCatchAverageYears),
                          aveYears = as.numeric(ave.years),
                          forecastYear = as.numeric(c(rep(0,fit$data$noYears),seq(1,nYears,length=nYears))),
                          FModel = as.numeric(FModel),
                          target = as.numeric(target),
                          selectivity = as.numeric(customSel),
                          recModel = as.numeric(recModel),
                          logRecruitmentMedian = as.numeric(logRecruitmentMedian),
                          logRecruitmentVar = as.numeric(logRecruitmentVar),
                          fsdTimeScaleModel = as.numeric(fsdTimeScaleModel),
                          simFlag = c(0,0),
                          uniroot = as.numeric(useUniroot),
                          hcrConf = hcrConf,
                          hcrCurrentSSB = hcrCurrentSSB)
    
    if(any(!is.na(findMSY))){
        args$map$logFScaleMSY <- NULL
        return(args)
    }
    ## Create forecast object
    if(!is.null(nosim) && nosim > 0){
        args$type <- "Fun"
        args$integrate <- NULL
        args$random <- NULL
    }

    ## Done with initial work chunk
    incpb()
    
    obj <- do.call(TMB::MakeADFun, args)
    ## Done with initial MakeADFun
    incpb()

    if(returnObj)
        return(obj)
    
    if(!is.null(nosim) && nosim > 0){
        if(year.base==max(fit$data$years)){
            est <- fit$sdrep$estY
            cov <- fit$sdrep$covY
        }
        if(year.base==(max(fit$data$years)-1)){
            stop("year.base before last assessment year is not implemented yet")
            est <- fit$sdrep$estYm1
            cov <- fit$sdrep$covYm1
        }
        names(est) <- gsub("(^.*[lL]ast)(Log[NF]$)","\\2",names(est))
        i0 <- which(fit$data$year == year.base)
        plMap <- pl
        map <- fit$obj$env$map
        with.map <- intersect(names(plMap), names(map))
        applyMap <- function(par.name) {
            tapply(plMap[[par.name]], map[[par.name]], mean)
        }
        plMap[with.map] <- sapply(with.map, applyMap, simplify = FALSE)
        sniii <- 1
        doSim <- function(){
            sim0 <- est 
            if(resampleFirst)
                sim0 <- rmvnorm(1, mu=est, Sigma=cov)
            estList0 <- split(as.vector(sim0), names(est))
            p <- unlist(plMap)
            names(p) <- rep(names(plMap), times = sapply(plMap,length))
            ## Only works when year.base is last assessment year
            indxN <- matrix(which(names(p) %in% "logN"),nrow=length(estList0$LogN))[,i0]
            indxF <- matrix(which(names(p) %in% "logF"),nrow=length(estList0$LogF))[,i0]
            p[indxN] <- estList0$LogN
            p[indxF] <- estList0$LogF
            v <- obj$simulate(par = p)
            sniii <<- sniii+1
            incpb()
            return(v)
        }
        simvals <- replicate(nosim, doSim(), simplify = FALSE)
        simlist <- vector("list",length(FModel) + 1)
        for(i in 0:(length(FModel))){
            y<-year.base+i
            ii <- i0 + i
            simlist[[i+1]] <- list(sim = do.call("rbind",lapply(simvals,function(x) c(x$logN[,ii], x$logF[,ii]))),
                                   fbar = sapply(simvals,function(x) exp(x$logfbar[ii])),
                                   catch = sapply(simvals,function(x) exp(x$logCatch[ii])),
                                   ssb = sapply(simvals,function(x) exp(x$logssb[ii])),
                                   rec = sapply(simvals,function(x) exp(x$logN[1,ii])),
                                   cwF = rep(NA, nosim),
                                   catchatage = do.call("cbind",lapply(simvals,function(x) exp(x$logCatchAge[,ii]))),
                                   land = sapply(simvals,function(x) exp(x$logLand[ii])),
                                   fbarL = sapply(simvals,function(x) exp(x$logfbarL[ii])),
                                   tsb = sapply(simvals,function(x) exp(x$logtsb[ii])),
                                   year=y)
            rownames(simlist[[i+1]]$catchatage) <- seq(fit$conf$minAge,fit$conf$maxAge,1)
        }
        attr(simlist, "fit")<-fit
        ## Similar to stockassessment::forecast
        collect <- function(x){
            est <- estimate(x)
            quan <- quantile(x, c(.025,.975), na.rm = TRUE)
            c(Estimate=est, low=quan[1], high=quan[2])
        }
        fbar <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$fbar))),3)
        fbarL <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$fbarL))),3)  
        rec <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$rec))))
        ssb <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$ssb))))
        tsb <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$tsb))))
        catch <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$catch))))
        land <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$land))))  
        caytable<-round(do.call(rbind, lapply(simlist, function(xx)apply(xx$catchatage,1,collect)))) 
        tab <- cbind(fbar,rec,ssb,catch)
        if(splitLD){
            tab<-cbind(tab,fbarL,fbar-fbarL,land,catch-land)
        }
        if(addTSB){
            tab<-cbind(tab,tsb)
        }
        ## if(!missing(customWeights)) tab <- cbind(tab,cwF=round(do.call(rbind, lapply(simlist, function(xx)collect(xx$cwF))),3))
        rownames(tab) <- unlist(lapply(simlist, function(xx)xx$year))
        nam <- c("Estimate","low","high")
        basename<-c("fbar:","rec:","ssb:","catch:")
        if(splitLD){
            basename<-c(basename,"fbarL:","fbarD:","Land:","Discard:")    
        }
        if(addTSB){
            basename<-c(basename,"tsb:")    
        }
        ## if(!missing(customWeights)){
        ##     basename<-c(basename,"cwF:")    
        ## }
        colnames(tab)<-paste0(rep(basename, each=length(nam)), nam)
        
        attr(simlist, "tab")<-tab
        shorttab<-t(tab[,grep("Estimate",colnames(tab))])
        rownames(shorttab)<-sub(":Estimate","",paste0(label,if(!is.null(label))":",rownames(shorttab)))
        attr(simlist, "shorttab")<-shorttab
        attr(simlist, "label") <- label    
        attr(simlist, "caytable")<-caytable
        class(simlist) <- "samforecast"
        ## Done with reporting
        incpb()
        close(pb)
        return(simlist)
    }else{
        obj$fn(fit$opt$par)
        ## Get results
        sdr <- TMB::sdreport(obj, fit$opt$par, svd_solve(fit$sdrep$cov.fixed),
                             bias.correct= biasCorrect,
                             skip.delta.method = biasCorrect,
                             bias.correct.control = list(sd = TRUE,
                                                         split = lapply(obj$env$ADreportIndex()[c("logfbar","logssb","logR","logCatch","logtsb","logLagR","logLand","logDis","loglandfbar","logdisfbar")[c(TRUE,TRUE,TRUE,TRUE,addTSB,lagR,splitLD,splitLD,splitLD,splitLD)]], utils::tail, n = nYears + 1 + (fit$data$noYears-1) * as.numeric(returnAllYears))
                                                         )
                             )
        ssdr <- summary(sdr)

        simlist <- list()
        for(i in 0:(length(FModel))){
            y<-year.base+i 
            simlist[[i+1]] <- list(sim=NA, fbar=NA, catch=NA, ssb=NA, rec=NA,
                                   cwF=NA, catchatage=NA, land=NA, fbarL=NA, tsb=NA, year=y)
        }

        attr(simlist, "fit")<-fit

        toCI <- function(x, trans = exp){
            trans(x %*% matrix(c(1,0,1,-2,1,2), nrow = 2, ncol =3))
        }

        indx <- 1:2
        if(biasCorrect)
            indx <- 3:4
        
        fbar <- toCI(ssdr[rownames(ssdr) %in% "logfbar",indx])
        rec <- toCI(ssdr[rownames(ssdr) %in% "logR",indx])
        if(lagR)
            rec <- toCI(ssdr[rownames(ssdr) %in% "logLagR",indx])
        ssb <- toCI(ssdr[rownames(ssdr) %in% "logssb",indx])
        catch <- toCI(ssdr[rownames(ssdr) %in% "logCatch",indx])
        ## caytable<-round(do.call(rbind, lapply(simlist, function(xx)apply(xx$catchatage,1,collect)))) 
        tab <- cbind(fbar,rec,ssb,catch)
        if(splitLD){
            land <- toCI(ssdr[rownames(ssdr) %in% "logLand",indx])
            fbarL <- toCI(ssdr[rownames(ssdr) %in% "loglandfbar",indx])
            dis <- toCI(ssdr[rownames(ssdr) %in% "logDis",indx])
            fbarD <- toCI(ssdr[rownames(ssdr) %in% "logdisfbar",indx])
            tab<-cbind(tab,fbarL,fbarD,land,dis)
        }
        if(addTSB){
            tsb <- toCI(ssdr[rownames(ssdr) %in% "logtsb",indx])
            tab<-cbind(tab,tsb)
        }
        
        ## Handle cwF !!
        
        ## Row and column names on table
        futureYears <- unlist(lapply(simlist, function(xx)xx$year))
        rownames(tab) <- sort(unique(c(fit$data$years, futureYears)))
        fullTable <- tab
        tab <- tab[as.numeric(rownames(tab)) %in% futureYears, , drop = FALSE]
        
        nam <- c("Estimate","low","high")
        basename<-c("fbar:","rec:","ssb:","catch:")
        if(splitLD){
            basename<-c(basename,"fbarL:","fbarD:","Land:","Discard:")    
        }
        if(addTSB){
            basename<-c(basename,"tsb:")    
        }
        ## if(!missing(customWeights)){
        ##     basename<-c(basename,"cwF:")    
        ## }

        colnames(tab)<-paste0(rep(basename, each=length(nam)), nam)
        if(returnAllYears || !biasCorrect)
            colnames(fullTable) <- colnames(tab)
        attr(simlist, "tab")<-tab
        shorttab<-t(tab[,grep("Estimate",colnames(tab))])
        rownames(shorttab)<-sub(":Estimate","",paste0(label,if(!is.null(label))":",rownames(shorttab)))
        attr(simlist, "shorttab")<-shorttab
        attr(simlist, "label") <- label    
        ## attr(simlist, "caytable")<-caytable
        attr(simlist, "fulltab") <- fullTable
        class(simlist) <- "samforecast"
        return(simlist)    
    }
}
