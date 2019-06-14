##' @export
forecast <- function(fit, ...){
    UseMethod("forecast")
}

##' @export
forecast.sam <- function(fit,
                     fscale = NULL,
                     catchval = NULL,
                     ## catchval.exact=NULL,
                     fval = NULL,
                     nextssb = NULL,
                     landval = NULL,
                     ## cwF=NULL,
                     nosim = NULL,
                     year.base = max(fit$data$years),
                     ave.years = max(fit$data$years)+(-4:0),
                     rec.years = max(fit$data$years)+(-9:0),
                     label = NULL,
                     overwriteSelYears = NULL,
                     deterministic = FALSE,
                     processNoiseF = TRUE,
                     ## customWeights = NULL,
                     customSel = NULL,
                     lagR = FALSE,
                     splitLD = FALSE,
                     addTSB = FALSE,
                     biasCorrect = TRUE,
                     returnAllYears = FALSE
                     ){

    ## Handle year.base < max(fit$data$years)
    if(year.base > max(fit$data$years)){
        stop("")
    }else if(year.base < max(fit$data$years)){
        warning("year.base is ignored for now")
    }
    
    ## Checks
    if(deterministic && length(fscale) > 0 && any(!is.na(fscale)) && is.null(customSel))
        warning("Forecasted F values may depend on the last estimated F vector and can therefore have uncertainty. Provide a custom selectivity to circumvent this.")

    
    
    
    ## Get number of forecast years
    lengthVec <- c(length(fscale),
                   length(catchval),
                   length(fval),length(nextssb),length(landval))
    if(any(lengthVec > 0 & lengthVec < max(lengthVec)))
        stop("...")    
    nYears <- max(length(fscale),length(catchval),length(fval),length(nextssb),length(landval))

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
    tab <- rbind(fscale,fval,catchval,nextssb,landval)
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
        logRecruitmentVar <- var(log(recpool))
    }

    
    ## Get F process time scale model
    ## By default, scale as random walk
    fsdTimeScaleModel <- rep(0,nYears)
    if(deterministic){ ## 'Zero' variance of F process
        fsdTimeScaleModel <- rep(2,nYears)
    }else if(!processNoiseF){ ## Constant variance of F process
        fsdTimeScaleModel <- rep(1,nYears)
    }
    ## When F model is used, fsdTimeScaleModel should be 1
    fsdTimeScaleModel[FModel == 0] <- 1

    
    ## Convert average years to indices
    ave.years <- match(ave.years, fit$data$years) - 1
    if(any(is.na(ave.years)))
        stop("ave.years has years without data.")
    
    ## Prepare forecast
    obj0 <- fit$obj
    args <- as.list(obj0$env)[formalArgs(TMB::MakeADFun)[formalArgs(TMB::MakeADFun) != "..."]]
    pl <- fit$pl
    pl$logF <- cbind(pl$logF,matrix(pl$logF[,ncol(pl$logF)],nrow(pl$logF),nYears))
    pl$logN <- cbind(pl$logN,matrix(pl$logN[,ncol(pl$logN)],nrow(pl$logN),nYears))
    args$parameters <- pl
    args$random <- unique(names(obj0$env$par[obj0$env$random]))
    args$data$forecast <- list(nYears = as.numeric(nYears),
                               aveYears = as.numeric(ave.years),
                               forecastYear = as.numeric(c(rep(0,fit$data$noYears),seq(1,nYears,length=nYears))),
                               FModel = as.numeric(FModel),
                               target = as.numeric(target),
                               selectivity = as.numeric(customSel),
                               recModel = as.numeric(recModel),
                               logRecruitmentMedian = as.numeric(logRecruitmentMedian),
                               logRecruitmentVar = as.numeric(logRecruitmentVar),
                               fsdTimeScaleModel = as.numeric(fsdTimeScaleModel),
                               simFlag = 0)

    ## Create forecast object
    obj <- do.call(TMB::MakeADFun, args)
    obj$fn(fit$opt$par)

    ## Get results
    sdr <- TMB::sdreport(obj, fit$opt$par, fit$opt$he,
                         bias.correct= biasCorrect,
                         skip.delta.method = biasCorrect,
                         bias.correct.control = list(sd = TRUE,
                                                     split = lapply(obj$env$ADreportIndex()[c("logfbar","logssb","logR","logCatch","logtsb","logLagR","logLand","logDis","loglandfbar","logdisfbar")[c(TRUE,TRUE,TRUE,TRUE,addTSB,lagR,splitLD,splitLD,splitLD,splitLD)]], tail, n = nYears + 1 + (fit$data$noYears-1) * as.numeric(returnAllYears))
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
    
    nam <- c("median","low","high")
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
    shorttab<-t(tab[,grep("median",colnames(tab))])
    rownames(shorttab)<-sub(":median","",paste0(label,if(!is.null(label))":",rownames(shorttab)))
    attr(simlist, "shorttab")<-shorttab
    attr(simlist, "label") <- label    
    ## attr(simlist, "caytable")<-caytable
    attr(simlist, "fulltab") <- fullTable
    class(simlist) <- "samforecast"
    simlist    
}
