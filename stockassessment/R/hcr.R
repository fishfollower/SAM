

softmin <- function(a,b){
    0.5 * (a + b - sqrt(1e-4 + (a-b)^2))
}
softmax <- function(a,b){
    0.5 * (a + b + sqrt(1e-4 + (a-b)^2))
}

## hcr <- function(SSB, Ftarget, Flim, Flow, Blim, Blow, Btrigger){
##     if(SSB < Blow)
##         return(Flow);
##     softmin(Ftarget, softmax(Flim, Flim + (SSB-Blim) * (Ftarget - Flim) / (Btrigger-Blim)))
## }




hcr <- function(fit,
                nYears = 20,
                Ftarget,
                Btrigger,
                Flim = 0,
                Blim = 0,
                Flow = 0,
                Blow = 0,
                ave.years = max(fit$data$years)+(-4:0),
                rec.years = c(), #max(fit$data$years)+(-9:0),
                nCatchAverageYears = 1,
                customSel = NULL,
                overwriteSelYears = NULL,
                lagR = FALSE,
                splitLD = FALSE,
                addTSB = FALSE,
                biasCorrect = TRUE,
                returnAllYears = FALSE,                    
                ...
                ){

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

    ave.years <- match(ave.years, fit$data$years) - 1
    if(any(is.na(ave.years)))
        stop("ave.years has years without data.")

    obj0 <- fit$obj
    obj0$fn(fit$opt$par)
    obj0$gr(fit$opt$par)
   
    arg <- as.list(obj0$env)[methods::formalArgs(TMB::MakeADFun)[methods::formalArgs(TMB::MakeADFun) != "..."]]
    pl <- fit$pl
    pl$logF <- cbind(pl$logF,matrix(pl$logF[,ncol(pl$logF)],nrow(pl$logF),nYears))
    pl$logN <- cbind(pl$logN,matrix(pl$logN[,ncol(pl$logN)],nrow(pl$logN),nYears))
    arg$parameters <- pl
    arg$random <- unique(names(obj0$env$par[obj0$env$random]))
    arg$data$simFlag <- c(1,1)
    arg$data$forecast <- list(nYears = as.numeric(nYears),
                               nCatchAverageYears = as.numeric(nCatchAverageYears),
                               aveYears = as.numeric(ave.years),
                               forecastYear = as.numeric(c(rep(0,fit$data$noYears),seq(1,nYears,length=nYears))),
                               FModel = as.numeric(rep(7, nYears)),
                               target = as.numeric(rep(0, nYears)),
                               selectivity = as.numeric(customSel),
                               recModel = as.numeric(recModel),
                               logRecruitmentMedian = as.numeric(logRecruitmentMedian),
                               logRecruitmentVar = as.numeric(logRecruitmentVar),
                               fsdTimeScaleModel = as.numeric(rep(1, nYears)),
                               simFlag = c(0,0),
                               uniroot = as.numeric(FALSE),
                               hcrConf = c(Ftarget, Flim, Flow, Blim, Blow, Btrigger)
                               )
    
    obj <- do.call(TMB::MakeADFun, arg)

    s1 <- replicate(1000,obj$simulate(), simplify = FALSE)
    fbar <- lapply(s1,function(x)exp(x$logfbar))
    ssb <- lapply(s1,function(x)exp(x$logssb))
    catch <- lapply(s1,function(x)exp(x$logCatch))
    
    list(fbar = fbar,
         ssb = ssb,
         catch = catch)    
}
