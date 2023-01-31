##' rmvnorm helper function to draw multivariate normal samples
##' @param n the number of samples.
##' @param mu the mean vector.
##' @param Sigma a positive-definite symmetric matrix specifying the covariance matrix.
##' @details Generates samples via the Cholesky decomposition, which is less platform dependent than eigenvalue decomposition.
##' @return If n = 1 a vector of the same length as mu, otherwise an n by length(mu) matrix with one sample in each row.
##' @export
rmvnorm <- function(n = 1, mu, Sigma){
    p <- length(mu)
    if(!all(dim(Sigma) == c(p, p))){
        stop("incompatible arguments")
    }
    idx <- diag(Sigma) > .Machine$double.xmin 
    L <- matrix(0,p,p)
    if(any(idx)){
        L[idx,idx] <- chol(Sigma[idx,idx])
    }
    X <- matrix(rnorm(p * n), n)
    X <- drop(mu) + t(X%*%L)
    if(n == 1){
        drop(X)
    }else{
        t(X)
    }
}

##' forecast function to do shortterm
##' @param fit an assessment object of type sam, as returned from the function sam.fit
##' @param fscale a vector of f-scales. See details.  
##' @param catchval a vector of target catches. See details.
##' @param catchval.exact a vector of target catches which will be met without noise. See details.   
##' @param fval a vector of target f values. See details.
##' @param nextssb a vector target SSB values the following year. See details
##' @param landval a vector of target catches. See details.   
##' @param cwF a vector target custom weighted F values. customWeights must also be specified
##' @param nosim number of simulations default is 1000
##' @param year.base starting year default last year in assessment. Currently it is only supported to use last assessment year or the year before  
##' @param ave.years vector of years to average for weights, maturity, M and such  
##' @param rec.years vector of years to use to resample recruitment from
##' @param label optional label to appear in short table
##' @param overwriteSelYears if a vector of years is specified, then the average selectivity of those years is used (not recommended)
##' @param deterministic option to turn all process noise off (not recommended, as it will likely cause bias)
##' @param processNoiseF option to turn off process noise in F
##' @param customWeights a vector of same length as number of age groups giving custom weights (currently only used for weighted average of F calculation)
##' @param customSel supply a custom selection vector that will then be used as fixed selection in all years after the final assessment year (not recommended)
##' @param lagR if the second youngest age should be reported as recruits
##' @param splitLD if TRUE the result is split in landing and discards
##' @param addTSB if TRUE the total stock biomass (TSB) is added
##' @param useSWmodel if TRUE the catch mean weight predicted from the assessment model is used (can only be used for configurations supporting this)
##' @param useCWmodel if TRUE the catch mean weight predicted from the assessment model is used (can only be used for configurations supporting this)
##' @param useMOmodel if TRUE the proportion mature predicted from the assessment model is used (can only be used for configurations supporting this)
##' @param useNMmodel if TRUE the natural mortality predicted from the assessment model is used (can only be used for configurations supporting this)
##' @param savesim save the individual simulations 
##' @param cf.cv.keep.cv exotic option 
##' @param cf.cv.keep.fv exotic option
##' @param cf.keep.fv.offset exotic option
##' @param estimate the summary function used (typically mean or median) 
##' @details There are three ways to specify a scenario. If e.g. four F values are specified (e.g. fval=c(.1,.2,.3,4)), then the first value is used in the last assessment year (base.year), and the three following in the three following years. Alternatively F's can be specified by a scale, or a target catch. Only one option can be used per year. So for instance to set a catch in the first year and an F-scale in the following one would write catchval=c(10000,NA,NA,NA), fscale=c(NA,1,1,1). The length of the vector specifies how many years forward the scenarios run. 
##' @return an object of type samforecast
##' @importFrom stats median uniroot quantile plogis
##' @importFrom Matrix bdiag
##' @importFrom TMB sdreport MakeADFun
##' @export
forecast <- function(fit,
                     fscale=NULL,
                     catchval=NULL,
                     catchval.exact=NULL,
                     fval=NULL,
                     nextssb=NULL,
                     landval=NULL,
                     cwF=NULL,
                     nosim=1000,
                     year.base=max(fit$data$years),
                     ave.years=max(fit$data$years)+(-4:0),
                     rec.years=max(fit$data$years)+(-9:0),
                     label=NULL, overwriteSelYears=NULL, deterministic=FALSE, processNoiseF=TRUE,  customWeights=NULL, customSel=NULL, lagR=FALSE, splitLD=FALSE, addTSB=FALSE, useSWmodel=(fit$conf$stockWeightModel>=1), useCWmodel=(fit$conf$catchWeightModel>=1), useMOmodel=(fit$conf$matureModel>=1), useNMmodel=(fit$conf$mortalityModel>=1), savesim=FALSE, cf.cv.keep.cv=matrix(NA, ncol=2*sum(fit$data$fleetTypes==0), nrow=length(catchval)), cf.cv.keep.fv=matrix(NA, ncol=2*sum(fit$data$fleetTypes==0), nrow=length(catchval)), cf.keep.fv.offset=matrix(0, ncol=sum(fit$data$fleetTypes==0), nrow=length(catchval)), estimate=median){
    ## if(sum(fit$data$fleetTypes==0) > 1)
    ##     stop("Forecast for multi fleet models not implemented yet")
    ##estimateLabel <- deparse1(substitute(estimate))
    estimateLabel <- paste(deparse(substitute(estimate), 500L), collapse = " ")
    idxN <- 1:nrow(fit$rep$nvar)
    
    idxF <- 1:nrow(fit$rep$fvar)+nrow(fit$rep$nvar)

    catchFleets <- which(fit$data$fleetType==0)
    noCatchFleets <- length(catchFleets)  
    idxFbyFleet <- lapply(catchFleets, function(i){x<-fit$conf$keyLogFsta[i, , drop = FALSE]; unique(x[x > -0.5]) + nrow(fit$rep$nvar) + 1})
    resample <- function(x, ...){
        if(deterministic){
            ret <- mean(x)
        }else{
            ret <- x[sample.int(length(x), ...)]
        }
        return(ret)
    }
    ns<-max(length(fscale), length(catchval), length(catchval.exact), length(fval), length(nextssb), length(cwF))    
    if(missing(fscale)&missing(fval)&missing(catchval)&missing(catchval.exact)&missing(nextssb)&missing(cwF))stop("No scenario is specified")    
    if(missing(fscale)) fscale <- rep(NA,ns)
    if(missing(fval)) fval <- rep(NA,ns)
    if(missing(catchval)) catchval <- rep(NA,ns)
    if(missing(catchval.exact)) catchval.exact <- rep(NA,ns)    
    if(missing(nextssb)) nextssb <-rep(NA,ns)
    if(missing(landval)) landval <-rep(NA,ns)  
    if(missing(cwF)) cwF <-rep(NA,ns)  

    if(!all(rowSums(!is.na(cbind(fscale, catchval, catchval.exact, fval, nextssb, landval, cwF)))==1)){
        stop("For each forecast year exactly one of fscale, catchval or fval must be specified (all others must be set to NA)")
    }

    if(!is.null(overwriteSelYears)&(noCatchFleets!=1))stop("overwriteSelYears option can only be used for single fleet models")
    
    if(!is.null(overwriteSelYears)){
        fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)  
        Ftab <- faytable(fit)
        fixedsel <- colMeans(Ftab[as.integer(rownames(Ftab))%in%overwriteSelYears,,drop=FALSE])
        fixedsel <- fixedsel/mean(fixedsel[fromto[1]:fromto[2]])
    }

    if(!is.null(customSel)&(noCatchFleets!=1))stop("customSel can only be used for single fleet models")
    
    if(!is.null(customSel)){
        fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)  
        customSel <- customSel/mean(customSel[fromto[1]:fromto[2]])
    }
    
    odat<-fit$obj$env$data
    opar<-fit$obj$env$parList(par = fit$obj$env$last.par.best)
    omap<-fit$obj$env$map
    omap[names(which(lapply(opar,length)==0))]<-NULL
    if(useSWmodel & (fit$conf$stockWeightModel==0)){
        stop("stockWeightModel cannot be used for forecast when it was not part of the fitted model")
    }
    if(useSWmodel){
        rnSW<-1:(nrow(opar$logSW)+ns)+as.integer(rownames(odat$stockMeanWeight)[1])-1  
        opar$logSW<-matrix(0,nrow=nrow(opar$logSW)+ns, ncol=ncol(opar$logSW))
        oran<-unique(names(fit$obj$env$par[fit$obj$env$random]))
        objsw <- MakeADFun(odat, opar, random = oran, DLL = "stockassessment", map=omap)
        sdrep<- sdreport(objsw, par.fixed=fit$opt$par, ignore.parm.uncertainty=TRUE)
        idx<-names(sdrep$value)=="logSW"
        simLogSw <-rmvnorm(nosim, sdrep$value[idx], sdrep$cov[idx,idx])
    }
    if(useCWmodel & (fit$conf$catchWeightModel==0)){
        stop("catchWeightModel cannot be used for forecast when it was not part of the fitted model")
    }
    if(useCWmodel){
        rnCW<-1:(nrow(opar$logCW)+ns)+as.integer(rownames(odat$catchMeanWeight)[1])-1      
        ## opar$logCW<-matrix(0,nrow=nrow(opar$logCW)+ns, ncol=ncol(opar$logCW))
        opar$logCW <- array(0, dim = dim(opar$logCW) + c(ns,0,0))
        oran<-unique(names(fit$obj$env$par[fit$obj$env$random]))
        objcw <- MakeADFun(odat, opar, random = oran, DLL = "stockassessment", map=omap)
        objcw$fn(fit$opt$par)
        sdrep<- sdreport(objcw, par.fixed=fit$opt$par, ignore.parm.uncertainty=TRUE,)
        idx<-names(sdrep$value)=="logCW"
        simLogCW <-rmvnorm(nosim, sdrep$value[idx], sdrep$cov[idx,idx])
    }
    if(useNMmodel & (fit$conf$mortalityModel==0)){
        stop("mortalityModel cannot be used for forecast when it was not part of the fitted model")
    }
    if(useNMmodel){
        rnNM<-1:(nrow(opar$logNM)+ns)+as.integer(rownames(odat$natMor)[1])-1            
        opar$logNM<-matrix(0,nrow=nrow(opar$logNM)+ns, ncol=ncol(opar$logNM))
        oran<-unique(names(fit$obj$env$par[fit$obj$env$random]))
        objnm <- MakeADFun(odat, opar, random = oran, DLL = "stockassessment", map=omap)
        sdrep<- sdreport(objnm, par.fixed=fit$opt$par, ignore.parm.uncertainty=TRUE)
        idx<-names(sdrep$value)=="logNM"
        simLogNM <-rmvnorm(nosim, sdrep$value[idx], sdrep$cov[idx,idx])
    }
    if(useMOmodel & (fit$conf$matureModel==0)){
        stop("matureModel cannot be used for forecast when it was not part of the fitted model")
    }
    if(useMOmodel){
        rnMO<-1:(nrow(opar$logitMO)+ns)+as.integer(rownames(odat$propMat)[1])-1                  
        opar$logitMO<-matrix(0,nrow=nrow(opar$logitMO)+ns, ncol=ncol(opar$logitMO))
        oran<-unique(names(fit$obj$env$par[fit$obj$env$random]))
        objmo <- MakeADFun(odat, opar, random = oran, DLL = "stockassessment", map=omap)
        sdrep<- sdreport(objmo, par.fixed=fit$opt$par, ignore.parm.uncertainty=TRUE)
        idx<-names(sdrep$value)=="logitMO"
        simLogitMO <-rmvnorm(nosim, sdrep$value[idx], sdrep$cov[idx,idx])
    }
    
    getF <- function(x, fleet=which(fit$data$fleetTypes==0), allowSelOverwrite = FALSE){
        getfleet <- function(f){
            idx <- fit$conf$keyLogFsta[f,]+2    
            ret <- c(NA,exp(x[idxF]))[idx]
            ret[is.na(ret)] <- 0
            ret
        }
        FF <- lapply(fleet,getfleet)
        ret <- Reduce("+",FF) 
        if(allowSelOverwrite){
            if(!is.null(overwriteSelYears)){
                fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)    
                thisfbar<-mean(ret[fromto[1]:fromto[2]])
                ret<-fixedsel*thisfbar
            }
            if(!is.null(customSel)){
                fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)    
                thisfbar<-mean(ret[fromto[1]:fromto[2]])
                ret<-customSel*thisfbar
            }
            FF[]<-lapply(fleet, function(x)ret)
        }
        attr(ret,"byFleet") <- do.call(cbind,FF)
        ret
    }

    fbar<-function(x){
        fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)  
        mean(getF(x)[fromto[1]:fromto[2]])
    }
    
    fbarFrac<-function(x, lf){
        fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)  
        mean((lf*getF(x))[fromto[1]:fromto[2]])
    }    

    getCWF<-function(x,w){
        sum(getF(x)*w)
    }    

    fbarbyfleet <- function(x) {
        fromto <- fit$conf$fbarRange - (fit$conf$minAge - 1)
        f<-function(fleet)mean(getF(x, fleet)[fromto[1]:fromto[2]])
        sapply(which(fit$data$fleetTypes == 0),f)
    }
    
    getN <- function(x){
        idx <- fit$conf$keyLogFsta[1,]+1
        nsize <- length(idx)        
        ret <- exp(x[1:nsize])
        ret
    }

    getState <- function(N,F){
        fbf <- fbf <- attr(getF(est),"byFleet")
        ka <- which(fit$conf$keyLogFsta>(-1),TRUE)
        Findx <- apply(ka, 1, function(j) fit$conf$keyLogFsta[j[1],j[2]])
        F <- apply(ka, 1, function(j) fbf[j[2],j[1]])
        F_unique_ordered <- unname(unlist(lapply(split(F,Findx),unique)))
        names(F_unique_ordered) <- rep("lastLogF",length(F_unique_ordered))
        ## k <- fit$conf$keyLogFsta[1,]
        ## F <- F[k>=0]
        ## k <- !duplicated(k[k>=0]) 
        x <- log(c(N,F_unique_ordered))
        x
    }    


    getProcessVar <- function(fit, set.rec.var.zero=TRUE){
        varList <- fit$rep[c("nvar","fvar")]
        if(processNoiseF==FALSE){
            varList$fvar<-varList$fvar*0
        }
        cov <- bdiag(varList)
        if(set.rec.var.zero){
            cov[1,1] <- 0
        }
        as.matrix(cov)
    }
    
    step <- function(x, nm, recpool, scale, inyear=FALSE){
        F <- getF(x, allowSelOverwrite=!inyear)
        N <- getN(x)
        if(!inyear){
            Z <- F+nm
            n <- length(N)
            N <- c(resample(recpool,1),N[-n]*exp(-Z[-n])+c(rep(0,n-2),N[n]*exp(-Z[n])))
            ## Recruitment does not use model prediction, so no effect of logNMeanAssumption
            if(fit$conf$logNMeanAssumption[2] == 1){ # Mean on natural scale
                N <- N + c(0,-0.5 * exp(2.0 * fit$pl$logSdLogN[conf$keyVarLogN[-1]]))
            }else if(fit$conf$logNMeanAssumption[2] == 2){ # Mode on natural scale
                N <- N + c(0,exp(2.0 * fit$pl$logSdLogN[fit$conf$keyVarLogN[-1]]))                
            }
        }
        xx <- rep(NA,length=length(x))
        xx[idxN] <- log(N)
                                        #if(inyear==FALSE)browser()
        if((noCatchFleets==1)){
          kk<-fit$conf$keyLogFsta[1,]+1
          xx[idxF] <- log(tapply(F[kk>0],kk[kk>0],mean))[idxF-min(idxF)+1]+log(scale)
        }else{
          xx[idxF] <- x[idxF]+log(scale)
        }
        return(xx)
    }

    scaleF <- function(x, scale){
        x[idxF] <- x[idxF]+log(scale)
        return(x)
    }
    
    scaleFbyFleet <- function(x, scale){
        f <- function(i){
            x[idxFbyFleet[[i]]] <<- x[idxFbyFleet[[i]]]+log(scale[i])
        }
        invisible(sapply(1:noCatchFleets,f))
        return(x)
    }
    
    sel<-function(x){
        getF(x)/fbar(x)
    }

    catch <- function(x, nm, cw){
        F <- getF(x)
        FF <- attr(F,"byFleet")
        Z <- F+nm
        N <- getN(x)
        C <- FF/Z*(1-exp(-Z))*N
        TCW <- sum(C*cw)
        attr(TCW, "byFleet") <- colSums(C*cw)
        return(TCW)
    }

    catchFrac <- function(x, nm, w, frac){
        F <- getF(x)
        Z <- F+nm
        N <- getN(x)
        C <- F/Z*(1-exp(-Z))*N
        return(sum(frac*w*C))
    }
    
    catchatage <- function(x, nm){
        F <- getF(x)
        Z <- F+nm
        N <- getN(x)
        C <- F/Z*(1-exp(-Z))*N
        return(C)
    }

    ssb <- function(x, nm, sw, mo, pm, pf){
        F <- getF(x)
        FF <- attr(F,"byFleet")
        ZZ <- pm*nm+rowSums(pf*FF)
        N <- getN(x)*exp(-ZZ)
        return(sum(N*mo*sw))
    }

    tsb <- function(x, sw){
        F <- getF(x)
        N <- getN(x)
        return(sum(N*sw))
    }

    rectab<-rectable(fit)
    recpool<-rectab[rownames(rectab)%in%rec.years,1]

                                        # Get final state
    if(year.base==max(fit$data$years)){
        est <- fit$sdrep$estY
        cov <- fit$sdrep$covY
    }
    if(year.base==(max(fit$data$years)-1)){  
        est <- fit$sdrep$estYm1
        cov <- fit$sdrep$covYm1
    }
    if(year.base<(max(fit$data$years)-1)){
        stop("State not saved, so cannot proceed from this year")
    }
    if(deterministic)cov<-cov*0
    sim<-rmvnorm(nosim, mu=est, Sigma=cov)

    if(is.null(overwriteSelYears) & is.null(customSel))  
        if(!isTRUE(all.equal(est,getState(getN(est),getF(est)))))
            stop("Sorry somthing is wrong here (check code for getN, getF, and getState)")  
    doAve <- function(x){
        if(length(dim(x))==2){
            ret <- colMeans(x[rownames(x)%in%ave.years,,drop=FALSE], na.rm=TRUE)
        }
        if(length(dim(x))==3){
            ret <- apply(x[rownames(x)%in%ave.years,,,drop=FALSE],c(2,3),mean, na.rm=TRUE)
        }
        ret
    }

    ave.sw <- doAve(fit$data$stockMeanWeight)
    ave.cw <- doAve(fit$data$catchMeanWeight)
    ave.mo <- doAve(fit$data$propMat)
    ave.nm <- doAve(fit$data$natMor)
    ave.lf <- doAve(fit$data$landFrac)
    ave.lw <- doAve(fit$data$landMeanWeight)
    ave.pm <- doAve(fit$data$propM)
    ave.pf <- doAve(fit$data$propF)
    getThisOrAve <- function(x,y, ave){
        if(y %in% rownames(x)){
            if(length(dim(x))==2){  
                ret <- x[which(rownames(x)==y),]
            }else{
                ret <- x[which(rownames(x)==y),,]
            }
            ret[is.na(ret)]<-ave[is.na(ret)]
        }else{
            ret <- ave
        }
        ret
    }
    procVar<-getProcessVar(fit)  
    simlist<-list()
    rs <- NULL  
    for(i in 0:(length(fscale)-1)){    
        y<-year.base+i  
        if(!is.null(rs)){
            assign(".Random.seed", rs, envir = .GlobalEnv)
            rs <- NULL    
        }
        sw<-getThisOrAve(fit$data$stockMeanWeight, y, ave.sw)
        cw<-getThisOrAve(fit$data$catchMeanWeight, y, ave.cw)
        mo<-getThisOrAve(fit$data$propMat, y, ave.mo)
        nm<-getThisOrAve(fit$data$natMor, y, ave.nm)
        lf<-getThisOrAve(fit$data$landFrac, y, ave.lf)    
        lw<-getThisOrAve(fit$data$landMeanWeight, y, ave.lw)
        pm<-getThisOrAve(fit$data$propM, y, ave.pm)
        pf<-getThisOrAve(fit$data$propF, y, ave.pf)
        if(useNMmodel){
            sim<-t(sapply(1:nrow(sim),
                          function(s){
                              yidx<-which(rnNM==as.integer(y))   
                              thisnm<-matrix(exp(simLogNM[s,]),nrow=nrow(opar$logNM))[yidx,]
                              step(sim[s,], nm=thisnm, recpool=recpool, scale=1, inyear=(i==0))
                          }
                          ))      
        }else{
            sim <- t(apply(sim, 1, function(s)step(s, nm=nm, recpool=recpool, scale=1, inyear=(i==0)))) 
        }
        if(i!=0){  
            cof <- coef(fit)
            if(deterministic){procVar<-procVar*0}
            if(!is.null(overwriteSelYears)){nn<-length(fit$conf$keyLogFsta[1,]); procVar[-c(1:nn),-c(1:nn)] <- 0}
            if(!is.null(customSel)){nn<-length(fit$conf$keyLogFsta[1,]); procVar[-c(1:nn),-c(1:nn)] <- 0}
            if(all(fit$conf$corFlag <3)){
                sim <- sim + rmvnorm(nosim, mu=rep(0,nrow(procVar)), Sigma=procVar)
            }else if(length(fit$conf$corFlag)==1 && fit$conf$corFlag[1] ==3){
                k<-fit$conf$keyLogFsta[1,]
                k <- unique(k[k >= 0] + 1)
                nF <- length(k)
                simTmp <- rmvnorm(nosim, mu=rep(0,nrow(procVar)), Sigma=procVar)
                rhoF = 2/(1 + exp(-2*cof[names(cof)=="sepFlogitRho"])) -1
                sepFalpha = cof[names(cof)=="sepFalpha"]
                for(jj in 1:nosim){
                                        #Calculates previous U and V.
                    V = mean(sim[jj,(dim(sim)[2]-nF+1):(dim(sim)[2])]) 
                    U = sim[jj,(dim(sim)[2]-nF+1):(dim(sim)[2]-1)] -V   
                    
                                        #Extract new contributions to U and V
                    Uepsilon = simTmp[jj,(dim(simTmp)[2]-nF+1):(dim(simTmp)[2]-1)]
                    Vepsilon = simTmp[jj,dim(simTmp)[2]]
                    
                                        #Calculates updated U and V
                    Unew = rhoF[1]*U + Uepsilon+ sepFalpha[1:(length(sepFalpha)-1)] 
                    Vnew = rhoF[2]*V + Vepsilon+ sepFalpha[length(sepFalpha)]
                    
                                        #Calculate F based on U and V
                    sim[jj,(dim(simTmp)[2]-nF+1):(dim(simTmp)[2]-1)] = Unew + Vnew
                    sim[jj,dim(simTmp)[2]] =  -sum(Unew) + Vnew #subtract sum(Unew) since U sum to 0 over all ages
                    
                                        #Update process for N
                    sim[jj,1:(dim(simTmp)[2]-nF)] = sim[jj,1:(dim(simTmp)[2]-nF)]   + simTmp[jj,1:(dim(simTmp)[2]-nF)] #Simulated N is not affected by the separable structure of F
                }
            }else{
                stop("forcast not implemented for the given corflag")
            }
        }
        
        
        if(!is.na(fscale[i+1])){
            sim<-t(apply(sim, 1, scaleF, scale=fscale[i+1]))    
        }

        if(!is.na(fval[i+1])){
            curfbar<-estimate(apply(sim, 1, fbar))
            adj<-fval[i+1]/curfbar
            sim<-t(apply(sim, 1, scaleF, scale=adj))    
        }

        if(!is.na(cwF[i+1])){
            if(missing(customWeights))stop("customWeights must be supplied when using the cwF option")  
            curcwF<-estimate(apply(sim, 1, getCWF, w=customWeights))
            adj<-cwF[i+1]/curcwF
            sim<-t(apply(sim, 1, scaleF, scale=adj))    
        }

        if(!is.na(catchval[i+1])){
            simtmp<-NA
            fun<-function(s){
                simtmp<<-t(apply(sim, 1, scaleF, scale=s))
                if(useCWmodel | useNMmodel){  
                    simcat<-sapply(1:nrow(simtmp), function(k){
                        if(useCWmodel){      
                            thisy <- which(rnCW==y);
                            thiscw <-matrix(exp(simLogCW[k,]),nrow=nrow(opar$logCW))[thisy,];
                        }else{
                            thiscw <- cw
                        }  
                        if(useNMmodel){      
                            thisy <- which(rnNM==y);
                            thisnm <-matrix(exp(simLogNM[k,]),nrow=nrow(opar$logNM))[thisy,];
                        }else{
                            thisnm <- nm
                        }  
                        catch(simtmp[k,], nm=thisnm, cw=thiscw)}
                        )
                }else{ 
                    simcat<-apply(simtmp, 1, catch, nm=nm, cw=cw)
                }
                return(catchval[i+1]-estimate(simcat))
            }
            ff <- uniroot(fun, c(0,100))$root
            sim <- simtmp
        }

        ## <<<<<<< HEAD !!!Conflict block 1; From master - NOT FULLY MERGED!!!
        if(!is.na(catchval.exact[i+1])){
            simtmp<-sim
            funk<-function(k){
                if(useCWmodel){
                    thisy <- which(rnCW==y)
                    thiscw <-matrix(exp(simLogCW[k,]),nrow=nrow(opar$logCW))[thisy,]
                }else{
                    thiscw=cw
                }
                if(useNMmodel){
                    thisy <- which(rnNM==y)
                    thisnm <-matrix(exp(simLogNM[k,]),nrow=nrow(opar$logNM))[thisy,]
                }else{
                    thisnm=nm
                }
                one<-function(s){
                    simtmp[k,]<<-scaleF(sim[k,],s)  
                    simcat<-catch(simtmp[k,], nm=thisnm, cw=thiscw)
                    return(catchval.exact[i+1]-simcat)
                }
                ff <- uniroot(one, c(0,100))$root
            }
            dd <- sapply(1:nrow(sim),funk)
            sim <- simtmp
        }
        
        if(!is.na(landval[i+1])){
            simtmp<-NA
            fun<-function(s){
                simtmp<<-t(apply(sim, 1, scaleF, scale=s))
                if(useNMmodel){
                    simcat<-sapply(1:nrow(simtmp), function(k){      
                        thisy <- which(rnNM==y);
                        thisnm <-matrix(exp(simLogNM[k,]),nrow=nrow(opar$logNM))[thisy,];
                        catchFrac(simtmp[k,], nm=thisnm, w=lw, frac=lf)
                    }
                    )
                }else{
                    simcat<-apply(simtmp, 1, catchFrac, nm=nm, w=lw, frac=lf)
                }
                return(landval[i+1]-estimate(simcat))
            }
            ff <- uniroot(fun, c(0,100))$root
            sim <- simtmp
        }

        if(!is.na(nextssb[i+1])){
            if((sum(pm)!=0) | (sum(pf)!=0))stop("nextssb option only available if SSB is calculated in the beginning of the year")  
            simtmp<-NA
            rs<-.Random.seed
            fun<-function(s){
                assign(".Random.seed", rs, envir = .GlobalEnv)
                simtmp<<-t(apply(sim, 1, scaleF, scale=s))

                if(useNMmodel){
                    simsim<-t(sapply(1:nrow(simtmp),
                                     function(s){
                                         yidx<-which(rnNM==as.integer(y))   
                                         thisnm<-matrix(exp(simLogNM[s,]),nrow=nrow(opar$logNM))[yidx,]
                                         step(simtmp[s,], nm=thisnm, recpool=recpool, scale=1, inyear=(i==0))
                                     }
                                     ))      
                }else{
                    simsim <- t(apply(simtmp, 1, function(s)step(s, nm=nm, recpool=recpool, scale=1, inyear=(i==0)))) 
                }
                if(i!=0){
                    if(deterministic)procVar<-procVar*0  
                    simsim <- simsim + rmvnorm(nosim, mu=rep(0,nrow(procVar)), Sigma=procVar)
                }
                if(useSWmodel | useMOmodel | useNMmodel){
                    ssbswmoi<-function(k){
                        if(useSWmodel){    
                            thisy <- which(rnSW==y)
                            thissw <-matrix(exp(simLogSw[k,]),nrow=nrow(opar$logSW))[thisy,]
                        }else{
                            thissw <-sw
                        }
                        if(useMOmodel){    
                            thisy <- which(rnMO==y)
                            thismo <-matrix(plogis(simLogitMO[k,]),nrow=nrow(opar$logitMO))[thisy,]
                        }else{
                            thismo <-mo
                        }
                        if(useNMmodel){    
                            thisy <- which(rnNM==y)
                            thisnm <-matrix(exp(simLogNM[k,]),nrow=nrow(opar$logNM))[thisy,]
                        }else{
                            thisnm <-nm
                        }
                        ssb(simsim[k,],nm=thisnm,sw=thissw,mo=thismo,pm=pm,pf=pf)
                    }  
                    simnextssb <- sapply(1:nrow(simsim), ssbswmoi)
                }else{
                    simnextssb <- apply(simsim, 1, ssb, nm=nm, sw=sw, mo=mo, pm=pm, pf=pf)
                }
                return(nextssb[i+1]-estimate(simnextssb))
            }
            ff <- uniroot(fun, c(0,100))$root
            sim <- simtmp
        }
        ## ======= !!!From multi -- NOT FULLY MERGED!!!
        if(any(!is.na(cf.cv.keep.cv[i+1,]))){
            cfcv <- cf.cv.keep.cv[i+1,]  
            cf <- cfcv[1:noCatchFleets]
            cv <- cfcv[1:noCatchFleets+noCatchFleets]
            cfcvtcv<-c(cfcv,catchval[i+1])
            ii <- which(apply(rbind(cv,cf),2,function(x)any(!is.na(x))))
            theta <- rep(1,length(ii)+1)
            if(length(theta)>noCatchFleets)stop("Over-specified in cf.cv.keep.cv")
            lsfun <- function(th_in) {
                th<-exp(th_in)
                s <- rep(NA, noCatchFleets)
                s[ii] <- th[1:length(ii)]
                s[-ii] <- th[length(ii) + 1]
                simtmp <<- t(apply(sim, 1, scaleFbyFleet, scale = s))
                cvfun <- function(x) {
                    tcv <- catch(x, nm = nm, cw = cw)
                    cv <- attr(tcv, "byFleet")
                    return(cv)
                }
                simcat <- apply(simtmp, 1, cvfun)
                medcv <- apply(simcat, 1, estimate)
                med <- c(medcv/sum(medcv), medcv, estimate(apply(simcat, 2, sum)))
                return(-sum(dnorm(log(cfcvtcv), log(med), sd=1 ,log=TRUE), na.rm = TRUE))
            }
            ff <- nlminb(log(theta), lsfun)
            sim <- simtmp
        }

        if(any(!is.na(cf.cv.keep.fv[i+1,]))){
            cfcv <- cf.cv.keep.fv[i+1,]  
            cf <- cfcv[1:noCatchFleets]
            cv <- cfcv[1:noCatchFleets+noCatchFleets]
            cfcvtfv<-c(cfcv,fval[i+1])
            ii <- which(apply(rbind(cv,cf),2,function(x)any(!is.na(x))))
            theta <- rep(1,length(ii)+1)
            if(length(theta)>noCatchFleets)stop("Over-specified in cf.cv.keep.fv")
            lsfun <- function(th_in){
                th <- th_in
                s <- rep(NA,noCatchFleets)
                s[ii] <- th[1:length(ii)]
                s[-ii] <- th[length(ii)+1]
                simtmp <<- t(apply(sim, 1, scaleFbyFleet, scale=s))
                cvfun <- function(x){
                    tcv <- catch(x,nm=nm,cw=cw)
                    cv <- attr(tcv,"byFleet")
                    return(cv)
                }
                simcat <- apply(simtmp, 1, cvfun)
                simfbar <- apply(simtmp, 1, fbar)
                medcv <- apply(simcat,1,estimate)
                med <- c((medcv-cf.keep.fv.offset[i+1,])/sum(medcv),medcv,estimate(simfbar))
                return(sum(((cfcvtfv-med)/cfcvtfv)^2, na.rm=TRUE))
            }
            ff <- nlminb(theta,lsfun, lower=0.001, upper=1000)
            ## >>>>>>> multi  !!!End of conflict block 1 -- NOT FULLY MERGED!!!
            sim <- simtmp
        }
        
        fbarsim <- apply(sim, 1, fbar)
        ## <<<<<<< HEAD  !!!Conflict block 2; From master -- NOT FULLY MERGED!!!
        fbarLsim <- apply(sim, 1, fbarFrac, lf=lf)
        if(useCWmodel | useNMmodel){
            catchcwi<-function(k){
                if(useCWmodel){      
                    thisy <- which(rnCW==y)
                    thiscw <-matrix(exp(simLogCW[k,]),nrow=nrow(opar$logCW))[thisy,]
                }else{
                    thiscw <- cw
                }
                if(useNMmodel){      
                    thisy <- which(rnNM==y)
                    thisnm <-matrix(exp(simLogNM[k,]),nrow=nrow(opar$logNM))[thisy,]
                }else{
                    thisnm <- nm
                }
                catch(sim[k,], nm=thisnm, cw=thiscw)
            }
            catchsim <- sapply(1:nrow(sim), catchcwi)
        }else{
            catchsim <- apply(sim, 1, catch, nm=nm, cw=cw)
        }

        if(useNMmodel){
            landsim<-sapply(1:nrow(sim), function(k){      
                thisy <- which(rnNM==y);
                thisnm <-matrix(exp(simLogNM[k,]),nrow=nrow(opar$logNM))[thisy,];
                catchFrac(sim[k,], nm=thisnm, w=lw, frac=lf)
            }
            )
        }else{
            landsim<-apply(sim, 1, catchFrac, nm=nm, w=lw, frac=lf)
        }
        catchatagesim <- apply(sim, 1, catchatage, nm=nm)
        if(useSWmodel | useMOmodel | useNMmodel){
            ssbswmoi<-function(k){
                if(useSWmodel){    
                    thisy <- which(rnSW==y)
                    thissw <-matrix(exp(simLogSw[k,]),nrow=nrow(opar$logSW))[thisy,]
                }else{
                    thissw <-sw
                }
                if(useMOmodel){    
                    thisy <- which(rnMO==y)
                    thismo <-matrix(plogis(simLogitMO[k,]),nrow=nrow(opar$logitMO))[thisy,]
                }else{
                    thismo <-mo
                }
                if(useNMmodel){    
                    thisy <- which(rnNM==y)
                    thisnm <-matrix(exp(simLogNM[k,]),nrow=nrow(opar$logNM))[thisy,]
                }else{
                    thisnm <-nm
                }
                ssb(sim[k,],nm=thisnm,sw=thissw,mo=thismo,pm=pm,pf=pf)
            }  
            ssbsim <- sapply(1:nrow(sim), ssbswmoi)
        }else{
            ssbsim <- apply(sim, 1, ssb, nm=nm, sw=sw, mo=mo, pm=pm, pf=pf)
        }
        tsbsim <- apply(sim, 1, tsb, sw=sw)
        if(lagR){
            recsim <- exp(sim[,2])
        }else{
            recsim <- exp(sim[,1])
        }
        cwFsim <- rep(NA,nrow(sim))
        if(!missing(customWeights)){
            cwFsim <- apply(sim, 1, getCWF, w=customWeights)
        }
        ## simlist[[i+1]] <- list(sim=sim, fbar=fbarsim, catch=catchsim, ssb=ssbsim, rec=recsim,
        ##                        cwF=cwFsim, catchatage=catchatagesim, land=landsim, fbarL=fbarLsim, tsb=tsbsim, year=y)
        ## =======   !!!From multi -- NOT FULLY MERGED!!!
        fbarbyfleetsim <- apply(sim, 1, fbarbyfleet)
        ## catchsim <- apply(sim, 1, catch, nm=nm, cw=cw)
        ## ssbsim <- apply(sim, 1, ssb, nm=nm, sw=sw, mo=mo, pm=pm, pf=pf)
        ## recsim <- exp(sim[,1])
        catchbysim <- apply(sim, 1, function(x)attr(catch(x, nm=nm, cw=cw), "byFleet"))
        simlist[[i+1]] <- list(sim=sim, fbar=fbarsim, fbarbyfleet = fbarbyfleetsim, catch=catchsim, ssb=ssbsim, rec=recsim, cwF=cwFsim, catchatage=catchatagesim, land=landsim, fbarL=fbarLsim, tsb=tsbsim, catchby=catchbysim,  year=y)
        ## >>>>>>> multi   !!!End of conflict block 2 -- NOT FULLY MERGED!!!
    }
    attr(simlist, "fit")<-fit

    collect <- function(x){
        quan <- quantile(x, c(.50,.025,.975))
        c(Estimate=estimate(x), low=quan[2], high=quan[3])
    }
    
    fbar <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$fbar))),3)
    fbarL <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$fbarL))),3)  
    rec <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$rec))))
    ssb <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$ssb))))
    tsb <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$tsb))))
    catch <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$catch))))
    ## <<<<<<< HEAD   !!!Conflict block 3; From master -- NOT FULLY MERGED!!!
    land <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$land))))  
    caytable<-round(do.call(rbind, lapply(simlist, function(xx)apply(xx$catchatage,1,collect)))) 
    tab <- cbind(fbar,rec,ssb,catch)
    if(splitLD){
        tab<-cbind(tab,fbarL,fbar-fbarL,land,catch-land)
    }
    if(addTSB){
        tab<-cbind(tab,tsb)
    }
    if(!missing(customWeights)) tab <- cbind(tab,cwF=round(do.call(rbind, lapply(simlist, function(xx)collect(xx$cwF))),3))
    rownames(tab) <- unlist(lapply(simlist, function(xx)xx$year))
    nam <- c(estimateLabel,"low","high")
    basename<-c("fbar:","rec:","ssb:","catch:")
    if(splitLD){
        basename<-c(basename,"fbarL:","fbarD:","Land:","Discard:")    
    }
    if(addTSB){
        basename<-c(basename,"tsb:")    
    }
    if(!missing(customWeights)){
        basename<-c(basename,"cwF:")    
    }
    colnames(tab)<-paste0(rep(basename, each=length(nam)), nam)
    attr(simlist, "tab")<-tab
    shorttab<-t(tab[,grep(estimateLabel,colnames(tab))])
    rownames(shorttab)<-sub(paste0(":",estimateLabel),"",paste0(label,if(!is.null(label))":",rownames(shorttab)))
    attr(simlist, "shorttab")<-shorttab
    attr(simlist, "label") <- label
    attr(simlist, "caytable")<-caytable    
    ## =======   !!!From multi -- NOT FULLY MERGED!!!
    if(sum(fit$data$fleetTypes==0)==1){
        catchby <- catch
    }else{
        catchby <- round(do.call(rbind, lapply(simlist,function(xx)as.vector(apply(xx$catchby,1,collect)))))
    }
    if(sum(fit$data$fleetTypes == 0) == 1) {
        fbarby <- fbar
    }
    else{
        fbarby <- round(do.call(rbind, lapply(simlist, function(xx) as.vector(apply(xx$fbarbyfleet, 1, collect)))),3)
    }
    ##   tab <- cbind(fbar, rec, ssb, catch)
    ##   rownames(tab) <- unlist(lapply(simlist, function(xx)xx$year))
      rownames(catchby) <- rownames(tab)
      rownames(fbarby) <- rownames(tab)
      nam <- c(estimateLabel,"low","high")
    ##   colnames(tab) <- paste0(rep(c("fbar:","rec:","ssb:","catch:"), each=length(nam)), nam)
      colnames(catchby) <- paste0(rep(paste0("F",1:sum(fit$data$fleetTypes==0),":"), each=length(nam)), nam)
      colnames(fbarby) <- paste0(rep(paste0("F", 1:sum(fit$data$fleetTypes == 0), ":"), each = length(nam)), nam)
    ##   attr(simlist, "tab") <- tab
      attr(simlist, "catchby") <- catchby
      attr(simlist, "fbarby") <- fbarby
    ##   shorttab <- t(tab[,grep("Estimate",colnames(tab))])
    ##   rownames(shorttab) <- sub(":Estimate","",paste0(label,if(!is.null(label))":",rownames(shorttab)))
    ##   attr(simlist, "shorttab") <- shorttab
    ##   attr(simlist, "label") <- label  
    ## >>>>>>> multi   !!!End of conflict block 3 -- NOT FULLY MERGED!!!
    class(simlist) <- "samforecast"
    attr(simlist,"estimateLabel") <- estimateLabel
    if(!savesim){  
        simlistsmall<-lapply(simlist, function(x)list(year=x$year))
        attributes(simlistsmall)<-attributes(simlist)
        return(simlistsmall)
    }else{
        return(simlist)
    }
}

