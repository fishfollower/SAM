##' small helper function
##' @param min from 
##' @param max to 
##' @details internal
setSeq<-function(min,max){
    if(min==max){
        ret <- 1
    }else{
        ret <- c(1:(max-min),max-min)
    }
    return(ret)
}

##' small helper function
##' @param x vector if indices  
##' @details internal
setS<-function(x){
    setSeq(1,length(x))
}

##' Setup basic minimal configuration for sam assessment
##' @param dat sam data object
##' @param level 1 or 2 (1 most basic configuration, 2 configuration with AR correlation structure on surveys)
##' @return a list containing the elements needed to configure a sam model (e.g. minAge, maxAge, maxAgePlusGroup, keyLogFsta, ...). 
##' @details The configuration returned by defcon is intended as a help to set up a syntactically correct configuration for the sam model. The dimensions are set from the data (years, age-classes, and fleet types available). The configuration is intended to be fairly simplistic in the hope that the model configured will at least converge (not guaranteed). Most importantly: No model validation has been performed, so it should not be assumed that the returned model configuration will result in a sensible assessment of the stock. The actual model configuration is the responsibility of the user.     
##' @export
defcon<-function(dat, level=1){
    fleetTypes <- dat$fleetTypes    
    ##ages <- do.call(rbind,tapply(dat$aux[,3], INDEX=dat$aux[,2], FUN=range))
    ages <- cbind(dat$minAgePerFleet, dat$maxAgePerFleet)    
    ages[fleetTypes%in%c(3,5,6,80,90,92),] <- NA
    minAge <- min(ages, na.rm=TRUE)
    maxAge <- max(ages, na.rm=TRUE)
    ages[is.na(ages)] <- minAge
    nAges <- maxAge-minAge+1
    nFleets <- nrow(ages)
    ret <- list()
    ret$minAge <- minAge
    ret$maxAge <- maxAge  
    ret$maxAgePlusGroup <- as.integer(ages[,2]==max(ages[,2], na.rm=TRUE))    
    x <- xOff <- matrix(0, nrow=nFleets, ncol=nAges)
    lastMax <- 0
    for(i in 1:nrow(x)){
        if(fleetTypes[i]==0){
            aa <- ages[i,1]:ages[i,2]
            aa<-aa[tapply(dat$logobs[dat$aux[,2]==i], INDEX=dat$aux[,3][dat$aux[,2]==i], function(x)!all(is.na(x)))]
            x[i,aa-minAge+1] <- setS(aa)+lastMax
            lastMax <- max(x)
        }
    }  
    ret$keyLogFsta <- x - 1
    ret$keyLogFmu <- xOff - 1
    ret$keyLogFrho <- xOff - 1

    ret$corFlag <- rep(2,length(which(fleetTypes==0)))
    x <- matrix(0, nrow=nFleets, ncol=nAges)
    lastMax <- 0
    for(i in 1:nrow(x)){
        if(fleetTypes[i]%in%c(1,2,3,6)){
            x[i,(ages[i,1]-minAge+1):(ages[i,2]-minAge+1)] <- setSeq(ages[i,1],ages[i,2])+lastMax
            lastMax <- max(x)
        }
    }  
    ret$keyLogFpar <- x - 1 
    ret$keyQpow <- matrix(-1, nrow=nFleets, ncol=nAges)
    x<-matrix(0, nrow=nFleets, ncol=nAges)
    lastMax <- 0
    for(i in 1:nrow(x)){
        if(fleetTypes[i]==0){
            x[i,(ages[i,1]-minAge+1):(ages[i,2]-minAge+1)]<-lastMax+1
            lastMax <- max(x)
        }
    }  
    ret$keyVarF <- x - 1

    ret$keyVarLogN <- c(1,rep(2,nAges-1)) - 1 

    ret$keyVarLogP <- numeric(0)
    if(any(fleetTypes==6))
        ret$keyVarLogP<- seq(1,length(which(fleetTypes==6))-1)-1

    x <- matrix(0, nrow=nFleets, ncol=nAges)
    lastMax <- 0
    for(i in 1:nrow(x)){
        if(fleetTypes[i]%in%c(0,1,2,3,6,81)){
            x[i,(ages[i,1]-minAge+1):(ages[i,2]-minAge+1)] <- lastMax+1
            lastMax <- max(x)
        }else if(fleetTypes[i]%in%c(80,90,92)){
            x[i,] <- 0
            x[i,1] <- lastMax+1
            lastMax <- max(x)
        }
                                        #    if(fleetTypes[i]==6 & i == max(which(fleetTypes==6)))
                                        #      x[i,(ages[i,1]-minAge+1):(ages[i,2]-minAge+1)] <- 0
    }  
    ret$keyVarObs <- x - 1

    ret$obsCorStruct <- factor(rep("ID",nFleets),levels=c("ID","AR","US"))
    if(level==2){
        ret$obsCorStruct[fleetTypes==2]<-"AR"
    }
    ret$obsCorStruct[fleetTypes==7] <- NA
    ret$keyCorObs <- matrix(-1, nrow=nFleets, ncol=nAges-1)
    colnames(ret$keyCorObs)<-paste(minAge:(maxAge-1),(minAge+1):maxAge,sep="-")
    nextpar<-0
    for(i in 1:nrow(x)){
        if(fleetTypes[i]!=7){
            if(ages[i,1]<ages[i,2]){
                if((level==2)&(fleetTypes[i]==2)){
                    ret$keyCorObs[i,(ages[i,1]-minAge+1):(ages[i,2]-minAge)]<-nextpar
                    nextpar<-nextpar+1
                }else{
                    ret$keyCorObs[i,(ages[i,1]-minAge+1):(ages[i,2]-minAge)]<-NA
                }           
            }
        }
    }
    
    ret$stockRecruitmentModelCode <- 0
    ret$noScaledYears <- 0
    ret$keyScaledYears <- numeric(0)
    ret$keyParScaledYA <- array(0,c(0,0))

    cs <- colSums(apply(dat$catchMeanWeight,1:2,median,na.rm=T), na.rm = TRUE)
    ii <- which(dat$fleetTypes==0)
    tc <- tapply(dat$logobs[dat$aux[,2]%in%ii], INDEX=dat$aux[,3][dat$aux[,2]%in%ii], function(x)sum(x,na.rm=TRUE))
    tc <- tc*cs[names(cs)%in%names(tc)]
    pp <- tc/sum(tc)
    ret$fbarRange <- c(min(which(cumsum(pp)>=0.25)), length(pp)-min(which(cumsum(rev(pp))>=0.25))+1)+(minAge-1)
    ret$fbarRange <- ifelse(is.finite(ret$fbarRange), ret$fbarRange, c(ret$minAge,ret$maxAge))
    ret$keyBiomassTreat <- ifelse(dat$fleetTypes==3, 0, -1)
    ret$obsLikelihoodFlag <- factor(rep("LN",nFleets),levels=c("LN","ALN"))
    ret$fixVarToWeight <- rep(0,nFleets)
    ret$fracMixF <- 0
    ret$fracMixN <- rep(0,nAges)
    ret$fracMixObs <- rep(0,nFleets)
    ret$constRecBreaks <- numeric(0)

    ret$predVarObsLink <- matrix(NA, nrow=nFleets, ncol=maxAge -minAge+1)
    for(i in 1:nrow(x)){
        if(ages[i,1]<ages[i,2]){
            ret$predVarObsLink[i,(ages[i,1]-minAge+1):(ages[i,2]-minAge+1)]<- -1
        }
    }
    
    ## ret$hockeyStickCurve <- 20
    ret$stockWeightModel <- 0
    ret$keyStockWeightMean <- rep(NA_integer_,nAges)
    ret$keyStockWeightObsVar <- rep(NA_integer_,nAges)  
    ret$catchWeightModel <- 0
    ret$keyCatchWeightMean <- matrix(NA_integer_,sum(fleetTypes==0),nAges)
    ret$keyCatchWeightObsVar <- matrix(NA_integer_,sum(fleetTypes==0),nAges)  
    ret$matureModel <- 0
    ret$keyMatureMean <- rep(NA_integer_,nAges)
    ret$mortalityModel <- 0
    ret$keyMortalityMean <- rep(NA_integer_,nAges)
    ret$keyMortalityObsVar <- rep(NA_integer_,nAges)  
    ret$keyXtraSd<-matrix(NA_integer_, nrow=0, ncol=4)
    ret$logNMeanAssumption <- c(0,0)
    ret$initState <- 0
    ret$recruitmentAutocorrelation <- 0
    ret$keyLogFseason <- xOff - 1
    ret$seasonTimes <- c(0,1)
    ret$isFishingSeason <- 1
    ret$seasonFirstYear <- -Inf
    ret$seasonFixedEffect <- 0
    return(ret) 
}

##' Saves a model configuration list to a file  
##' @param x sam configuration list as returned from defcon or loadConf
##' @param file the file to save the configuration to
##' @param overwrite logical if an existing file should be overwritten (FALSE by default) 
##' @details function useful for saving a model configuration. A saved configuration can be read back in via the loadConf function 
##' @export
saveConf <- function(x, file="", overwrite=FALSE){
    writeConf <- function(x,...) UseMethod("writeConf")

    writeConf.default <- function(x,...){
        stop("Unimplemented class in writeConf")
    }

    writeConf.integer <- function(x,...){
        cat("\n",x,"\n",...)
    }

    writeConf.numeric <- function(x,...){
        cat("\n",x,"\n",...)
    }

    writeConf.matrix <- function(x,...){
        if(nrow(x)>0){
            cat(capture.output(prmatrix(x, rowlab=rep("", nrow(x)), collab=rep("   ",ncol(x)))), sep="\n", ...)
        }else{
            cat("\n", ...)
        }
    }

    writeConf.factor <- function(x,...){
        cat(" | Possible values are:", paste0('\"',levels(x),'\"'), ...)
        cat("\n", paste0('\"',x,'\"'),"\n", ...)
    }
    
    if(file.exists(file) & !overwrite){
        cat("Notice: Did not overwrite exsisting file\n")
    }else{
                                        # Intro txt   
        cat(paste0("# Configuration saved: ",date()), file=file)
        cat("\n#\n# Where a matrix is specified rows corresponds to fleets and columns to ages.\n" , file=file, append=TRUE)
        cat("# Same number indicates same parameter used\n" , file=file, append=TRUE)
        cat("# Numbers (integers) starts from zero and must be consecutive\n" , file=file, append=TRUE)
        cat("# Negative numbers indicate that the parameter is not included in the model\n#" , file=file, append=TRUE)
                                        #

        txt<-list()
        txt$minAge <- "The minimium age class in the assessment"
        txt$maxAge <- "The maximum age class in the assessment"
        txt$maxAgePlusGroup <- "Is last age group considered a plus group for each fleet (1 yes, or 0 no)." 
        txt$keyLogFsta <- "Coupling of the fishing mortality states processes for each age (normally only \n# the first row (= fleet) is used). \n# Sequential numbers indicate that the fishing mortality is estimated individually \n# for those ages; if the same number is used for two or more ages, F is bound for \n# those ages (assumed to be the same). Binding fully selected ages will result in a \n# flat selection pattern for those ages."
        txt$corFlag <-"Correlation of fishing mortality across ages (0 independent, 1 compound symmetry, \n# 2 AR(1), 3 separable AR(1). \n# 0: independent means there is no correlation between F across age \n# 1: compound symmetry means that all ages are equally correlated; \n# 2: AR(1) first order autoregressive - similar ages are more highly correlated than \n# ages that are further apart, so similar ages have similar F patterns over time. \n# if the estimated correlation is high, then the F pattern over time for each age \n# varies in a similar way. E.g if almost one, then they are parallel (like a \n# separable model) and if almost zero then they are independent. \n# 3: Separable AR - Included for historic reasons . . .  more later"
        txt$keyLogFpar <- "Coupling of the survey catchability parameters (nomally first row is \n# not used, as that is covered by fishing mortality)."
        txt$keyQpow <- "Density dependent catchability power parameters (if any)."
        txt$keyVarF <- "Coupling of process variance parameters for log(F)-process (Fishing mortality \n# normally applies to the first (fishing) fleet; therefore only first row is used)"
        txt$keyVarLogN <- "Coupling of the recruitment and survival process variance parameters for the \n# log(N)-process at the different ages. It is advisable to have at least the first age \n# class (recruitment) separate, because recruitment is a different process than \n# survival."
        txt$keyVarObs <- "Coupling of the variance parameters for the observations. \n# First row refers to the coupling of the variance parameters for the catch data \n# observations by age \n# Second and further rows refers to coupling of the variance parameters for the \n# index data observations by age"
        txt$obsCorStruct <- "Covariance structure for each fleet (\"ID\" independent, \"AR\" AR(1), or \"US\" for unstructured)."
        txt$keyCorObs <- paste0("Coupling of correlation parameters can only be specified if the AR(1) structure is chosen above.",
                                "\n# NA's indicate where correlation parameters can be specified (-1 where they cannot).",
                                paste0("\n#",paste0(colnames(x$keyCorObs), collapse=" ")))
        txt$stockRecruitmentModelCode <- "Stock recruitment code (0 for plain random walk, 1 for Ricker, 2 for Beverton-Holt, 3 piece-wise constant, 61 for segmented regression/hockey stick, 62 for AR(1), 63 for bent hyperbola / smooth hockey stick, 64 for power function with degree < 1, 65 for power function with degree > 1, 66 for Shepher, 67 for Deriso, 68 for Saila-Lorda, 69 for sigmoidal Beverton-Holt, 90 for CMP spline, 91 for more flexible spline, and 92 for most flexible spline)."
        txt$noScaledYears <- "Number of years where catch scaling is applied."
        txt$keyScaledYears <- "A vector of the years where catch scaling is applied."
        txt$keyParScaledYA <- "A matrix specifying the couplings of scale parameters (nrow = no scaled years, ncols = no ages)."
        txt$fbarRange <- "lowest and higest age included in Fbar"
        txt$keyBiomassTreat <- "To be defined only if a biomass survey is used (0 SSB index, 1 catch index, 2 FSB index, 3 total catch, 4 total landings, 5 TSB index,  6 TSN index, and 10 Fbar idx)."
        txt$obsLikelihoodFlag <- "Option for observational likelihood"
        txt$fixVarToWeight <- "If weight attribute is supplied for observations this option sets the treatment (0 relative weight, 1 fix variance to weight). Can be specified fleetwise."
        txt$fracMixF <- "The fraction of t(3) distribution used in logF increment distribution" 
        txt$fracMixN <- "The fraction of t(3) distribution used in logN increment distribution (for each age group)"
        txt$fracMixObs <- "A vector with same length as number of fleets, where each element is the fraction of t(3) distribution used in the distribution of that fleet"
        txt$constRecBreaks <- "For stock-recruitment code 3: Vector of break years between which recruitment is at constant level. The break year is included in the left interval. For spline stock-recruitment: Vector of log-ssb knots. (This option is only used in combination with stock-recruitment code 3, 90-92, and 290)"
        txt$predVarObsLink <- "Coupling of parameters used in a prediction-variance link for observations."
        txt$stockWeightModel <- "Integer code describing the treatment of stock weights in the model (0 use as known, 1 use as observations to inform stock weight process (GMRF with cohort and within year correlations)), 2 to add extra correlation to plusgroup"
        txt$keyStockWeightMean <- "Coupling of stock-weight process mean parameters (not used if stockWeightModel==0)"
        txt$keyStockWeightObsVar <- "Coupling of stock-weight observation variance parameters (not used if stockWeightModel==0)"
        txt$catchWeightModel <- "Integer code describing the treatment of catch weights in the model (0 use as known, 1 use as observations to inform catch weight process (GMRF with cohort and within year correlations)), 2 to add extra correlation to plusgroup"
        txt$keyCatchWeightMean <- "Coupling of catch-weight process mean parameters (not used if catchWeightModel==0)"
        txt$keyCatchWeightObsVar <- "Coupling of catch-weight observation variance parameters (not used if catchWeightModel==0)"
        txt$matureModel <- "Integer code describing the treatment of proportion mature in the model (0 use as known, 1 use as observations to inform proportion mature process (GMRF with cohort and within year correlations on logit(proportion mature))), 2 to add extra correlation to plusgroup"
        txt$keyMatureMean <- "Coupling of mature process mean parameters (not used if matureModel==0)"
        txt$mortalityModel <- "Integer code describing the treatment of natural mortality in the model (0 use as known, 1 use as observations to inform natural mortality process (GMRF with cohort and within year correlations)), 2 to add extra correlation to plusgroup"
        txt$MortalityMean <- "Coupling of natural mortality process mean parameters (not used if mortalityModel==0)"
        txt$keyMortalityObsVar <- "Coupling of natural mortality observation variance parameters (not used if mortalityModel==0)"
        txt$keyXtraSd<-"An integer matrix with 4 columns (fleet year age coupling), which allows additional uncertainty to be estimated for the specified observations"
        txt$logNMeanAssumption <- "Flags indicating what the population model should correspond to. 0: Median, 1: Mean, 2: Mode. Two values are are given to differentiate recruitment and other ages."
        txt$initState <- "Flag indicating whether initial parameters should be added for the latent processes."
        txt$recruitmentAutocorrelation <- "Number of auto-correlation parameters for recruitment. The auto-regressive process is forced to be stationary with real characteristic roots."
        nam<-names(x)
        dummy<-lapply(1:length(nam), function(i){
            cat('\n$', file=file, append=TRUE)
            cat(nam[i], file=file, append=TRUE)
            cat('\n#', txt[[nam[i]]], file=file, append=TRUE)
            writeConf(x[[i]], file=file, append=TRUE)
        }
        )
    }
}

##' Loads a model configuration from a file  
##' @param dat sam data list as returned from the function setup.sam.data
##' @param file the file to read the configuration from
##' @param patch logical if TRUE missing entries will be automatically filled with default values   
##' @details function useful loading a model configuration. Such a configuration can be saved via the saveConf function
##' @importFrom utils capture.output
##' @export
loadConf <- function(dat, file, patch=TRUE){
    dconf <- defcon(dat)
    confWithName<-lapply(1:length(dconf), function(i){x<-dconf[[i]]; attr(x,"nam")<-names(dconf)[i]; x})
    lin <- c(readLines(file),"$end")
    lin<-lin[-grep("^#",lin)]
    keyIdx <- grep("^\\$",lin)
    getIdx <- function(nam){
        idx1<-grep(paste0("^\\$",nam, "( |$)"),lin)+1
        idx2<-min(keyIdx[keyIdx>(idx1-1)])-1
        ret <- NULL
        if(idx1<=idx2){
            ret <- idx1:idx2
        }
        ret
    }
    readConf <- function(x) UseMethod("readConf")

    readConf.default <- function(x){
        stop("Unimplemented class in readConf")
    }
    readConf.numeric <- function(x){
        nam <- attr(x,"nam")
        scan(textConnection(lin[getIdx(nam)]), quiet=TRUE)
    }
    readConf.integer <- function(x){
        nam <- attr(x,"nam")
        scan(textConnection(lin[getIdx(nam)]), quiet=TRUE)
    }
    readConf.matrix <- function(x){
        nam <- attr(x,"nam")
        x <- try(read.table(text=lin[getIdx(nam)], header=FALSE), silent = TRUE)
        if(inherits(x, "try-error")){
            return(matrix(NA_real_, nrow=0, ncol=0))
        }else{
            return(as.matrix(x))
        }
    }
    readConf.factor <- function(x){
        nam <- attr(x,"nam")
        factor(scan(textConnection(lin[getIdx(nam)]), what="character", quiet="TRUE"), levels=levels(x))
    }
    isOK <- sapply(names(dconf), function(n)length(grep(paste0("^\\$",n, "( |$)"),lin))==1)
    if(!all(isOK) & !patch)stop("The configuration file is not compatible with model version. Consider running with patch=TRUE")
    fun<-function(i){
        if(isOK[i]){
            readConf(confWithName[[i]])
        }else{
            dconf[[i]]
        }
    }
    conf <- lapply(1:length(confWithName), fun)
    names(conf) <- names(dconf)
    conf
}
