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
##' @return a list containing the elements needed to configure a sam model (e.g. minAge, maxAge, maxAgePlusGroup, keyLogFsta, ...). 
##' @details The configuration returned by defcon is intended as a help to set up a syntactically correct configuration for the sam model. The dimensions are set from the data (years, age-classes, and fleet types available). The configuration is intended to be fairly simplistic in the hope that the model configured will at least converge (not guaranteed). Most importantly: No model validation has been performed, so it should not be assumed that the returned model configuration will result in a sensible assessment of the stock. The actual model configuration is the responsibility of the user.     
##' @export
defcon<-function(dat){
  fleetTypes <- dat$fleetTypes
  ages <- do.call(rbind,tapply(dat$aux[,3], INDEX=dat$aux[,2], FUN=range))
  ages[fleetTypes%in%c(3,5),] <- NA
  minAge <- min(ages, na.rm=TRUE)
  maxAge <- max(ages, na.rm=TRUE)
  ages[is.na(ages)] <- minAge
  nAges <- maxAge-minAge+1
  nFleets <- nrow(ages)
  ret <- list()
  ret$minAge <- minAge
  ret$maxAge <- maxAge  
  ret$maxAgePlusGroup <- as.integer(ages[,2]==max(ages[,2], na.rm=TRUE))
  x <- matrix(0, nrow=nFleets, ncol=nAges)
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
  ret$corFlag <- 2
  x <- matrix(0, nrow=nFleets, ncol=nAges)
  lastMax <- 0
  for(i in 1:nrow(x)){
    if(fleetTypes[i]%in%c(1,2,3)){
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
  x <- matrix(0, nrow=nFleets, ncol=nAges)
  lastMax <- 0
  for(i in 1:nrow(x)){
    if(fleetTypes[i]%in%c(0,1,2,3)){
      x[i,(ages[i,1]-minAge+1):(ages[i,2]-minAge+1)] <- lastMax+1
      lastMax <- max(x)
    }
  }  
  ret$keyVarObs <- x - 1
  ret$obsCorStruct <- factor(rep("ID",nFleets),levels=c("ID","AR","US"))
  ret$keyCorObs <- matrix(-1, nrow=nFleets, ncol=nAges-1)
  colnames(ret$keyCorObs)<-paste(minAge:(maxAge-1),(minAge+1):maxAge,sep="-")
  for(i in 1:nrow(x)){
      if(ages[i,1]<ages[i,2]){
        ret$keyCorObs[i,(ages[i,1]-minAge+1):(ages[i,2]-minAge)]<-NA
      }
  }
    
  ret$stockRecruitmentModelCode <- 0
  ret$noScaledYears <- 0
  ret$keyScaledYears <- numeric(0)
  ret$keyParScaledYA <- array(0,c(0,0))

  cs <- colSums(dat$catchMeanWeight)
  ii <- min(which(dat$fleetTypes==0))
  tc <- tapply(dat$logobs[dat$aux[,2]==ii], INDEX=dat$aux[,3][dat$aux[,2]==ii], function(x)sum(x,na.rm=TRUE))
  tc <- tc*cs[names(cs)%in%names(tc)]
  pp <- tc/sum(tc)
  ret$fbarRange <- c(min(which(cumsum(pp)>=0.25)), length(pp)-min(which(cumsum(rev(pp))>=0.25))+1)+(minAge-1)
  ret$keyBiomassTreat <- ifelse(dat$fleetTypes==3, 0, -1)
  ret$obsLikelihoodFlag <- factor(rep("LN",nFleets),levels=c("LN","ALN"))
  ret$fixVarToWeight <- 0
  ret$fracMixF <- 0
  ret$fracMixN <- 0
  ret$fracMixObs <- rep(0,nFleets)
  ret$constRecBreaks <- numeric(0)
  ret$meanVarObsLink <- matrix(-1, nrow=nFleets, ncol=maxAge -minAge+1)
  ret$meanVarFprocLink <- matrix(-1, nrow=nFleets, ncol=maxAge -minAge+1)
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
    cat("# Numbers (integers) starts from zero and must be consecutive\n#" , file=file, append=TRUE)
    #

    txt<-list()
    txt$minAge <- "The minimium age class in the assessment"
    txt$maxAge <- "The maximum age class in the assessment"
    txt$maxAgePlusGroup <- "Is last age group considered a plus group for each fleet (1 yes, or 0 no)." 
    txt$keyLogFsta <- "Coupling of the fishing mortality states (nomally only first row is used)."
    txt$corFlag <- "Correlation of fishing mortality across ages (0 independent, 1 compound symmetry, 2 AR(1), 3 separable AR(1)."
    txt$keyLogFpar <- "Coupling of the survey catchability parameters (nomally first row is not used, as that is covered by fishing mortality)."
    txt$keyQpow <- "Density dependent catchability power parameters (if any)."
    txt$keyVarF <- "Coupling of process variance parameters for log(F)-process (nomally only first row is used)"
    txt$keyVarLogN <- "Coupling of process variance parameters for log(N)-process"
    txt$keyVarObs <- "Coupling of the variance parameters for the observations."
    txt$obsCorStruct <- "Covariance structure for each fleet (\"ID\" independent, \"AR\" AR(1), or \"US\" for unstructured)."
    txt$keyCorObs <- paste0("Coupling of correlation parameters can only be specified if the AR(1) structure is chosen above.",
                            "\n# NA's indicate where correlation parameters can be specified (-1 where they cannot).",
                            paste0("\n#",paste0(colnames(x$keyCorObs), collapse=" ")))
    txt$stockRecruitmentModelCode <- "Stock recruitment code (0 for plain random walk, 1 for Ricker, 2 for Beverton-Holt, and 3 piece-wise constant)."
    txt$noScaledYears <- "Number of years where catch scaling is applied."
    txt$keyScaledYears <- "A vector of the years where catch scaling is applied."
    txt$keyParScaledYA <- "A matrix specifying the couplings of scale parameters (nrow = no scaled years, ncols = no ages)."
    txt$fbarRange <- "lowest and higest age included in Fbar"
    txt$keyBiomassTreat <- "To be defined only if a biomass survey is used (0 SSB index, 1 catch index, 2 FSB index, 3 total catch, 4 total landings and 5 TSB index)."
    txt$obsLikelihoodFlag <- "Option for observational likelihood"
    txt$fixVarToWeight <- "If weight attribute is supplied for observations this option sets the treatment (0 relative weight, 1 fix variance to weight)."
    txt$fracMixF <- "The fraction of t(3) distribution used in logF increment distribution" 
    txt$fracMixN <- "The fraction of t(3) distribution used in logN increment distribution"
    txt$fracMixObs <- "A vector with same length as number of fleets, where each element is the fraction of t(3) distribution used in the distribution of that fleet"
    txt$constRecBreaks <- "Vector of break years between which recruitment is at constant level. The break year is included in the left interval. (This option is only used in combination with stock-recruitment code 3)"
    txt$meanVarObsLink <- "Coupling of parameters used in a mean-variance link for observations."
    txt$meanVarFprocLink <- "Coupling of parameters used in a mean-variance link for log(F) process increments"
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
    idx2<-min(keyIdx[keyIdx>idx1])-1
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
