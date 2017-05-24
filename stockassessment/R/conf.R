##' small helper function
##' @param min from 
##' @param max to 
##' @details ...
setSeq<-function(min,max){
  if(min==max){
    ret <- 1
  }else{
    ret <- c(1:(max-min),max-min)
  }
  return(ret)
}

##' Setup basic minimal configuration for sam assessment
##' @param dat sam data object 
##' @return a list containing 
##' \item{minAge}{}
##' \item{maxAge}{}
##' \item{maxAgePlusGroup}{}
##' \item{keyLogFsta}{}
##' \item{keyLogFpar}{}
##' \item{keyQpow}{}
##' \item{keyVarF}{}
##' \item{keyVarLogN}{}
##' \item{keyVarObs}{}
##' \item{obsCorStruct}{}
##' \item{keyCorObs}{}
##' \item{corFlag}{set to zero as a placeholder here.}
##' \item{stockRecruitmentModelCode}{set to zero as a placeholder here.}
##' \item{noScaledYears}{set to zero as a placeholder here.}
##' \item{keyScaledYears}{a scalar set to zero as a placeholder here.}
##' \item{keyParScaledYA}{an array set to zero as a placeholder here.}
##' \item{fbarRange}{}
##' \item{obsLikelihoodFlag}{A vector with an element for each fleet. "LN" means, and "ALN" means}
##' @details ...
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
  ret$maxAgePlusGroup <- 1
  x <- matrix(0, nrow=nFleets, ncol=nAges)
  lastMax <- 0
  for(i in 1:nrow(x)){
    if(fleetTypes[i]==0){
      x[i,(ages[i,1]-minAge+1):(ages[i,2]-minAge+1)] <- setSeq(ages[i,1],ages[i,2])+lastMax
      lastMax <- max(x)
    }
  }  
  ret$keyLogFsta <- x - 1
  ret$corFlag <- 0
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
  ret$fbarRange <- ages[fleetTypes==0,]
  ret$keyBiomassTreat <- ifelse(dat$fleetTypes==3, 0, -1)
  ret$obsLikelihoodFlag <- factor(rep("LN",nFleets),levels=c("LN","ALN"))
  ret$fixVarToWeight <- 0
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
      cat(capture.output(prmatrix(x, rowlab=rep("", nrow(x)), collab=rep("",ncol(x)))), sep="\n", ...)
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
    cat(paste0("# Configuration saved: ",date()), file=file)  
    nam<-names(x)
    dummy<-lapply(1:length(nam), function(i){
        cat('\n$', file=file, append=TRUE)
        cat(nam[i], file=file, append=TRUE)
        writeConf(x[[i]], file=file, append=TRUE)
      }
    )
  }
}

##' Loads a model configuration from a file  
##' @param dat sam data list as returned from the function setup.sam.data
##' @param file the file to read the configuration from
##' @details function useful loading a model configuration. Such a configuration can be saved via the saveConf function
##' @importFrom utils capture.output
##' @export
loadConf <- function(dat, file){
  dconf <- defcon(dat)
  confWithName<-lapply(1:length(dconf), function(i){x<-dconf[[i]]; attr(x,"nam")<-names(dconf)[i]; x})
  lin <- c(readLines(file),"$end")
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
  #lapply(confWithName, readConf) does not work?
  conf <- lapply(1:length(confWithName), function(i)readConf(confWithName[[i]]))
  names(conf) <- names(dconf)
  conf
}
