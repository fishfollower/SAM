##' Write ICES/CEFAS data file from matrix 
##' @param x a matrix where rownames are taken as years and colnames are taken as ages
##' @param fileout file name or connection
##' @param ... Arguments to be passed to write 
##' @details
##' 
##' Takes the data and writes them in the ICES/CEFAS format. It is assumed that rows represent consecutive years and cols consecutive ages 
##' 
##' @export
write.ices <- function(x, fileout, ...){
    top <- paste0(fileout, " auto written\n1 2\n", paste0(range(as.integer(rownames(x))), collapse="  "),"\n",paste0(range(as.integer(colnames(x))), collapse="  "), "\n1")    
    write(top,fileout,...)
    write(t(x),fileout,ncolumns=ncol(x),append=TRUE,sep="  \t",...)
}

##' Extract a fleet from a fitted object 
##' @param fit A fitted object as returned from sam.fit
##' @param fleet The number of the fleet 
getFleet <- function(fit, fleet){
  fidx <- fit$data$aux[,"fleet"]==fleet
  aux <- fit$data$aux[fidx,]
  logobs <- fit$data$logobs[fidx ]
  .goget <- function(y, a) {
    ret <- exp(logobs[aux[, "year"] == y & aux[, "age"] == a])
    ifelse(length(ret) == 0, 0, ret)
  }
  yr <- min(aux[,"year"]):max(aux[,"year"])
  ar <- min(aux[,"age"]):max(aux[,"age"])
  tmp <- outer(yr, ar, Vectorize(.goget))
  dimnames(tmp)[[1]] <- yr
  dimnames(tmp)[[2]] <- ar
  return(tmp)
}

##' Write surveys in ICES/CEFAS data file from a model object  
##' @param fit A fitted object as returned from sam.fit
##' @param fileout file name or connection
##' @param ... Arguments to be passed to write 
##' @details
##' 
##' Takes the survey data from the fitted object and writes them in the ICES/CEFAS format. 
##' 
##' @export
write.surveys <- function(fit,fileout,...){
  sidx <- which(fit$data$fleetTypes%in%c(2,3,4))
  top <- paste0(fileout," auto written\n", 100+length(sidx))
  write(top,fileout,...)
  for(s in sidx){
    write(paste0(attr(fit$data, "fleetNames")[s]), fileout, append=TRUE, ...)
    S <- getFleet(fit,s)
    yr <- range(as.integer(rownames(S)))
    ar <- range(as.integer(colnames(S)))
    write(paste0(yr[1], " ", yr[2]), fileout, append=TRUE, ...)
    st <- fit$data$sampleTimes[s]
    write(paste0(1," ", 1, " ", st, " ", st), fileout, append=TRUE, ...)
    write(paste0(ar[1], " ", ar[2]), fileout, append=TRUE, ...)
    x <- cbind(1,S)
    write(t(x),fileout,ncolumns=ncol(x),append=TRUE,sep="  \t",...)
  }
}


##' Write all data files from a model object in the usual ICES/CEFAS data files  
##' @param fit A fitted object as returned from sam.fit
##' @param dir Directory where the files are written  
##' @details
##' 
##' Write all data files from a model object in the usual ICES/CEFAS data files  
##' 
##' @export
write.data.files<-function(fit, dir="."){
  od <- setwd(dir)
  write.ices(fit$data$catchMeanWeight, "cw.dat")
  write.ices(fit$data$disMeanWeight, "dw.dat")
  write.ices(fit$data$landMeanWeight, "lw.dat")
  write.ices(fit$data$landFrac, "lf.dat")  
  write.ices(fit$data$propMat, "mo.dat")    
  write.ices(fit$data$stockMeanWeight, "sw.dat")
  write.ices(fit$data$propF, "pf.dat")
  write.ices(fit$data$propM, "pm.dat")
  write.ices(fit$data$natMor, "nm.dat")
  write.ices(getFleet(fit,1), "cn.dat")
  write.surveys(fit, "survey.dat")
  setwd(od)
}
