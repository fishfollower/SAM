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

##' Extract a fleet observed or predicted value from a fitted object 
##' @param fit A fitted object as returned from sam.fit
##' @param fleet The fleet number
##' @param pred Should it be predicted value, default is observed
##' @details Extract for example the observed or predicted catch at age of fleet "fleet"
##' @return A matrix of observed or predicted values for fleet "fleet"
##' @export
getFleet <- function(fit, fleet, pred="FALSE"){
  fidx <- fit$data$aux[,"fleet"]==fleet
  aux <- fit$data$aux[fidx,]
  if(pred) logout <- fit$rep$predObs[fidx] else logout <-fit$data$logobs[fidx]
  .goget <- function(y, a) {
    ret <- exp(logout[aux[, "year"] == y & aux[, "age"] == a])
    ifelse(length(ret) == 0, 0, ret)
  }
  yr <- min(aux[,"year"]):max(aux[,"year"])
  ar <- min(aux[,"age"]):max(aux[,"age"])
  tmp <- outer(yr, ar, Vectorize(.goget))
  dimnames(tmp)[[1]] <- yr
  dimnames(tmp)[[2]] <- ar
  return(tmp)
}

##' Extract a list of catch fleets
##'
##' @param fit A fitted object as returned from sam.fit
##' @param pred Should it be predicted value, default is observed
##' @return A list of matrices of observed or predicted values for catch fleets
##' @export
getResidualFleets <- function(fit, pred="FALSE"){
    ii <- which(fit$data$fleetTypes == 0)
    r <- lapply(ii, getFleet, fit=fit, pred=pred)
    names(r) <- attr(fit$data,"fleetNames")[ii]
    return(r)
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

##' Write all data files from a list as created by 'setup.sam.data'  
##' @param dat A list as created by 'setup.sam.data'
##' @param dir Directory where the files are written  
##' @details
##' 
##' Write all data files from a list as created by 'setup.sam.data'
##' 
##' @export
write.data.files<-function(dat, dir="."){
  od <- setwd(dir)
  write.ices(dat$catchMeanWeight, "cw.dat")
  write.ices(dat$disMeanWeight, "dw.dat")
  write.ices(dat$landMeanWeight, "lw.dat")
  write.ices(dat$landFrac, "lf.dat")  
  write.ices(dat$propMat, "mo.dat")    
  write.ices(dat$stockMeanWeight, "sw.dat")
  write.ices(dat$propF, "pf.dat")
  write.ices(dat$propM, "pm.dat")
  write.ices(dat$natMor, "nm.dat")
  fit <- list(data=dat)
  cn <- getResidualFleets(fit)
  if(length(cn) == 1){
      write.ices(cn[[1]], "cn.dat")
  }else{
      sapply(seq_along(cn), function(i) write.ices(cn[[i]], sprintf("cn_%s.dat",LETTERS[i])))
  }
  write.surveys(fit, "survey.dat")
  setwd(od)
}

##' Read all standard data SAM files and return a list as created by 'setup.sam.data'
##' @param dir Directory to read from
##' @return list (as created by 'setup.sam.data')
##' @details
##'
##' Read all standard SAM data files
##'
##' @export
read.data.files<-function(dir="."){
    od <- setwd(dir); on.exit(setwd(od));

    cn<-read.ices("cn.dat")
    cw<-read.ices("cw.dat")
    dw<-read.ices("dw.dat")
    lw<-read.ices("lw.dat")
    mo<-read.ices("mo.dat")
    nm<-read.ices("nm.dat")
    pf<-read.ices("pf.dat")
    pm<-read.ices("pm.dat")
    sw<-read.ices("sw.dat")
    lf<-read.ices("lf.dat")
    surveys<-read.ices("survey.dat")

    dat<-setup.sam.data(surveys=surveys,
                    residual.fleets=cn,
                    prop.mature=mo,
                    stock.mean.weight=sw,
                    catch.mean.weight=cw,
                    dis.mean.weight=dw,
                    land.mean.weight=lw,
                    prop.f=pf,
                    prop.m=pm,
                    natural.mortality=nm,
                    land.frac=lf)
    dat
}
