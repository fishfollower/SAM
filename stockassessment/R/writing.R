##' Write ICES/CEFAS data file from matrix 
##' @param x a matrix where rownames are taken as years and colnames are taken as ages
##' @param fileout file name or connection
##' @param writeToOne Write multi fleet data to one file if data is equal for all fleets
##' @param ... Arguments to be passed to write 
##' @details
##' 
##' Takes the data and writes them in the ICES/CEFAS format. It is assumed that rows represent consecutive years and cols consecutive ages 
##' 
##' @export
write.ices <- function(x, fileout, writeToOne = TRUE, ...){
    writeOne <- function(x, fileout, ...){
        top <- paste0(fileout, " auto written\n1 2\n", paste0(range(as.integer(rownames(x))), collapse="  "),"\n",paste0(range(as.integer(colnames(x))), collapse="  "), "\n1")    
        write(top,fileout,...)
        write(t(x),fileout,ncolumns=ncol(x),append=TRUE,sep="  \t",...)
        if(any(!(names(attributes(x)) %in% c("dim","dimnames")))){
            a <- attributes(x)
            customA <- a[!(names(attributes(x)) %in% c("dim","dimnames"))]
            invisible(sapply(seq_along(customA), function(i){
                write(paste0(sprintf("#'@%s\n", names(customA)[i]),paste(paste("##'",deparse(customA[[i]])),collapse="\n")), fileout, append = TRUE,...)
            }))
        }
    }
    if(is.matrix(x)){
        writeOne(x,fileout, ...)
        return(NULL);
    }

    if(is.array(x) && length(dim(x))==3){
        ## Convert to list of matrices
        x <- lapply(split(x,rep(seq_len(dim(x)[3]), each = prod(dim(x)[1:2]))), matrix, nrow=dim(x)[1], ncol = dim(x)[2], dimnames = dimnames(x)[1:2])
    }

    if(!is.list(x))
        stop("x must be a matrix, and array or a list of matrices")

    if(writeToOne && all(sapply(x, function(y) isTRUE(all.equal(y,x[[1]])))))
        x <- x[1]
    
    if(length(x) == 1){
        writeOne(x[[1]], fileout)
    }else{
        f2 <- sub("\\.dat","_%05d.dat",fileout)
        sapply(seq_along(x), function(i) writeOne(x[[i]],sprintf(f2,i)))
    }    
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

getResidualSumFleets <- function(fit, pred="FALSE"){
    ii <- which(fit$data$fleetTypes == 7)
    r <- lapply(ii, getFleet, fit=fit, pred=pred)
    names(r) <- attr(fit$data,"fleetNames")[ii]
    for(i in seq_along(ii))
        attr(r[[i]],"sumof") <- which(fit$data$sumKey[ii[i],]>0)
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
##' @param writeToOne Write multi fleet data to one file if data is equal for all fleets
##' @details
##' 
##' Write all data files from a list as created by 'setup.sam.data'
##' 
##' @export
write.data.files<-function(dat, dir=".", writeToOne = TRUE, ...){
    od <- setwd(dir); on.exit(setwd(od));
    write.ices(dat$catchMeanWeight, "cw.dat", writeToOne=writeToOne, ...)
    write.ices(dat$disMeanWeight, "dw.dat", writeToOne=writeToOne,...)
    write.ices(dat$landMeanWeight, "lw.dat", writeToOne=writeToOne,...)
    write.ices(dat$landFrac, "lf.dat", writeToOne=writeToOne,...)  
    write.ices(dat$propMat, "mo.dat", writeToOne=writeToOne,...)    
    write.ices(dat$stockMeanWeight, "sw.dat", writeToOne=writeToOne,...)
    write.ices(dat$propF, "pf.dat", writeToOne=writeToOne,...)
    write.ices(dat$propM, "pm.dat", writeToOne=writeToOne,...)
    write.ices(dat$natMor, "nm.dat", writeToOne=writeToOne,...)
    fit <- list(data=dat)
    write.ices(getResidualFleets(fit),"cn.dat", writeToOne=writeToOne,...)
    write.surveys(fit, "survey.dat", ...)
    sumFleets <- getResidualSumFleets(fit)
    if(length(sumFleets) > 0)
        write.ices(sumFleets,"cn_sum.dat", writeToOne=writeToOne,...)
    setwd(od)
}
