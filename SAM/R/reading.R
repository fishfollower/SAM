##' Function to supress incomplete final line warning
##' @param ...
##' @details ...
read.table.nowarn<-function(...){
  tryCatch.W.E <- function(expr)
  {
    W <- NULL
    w.handler <- function(w){ # warning handler
      if(!grepl('incomplete final line',w))W<<-w
      invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),warning = w.handler),warning = W)
  }
  lis<-tryCatch.W.E(read.table(...))
  if(!is.null(lis$warning))warning(lis$warning)
  lis$value
}

##' Function to test if x is ...
##' @param x
##' @details ...
is.whole.positive.number <- function(x, tol = .Machine$double.eps^0.5){
  (abs(x - round(x)) < tol)&(x>=0)
}

##' Function to read ices survey format
##' @param filen
##' @details ...
read.surveys<-function(filen){
  # Function to read ices survey file 

  if(!file.exists(filen)){
    stop(paste("File",filen, "does not exsist"))
  }

  lin<-readLines(filen,warn=FALSE)[-c(1:2)]
  empty<-which(lapply(lapply(strsplit(lin, split='[[:space:]]+'), 
               paste, collapse=''), nchar)==0)
  if(length(empty)>0){
    lin<-lin[-empty]
  }
  lin<-sub("^\\s+", "",lin)
  idx1<-grep('^[A-Z#]', lin, ignore.case=TRUE)
  idx2<-c(idx1[-1]-1,length(lin))
  names<-lin[idx1]
  years<-matrix(as.numeric(unlist(strsplit(lin[idx1+1], '[[:space:]]+'))), ncol=2, byrow=TRUE)
  times<-matrix(as.numeric(unlist(strsplit(lin[idx1+2], '[[:space:]]+'))), ncol=4, byrow=TRUE)[,3:4,drop=FALSE]
  ages<-matrix(as.numeric(unlist(lapply(strsplit(lin[idx1+3], '[[:space:]]+'), function(x)x[1:2]))), ncol=2, byrow=TRUE)
  for(i in 1:length(names)){
    # Check years 
    if(!is.whole.positive.number(years[i,1])){
      stop(paste("In file",filen, ": Minimum year is expected to be a positive integer number for fleet number",i))
    }
    if(!is.whole.positive.number(years[i,2])){
      stop(paste("In file",filen, ": Maximum year is expected to be a positive integer number for fleet number",i))
    }
    if(years[i,1]>years[i,2]){
      stop(paste("In file",filen, ": Maximum year is expected to be greater than minimum year for fleet number",i))
    }
    # Check ages 
    if(!is.whole.positive.number(ages[i,1])){
      stop(paste("In file",filen, ": Minimum age is expected to be a positive integer number for fleet number",i))
    }
    if(!is.whole.positive.number(ages[i,2])){
      stop(paste("In file",filen, ": Maximum age is expected to be a positive integer number for fleet number",i))
    }    
    if(ages[i,1]>ages[i,2]){
      stop(paste("In file",filen, ": Maximum age is expected to be greater than minimum age for fleet number",i))
    }
    # Check times
    if((times[i,1]<0)|(times[i,1]>1)){
      stop(paste("In file",filen, ": Minimum survey time is expected to be within [0,1] for fleet number",i))
    } 
    if((times[i,2]<0)|(times[i,2]>1)){
      stop(paste("In file",filen, ": Maximum survey time is expected to be within [0,1] for fleet number",i))
    } 
    if(times[i,2]<times[i,1]){
      stop(paste("In file",filen, ": Maximum survey time is expected greater than minimum survey time for fleet number",i))
    } 
  }
  onemat<-function(i){
    lin.local<-gsub('^[[:blank:]]*','',lin[(idx1[i]+4):idx2[i]])
    nr<-idx2[i]-idx1[i]-3
    ret<-matrix(as.numeric(unlist((strsplit(lin.local,'[[:space:]]+')))),nrow=nr, byrow=TRUE)   #[,1:(2+ages[i,2]-ages[i,1]),drop=FALSE]
    if(nrow(ret)!=(years[i,2]-years[i,1]+1)){
      stop(paste("In file",filen, ": Year range specified does not match number of rows for survey fleet number",i))
    } 
    if((ncol(ret)-1)<(ages[i,2]-ages[i,1]+1)){
      stop(paste("In file",filen, ": Fewer columns than indicated by age range for survey fleet number",i))
    } 
    if(!is.numeric(ret)){
      stop(paste("In file",filen, ": Non numeric data values detected for survey fleet number",i))
    }
    ret<-as.matrix(ret[,-1]/ret[,1])
    rownames(ret)<-years[i,1]:years[i,2]
    ret<-ret[,1:length(ages[i,1]:ages[i,2])]
    colnames(ret)<-ages[i,1]:ages[i,2]
    attr(ret,'time')<-times[i,]
    ret[ret<0]<-NA
    ret
  }
  obs<-lapply(1:length(names),onemat)  
  names(obs)<-names
  obs
}

##' Function to read ICES/CEFAS data files and validate if input makes sense  
##' @param filen The filename
##' @importFrom TMB MakeADFun sdreport
##' @details
##' First two lines are ignored and can be used for comments. 
##' Can read formats 1 full, 2 row, 3 scalar, and 5 column
##'
##' Tests: 
##' Formatcode is valid, years and ages are pos. integers 
##' minimum <= maximum for years and ages
##' number of rows and coulmns match year and age ranges
##' data contains only numbers.  
##' 
##' Returns: A validated data matrix.
##' @export
read.ices<-function(filen){

  if(!file.exists(filen)){
    stop(paste("File",filen, "does not exsist"))
  }

  if(grepl("[0-9]", scan(filen, skip=2, n=1, quiet=TRUE, what=""))){ # is not a survey file 
    
    head<-scan(filen, skip=2, n=5, quiet=TRUE)
    minY<-head[1]
    maxY<-head[2]
    minA<-head[3]
    maxA<-head[4]
    datatype<-head[5]
  
    if(!is.whole.positive.number(minY)){
      stop(paste("In file",filen, ": Minimum year is expected to be a positive integer number"))
    }
    if(!is.whole.positive.number(maxY)){
      stop(paste("In file",filen, ": Maximum year is expected to be a positive integer number"))
    }
    if(!is.whole.positive.number(minA)){
      stop(paste("In file",filen, ": Minimum age is expected to be a positive integer number"))
    }
    if(!is.whole.positive.number(maxA)){
      stop(paste("In file",filen, ": Maximum age is expected to be a positive integer number"))
    }
  
    if(!(datatype%in%c(1,2,3,5))){
      stop(paste("In file",filen, ": Datatype code is expected to be one of the numbers 1, 2, 3, or 5"))
    }
  
    if(minY>maxY){
      stop(paste("In file",filen, ": Minimum year is expected to be less than maximum year"))
    }
    if(minA>maxA){
      stop(paste("In file",filen, ": Minimum age is expected to be less than maximum age"))
    }
  
    C<-as.matrix(read.table.nowarn(filen, skip=5, header=FALSE))
  
    if(datatype==1){
      if((maxY-minY+1)!=nrow(C)){
        stop(paste("In file",filen, ": Number of rows does not match the year range given"))
      }
      if((maxA-minA+1)>ncol(C)){
        stop(paste("In file",filen, ": Fewer columns than the age range given"))
      }
    }
  
    if(datatype==2){
      C<-as.matrix(read.table.nowarn(filen, skip=5, header=FALSE))
      if(1!=nrow(C)){
        stop(paste("In file",filen, ": For datatype 2 only one row of data is expected"))
      }
      if((maxA-minA+1)>ncol(C)){
        stop(paste("In file",filen, ": Fewer columns than the age range given"))
      }
      C<-C[rep(1,maxY-minY+1),]
    }
  
    if(datatype==3){
      C<-as.matrix(read.table.nowarn(filen, skip=5, header=FALSE))
      if(1!=nrow(C)){
        stop(paste("In file",filen, ": For datatype 3 only one row of data is expected"))
      }
      if(1!=ncol(C)){
        stop(paste("In file",filen, ": For datatype 3 only one column of data is expected"))
      }
      C<-C[rep(1,maxY-minY+1),rep(1,maxA-minA+1)]
    }
  
    if(datatype==5){
      C<-as.matrix(read.table.nowarn(filen, skip=5, header=FALSE))
      if((maxY-minY+1)!=nrow(C)){
        stop(paste("In file",filen, ": Number of rows does not match the year range given"))
      }
      if(1!=ncol(C)){
        stop(paste("In file",filen, ": For datatype 5 only one column of data is expected"))
      }
      C<-C[,rep(1,maxA-minA+1)]
    }
    rownames(C)<-minY:maxY
    C<-C[,1:length(minA:maxA)]
    colnames(C)<-minA:maxA
  
    if(!is.numeric(C)){
        stop(paste("In file",filen, ": Non numeric data values detected (could for instance be comma used as decimal operator)"))
    }
    return(C)
  }else{
    return(read.surveys(filen))  
  }
}
