##' Table helper 
##' @param fit returned object from sam.fit
##' @param quoted name of what to extract
##' @param x rownames of table
##' @param trans function to be applied
##' @details ...
.tableit <-function (fit, what, x=fit$data$years, trans=function(x)x){
   idx<-names(fit$sdrep$value)==what
   y<-fit$sdrep$value[idx]
   ci<-y+fit$sdrep$sd[idx]%o%c(-2,2)
   ret<-trans(cbind(y,ci))
   rownames(ret)<-x
   colnames(ret)<-c("Estimate","Low","Hig")
   return(ret)
}

##' SSB table 
##' @param  fit ...
##' @details ...
##' @export
ssbtable<-function(fit){
   ret<-.tableit(fit, "logssb", trans=exp) 
   return(ret)
}

##' TSB table 
##' @param  fit ... 
##' @details ...
##' @export
tsbtable<-function(fit){
   ret<-.tableit(fit, "logtsb", trans=exp) 
   return(ret)
}

##' Fbar table 
##' @param  fit ... 
##' @details ...
##' @export
fbartable<-function(fit){
   ret<-.tableit("logfbar", trans=exp) 
   return(ret)
}

##' Recruit table 
##' @param  fit ... 
##' @details ...
##' @export
rectable<-function(fit){
   ret<-.tableit(fit, "logR", trans=exp) 
   return(ret)
}

