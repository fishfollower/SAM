##' Table helper 
##' @param  fit ... 
##' @details ...
.tableit <-function (what, x=fit$data$years, ylab=what, trans=function(x)x ,...){
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
ssbtable<-function(fit,...){
   ret<-.tableit("logssb", trans=exp, ...) 
   return(ret)
}

##' Fbar table 
##' @param  fit ... 
##' @details ...
##' @export
fbartable<-function(fit,...){
   ret<-.tableit("logfbar", trans=exp, ...) 
   return(ret)
}

##' Recruit table 
##' @param  fit ... 
##' @details ...
##' @export
rectable<-function(fit,...){
   ret<-.tableit("logR", trans=exp, ...) 
   return(ret)
}

