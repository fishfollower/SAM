##' SSB table 
##' @param  fit ... 
##' @details ...
##' @export
ssbtable<-function(fit,...){
   idx<-names(fit$sdrep$value)=="logssb"
   y<-fit$sdrep$value[idx]
   ci<-y+fit$sdrep$sd[idx]%o%c(-2,2)
   ret<-cbind(y,ci)
   rownames(ret)<-fit$data$years
   colnames(ret)<-c("Estimate","Low","Hig")
   return(exp(ret))
}
