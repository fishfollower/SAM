# Script to validate data to match the state-space assessment model 
# 
# Anders Nielsen <anders@nielsensweb.org> Jul. 2012

getsourcelist<-function(defaults=c("common.R", "dataplot.R", "datascript.R", 
                                   "datavalidator.R", "forecast.R", "leaveout.R", 
                                   "model.R", "plotscript.R", "residuals.R", "retro.R",
                                   "sourcevalidator.R")){
  ret<-unique(c(defaults,dir(pattern="*.R")))
  return(ret)
}  


test.sourcenames<-function(path, sourcelist=getsourcelist()){
  # Function to test if specified filenames exist  
  owd<-getwd()
  setwd(path)
  ret<-file.exists(sourcelist)
  setwd(owd)
  return(data.frame(file=sourcelist,exists=ret))  
}

check.one.source<-function(filen){
    dummy<-parse(filen)
}

check.all.source<-function(path,filen=getsourcelist()){  
  owd<-getwd()
  setwd(path)
  out<-test.sourcenames('.',filen)
  ret<-list()
  self<-rep(NA,length(filen))
  message<-rep("",length(filen))
  wrap<-function(expr){
    tryCatch(expr,error=function(e){E<<-e})
  }
  for(i in 1:length(filen)){
    if(out[i,2]){
      fil<-filen[i]
      E<-NULL
      dummy<-wrap(parse(fil))
      if(is.null(E)){
        self[i]<-TRUE
      }else{
        self[i]<-FALSE
        message[i]<-E$message
      }
    }
  }
  out$self<-self
  out$message<-message
  setwd(owd)
  return(out)
}
