library(stockassessment)
source("src/common.R")
# load what has been saved
setwd("run")
for(f in dir(pattern="RData"))load(f) 
setwd("..")

basefit<-NULL
if(file.exists("baserun/model.RData")){
  local({load("baserun/model.RData"); basefit<<-fit})
  if(abs(logLik(basefit)-logLik(fit))<1.0e-9)basefit<<-NULL
}


plotcounter<-1
tit.list<-list()

setcap<-function(title="", caption=""){   
 tit.list[length(tit.list)+1]<<-paste("# Title",plotcounter,"\n")
 tit.list[length(tit.list)+1]<<-paste(title,"\n")
 tit.list[length(tit.list)+1]<<-paste("# Caption",plotcounter,"\n")
 tit.list[length(tit.list)+1]<<-paste(caption,"\n")
 plotcounter<<-plotcounter+1 
}


############################## plots ##############################
plots<-function(){
  par(cex.lab=1, cex.axis=1, mar=c(5,5,1,1))
    
  if(exists("fit")){
    fits <- c(fit,basefit)
    
    ssbplot(fits, addCI=TRUE)
    stampit()
    setcap("Spawning stock biomass", "Spawning stock biomass. 
            Estimates from the current run and point wise 95% confidence 
            intervals are shown by black line and shaded area.")
    
    fbarplot(fits, addCI=TRUE)
    stampit()

    recplot(fits, addCI=TRUE, las=0)

    stampit()

    catchplot(fits, addCI=TRUE)
    stampit()
    
    plot(ypr(fit, aveYears = 20))
    stampit()
    
    obscorrplot(fit)
    stampit()

    
    for(f in 1:fit$data$noFleets){
      fitplot(fit, fleets=f)
      stampit()
    }
    
    

    #Q<-fit$pl$logFpar
    #Qsd<-fit$plsd$logFpar
    #key<-fit$conf$keyLogFpar
    #fun<-function(x)if(x<0){NA}else{Q[x+1]}
    #FF<-Vectorize(fun)
    #ages<-fit$conf$minAge:fit$conf$maxAge
    #matplot(ages, exp(t(matrix(FF(key), nrow=5))), type="l", lwd=5, lty="solid", xlab="Ages", ylab="Q")
    #legend("topright", lwd=5, col=2:5, legend=attr(fit$data, "fleetNames")[2:5])
    #stampit()

  }  
  
  if(exists("RES")){  
    plot(RES)
    par(mfrow=c(1,1))
    stampit()
  
  }
  
   
  if(exists("RESP")){  
    plot(RESP)
    par(mfrow=c(1,1))
    stampit()
  } 
 
  
  if(exists("LO")){  
    ssbplot(LO)
    stampit()
    
    fbarplot(LO)
    stampit()

    recplot(LO)
    stampit()

    catchplot(LO)
    stampit()
    
  } 
  
  if(exists("RETRO")){  
    ssbplot(RETRO, las=0, drop=1)
    stampit()
    
    fbarplot(RETRO, las=0, drop=1)
    stampit()
    
    recplot(RETRO, las=0, drop=1)
    stampit()

    catchplot(RETRO)
    stampit()
  } 
  if(exists("FC")){  
    lapply(FC, function(f){plot(f); title(attr(f,"label"), outer=TRUE, line=-1); stampit()})
  }  
  
}


setwd('res')
file.remove(dir(pattern='png$'))
stamp<-gsub('-[[:digit:]]{4}$','',gsub(':','.',gsub(' ','-',gsub('^[[:alpha:]]{3} ','',date()))))
png(filename = paste(stamp,"_%03d.png", sep=''), width = 480, height = 480,
    units = "px", pointsize = 10, bg = "white")
  plots()    
dev.off()

writeLines(unlist(tit.list),'titles.cfg') 

png(filename = paste("big_",stamp,"_%03d.png", sep=''), width = 1200, height = 1200, 
    units = "px", pointsize = 20, bg = "white")
  plots()    
dev.off()

#pdf(onefile=FALSE, width = 8, height = 8)
#  plots()    
#dev.off()

file.remove(dir(pattern='html$'))

tsb<-tsbtable(fit)
colnames(tsb)<-c("TSB","Low", "High")
tab.summary <- cbind(summary(fit), tsb)
xtab(tab.summary, caption=paste('Table 1. Estimated recruitment, spawning stock biomass (SSB), 
     and average fishing mortality','.',sep=''), cornername='Year', 
     file=paste(stamp,'_tab1.html',sep=''), dec=c(0,0,0,0,0,0,3,3,3,0,0,0))

ftab <- faytable(fit)
xtab(ftab, caption=paste('Table 2. Estimated fishing mortality at age','.',sep=''), cornername='Year \ Age', 
     file=paste(stamp,'_tab2.html',sep=''), dec=rep(3,ncol(ftab)))

ntab <- ntable(fit)
xtab(ntab, caption=paste('Table 3. Estimated stock numbers at age','.',sep=''), cornername='Year \ Age', 
     file=paste(stamp,'_tab3.html',sep=''), dec=rep(0,ncol(ntab)))

ptab <- partable(fit)
xtab(ptab, caption=paste('Table 4. Table of model parameters','.',sep=''), cornername='Parameter name', 
     file=paste(stamp,'_tab4.html',sep=''), dec=rep(3,ncol(ptab)))

mtab <- modeltable(c(Current=fit, base=basefit))
mdec <- c(2,0,2,6)[1:ncol(mtab)]
xtab(mtab, caption=paste('Table 5. Model fitting','.',sep=''), cornername='Model', 
     file=paste(stamp,'_tab5.html',sep=''), dec=mdec)
     
sdState<-function(fit, y=max(fit$data$years)-1:0){
  idx <- names(fit$sdrep$value) == "logR"
  sdLogR<-fit$sdrep$sd[idx][fit$data$years%in%y]
  idx <- names(fit$sdrep$value) == "logssb"
  sdLogSSB<-fit$sdrep$sd[idx][fit$data$years%in%y]
  idx <- names(fit$sdrep$value) == "logfbar"
  sdLogF<-fit$sdrep$sd[idx][fit$data$years%in%y]
  ret<-cbind(sdLogR, sdLogSSB, sdLogF)
  rownames(ret)<-y
  colnames(ret)<-c("sd(log(R))", "sd(log(SSB))", "sd(log(Fbar))")
  return(ret)
}

sdtab<-sdState(fit)
xtab(sdtab, caption=paste('Table 6. Table of selected sd','.',sep=''), cornername='Year', 
     file=paste(stamp,'_tab6.html',sep=''), dec=rep(3,ncol(sdtab)))


if(exists("FC")){  
    ii<-0
    lapply(FC, function(f){
       ii<<-ii+1;
       tf<-attr(f,"tab");
       dec<-c(3,3,3,rep(0,ncol(tf)-3));
       xtab(tf, caption=paste0('Forecast table ',ii,'. ', attr(f,"label"),'.'), 
       cornername='Year', file=paste(stamp,'_tabX',ii,'.html',sep=''), dec=dec);      
       })
}  

setwd("..") 

