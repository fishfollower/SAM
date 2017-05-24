source('src/datavalidator.R')

dummyplot<-function(text='An error was found in the data, please check error logs'){
  plot(c(0,1),c(0,1),axes=FALSE, xlab='', ylab='' , type='n')
  box()
  text(.5,.5,labels=text)
}

plotit<-function(fil){
  key=c(cn='Catch in numbers', survey='Survey index', nm='Natural mortality', mo='Maturity ogive', 
        sw='Stock mean weight', cw='Catch mean weight', lw='Land mean weight', dw='Discard mean weight',
        lf='Fraction landed', pf='F before spawning', pm='M before spawning')
  wrap<-function(expr){
    tryCatch(expr,error=function(e){E<<-e})
  }
  if(file.exists(fil)){
    basename<-sub('\\..*$','',fil)
    png(filename = paste(basename,".png", sep=""), width = 600, height = 600, units = "px", pointsize = 10, bg = "white")    
      E<-NULL
      if(fil!='survey.dat'){
        x<-wrap(read.ices(fil))
        if(is.null(E)){
          if(all(apply(x,2,diff)==0)){
            plot(as.numeric(colnames(x)), x[1,], xlab='Age', ylab=key[basename], type='b', 
                 lwd=2, cex=2, col='blue')
          }else{
            matplot(rownames(x),x,type='l', xlab='Year', ylab=key[basename], lwd=2)
            legend('topleft',bty='n', legend=paste('Age',colnames(x)), lty=1:5,col=1:6, ncol=5, lwd=2)      
          }
        }else{dummyplot()}
      }else{
        x<-wrap(read.ices(fil))
        if(is.null(E)){
          div<-list(c(1,1),c(2,1),c(3,1),c(2,2),c(3,2),c(3,2),c(3,3),c(3,3),c(3,3),c(4,3),c(4,3),c(4,3))
          op<-par(mfrow=div[[length(x)]], mar=c(4,4,1,2), mgp=c(2,1,0))
          fun<-function(i){
            xx<-x[[i]]
            matplot(rownames(xx),xx,type='l', xlab='Year', ylab=key[basename], main=names(x)[i],lwd=2)
            legend('topleft',bty='n', legend=paste('Age',colnames(xx)), lty=1:5,col=1:6, ncol=5, lwd=2)      
          } 
          lapply(1:length(x),fun)->dummy
        }else{dummyplot()}
      }
    dev.off() 
  }
}

owd<-getwd()
setwd('data')
plotit('cn.dat')
plotit('cw.dat')
plotit('dw.dat')
plotit('lw.dat')
plotit('sw.dat')
plotit('mo.dat')
plotit('nm.dat')
plotit('survey.dat')
plotit('lf.dat')
plotit('pf.dat')
plotit('pm.dat')
setwd(owd)


