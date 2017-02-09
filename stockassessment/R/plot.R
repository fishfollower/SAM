##' Plot helper
##' @param fit the fitted object from sam.fit
##' @param what quoted name of object to extract
##' @param x x-alues
##' @param ylab label on y-axis
##' @param xlab label on x-axis
##' @param ex extra y's to make room for
##' @param trans function to transform values by
##' @param add logical, plotting is to be added on existing plot
##' @param ci logical, confidence intervals should be plotted
##' @param cicol color to plot the confidence polygon
##' @param addCI logicial plot CI for added fits 
##' @param drop number of years to be left unplotted at the end.
##' @param ... extra arguments transferred to plot
##' @importFrom graphics plot polygon grid lines
##' @importFrom grDevices gray
##' @details The basic plotting used bu many of the plotting functions (e.g. ssbplot, fbarplot ...) 
.plotit <-function (fit, what, x=fit$data$years, ylab=what, xlab="Years", ex=numeric(0), trans=function(x)x, add=FALSE, ci=TRUE, cicol=gray(.5,alpha=.5), addCI=FALSE, drop=0, ...){
  if(class(fit)!="samset"){ 
    idx <- names(fit$sdrep$value)==what
    y <- fit$sdrep$value[idx]
    lowhig <- y+fit$sdrep$sd[idx]%o%c(-2,2)
    didx <- 1:(length(x)-drop)
    xr<-range(x)
    x<-x[didx]
    y<-y[didx]
    lowhig<-lowhig[didx,]
    if(add){
      lines(x, trans(y), lwd=3,...)
    }else{
      plot(x, trans(y), xlab=xlab, ylab=ylab, type="n", lwd=3, xlim=xr, ylim=range(c(trans(lowhig),0,ex)), las=1,...)
      grid(col="black")
      lines(x, trans(y), lwd=3, ...)
    }
    if(ci){
      polygon(c(x,rev(x)), y = c(trans(lowhig[,1]),rev(trans(lowhig[,2]))), border = gray(.5,alpha=.5), col = cicol)
    }
  }else{
    colSet=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")  
    idxfrom <- 1
    if(is.null(attr(fit,"fit"))){
      attr(fit,"fit") <- fit[[1]]
      idxfrom <- 2
    }
    if(is.null(x))x=attr(fit,"fit")$data$years
    .plotit(attr(fit,"fit"), what=what, x=x, ylab=ylab, xlab=xlab, ex=ex, trans=trans, add=add, ci=ci, cicol=cicol, drop=drop,...)
    if(addCI){
      d<-lapply(idxfrom:length(fit), function(i).plotit(fit[[i]], what=what, trans=trans, add=TRUE, col=colSet[(i-1)%%length(colSet)+1], cicol=paste0(colSet[(i-1)%%length(colSet)+1],"80"), drop=drop, ...)) 
    }else{
      d<-lapply(idxfrom:length(fit), function(i).plotit(fit[[i]], what=what, trans=trans, add=TRUE, ci=FALSE, col=colSet[(i-1)%%length(colSet)+1], drop=drop, ...))
    }
  }
}
   
##' Plot by one or two  
##' @param x numeric vector
##' @param y numeric vector
##' @param z numeric vector
##' @param by vector or two column matrix to create sub sets from  
##' @param bubblescale scaling of bubble size
##' @param x.common logical: use same x-axis for all plots
##' @param y.common logical: use same y-axis for all plots
##' @param xlab normal graphical parameter
##' @param ylab normal graphical parameter
##' @param xlim normal graphical parameter
##' @param ylim normal graphical parameter
##' @param axes normal graphical parameter
##' @param ... additional graphical parameters 
##' @importFrom graphics plot layout axis box mtext legend
##' @importFrom grDevices gray rgb
##' @details Function used for splitting plots e.g. used to plot residuals 
##' @export
##' @examples
##' exdat<-expand.grid(age=1:5, year=1950:2016, fleet=1:3)
##' exdat$perfectres<-rnorm(nrow(exdat))
##' attach(exdat)
##' par(ask=FALSE)
##' plotby(year,age,perfectres, by=fleet)
##' detach(exdat)
plotby <-function(x=NULL, y=NULL, z=NULL, by=NULL, bubblescale=1, x.common=TRUE, y.common=TRUE, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, axes=TRUE, ...){
  if(is.null(by) | length(unique(by))==1){
    # 0
    if(length(x)==0 & length(y)==0 & length(z)==0){
      plot(NA, xlim=0:1, ylim=0:1, xlab="", ylab="", axes=FALSE, ...)
    }
    # 1
    if(length(x)>0 & length(y)==0 & length(z)==0){
      if(missing(ylab)) ylab <- deparse(substitute(x))
      if(missing(ylim)) ylim <- c(min(x)-1, max(x)+1)      
      plot(x, ylab=ylab, ylim=ylim, axes=axes, ...)
    }
    if(length(x)==0 & length(y)>0 & length(z)==0){
      if(missing(ylab)) ylab <- deparse(substitute(y))
      if(missing(ylim)) ylim <- c(min(y)-1, max(y)+1)      
      plot(y, ylab=ylab, ylim=ylim, axes=axes, ...)
    }
    if(length(x)==0 & length(y)==0 & length(z)>0){
      if(missing(ylab)) ylab <- deparse(substitute(z))
      if(missing(ylim)) ylim <- c(min(z)-1, max(z)+1)      
      plot(z, ylab=ylab, ylim=ylim, axes=axes, ...)
    }
    # 2
    if(length(x)>0 & length(y)>0 & length(z)==0){
      if(missing(xlab)) xlab <- deparse(substitute(x))
      if(missing(ylab)) ylab <- deparse(substitute(y))
      if(missing(xlim)) xlim <- c(min(x)-1, max(x)+1)
      if(missing(ylim)) ylim <- c(min(y)-1, max(y)+1)
      plot(x, y, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, axes=axes, ...)
    }
    if(length(x)>0 & length(y)==0 & length(z)>0){
      if(missing(xlab)) xlab <- deparse(substitute(x))
      if(missing(ylab)) ylab <- deparse(substitute(z))
      if(missing(xlim)) xlim <- c(min(x)-1, max(x)+1)
      if(missing(ylim)) ylim <- c(min(z)-1, max(z)+1)
      plot(x, z, xlab=xlab,  ylab=ylab, xlim=xlim, ylim=ylim, axes=axes, ...)
    }
    if(length(x)==0 & length(y)>0 & length(z)>0){
      if(missing(xlab)) xlab <- deparse(substitute(y))
      if(missing(ylab)) ylab <- deparse(substitute(z))
      if(missing(xlim)) xlim <- c(min(y)-1, max(y)+1)
      if(missing(ylim)) ylim <- c(min(z)-1, max(z)+1)
      plot(y, x, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, axes=axes, ...)
    }
    # 3
    if(length(x)>0 & length(y)>0 & length(z)>0){
      if(missing(xlab)) xlab <- deparse(substitute(x))
      if(missing(ylab)) ylab <- deparse(substitute(y))
      if(missing(xlim)) xlim <- c(min(x)-1, max(x)+1)
      if(missing(ylim)) ylim <- c(min(y)-1, max(y)+1)
      cex=sqrt(abs(z))/max(sqrt(abs(z)))*5*bubblescale
      plot(x, y, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, type="n", axes=axes, ...)
      neg <- z<0
      points(x[neg],y[neg], cex=cex[neg], col=rgb(1, 0, 0, alpha=.5), pch=19, ...)
      points(x[!neg],y[!neg], cex=cex[!neg], col=rgb(0, 0, 1, alpha=.5), pch=19, ...)
    }
  }
  if(is.vector(by)|is.factor(by)){
    uby <- unique(as.vector(by))
    if(length(uby)==3){
      div<-c(3,1)
    }else{    
      div<-rep(ceiling(sqrt(length(uby))),2)
      if(div[1]*(div[2]-1)>=length(uby))div[2] <- div[2]-1
    }
    laym<-matrix(1:(div[1]*div[2]),nrow=div[1], ncol=div[2])
    layout(laym)
    if(missing(xlab)) xlab <- deparse(substitute(x))
    if(missing(ylab)) ylab <- deparse(substitute(y))
    if(x.common) xlim <- c(min(x)-1, max(x)+1)
    if(y.common) ylim <- c(min(y)-1, max(y)+1)
    if(x.common) op1<-par(oma=c(par("mar")[1], par("oma")[2], par("mar")[3], par("oma")[4]),
                          mar=c(.2,par("mar")[2],.2,par("mar")[4]))
    if(y.common) op2<-par(oma=c(par("oma")[1], par("mar")[2], par("oma")[3], par("mar")[4]),
                          mar=c(par("mar")[1],.2,par("mar")[3],.2))
    
    .fun<-function(i){
      b<-uby[i]
      plotby(x[by==b], y[by==b], z[by==b], xlab=ifelse(x.common,"",xlab), ylab=ifelse(y.common,"",ylab),
             xlim=xlim, ylim=ylim, bubblescale=bubblescale, axes=FALSE, ...)
      legend("top", bty="n", legend=uby[i], text.col=gray(.5))
      if(!x.common)axis(1)
      if(!y.common)axis(2)
      if(x.common&(row(laym)[i]==div[1]))axis(1)
      if(x.common&(i==length(uby)))axis(1)
      if(y.common&(col(laym)[i]==1))axis(2)
      box()
    }
    d <- lapply(1:length(uby), .fun)
    if(x.common)mtext(xlab, side = 1, line = 3, outer = TRUE, ...)
    if(y.common)mtext(ylab, side = 2, line = 2.5, outer = TRUE, ...)
    if(x.common)par(op1)
    if(y.common&!x.common)par(op2)
  }  
  if(is.matrix(by)){      
    uby1 <- unique(by[,1])
    uby2 <- unique(by[,2])
    div<-c(length(uby1), length(uby2))
    laym<-matrix(1:(div[1]*div[2]),nrow=div[1], ncol=div[2])
    layout(laym)
    if(missing(xlab)) xlab <- deparse(substitute(x))
    if(missing(ylab)) ylab <- deparse(substitute(y))
    if(x.common) xlim <- c(min(x)-1, max(x)+1)
    if(y.common) ylim <- c(min(y)-1, max(y)+1)
    if(x.common) op1<-par(oma=c(par("mar")[1], par("oma")[2], par("mar")[3], par("oma")[4]),
                          mar=c(.2,par("mar")[2],.2,par("mar")[4]))
    if(y.common) op2<-par(oma=c(par("oma")[1], par("mar")[2], par("oma")[3], par("mar")[4]),
                          mar=c(par("mar")[1],.2,par("mar")[3],.2))
    
    fun<-function(i,j){
      b<-(by[,1]==uby1[i])&(by[,2]==uby2[j])
      plotby(x[b], y[b], z[b], xlab=ifelse(x.common,"",xlab), ylab=ifelse(y.common,"",ylab),
             xlim=xlim, ylim=ylim, bubblescale=bubblescale, axes=FALSE, ...)
      legend("top", bty="n", legend=paste(uby1[i],":",uby2[j], sep=""), text.col=gray(.5))
      if(!x.common)axis(1)
      if(!y.common)axis(2)
      if(x.common&(i==div[1]))axis(1)
      if(y.common&(j==1))axis(2)
      box()
    }
    for(j in 1:length(uby2)){for(i in 1:length(uby1))fun(i,j)}
    if(x.common)mtext(xlab, side = 1, line = 3, outer = TRUE, ...)
    if(y.common)mtext(ylab, side = 2, line = 2.5, outer = TRUE, ...)
    if(x.common)par(op1)
    if(y.common&!x.common)par(op2)
  }  
}

##' SAM Fbar plot 
##' @param fit the object returned from sam.fit 
##' @param partial true if included partial F's are to be plotted
##' @param drop number of years to be left unplotted at the end. Default (NULL) is to not show years at the end with no catch information  
##' @param ... extra arguments transferred to plot including the following: \cr
##' \code{add} logical, plotting is to be added on existing plot \cr
##' \code{ci} logical, confidence intervals should be plotted \cr
##' \code{cicol} color to plot the confidence polygon
##' @importFrom graphics matplot
##' @details Plot the defined fbar. 
##' @export
fbarplot<-function(fit,partial=(class(fit)=="sam"), drop=NULL,...){
  if(class(fit)=="sam"){
    fitlocal <- fit
  }
  if(class(fit)=="samset"){
    if(!is.null(attr(fit,"fit"))){
      fitlocal <- attr(fit,"fit")
    }else{
      fitlocal <- fit[[1]]
    }
  }
  if(is.null(drop)){
    drop=max(fitlocal$data$aux[,"year"])-max(fitlocal$data$aux[fitlocal$data$aux[,"fleet"]==1,"year"])
  }
  
  fbarRange<-fitlocal$conf$fbarRange
  fbarlab=substitute(bar(F)[X-Y],list(X=fbarRange[1],Y=fbarRange[2]))
  fmat<-fitlocal$pl$logF[fitlocal$conf$keyLogFsta[1,]+1,]
  idx<-which(fitlocal$conf$minAge:fitlocal$conf$maxAge %in% fbarRange[1]:fbarRange[2])
  exx <- if(partial){exp(fmat[idx,])}else{numeric(0)}
  .plotit(fit, "logfbar", ylab=fbarlab, trans=exp, ex=exx, drop=drop, ...)
  if(partial){
    matplot(fitlocal$data$years, t(exp(fmat[idx,])), add=TRUE, type="b", col="lightblue", pch=as.character(fbarRange[1]:fbarRange[2]))
  }  
}

##' SAM SSB plot 
##' @param fit the object returned from sam.fit
##' @param ... extra arguments transferred to plot including the following: \cr
##' \code{add} logical, plotting is to be added on existing plot \cr
##' \code{ci} logical, confidence intervals should be plotted \cr
##' \code{cicol} color to plot the confidence polygon
##' @details Plot of spawning stock biomass 
##' @export
ssbplot<-function(fit, ...){
  .plotit(fit, "logssb", ylab="SSB", trans=exp,...)
}

##' SAM TSB plot 
##' @param fit the object returned from sam.fit
##' @param ... extra arguments transferred to plot including the following: \cr
##' \code{add} logical, plotting is to be added on existing plot \cr
##' \code{ci} logical, confidence intervals should be plotted \cr
##' \code{cicol} color to plot the confidence polygon
##' @details Plot of total stock biomass
##' @export
tsbplot<-function(fit, ...){
  .plotit(fit, "logtsb", ylab="TSB", trans=exp,...)
}

##' SAM Recruits plot 
##' @param fit the object returned from sam.fit
##' @param ... extra arguments transferred to plot including the following: \cr
##' \code{add} logical, plotting is to be added on existing plot \cr
##' \code{ci} logical, confidence intervals should be plotted \cr
##' \code{cicol} color to plot the confidence polygon
##' @details Plot of numbers of recruits (youngest age class)
##' @export
recplot<-function(fit,...){
   lab<-paste("Recruits (age ", fit$conf$minAge, ")", sep="")
  .plotit(fit, "logR", ylab=lab, trans=exp,...)
}

##' SAM catch plot 
##' @param fit the object returned from sam.fit
##' @param obs.show if observations are to be shown also
##' @param drop number of years to be left unplotted at the end. Default (NULL) is to not show years at the end with no catch information 
##' @param ... extra arguments transferred to plot including the following: \cr
##' \code{add} logical, plotting is to be added on existing plot \cr
##' \code{ci} logical, confidence intervals should be plotted \cr
##' \code{cicol} color to plot the confidence polygon
##' @details Plot of estimated (and optionally observed) total catch in weight  
##' @importFrom graphics points
##' @export
catchplot<-function(fit, obs.show=TRUE, drop=NULL,...){
  if(class(fit)=="sam"){
    fitlocal <- fit
  }
  if(class(fit)=="samset"){
    if(!is.null(attr(fit,"fit"))){
      fitlocal <- attr(fit,"fit")
    }else{
      fitlocal <- fit[[1]]
    }
  }
  if(is.null(drop)){
    drop=max(fitlocal$data$aux[,"year"])-max(fitlocal$data$aux[fitlocal$data$aux[,"fleet"]==1,"year"])
  }
  CW <- fitlocal$data$catchMeanWeight
  x <- as.numeric(rownames(CW))
  .plotit(fit, "logCatch", ylab="Catch", trans=exp, drop=drop,...)
  if(obs.show){
    aux <- fitlocal$data$aux
    logobs <- fitlocal$data$logobs
    .goget <- function(y,a){
        ret <- exp(logobs[aux[,"fleet"]==1 & aux[,"year"]==y & aux[,"age"]==a])
        ifelse(length(ret)==0,0,ret)
    }
    points(x, rowSums(outer(rownames(CW), colnames(CW), Vectorize(.goget))*CW, na.rm=TRUE), pch=4, lwd=2, cex=1.2)
  }  
}


##' SAM parameter plot 
##' @param fit the object returned from sam.fit
##' @param cor.report.limit correlations with absolute value > this number is reported in the plot 
##' @param ... extra arguments transferred to plot
##' @details Plot of all estimated model parameters (fixed effects). Shown with confidence interval.  
##' @export
##' @importFrom stats cov2cor
parplot<-function(fit, cor.report.limit=0.95, ...){
  if(class(fit)=="sam"){
    fit <- list(fit)
    class(fit) <- "samset"
  }
  if(!is.null(attr(fit,"fit"))){
    fit <- c(list(attr(fit,"fit")), fit)
  }
  param <- lapply(fit, coef) 
  nam <- names(param[[1]])
  dup <- duplicated(nam)
  namadd <- rep(0,length(nam))
  for(i in 2:length(dup)){
    if(dup[i])namadd[i] <- namadd[i-1]+1
  }
  nam <- paste(nam, namadd, sep="_")
  corrs <- cov2cor(attr(param[[1]], "cov"))-diag(length(param[[1]]))
  rownames(corrs)<-nam
  colnames(corrs)<-nam
  higcor <- lapply(1:nrow(corrs),function(i)corrs[i,][corrs[i,]>cor.report.limit]*100)
  names(higcor) <- nam
  lowcor<-lapply(1:nrow(corrs),function(i)corrs[i,][corrs[i,]<(-cor.report.limit)]*100)
  names(lowcor) <- nam
  for(i in 1:length(param)){
    m <- param[[i]]+t(c(-2,0,2)%o%attr(param[[i]],"sd"))  
    if(i==1){
      mat <- cbind(m,0)
      rownames(mat) <- nam
    }else{
      if(nrow(m)==length(nam)){
        rownames(m) <- nam
        mat<-rbind(mat,cbind(m,-i+1))
      }
    }
  }
  .plotapar<-function(name){
    sub <- mat[rownames(mat)==name,,drop=FALSE]
    xold <- sub[,4]
    if(nrow(sub)==1)sub<-rbind(cbind(sub[,1:3,drop=FALSE],-.5), cbind(sub[,1:3,drop=FALSE],.5))
    x <- sub[,4]
    y <- sub[,2]
    plot(x, y, xlim=c(min(x)-1,5), ylim=range(sub[,1:3]), ylab=name, xlab="", axes=FALSE, type="n",...)
    box()
    axis(2, las=1)
    lines(x, y, lwd=3,...)
    polygon(c(x,rev(x)), y = c(sub[,1],rev(sub[,3])), border = gray(.5,alpha=.5), col = gray(.5,alpha=.5))
    d <- sapply(1:length(xold), function(i)lines(xold[c(i,i)],c(sub[i,1],sub[i,3]), lty="dotted", lwd=.5))
    idx <- which(nam==name)

    if(length(higcor[[name]])!=0){
      legend("topright", legend=paste0(names(higcor[[name]]),": ",round(higcor[[name]]),"%"), bty="n", text.col="blue")
    }
    if(length(lowcor[[name]])!=0){
      legend("bottomright", legend=paste0(names(lowcor[[name]]),": ",round(lowcor[[name]]),"%"), bty="n", text.col="red")
    }
  }
  div <- rep(ceiling(sqrt(length(nam))),2)
  if(div[1]*(div[2]-1)>=length(nam))div[2] <- div[2]-1
  op <- par(mar=c(.2,par("mar")[2],.2,par("mar")[4]))
  laym <- matrix(1:(div[1]*div[2]),nrow=div[1], ncol=div[2])
  layout(laym)
  d <- sapply(nam, .plotapar)
  par(op)
}

##' Extract observation covariance matrices from a SAM fit
##' @param fit the object returned from sam.fit
##' @param corr if TRUE return correlation matrices rather than covariances
##' @return a list of matrices
##' @export
obscov<-function(fit, corr=FALSE){
    stopifnot(class(fit)=="sam")
    res<-fit$rep$obsCov
    for(i in 1:length(res)) rownames(res[[i]])<-fit$data$minAgePerFleet[i]:fit$data$maxAgePerFleet[i]
    if(corr) for(i in 1:length(res)) res[[i]]<-cov2cor(res[[i]])
    res
}

##' Plots the estimated correlation matrices by fleet.
##' @param fit the object returned from sam.fit
##' @param ... extra arguments to plot
##' @importFrom ellipse plotcorr 
##' @export
obscorrplot<-function(fit,...){
    ccolors <- c("#A50F15","#DE2D26","#FB6A4A","#FCAE91","#FEE5D9","white",
                 "#EFF3FF","#BDD7E7","#6BAED6","#3182BD","#08519C")
    x<-obscov(fit,TRUE)

    if(length(x)==3){
      div<-c(3,1)
    }else{    
      div<-rep(ceiling(sqrt(length(x))),2)
      if(div[1]*(div[2]-1)>=length(x))div[2] <- div[2]-1
    }
    laym<-matrix(1:(div[1]*div[2]),nrow=div[1], ncol=div[2])
    layout(laym)
    
    for(i in 1:length(x)) 
        plotcorr(x[[i]],col=ccolors[5*x[[i]]+6],mar=0.1+c(2,2,2,2),...)
    
}


##' Plots the classic yield per recruit
##' @param fit the object returned from sam.fit
##' @param Flimit Upper limit for Fbar
##' @param Fdelta increments on the Fbar axis 
##' @param aveYears Number of years back to use when calculating averages (selection, weights, ...)
##' @param ageLimit Oldest age used (should be high)
##' @importFrom graphics title
##' @export
yprplot<-function(fit, Flimit=2, Fdelta=0.01, aveYears=15, ageLimit=100){
  idxF <- fit$conf$keyLogFsta[1,]+1
  barAges <- do.call(":",as.list(fit$conf$fbarRange))+(1-fit$conf$minAge) 
  last.year.used=max(fit$data$years)
  idxno<-which(fit$data$years==last.year.used)
  #dim<-fit.current$stateDim
  #idxN<-1:ncol(stock.mean.weight) 
  F <- exp(fit$pl$logF[idxF,])
  
  sel<-function(){
    Sa<-rep(0,nrow(F))
    K<-0
    for(i in 0:(aveYears-1)){
      thisF<-F[,idxno-i]
      Sa<-Sa+thisF
      K<-K+fbartable(fit)[idxno-i]
    }
    return(Sa/K)
  }

  extend<-function(x,len=100){
    ret<-numeric(len)
    ret[1:length(x)]<-x
    ret[-c(1:length(x))]<-x[length(x)]
    ret
  }

  ave.sl<-sel()
  ave.sw<-colMeans(fit$data$stockMeanWeight[(idxno-aveYears+1):idxno,,drop=FALSE])
  ave.cw<-colMeans(fit$data$catchMeanWeight[(idxno-aveYears+1):(idxno-1),,drop=FALSE])
  ave.pm<-colMeans(fit$data$propMat[(idxno-aveYears+1):idxno,,drop=FALSE])
  ave.nm<-colMeans(fit$data$natMor[(idxno-aveYears+1):idxno,,drop=FALSE])
  ave.lf<-colMeans(fit$data$landFrac[(idxno-aveYears+1):(idxno-1),,drop=FALSE])
  ave.cw.land<-colMeans(fit$data$landMeanWeight[(idxno-aveYears+1):(idxno-1),,drop=FALSE])

  N<-numeric(ageLimit)
  N[1]<-1.0
  M<-extend(ave.nm)
  sw<-extend(ave.sw)
  cw<-extend(ave.cw.land)
  pm<-extend(ave.pm)
  lf<-extend(ave.lf)
 
  deltafirst <- 0.00001
  delta <- Fdelta
  scales<-c(0, deltafirst, seq(0.01, Flimit, by=delta))
  yields<-numeric(length(scales))
  ssbs<-numeric(length(scales))
  for(i in 1:length(scales)){
    scale<-scales[i]
    F<-extend(ave.sl*scale)
    Z<-M+F
    for(a in 2:length(N)){
      N[a]<-N[a-1]*exp(-Z[a-1])  
    }
    C<-F/Z*(1-exp(-Z))*N*lf  
    Y<-sum(C*cw)
    yields[i]<-Y
    ssbs[i]<-sum(N*pm*sw)
  }


  par(mar=c(5.1,4.1,4.1,5.1))
  fbarlab = substitute(bar(F)[X - Y], list(X = fit$conf$fbarRange[1], Y = fit$conf$fbarRange[2]))

  plot(scales, yields, type='l', xlab=fbarlab, ylab='Yield per recruit')

  fmaxidx<-which.max(yields)
  fmax<-scales[fmaxidx]
  deltaY<-diff(yields)
  f01idx<-which.min((deltaY/delta-0.1*deltaY[1]/deltafirst)^2)+1
  f01<-scales[f01idx]
  lines(c(fmax,fmax), c(par('usr')[1],yields[fmaxidx]), lwd=3, col='red')
  lines(c(f01,f01), c(par('usr')[1],yields[f01idx]), lwd=3, col='blue')

  ssbscale<-max(yields)/max(ssbs)

  lines(scales, ssbscale*ssbs, lty='dotted')
  ssbtick<-pretty(ssbs)
  ssbat<-ssbtick*ssbscale
  axis(4,at=ssbat, labels=ssbtick)
  mtext('SSB per recruit', side=4, line=2)

  f35spridx<-which.min((ssbs-0.35*ssbs[1])^2)+1
  f35<-scales[f35spridx]
  lines(c(f35,f35), c(par('usr')[1],ssbs[f35spridx]*ssbscale), lwd=3, col='green')

  title(eval(substitute(expression(F[max]==fmax~ ~ ~ ~ ~F[0.10]==f01~ ~ ~ ~ ~F[0.35*SPR]==f35), 
                        list(fmax=round(fmax,2), f01=round(f01,2), f35=round(f35,2)))))
}
