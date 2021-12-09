##' Plot helper
##' @param fit the fitted object from sam.fit of a set of such fits c(fit1,fit2)
##' @param what quoted name of object to extract
##' @param x x-alues
##' @param ylab label on y-axis
##' @param xlab label on x-axis
##' @param ex extra y's to make room for
##' @param trans function to transform values by
##' @param add logical, plotting is to be added on existing plot
##' @param ci logical, confidence intervals should be plotted
##' @param cicol color to plot the confidence polygon
##' @param addCI A logical vector indicating if confidence intervals should be plotted for the added fits.  
##' @param drop number of years to be left unplotted at the end.
##' @param unnamed.basename the name to assign an unnamed basefit 
##' @param xlim ...
##' @param ... extra arguments transferred to plot
##' @importFrom graphics plot polygon grid lines
##' @importFrom grDevices gray
##' @details The basic plotting used bu many of the plotting functions (e.g. ssbplot, fbarplot ...) 
plotit <-function (fit, what,...){
    UseMethod("plotit")
}
##' @rdname plotit
##' @method plotit sam
##' @export
plotit.sam <- function(fit, what, x=fit$data$years, ylab=what, xlab="Years", ex=numeric(0), trans=function(x)x, add=FALSE, ci=TRUE, cicol=gray(.5,alpha=.5),
                   addCI=NA, drop=0, unnamed.basename="current", xlim=NULL,ylim=NULL,ylimAdd=NA,...){
    idx <- names(fit$sdrep$value)==what
    y <- fit$sdrep$value[idx]
    lowhig <- y+fit$sdrep$sd[idx]%o%c(-2,2)
    didx <- 1:(length(x)-drop)
    if(missing(xlim)){
      xr <- range(x)
    }else{
      xr <- xlim
    }
    x<-x[didx]
    y<-y[didx]
    lowhig<-lowhig[didx,]
    if(add){
      lines(x, trans(y), lwd=3,...)
    }else{
        if(missing(ylim)){
            yr <- range(c(trans(lowhig),0,ex,ylimAdd), na.rm = TRUE)
        }else{
            yr <- ylim
        }
        plot(x, trans(y), xlab=xlab, ylab=ylab, type="n", lwd=3, xlim=xr, ylim=yr, las=1,...)
      grid(col="black")
      lines(x, trans(y), lwd=3, ...)
    }
    if(ci){
      polygon(c(x,rev(x)), y = c(trans(lowhig[,1]),rev(trans(lowhig[,2]))), border = gray(.5,alpha=.5), col = cicol)
      lines(x, trans(y), lwd=3, col=cicol)
      lines(x, trans(y), lwd=2, col="black", lty="dotted")
    }
}
##' @rdname plotit
##' @method plotit samset
##' @export
plotit.samset <- function(fit, what, x=fit$data$years, ylab=what, xlab="Years", ex=numeric(0), trans=function(x)x, add=FALSE, ci=TRUE, cicol=gray(.5,alpha=.5),
                   addCI=rep(FALSE,length(fit)), drop=0, unnamed.basename="current", xlim=NULL,...){
    if(is.logical(addCI) & (length(addCI)==1))addCI=rep(addCI,length(fit))
    #colSet <- c("#FF0000", "#FF7B00", "#FFF500", "#A8FF00", "#14FF00", "#00FFA3", "#00A4FF", "#5100FF", "#CC00FF", "#969696")
    colSet <- c("#332288"  , "#88CCEE"  , "#44AA99"  , "#117733"  , "#999933"  , "#DDCC77"  , "#661100"  , "#CC6677"  , "#882255"  , "#AA4499")  
    idxfrom <- 1
    leg<-names(fit)
    if(is.null(attr(fit,"fit"))){
      attr(fit,"fit") <- fit[[1]]
      idxfrom <- 2
    }else{
      leg<-c(unnamed.basename, leg)
    }
    if(is.null(x))x=attr(fit,"fit")$data$years
    if(missing(xlim)){
      xr <- range(x)
    }else{
      xr <- xlim
    }
    plotit(attr(fit,"fit"), what=what, x=x, ylab=ylab, xlab=xlab, ex=ex, trans=trans, add=add, ci=ci, cicol=cicol, drop=drop, xlim=xr,...)
      d<-lapply(idxfrom:length(fit), function(i)plotit(fit[[i]], what=what, trans=trans, add=TRUE, ci=addCI[i],
                                                        col=colSet[(i-1)%%length(colSet)+1], cicol=paste0(colSet[(i-1)%%length(colSet)+1],"80"), drop=drop, ...)) 
    if(!is.null(names(fit))){
        legend("bottom",legend=leg, lwd=3, col=c(par("col"),colSet[((idxfrom:length(fit))-1)%%length(colSet)+1]), ncol=3, bty="n")
        legend("bottom",legend=leg, lwd=2, col=rep("black", length(leg)), lty="dotted", ncol=3, bty="n")
    }
}
##' @rdname plotit
##' @method plotit samforecast
##' @export
plotit.samforecast <- function(fit, what, x=fit$data$years, ylab=what, xlab="Years", ex=numeric(0), trans=function(x)x, add=FALSE, ci=TRUE, cicol=gray(.5,alpha=.5),
                   addCI=NA, drop=0, unnamed.basename="current", xlim=NULL,ylim=NULL,...){
    xy <- unlist(lapply(fit, function(xx) xx$year))
    thisfit<-attr(fit,"fit")
    xr <- range(thisfit$data$years, xy)
    if(missing(ylim)){
        v1 <- tableit(fit, what = what, trans = trans)
        v2 <- tableit(thisfit, what = what, trans = trans)
        ylim <- range(v1,v2)
    }
    plotit(thisfit, what=what, ylab=ylab, xlab=xlab, ex=ex, trans=trans, add=add, ci=ci, cicol=cicol, drop=drop, xlim=xr,ylim=ylim,...)
}

##' @rdname plotit
##' @method plotit hcr
##' @export
plotit.hcr <- function(fit, what, x=fit$data$years, ylab=what, xlab="Years", ex=numeric(0), trans=function(x)x, add=FALSE, ci=TRUE, cicol=gray(.5,alpha=.5),
                   addCI=NA, drop=0, unnamed.basename="current", xlim=NULL,...){
    plotit.samforecast(fit$forecast, what = what, x = x, ylab = ylab, xlab = xlab, ex = ex, trans = trans, add = add, ci = ci, cicol = cicol, addCI = addCI, drop = drop, unnamed.basename = unnamed.basename, xlim = xlim, ...)
}


##' SAM add forecasts 
##' @param fit the object returned from sam.fit
##' @param what what to plot
##' @param dotcol color for dot
##' @param dotpch pch for dot
##' @param dotcex cex for dot
##' @param intervalcol color for interval
##' @param ... extra arguments not currently used
##' @details internal plotting fun
##' @importFrom graphics arrows
addforecast<-function(fit, what, dotcol="black", dotpch=19, dotcex=1.5, intervalcol=gray(.5,alpha=.5),...){
    UseMethod("addforecast")
}
##' @rdname addforecast
##' @method addforecast samforecast
##' @export
addforecast.samforecast <- function(fit, what, dotcol="black", dotpch=19, dotcex=1.5, intervalcol=gray(.5,alpha=.5),...){
    x <- attr(fit,"tab")
    y <- as.numeric(rownames(x))
    dummy <- sapply(1:length(y),
                    function(i){
                        xx<-c(x[i,paste(what,"low", sep=":")],x[i,paste(what,"high", sep=":")]);
                        units = par(c('usr', 'pin'))
                        xx_to_inches = with(units, pin[2L]/diff(usr[3:4]))
                        if(abs(xx_to_inches * diff(xx))>0.01){
                            arrows(y[i],xx[1],y[i],xx[2],lwd=3, col=intervalcol, angle=90, code=3, length=.1)
                        }
                    }
                    )
    points(y,x[,paste(what,"median", sep=":")], pch=dotpch, cex=dotcex, col=dotcol)
}

addforecast.hcr <- function(fit, what, dotcol="black", dotpch=19, dotcex=1.5, intervalcol=gray(.5,alpha=.5),...){
    addforecast(fit$forecast, what = what, dotcol = dotcol, dotpch = dotpch, dotcex = dotcex, intervalcol = intervalcol, ...)
}


##' Plot by one or two  
##' @param x numeric vector of points to be plotted
##' @param y numeric vector of points to be plotted
##' @param z numeric vector of points to be plotted
##' @param x.line numeric vector of points of line to be added
##' @param y.line numeric vector of points of line to be added
##' @param z.line numeric vector of points of line to be added
##' @param by vector or two column matrix to create sub sets from
##' @param bubblescale scaling of bubble size
##' @param x.common logical: use same x-axis for all plots
##' @param y.common logical: use same y-axis for all plots
##' @param z.common logical: use same z-axis for all plots
##' @param xlab normal graphical parameter
##' @param ylab normal graphical parameter
##' @param xlim normal graphical parameter
##' @param ylim normal graphical parameter
##' @param zmax internally used to scale bubbles similarly 
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
plotby <-function(x=NULL, y=NULL, z=NULL, x.line=NULL, y.line=NULL, z.line=NULL, by=NULL,  bubblescale=1, x.common=!is.null(x), y.common=!is.null(y), z.common=!is.null(z), xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, zmax=NULL, axes=TRUE, ...){
  if(is.null(by)){
    # 0
    if(length(x)==0 & length(y)==0 & length(z)==0){
      if(length(xlim)==0) xlim <- 0:1
      if(length(ylim)==0) ylim <- 0:1
      plot(NA, xlim=xlim, ylim=ylim, xlab="", ylab="", axes=FALSE, ...)
    }
    # 1
    if(length(x)>0 & length(y)==0 & length(z)==0){
      if(missing(ylab)) ylab <- deparse(substitute(x))
      if(missing(ylim)) ylim <- c(min(x)-1, max(x)+1)      
      plot(x, ylab=ylab, ylim=ylim, axes=axes, ...)
      if(length(x.line)>0){
        lines(x.line,...)    
      }
    }
    if(length(x)==0 & length(y)>0 & length(z)==0){
      if(missing(ylab)) ylab <- deparse(substitute(y))
      if(missing(ylim)) ylim <- c(min(y)-1, max(y)+1)      
      plot(y, ylab=ylab, ylim=ylim, axes=axes, ...)
      if(length(y.line)>0){
        lines(y.line, ...)    
      }
    }
    if(length(x)==0 & length(y)==0 & length(z)>0){
      if(missing(ylab)) ylab <- deparse(substitute(z))
      if(missing(ylim)) ylim <- c(min(z)-1, max(z)+1)      
      plot(z, ylab=ylab, ylim=ylim, axes=axes, ...)
      if(length(z.line)>0){
        lines(z.line, ...)    
      }
    }
    # 2
    if(length(x)>0 & length(y)>0 & length(z)==0){
      if(missing(xlab)) xlab <- deparse(substitute(x))
      if(missing(ylab)) ylab <- deparse(substitute(y))
      if(missing(xlim)) xlim <- c(min(x)-1, max(x)+1)
      if(missing(ylim)) ylim <- c(min(y)-1, max(y)+1)
      plot(x, y, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, axes=axes, ...)
      if(length(y.line)>0){
        o<-order(x)  
        lines(x[o], y.line[o],...)    
      }
    }
    if(length(x)>0 & length(y)==0 & length(z)>0){
      if(missing(xlab)) xlab <- deparse(substitute(x))
      if(missing(ylab)) ylab <- deparse(substitute(z))
      if(missing(xlim)) xlim <- c(min(x)-1, max(x)+1)
      if(missing(ylim)) ylim <- c(min(z)-1, max(z)+1)
      plot(x, z, xlab=xlab,  ylab=ylab, xlim=xlim, ylim=ylim, axes=axes, ...)
      if(length(z.line)>0){
        o<-order(x)  
        lines(x[o], z.line[o],...)    
      }
    }
    if(length(x)==0 & length(y)>0 & length(z)>0){
      if(missing(xlab)) xlab <- deparse(substitute(y))
      if(missing(ylab)) ylab <- deparse(substitute(z))
      if(missing(xlim)) xlim <- c(min(y)-1, max(y)+1)
      if(missing(ylim)) ylim <- c(min(z)-1, max(z)+1)
      plot(y, z, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, axes=axes, ...)
      if(length(z.line)>0){
        o<-order(y)  
        lines(y[o], z.line[o],...)    
      }
    }
    # 3
    if(length(x)>0 & length(y)>0 & length(z)>0){
      if(missing(xlab)) xlab <- deparse(substitute(x))
      if(missing(ylab)) ylab <- deparse(substitute(y))
      if(missing(xlim)) xlim <- c(min(x)-1, max(x)+1)
      if(missing(ylim)) ylim <- c(min(y)-1, max(y)+1)
      if(is.null(zmax)) zmax <- max(sqrt(abs(z)), na.rm=TRUE)
      cex=sqrt(abs(z))/zmax*5*bubblescale
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
    if(z.common) zmax <- max(sqrt(abs(z)), na.rm=TRUE)
    if(x.common) op1<-par(oma=c(par("mar")[1], par("oma")[2], par("mar")[3], par("oma")[4]),
                          mar=c(.2,par("mar")[2],.2,par("mar")[4]))
    if(y.common) op2<-par(oma=c(par("oma")[1], par("mar")[2], par("oma")[3], par("mar")[4]),
                          mar=c(par("mar")[1],.2,par("mar")[3],.2))
    
    .fun<-function(i){
      b<-uby[i]
      plotby(x[by==b], y[by==b], z[by==b], x.line[by==b], y.line[by==b], z.line[by==b], xlab=ifelse(x.common,"",xlab), ylab=ifelse(y.common,"",ylab),
             xlim=xlim, ylim=ylim, zmax=zmax, bubblescale=bubblescale, axes=FALSE, ...)
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
    par(mfrow=c(1,1))
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
    if(z.common) zmax <- max(sqrt(abs(z)), na.rm=TRUE)
    if(x.common) op1<-par(oma=c(par("mar")[1], par("oma")[2], par("mar")[3], par("oma")[4]),
                          mar=c(.2,par("mar")[2],.2,par("mar")[4]))
    if(y.common) op2<-par(oma=c(par("oma")[1], par("mar")[2], par("oma")[3], par("mar")[4]),
                          mar=c(par("mar")[1],.2,par("mar")[3],.2))
    
    fun<-function(i,j){
      b<-(by[,1]==uby1[i])&(by[,2]==uby2[j])
      plotby(x[b], y[b], z[b], x.line[b], y.line[b], z.line[b], xlab=ifelse(x.common,"",xlab), ylab=ifelse(y.common,"",ylab),
             xlim=xlim, ylim=ylim, zmax=zmax, bubblescale=bubblescale, axes=FALSE, ...)
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
    par(mfrow=c(1,1))
  }  
}

##' SAM Fbar plot 
##' @param fit the object returned from sam.fit 
##' @param partial true if included partial F's are to be plotted
##' @param drop number of years to be left unplotted at the end. Default (NULL) is to not show years at the end with no catch information
##' @param pcol color of partial lines
##' @param page partial ages to plot
##' @param ... extra arguments transferred to plot including the following: \cr
##' \code{add} logical, plotting is to be added on existing plot \cr
##' \code{ci} logical, confidence intervals should be plotted \cr
##' \code{cicol} color to plot the confidence polygon
##' @importFrom graphics matplot
##' @details Plot the defined fbar. 
##' @export
fbarplot<-function(fit,...){
    UseMethod("fbarplot")
}
##' @rdname fbarplot
##' @method fbarplot sam
##' @param plot true if fbar should be plotted
##' @export
fbarplot.sam <- function(fit,partial = TRUE, drop=NULL, pcol="lightblue", page=NULL, plot = TRUE,...){
     if(is.null(drop)){
        drop=max(fit$data$aux[,"year"])-max(fit$data$aux[fit$data$aux[,"fleet"]==1,"year"])
    }  
    fbarRange<-fit$conf$fbarRange
    fbarlab=substitute(bar(F)[X-Y],list(X=fbarRange[1],Y=fbarRange[2]))
    fmat<-t(faytable(fit))#fitlocal$pl$logF[fitlocal$conf$keyLogFsta[1,]+1,]
    if(is.null(page)){
        page<-fbarRange[1]:fbarRange[2]
    }
    idx <- which(fit$conf$minAge:fit$conf$maxAge %in% page)
    exx <- if(partial){fmat[idx,]}else{numeric(0)}
    if(plot){
        plotit(fit, "logfbar", ylab=fbarlab, trans=exp, ex=exx, drop=drop, ...)
        if(partial){
            idxx <- 1:(length(fit$data$years)-drop)
            matplot(fit$data$years[idxx], t(fmat[idx,idxx]), add=TRUE, type="b", col=pcol, pch=as.character(page))
        }
    }
    invisible(list(drop=drop,fbarRange=fbarRange,fbarlab=fbarlab,fmat=fmat,page=page,idx=idx,exx=exx))
}
##' @rdname fbarplot
##' @method fbarplot samset
##' @export
fbarplot.samset <- function(fit,partial = FALSE, drop=NULL, pcol="lightblue", page=NULL,...){
    if(!is.null(attr(fit,"fit"))){
        fitlocal <- attr(fit,"fit")
    }else{
        fitlocal <- fit[[1]]
    }
    tmp <- fbarplot(fitlocal,partial,drop,pcol,page,plot=FALSE,...)
    plotit(fit, "logfbar", ylab=tmp$fbarlab, trans=exp, ex=tmp$exx, drop=tmp$drop, ...)
    if(partial){
        idxx <- 1:(length(fitlocal$data$years)-tmp$drop)
        matplot(fitlocal$data$years[idxx], t(tmp$fmat[tmp$idx,idxx]), add=TRUE, type="b", col=pcol, pch=as.character(tmp$page))
    }
}
##' @rdname fbarplot
##' @method fbarplot samforecast
##' @export
fbarplot.samforecast <- function(fit,partial = FALSE, drop=NULL, pcol="lightblue", page=NULL,...){
    fitlocal <- attr(fit,"fit")
    tmp <- fbarplot(fitlocal,partial,drop,pcol,page,plot=FALSE,...)
    plotit(fit, "logfbar", ylab=tmp$fbarlab, trans=exp, ex=tmp$exx, drop=tmp$drop, ...)
    if(partial){
        idxx <- 1:(length(fitlocal$data$years)-tmp$drop)
        matplot(fitlocal$data$years[idxx], t(tmp$fmat[tmp$idx,idxx]), add=TRUE, type="b", col=pcol, pch=as.character(tmp$page))
    }
    addforecast(fit,"fbar")
}

##' @rdname fbarplot
##' @method fbarplot hcr
##' @export
fbarplot.hcr <- function(fit,partial = FALSE, drop=NULL, pcol="lightblue", page=NULL,...){
    fitlocal <- attr(fit,"fit")
    tmp <- fbarplot(fitlocal,partial,drop,pcol,page,plot=FALSE,...)
    plotit(fit, "logfbar", ylab=tmp$fbarlab, trans=exp, ex=tmp$exx, drop=tmp$drop, ...)
    if(partial){
        idxx <- 1:(length(fitlocal$data$years)-tmp$drop)
        matplot(fitlocal$data$years[idxx], t(tmp$fmat[tmp$idx,idxx]), add=TRUE, type="b", col=pcol, pch=as.character(tmp$page))
    }
    addforecast(fit$forecast,"fbar")
}



##' SAM F-selectivity plot 
##' @param fit An object returned from sam.fit 
##' @param cexAge cex variable giving the size of the age numbers 
##' @param ... extra arguments transferred to barplot and text
##' @importFrom graphics barplot text
##' @details Plots selectivity in F. 
##' @export
fselectivityplot<-function(fit, cexAge = 1,...){
  UseMethod("fselectivityplot")
}
##' @rdname fselectivityplot
##' @method fselectivityplot sam
##' @export
fselectivityplot.sam <- function(fit, cexAge = 1,...){
    fmat<- faytable(fit)
    fmat[is.na(fmat)] <- 0
    barplot(t(fmat/rowSums(fmat)),border=NA,space=c(0),xlab="Year", main = "Selectivity in F", ...)
    text(1,cumsum(t(fmat/rowSums(fmat))[,1]) - 0.5*t(fmat/rowSums(fmat))[,1] ,label=as.character(1:ncol(fmat)),adj=c(0.0,0.2), cex = cexAge)
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
    UseMethod("ssbplot")
}
##' @rdname ssbplot
##' @method ssbplot default
##' @export
ssbplot.default <- function(fit,...){
    plotit(fit, "logssb", ylab="SSB", trans=exp,...)
}
##' @rdname ssbplot
##' @method ssbplot samforecast
##' @export
ssbplot.samforecast <- function(fit,...){
    plotit(fit, "logssb", ylab="SSB", trans=exp,...)
    addforecast(fit,"ssb")
}

##' @rdname ssbplot
##' @method ssbplot hcr
##' @export
ssbplot.hcr <- function(fit,...){
    plotit(fit, "logssb", ylab="SSB", trans=exp,...)
    addforecast(fit,"ssb")
}


##' SAM life expectancy plot 
##' @param fit the object returned from sam.fit
##' @param atRecruit If true, show life expectancy given survival until minAge, otherwise show life expectancy at birth
##' @param ... extra arguments transferred to plot including the following: \cr
##' \code{add} logical, plotting is to be added on existing plot \cr
##' \code{ci} logical, confidence intervals should be plotted \cr
##' \code{cicol} color to plot the confidence polygon
##' @details Plot of life expectancy 
##' @export
lifeexpectancyplot<-function(fit, atRecruit = TRUE, ...){
    UseMethod("lifeexpectancyplot")
}
##' @rdname lifeexpectancyplot
##' @method lifeexpectancyplot default
##' @export
lifeexpectancyplot.default <- function(fit, atRecruit = TRUE, ylimAdd = fit$conf$maxAge, ...){
    if(atRecruit){
        plotit(fit, "logLifeExpectancyRec", ylab="Life expectancy at recruitment", xlab="Year", trans=exp, ylimAdd = ylimAdd, ...)
    }else{
        plotit(fit, "logLifeExpectancy", ylab="Life expectancy at birth", xlab="Year", trans=exp, ylimAdd = ylimAdd, ...)
    }
    abline(h = c(fit$conf$maxAge, fit$conf$minAge), col = "darkgrey",lwd=3, lty = 4)
}
##' @rdname lifeexpectancyplot
##' @method lifeexpectancyplot samforecast
##' @export
lifeexpectancyplot.samforecast <- function(fit, atRecruit = TRUE, ylimAdd = fit$conf$maxAge,...){
    if(atRecruit){
        plotit(fit, "logLifeExpectancyRec", ylab="Life expectancy at recruitment", xlab="Year", trans=exp, ylimAdd = ylimAdd, ...)
        addforecast(fit,"logLifeExpectancyRec")
    }else{
        plotit(fit, "logLifeExpectancy", ylab="Life expectancy at birth", xlab="Cohort", trans=exp, ylimAdd = ylimAdd, ...)
        addforecast(fit,"logLifeExpectancy")
    }
    abline(h = c(fit$conf$maxAge, fit$conf$minAge), col = "darkgrey",lwd=3, lty = 4)
}

##' @rdname lifeexpectancyplot
##' @method lifeexpectancyplot hcr
##' @export
lifeexpectancyplot.hcr <- function(fit, atRecruit = TRUE, ylimAdd = fit$conf$maxAge, ...){
    if(atRecruit){
        plotit(fit, "logLifeExpectancyRec", ylab="Life expectancy at recruitment", xlab="Year", trans=exp, ylimAdd = ylimAdd, ...)
        addforecast(fit,"logLifeExpectancyRec")
    }else{
        plotit(fit, "logLifeExpectancy", ylab="Life expectancy at birth", xlab="Cohort", trans=exp, ylimAdd = ylimAdd, ...)
        addforecast(fit,"logLifeExpectancy")
    }
    abline(h = c(fit$conf$maxAge, fit$conf$minAge), col = "darkgrey",lwd=3, lty = 4)
}


##' SAM years lost to fishing plot 
##' @param fit the object returned from sam.fit
##' @param cause Fisning, Other, or LifeExpectancy
##' ##' @param ... extra arguments transferred to plot including the following: \cr
##' \code{add} logical, plotting is to be added on existing plot \cr
##' \code{ci} logical, confidence intervals should be plotted \cr
##' \code{cicol} color to plot the confidence polygon
##' @details Plot of years lost to fishing
##' @export
yearslostplot<-function(fit,cause, ...){
    UseMethod("yearslostplot")
}
##' @rdname yearslostplot
##' @method yearslostplot default
##' @export
yearslostplot.default<-function(fit, cause=c("Fishing","Other","LifeExpectancy"), ...){
    cv <- match.arg(cause)
    if(cv == "Fishing"){
        what <- "logYLTF"
        lab <- sprintf("Life years lost to fishing between age %d and %d",fit$conf$minAge,fit$conf$maxAge)
    }else if(cv == "Other"){
        what <- "logYLTM"
        lab <- sprintf("Life years lost to other causes between age %d and %d",fit$conf$minAge,fit$conf$maxAge)
    }else{
        what <- "logYNL"
        lab <- sprintf("Temporary life expectancy between age %d and %d",fit$conf$minAge,fit$conf$maxAge)
    }
    plotit(fit, what, ylab=lab, xlab="Year", trans=exp,...)
}
##' @rdname yearslostplot
##' @method yearslostplot samforecast
##' @export
yearslostplot.samforecast <- function(fit,cause=c("Fishing","Other","LifeExpectancy"), ...){
   cause <- match.arg(cause)
    if(case == "Fishing"){
        what <- "logYLTF"
        lab <- sprintf("Life years lost to fishing between age %d and %d",fit$conf$minAge,fit$conf$maxAge)
    }else if(case == "Other"){
        what <- "logYLTM"
        lab <- sprintf("Life years lost to other causes between age %d and %d",fit$conf$minAge,fit$conf$maxAge)
    }else{
        what <- "logYNL"
        lab <- sprintf("Temporary life expectancy between age %d and %d",fit$conf$minAge,fit$conf$maxAge)
    }
    plotit(fit, what, ylab=lab, xlab="Year", trans=exp,...)
    addforecast(fit,what)
}

##' @rdname yearslostplot
##' @method yearslostplot hcr
##' @export
yearslostplot.hcr <- function(fit,cause=c("Fishing","Other","LifeExpectancy"), ...){
   cause <- match.arg(cause)
    if(case == "Fishing"){
        what <- "logYLTF"
        lab <- sprintf("Life years lost to fishing between age %d and %d",fit$conf$minAge,fit$conf$maxAge)
    }else if(case == "Other"){
        what <- "logYLTM"
        lab <- sprintf("Life years lost to other causes between age %d and %d",fit$conf$minAge,fit$conf$maxAge)
    }else{
        what <- "logYNL"
        lab <- sprintf("Temporary life expectancy between age %d and %d",fit$conf$minAge,fit$conf$maxAge)
    }
    plotit(fit, what, ylab=lab, xlab="Year", trans=exp,...)
    addforecast(fit,what)
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
    UseMethod("tsbplot")
}
##' @rdname tsbplot
##' @method tsbplot default
##' @export
tsbplot.default <- function(fit,...){
    plotit(fit, "logtsb", ylab="TSB", trans=exp,...)
}

##' SAM Recruits plot 
##' @param fit the object returned from sam.fit
##' @param lagR use the age after the youngest as R 
##' @param ... extra arguments transferred to plot including the following: \cr
##' \code{add} logical, plotting is to be added on existing plot \cr
##' \code{ci} logical, confidence intervals should be plotted \cr
##' \code{cicol} color to plot the confidence polygon
##' @details Plot of numbers of recruits (youngest age class)
##' @export
recplot<-function(fit, lagR=FALSE, ...){
    UseMethod("recplot")
}
##' @rdname recplot
##' @method recplot sam
##' @export
recplot.sam <- function(fit, lagR=FALSE, ...){
    if(!lagR){
      lab<-paste("Recruits (age ", fit$conf$minAge, ")", sep="")
      plotit(fit, "logR", ylab=lab, trans=exp,...)
    }else{
      lab<-paste("Recruits (age ", fit$conf$minAge+1, ")", sep="")
      plotit(fit, "logLagR", ylab=lab, trans=exp,...)
    }
}
##' @rdname recplot
##' @method recplot samset
##' @export
recplot.samset <- function(fit, lagR=FALSE, ...){
    if(!is.null(attr(fit,"fit"))){
        fitlocal <- attr(fit,"fit")
    }else{
        fitlocal <- fit[[1]]
    }
    if(!lagR){
      lab<-paste("Recruits (age ", fitlocal$conf$minAge, ")", sep="")
      plotit(fit, "logR", ylab=lab, trans=exp,...)
    }else{
      lab<-paste("Recruits (age ", fitlocal$conf$minAge+1, ")", sep="")
      plotit(fit, "logLagR", ylab=lab, trans=exp,...)
    }
}
##' @rdname recplot
##' @method recplot samforecast
##' @export
recplot.samforecast <- function(fit, lagR=FALSE, ...){
    fitlocal <- attr(fit,"fit")
    if(!lagR){
      lab<-paste("Recruits (age ", fitlocal$conf$minAge, ")", sep="")
      plotit(fit, "logR", ylab=lab, trans=exp,...)
      addforecast(fit, "rec")
    }else{
      lab<-paste("Recruits (age ", fitlocal$conf$minAge+1, ")", sep="")
      plotit(fit, "logLagR", ylab=lab, trans=exp,...)
      addforecast(fit, "rec")
    }
}

##' @rdname recplot
##' @method recplot hcr
##' @export
recplot.hcr <- function(fit, lagR=FALSE, ...){
    fitlocal <- attr(fit,"fit")
    if(!lagR){
      lab<-paste("Recruits (age ", fitlocal$conf$minAge, ")", sep="")
      plotit(fit, "logR", ylab=lab, trans=exp,...)
      addforecast(fit, "rec")
    }else{
      lab<-paste("Recruits (age ", fitlocal$conf$minAge+1, ")", sep="")
      plotit(fit, "logLagR", ylab=lab, trans=exp,...)
      addforecast(fit, "rec")
    }
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
    UseMethod("catchplot")
}
##' @rdname catchplot
##' @param plot true if catch should be plotted
##' @method catchplot sam
##' @export
catchplot.sam <- function(fit, obs.show=TRUE, drop=NULL,plot=TRUE,...){
    if(is.null(drop)){
        drop=max(fit$data$aux[,"year"])-max(fit$data$aux[fit$data$aux[,"fleet"]==1,"year"])
    }
    CW <- fit$data$catchMeanWeight
    CW <- CW[apply(!is.na(CW),1,all),]
    x <- as.numeric(rownames(CW))
    obs <- NULL
    if(plot)
        plotit(fit, "logCatch", ylab="Catch", trans=exp, drop=drop,...)
    if(obs.show){
        aux <- fit$data$aux
        logobs <- fit$data$logobs
        .goget <- function(y,a){
            ret <- exp(logobs[aux[,"fleet"]==1 & aux[,"year"]==y & aux[,"age"]==a])
            ifelse(length(ret)==0,0,ret)
        }
        if(plot)
            points(x, rowSums(outer(rownames(CW), colnames(CW), Vectorize(.goget))*CW, na.rm=TRUE), pch=4, lwd=2, cex=1.2)
        obs <- list(x=x,y=rowSums(outer(rownames(CW), colnames(CW), Vectorize(.goget))*CW, na.rm=TRUE))
    }
    invisible(list(drop=drop,obs=obs))
}
##' @rdname catchplot
##' @method catchplot samset
##' @export
catchplot.samset <- function(fit, obs.show=TRUE, drop=NULL,...){
    if(!is.null(attr(fit,"fit"))){
        fitlocal <- attr(fit,"fit")
    }else{
        fitlocal <- fit[[1]]
    }
    tmp <- catchplot(fitlocal,obs.show,drop,plot=FALSE,...)
    plotit(fit, "logCatch", ylab="Catch", trans=exp, drop=tmp$drop,...)
    if(obs.show){
        points(tmp$obs$x, tmp$obs$y, pch=4, lwd=2, cex=1.2)
    }
}
##' @rdname catchplot
##' @method catchplot samforecast
##' @export
catchplot.samforecast <- function(fit, obs.show=TRUE, drop=NULL,...){
    fitlocal <- attr(fit,"fit")
    tmp <- catchplot(fitlocal,obs.show,drop,plot=FALSE,...)
    plotit(fit, "logCatch", ylab="Catch", trans=exp, drop=tmp$drop,...)
    if(obs.show){
        points(tmp$obs$x, tmp$obs$y, pch=4, lwd=2, cex=1.2)
    }
    addforecast(fit, "catch")
}
##' @rdname catchplot
##' @method catchplot hcr
##' @export
catchplot.hcr <- function(fit, obs.show=TRUE, drop=NULL,...){
    fitlocal <- attr(fit,"fit")
    tmp <- catchplot(fitlocal,obs.show,drop,plot=FALSE,...)
    plotit(fit, "logCatch", ylab="Catch", trans=exp, drop=tmp$drop,...)
    if(obs.show){
        points(tmp$obs$x, tmp$obs$y, pch=4, lwd=2, cex=1.2)
    }
    addforecast(fit, "catch")
}

##' SAM parameter plot 
##' @param fit the object returned from sam.fit
##' @param cor.report.limit correlations with absolute value > this number is reported in the plot 
##' @param ... extra arguments transferred to plot
##' @details Plot of all estimated model parameters (fixed effects). Shown with confidence interval.  
##' @export
##' @importFrom stats cov2cor
parplot<-function(fit, cor.report.limit=0.95, ...){
    UseMethod("parplot")
}
##' @rdname parplot
##' @method parplot sam
##' @export
parplot.sam <- function(fit, cor.report.limit=0.95, ...){
    parplot(c(fit),cor.report.limit,...)
}
##' @rdname parplot
##' @method parplot samset
##' @export
parplot.samset <- function(fit, cor.report.limit=0.95, ...){
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
##' @param ... extra arguments not currently used
##' @return a list of matrices
##' @export
obscov<-function(fit, corr=FALSE,...){
    UseMethod("obscov")
}
##' @rdname obscov
##' @method obscov sam
##' @export
obscov.sam<-function(fit, corr=FALSE,...){
    res<-fit$rep$obsCov
    for(i in 1:length(res)) rownames(res[[i]])<-fit$data$minAgePerFleet[i]:(fit$data$minAgePerFleet[i]+nrow(res[[i]])-1)
    if(corr) for(i in 1:length(res)) if(any(is.na(res[[i]])))res[[i]][]<-NA else res[[i]]<-cov2cor(res[[i]])
    res
}
##' @rdname obscov
##' @method obscov samset
##' @export
obscov.samset <- function(fit, corr=FALSE,...){
    return(lapply(fit,obscov))
}

##' Plots the estimated correlation matrices by fleet.
##' @param fit the object returned from sam.fit
##' @param ... extra arguments to plot
##' @export
obscorrplot<-function(fit,...){
    UseMethod("obscorrplot")
}
##' @rdname obscorrplot
##' @method obscorrplot sam
##' @export
obscorrplot.sam <- function(fit,...){
    x <- obscov(fit,TRUE)
    fn <- attr(fit$data,"fleetNames")
    for(i in 1:length(x)){
        xx <- x[[i]]
        ages <- fit$data$minAgePerFleet[i]:fit$data$maxAgePerFleet[i]
        if( fit$conf$obsLikelihoodFlag[i]=="ALN" ) ages<-ages[ -length(ages) ]
        rownames(xx) <- ages
        colnames(xx) <- ages
        x[[i]] <- xx
    }
    corplotcommon(x,fn,...)
}

##' Plots the residual between-age correlation matrices by fleet.
##' @param res the object returned from residuals.sam
##' @param ... extra arguments to plot
##' @importFrom stats cor xtabs
##' @export
empirobscorrplot<-function(res,...){
    UseMethod("empirobscorrplot")
}
##' @rdname empirobscorrplot
##' @method empirobscorrplot samres
##' @export
empirobscorrplot.samres <- function(res,...){
    dat <- data.frame(resid=res$residual,age=res$age,year=res$year,fleet=res$fleet)
    fleets <- unique( dat$fleet)
    fn <- attr(res,"fleetNames")
    x <- list()
    for(i in 1:length(fleets)){
        tmp <- xtabs( resid ~ age + year, data=dat[dat$fleet==fleets[i],])
        xx <- cor( t(tmp) )
        x[[ length(x) + 1 ]] <- xx    
    }
    corplotcommon(x,fn,...)
}
##' Common function for plotting correlation matrices.
##' @param x a list of correlation matrices
##' @param fn a vector of fleet names
##' @param ... extra arguments to plotcorr
##' @importFrom ellipse plotcorr
##' @export
corplotcommon<-function(x,fn,...){
    op <- par(no.readonly=TRUE)
    ccolors <- c("#A50F15","#DE2D26","#FB6A4A","#FCAE91","#FEE5D9","white",
                 "#EFF3FF","#BDD7E7","#6BAED6","#3182BD","#08519C")

    if(length(x)==3){
      div <- c(3,1)
    }else{    
      div <- rep(ceiling(sqrt(length(x))),2)
      if(div[1]*(div[2]-1)>=length(x))div[2] <- div[2]-1
    }
    laym <- matrix(1:(div[1]*div[2]),nrow=div[1], ncol=div[2])
    layout(laym)

    for(i in 1:length(x)){
        xx <- x[[i]]
        if(!any(is.na(xx))){
          plotcorr(xx,col=ccolors[5*xx+6],mar=0.1+c(2,2,2,2), main=substr(fn[i], 1, 20),...)
        }
    }
    par(op)
}

##' Plots between-age correlations by fleet, either estimated or empirical using residuals.
##' @param x Either a sam fit as returned by sam.fit OR the object returned from residuals.sam
##' @param ... extra arguments to plot
##' @importFrom ellipse plotcorr 
##' @export
corplot<-function(x,...){
    UseMethod("corplot")
}
##' @rdname corplot
##' @method corplot sam
##' @export
corplot.sam <- function(x,...){
    obscorrplot(x,...)
}
##' @rdname corplot
##' @method corplot samres
##' @export
corplot.samres <- function(x,...){
    empirobscorrplot(x,...)
}
 
##' Plots the stock recruitment 
##' @param fit the object returned from sam.fit
##' @param ... extra arguments to plot
##' @importFrom graphics text
##' @export
srplot<-function(fit, ...){
    UseMethod("srplot")
}
##' @rdname srplot
##' @method srplot sam
##' @param textcol color of years on plot
##' @param years the plotting symbols are the years
##' @param linetype type for the plot (default line)
##' @param linecol color of lines between points
##' @param xlim bounds for x-axis
##' @param ylim bounds for y-axis
##' @param add false if a new plot should be created
##' @export
srplot.sam <- function(fit, textcol="red", years=TRUE, linetype="l", linecol="black", xlim, ylim, add=FALSE, ...){
  X <- summary(fit)
  n<-nrow(X)
  lag <- fit$conf$minAge
  idxR <- (lag+1):n
  idxS <- 1:(n-lag)
  R<-X[idxR,1]
  S<-X[idxS,4]
  Rnam<-colnames(X)[1]
  Snam<-colnames(X)[4]
  y<-rownames(X)
  if(add){
    lines(S,R)
  }else{
    if (missing(xlim)) xlim=range(0,S)
    if (missing(ylim)) ylim=range(0,R)
    plot(S,R, xlab=Snam, ylab=Rnam, type=linetype, col=linecol, xlim=xlim, ylim=ylim, ...)
  }
  if (years) text(S,R, labels=y[idxR], cex=.7, col=textcol )
}

##' Plots fit to data 
##' @param fit the object returned from sam.fit
##' @param log should the plot be against log-obs
##' @param fleets an integer vector of fleets to plot. Default is all of them
##' @param ... extra arguments to plot
##' @export
fitplot <- function(fit, log=TRUE, ...){
    UseMethod("fitplot")
}
##' @rdname fitplot
##' @method fitplot sam
##' @export
fitplot.sam <- function(fit, log=TRUE,fleets=unique(fit$data$aux[,"fleet"]), ...){  
  idx<-fit$data$aux[,"fleet"]%in%fleets  
  trans <- function(x)if(log){x}else{exp(x)}  
  p <- trans(fit$obj$report(c(fit$sdrep$par.fixed,fit$sdrep$par.random))$predObs[idx])
  o <- trans(fit$data$logobs[idx])
  aa <- fit$data$aux[idx,"age"]
  neg.age <- (aa < -1.0e-6)
  aa[neg.age] <- NA
  a <- paste0("a=",aa," ")
  f <- paste0(" f=",strtrim(attr(fit$data,"fleetNames")[fit$data$aux[idx,"fleet"]],50))
  Year <- fit$data$aux[idx,"year"]
  if(length(fleets)==1){
    myby <- paste(a, ":", f)
  }else{
    myby <- cbind(a,f)
  }
  plotby(Year, o, y.line=p, by=myby, y.common=FALSE, ylab="", ...)
}


##' plot survey catchabilities
##' @param qt An object of class 'samqtable' as returned from qtable
##' @param exp if true return on natural scale rather than log
##' @importFrom grDevices n2mfrow
##' @export
qtableplot<-function(qt, exp=FALSE){
    UseMethod("qtableplot")
}
##' plot survey catchabilities
##' @rdname qtableplot
##' @method qtableplot samqtable
##' @export
qtableplot.samqtable<-function(qt,exp=FALSE){
    sds<-attr(qt,"sd")
    hi <- qt + 2*sds
    lo <- qt - 2*sds
    ylabel <- ifelse(exp,"Q", "logQ")
    if(exp == TRUE) { qt<-exp(qt); hi<-exp(hi); lo<-exp(lo) }
    op<-par(mfrow=n2mfrow(nrow(qt)))
    on.exit(par(op))
    for(f in 1:nrow(qt)){
        yl <- range(rbind(lo[f,]-0.15*abs(lo[f,]),hi[f,]+0.15*abs(hi[f,])),na.rm=TRUE)
        plot(as.numeric(colnames(qt)),qt[f,],main=rownames(qt)[f],type="b",ylim=yl,ylab=ylabel,xlab="Age")
        arrows(as.numeric(colnames(qt)),lo[f,],y1=hi[f,],angle=90,code=3,length=0.1)
    }
}


##' SAM Data plot 
##' @param fit the object returned from sam.fit
##' @param col color to use for each fleet, default is two sequential colors \cr
##' @param fleet_type character vector giving the type of data per fleet. The default uses fit$data$fleetTypes as follows: \cr
##' \code{fit$data$fleetTypes==0} "Catch at age" \cr
##' \code{fit$data$fleetTypes==1} "Catch at age with effort" \cr
##' \code{fit$data$fleetTypes==2 or 6} "Index at age" \cr
##' \code{fit$data$fleetTypes==3} "Biomass or catch index" \cr
##' \code{fit$data$fleetTypes==5} "Tagging data" \cr
##' \code{fit$data$fleetTypes==7} "Sum of fleets"
##' @param fleet_names character vector giving fleet names. The default is given by attr(fit$data,"fleetNames")
##' @details Plot data available for the stock 
##' @export
dataplot<-function(fit, col=NULL, fleet_type=NULL, fleet_names=NULL){
  UseMethod("dataplot")
}
##' @rdname dataplot
##' @method dataplot sam
##' @export
dataplot.sam <- function(fit, col=NULL, fleet_type=NULL, fleet_names=NULL){
  years <- fit$data$years
  nf <- fit$data$noFleets
  for (k in 1:nf){ # Remove -1 for weight indices
    if(fit$data$minAgePerFleet[k]==-1) fit$data$minAgePerFleet[k]<- NA
    if(fit$data$maxAgePerFleet[k]==-1) fit$data$maxAgePerFleet[k]<- NA
  }
  yspace <- 2 # space between fleets
  noage <- length(min(fit$data$minAgePerFleet, na.rm = TRUE):max(fit$data$maxAgePerFleet, na.rm = TRUE))
  
  dat<-cbind(fit$data$aux, fit$data$logobs)
  
  if (missing(col)) col <- c("#67a9cf","#ef8a62")
  col=rep(col, length.out=nf)
  
  ynum <- noage*nf+(nf-1)*2
  ylab=seq((noage/2)-1, ynum, (noage+yspace))
  
  if (missing(fleet_names)) fleet_names <- attr(fit$data, "fleetNames")
  for (k in 1:nf){ # resize fleet name length if too long (>19 characters)
    fleet_names[k] <- abbreviate(fleet_names[k],minlength = 19, dot=TRUE, use.classes=FALSE)
  }
  
  if (missing(fleet_type)){
    fleet_type <- fit$data$fleetTypes
    for (i in 1:nf){
      if(fleet_type[i]==0) fleet_type[i]<-"Catch at age"
      if(fleet_type[i]==1) fleet_type[i]<-"Catch at age with effort"
      if(fleet_type[i]==2 || fleet_type[i]==6) fleet_type[i]<-"Index at age"
      if(fleet_type[i]==3) fleet_type[i]<-"Biomass or catch index"
      if(fleet_type[i]==5) fleet_type[i]<-"Tagging data"
      if(fleet_type[i]==7) fleet_type[i]<-"Sum of fleets"
    }
  }

  layout(matrix(c(rep(1,3),2), nrow=1))
  par(oma=c(2,10,2,0),mar=c(0,0,0,0), xpd=NA)
  plot(x=rep(years,(ynum+1) ), y= rep(0:ynum,length(years)), type="n", xlab="", ylab="", yaxt="n")
  axis(2, labels = fleet_names, at =ylab ,las=1, cex=0.8)
  x=0-min(fit$data$minAgePerFleet, na.rm = TRUE)
  for (i in 1:nf){
    if (!is.na(fit$data$minAgePerFleet[i])){
      for (a in min(fit$data$minAgePerFleet, na.rm = TRUE):max(fit$data$maxAgePerFleet, na.rm = TRUE)){
        lines(x=years,y=rep(x+a,length(years)), col="grey87")
      }  
      for (a in fit$data$minAgePerFleet[i]:fit$data$maxAgePerFleet[i]){
        yval=dat[which(dat[,2]==i & dat[,3]==a),4]
        for (k in 1:length(yval)) if(!is.na(yval[k]) & yval[k]!=0) yval[k]=x+a else  yval[k]=NA
        lines(x=dat[which(dat[,2]==i & dat[,3]==a),1], y=yval, lwd=1, col=col[i])
        points(x=dat[which(dat[,2]==i & dat[,3]==a),1], y=yval, lwd=1, col=col[i], pch=16)
        text(x=years[length(years)]+0.5, y=x+a,labels=a, cex=0.7)
      }
    } else {
      a=-1
      yval=dat[which(dat[,2]==i & dat[,3]==a),4]
      for (k in 1:length(yval)) if(!is.na(yval[k]) & yval[k]!=0) yval[k]=x+noage/2 else  yval[k]=NA
      lines(x=dat[which(dat[,2]==i & dat[,3]==a),1], y=yval, lwd=1, col=col[i])
      points(x=dat[which(dat[,2]==i & dat[,3]==a),1], y=yval, lwd=1, col=col[i], pch=16)
    }
    x=x+noage+yspace
  }
  plot(x=rep(1,(ynum+1) ), y= rep(0:ynum), type="n", xlab="", ylab="", yaxt="n", bty="n", xaxt="n")
  for (i in 1:nf){
    text(x=1, y=ylab[i], labels=fleet_type[i])
  }
  mtext(text = "Available data", side=3, line=0.5, at=3/8, outer = TRUE)
  mtext(text = "Data type", side=3, line=0.5, at=7/8, outer = TRUE)
}



##' Plots the sd of the log observations as estimated in SAM in increasing order
##' @param fit the object returned from sam.fit
##' @param barcol color for each fleet and age
##' @param marg margin for plot (mar in par())
##' @param ylim bounds for y-axis
##' @param ... extra arguments to plot
##' @importFrom graphics barplot
##' @export
sdplot<-function(fit, barcol=NULL, marg=NULL, ylim=NULL, ...){
  UseMethod("sdplot")
}
##' @rdname sdplot
##' @method sdplot sam
##' @export
sdplot.sam <- function(fit, barcol=NULL, marg=NULL, ylim=NULL, ...){
  cf <- fit$conf$keyVarObs
  fn <- attr(fit$data, "fleetNames")
  ages <- fit$conf$minAge:fit$conf$maxAge
  pt <- partable(fit)
  sd <- unname(exp(pt[grep("logSdLogObs",rownames(pt)),1]))
  v<-cf
  v[] <- c(NA,sd)[cf+2]
  res<-data.frame(fleet=fn[as.vector(row(v))],name=paste0(fn[as.vector(row(v))]," age ",ages[as.vector(col(v))]), sd=as.vector(v))
  res<-res[complete.cases(res),]
  o<-order(res$sd)
  res<-res[o,]
  if (missing(barcol)) barcol <- colors()[as.integer(as.factor(res$fleet))*10]
  if (missing(marg)) marg <- c(max(nchar(fn))*0.8,4,2,1)
  par(mar=marg)
  barplot(res$sd, names.arg=res$name,las=2, col=barcol, ylab="SD", ylim=ylim); box()
}
