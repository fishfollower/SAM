##' Plot helper
##' @param fit the fitted object from sam.fit
##' @param what quoted name of object to extract
##' @param x x-alues
##' @param ymab label on y-axis
##' @param xlab label on x-axis
##' @param ex extra y's to make room for
##' @param trans function to transform values by
##' @param Add if plotting is to be added on existing plot
##' @param ci if confidence intervals should be plottet
##' @param cicol color to plot the confidence polygon
##' @param drop number of years to be left unplotted at the end. 
##' @param ... extra arguments transferred to plot
##' @importFrom graphics plot polygon grid lines
##' @importFrom grDevices gray
##' @details ...
.plotit <-function (fit, what, x=fit$data$years, ylab=what, xlab="Years", ex=numeric(0), trans=function(x)x, add=FALSE, ci=TRUE, cicol=gray(.5,alpha=.5), drop=0,...){
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
    col10=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
    if(!is.null(attr(fit,"fit"))){
      fitlocal <- attr(fit,"fit")
    }else{
      fitlocal <- fit[[1]]
    }
    if(is.null(x))x=fitlocal$data$years
    .plotit(fitlocal, what=what, x=x, ylab=ylab, xlab=xlab, ex=ex, trans=trans, add=add, ci=ci, cicol=cicol, drop=drop,...)
    d<-lapply(1:length(fit), function(i).plotit(fit[[i]], what=what, trans=trans, add=TRUE, ci=FALSE, col=col10[(i-1)%%10+1], drop=drop, ...))
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
##' @details ...
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
##' @param ... extra arguments transferred to plot
##' @importFrom graphics matplot
##' @details ...
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
##' @param ... extra arguments transferred to plot
##' @details ...
##' @export
ssbplot<-function(fit, ...){
  .plotit(fit, "logssb", ylab="SSB", trans=exp,...)
}

##' SAM TSB plot 
##' @param fit the object returned from sam.fit
##' @param ... extra arguments transferred to plot
##' @details ...
##' @export
tsbplot<-function(fit, ...){
  .plotit(fit, "logtsb", ylab="TSB", trans=exp,...)
}

##' SAM Recruits plot 
##' @param fit the object returned from sam.fit
##' @param ... extra arguments transferred to plot
##' @details ...
##' @export
recplot<-function(fit,...){
   lab<-paste("Recruits (age ", fit$conf$minAge, ")", sep="")
  .plotit(fit, "logR", ylab=lab, trans=exp,...)
}

##' SAM catch plot 
##' @param fit the object returned from sam.fit
##' @param obs.show if observations are to be shown also
##' @param ... extra arguments transferred to plot
##' @param drop number of years to be left unplotted at the end. Default (NULL) is to not show years at the end with no catch information 
##' @details ...
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
##' @details ...
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
##' @author Casper
obscov<-function(fit, corr=FALSE){
    stopifnot(class(fit)=="sam")
    res<-fit$rep$obsCov
    for(i in 1:length(res)) rownames(res[[i]])<-fit$data$minAgePerFleet[i]:fit$data$maxAgePerFleet[i]
    if(corr) for(i in 1:length(res)) res[[i]]<-cov2cor(res[[i]])
    res
}

##' A generic function returning an ellipse or other outline of a
##' confidence region for two parameters.
##' 
##' Taken from the 'ellipse' package.
##' @param x An object. In the default method the parameter ‘x’ should be
##' a correlation between -1 and 1 or a square positive definite
##' matrix at least 2x2 in size. It will be treated as the
##' correlation or covariance of a multivariate normal
##' distribution.
##' @param scale If ‘x’ is a correlation matrix, then the standard deviations
##' of each parameter can be given in the scale parameter.  This
##' defaults to ‘c(1, 1)’, so no rescaling will be done.
##' @param centre The centre of the ellipse will be at this position.
##' @param level The confidence level of a pairwise confidence region.  The
##' default is 0.95, for a 95% region.  This is used to control
##' the size of the ellipse being plotted.  A vector of levels
##' may be used.
##' @param t The size of the ellipse may also be controlled by specifying
##' the value of a t-statistic on its boundary.  This defaults to
##' the appropriate value for the confidence region.
##' @param which This parameter selects which pair of variables from the
##' matrix will be plotted.  The default is the first 2.
##' @param npoints The number of points used in the ellipse.  Default is 100.
##' @param ... Descendant methods may require additional parameters.
##' @return An ‘npoints’ x ‘2’ matrix is returned with columns named according
##' to the row names of the matrix ‘x’ (default ‘'x'’ and ‘'y'’),
##' suitable for plotting.
ellipse<-function (x, scale = c(1, 1), centre = c(0, 0), level = 0.95, 
                   t = sqrt(qchisq(level, 2)), which = c(1, 2), npoints = 100, 
                   ...){
    names <- c("x", "y")
    if (is.matrix(x)) {
        xind <- which[1]
        yind <- which[2]
        r <- x[xind, yind]
        if (missing(scale)) {
            scale <- sqrt(c(x[xind, xind], x[yind, yind]))
            if (scale[1] > 0) 
                r <- r/scale[1]
            if (scale[2] > 0) 
                r <- r/scale[2]
        }
        if (!is.null(dimnames(x)[[1]])) 
            names <- dimnames(x)[[1]][c(xind, yind)]
    }
    else r <- x
    r <- min(max(r, -1), 1)
    d <- acos(r)
    a <- seq(0, 2 * pi, len = npoints)
    matrix(c(t * scale[1] * cos(a + d/2) + centre[1], t * scale[2] * 
        cos(a - d/2) + centre[2]), npoints, 2, dimnames = list(NULL, 
        names))
}

##' This function plots a correlation matrix using ellipse-shaped
##' glyphs for each entry.  The ellipse represents a level curve of
##' the density of a bivariate normal with the matching correlation.
##'
##' Taken from the 'ellipse' package.
##' @title 
##' @param corr A matrix containing entries between ‘-1’ and ‘1’ to be plotted as correlations.
##' @param outline Whether the ellipses should be outlined in the default colour.
##' @param col Which colour(s) to use to fill the ellipses.
##' @param numbers Whether to plot numerical correlations in place of ellipses.
##'           If numbers is ‘TRUE’, then the correlations will be rounded
##'           to a single decimal place and placed on the plot.
##' @param type Character. Plot ‘"full"’ matrix or just ‘"upper"’ or
##' ‘"lower"’ triangular part of it.
##' @param diag Logical. Plot diagonal elements or not.
##' @param bty for ‘plot’
##' @param axes for ‘plot’
##' @param xlab for ‘plot’
##' @param ylab for ‘plot’
##' @param asp for ‘plot’
##' @param cex.lab for ‘plot’
##' @param cex for ‘plot’
##' @param mar for ‘plot’
##' @param ... for ‘plot’
plotcorr<-function(corr, outline = TRUE, col = "grey", numbers = FALSE, 
    type = c("full", "lower", "upper"), diag = (type == "full"), 
    bty = "n", axes = FALSE, xlab = "", ylab = "", asp = 1, cex.lab = par("cex.lab"), 
    cex = 0.75 * par("cex"), mar = 0.1 + c(2, 2, 4, 2), ...) 
{
    savepar <- par(pty = "s", mar = mar)
    on.exit(par(savepar))
    if (is.null(corr)) 
        return(invisible())
    if ((!is.matrix(corr)) || (round(min(corr, na.rm = TRUE), 
        6) < -1) || (round(max(corr, na.rm = TRUE), 6) > 1)) 
        stop("Need a correlation matrix")
    plot.new()
    par(new = TRUE)
    rowdim <- dim(corr)[1]
    coldim <- dim(corr)[2]
    rowlabs <- dimnames(corr)[[1]]
    collabs <- dimnames(corr)[[2]]
    if (is.null(rowlabs)) 
        rowlabs <- 1:rowdim
    if (is.null(collabs)) 
        collabs <- 1:coldim
    rowlabs <- as.character(rowlabs)
    collabs <- as.character(collabs)
    col <- rep(col, length = length(corr))
    dim(col) <- dim(corr)
    type <- match.arg(type)
    cols <- 1:coldim
    rows <- 1:rowdim
    xshift <- 0
    yshift <- 0
    if (!diag) {
        if (type == "upper") {
            cols <- 2:coldim
            rows <- 1:(rowdim - 1)
            xshift <- 1
        }
        else if (type == "lower") {
            cols <- 1:(coldim - 1)
            rows <- 2:rowdim
            yshift <- -1
        }
    }
    maxdim <- max(length(rows), length(cols))
    plt <- par("plt")
    xlabwidth <- max(strwidth(rowlabs[rows], units = "figure", 
        cex = cex.lab))/(plt[2] - plt[1])
    xlabwidth <- xlabwidth * maxdim/(1 - xlabwidth)
    ylabwidth <- max(strwidth(collabs[cols], units = "figure", 
        cex = cex.lab))/(plt[4] - plt[3])
    ylabwidth <- ylabwidth * maxdim/(1 - ylabwidth)
    plot(c(-xlabwidth - 0.5, maxdim + 0.5), c(0.5, maxdim + 1 + 
        ylabwidth), type = "n", bty = bty, axes = axes, xlab = "", 
        ylab = "", asp = asp, cex.lab = cex.lab, ...)
    text(rep(0, length(rows)), length(rows):1, labels = rowlabs[rows], 
        adj = 1, cex = cex.lab)
    text(cols - xshift, rep(length(rows) + 1, length(cols)), 
        labels = collabs[cols], srt = 90, adj = 0, cex = cex.lab)
    mtext(xlab, 1, 0)
    mtext(ylab, 2, 0)
    mat <- diag(c(1, 1))
    plotcorrInternal <- function() {
        if (i == j && !diag) 
            return()
        if (!numbers) {
            mat[1, 2] <- corr[i, j]
            mat[2, 1] <- mat[1, 2]
            ell <- ellipse(mat, t = 0.43)
            ell[, 1] <- ell[, 1] + j - xshift
            ell[, 2] <- ell[, 2] + length(rows) + 1 - i - yshift
            polygon(ell, col = col[i, j])
            if (outline) 
                lines(ell)
        }
        else {
            text(j + 0.3 - xshift, length(rows) + 1 - i - yshift, 
                round(10 * corr[i, j], 0), adj = 1, cex = cex)
        }
    }
    for (i in 1:dim(corr)[1]) {
        for (j in 1:dim(corr)[2]) {
            if (type == "full") {
                plotcorrInternal()
            }
            else if (type == "lower" && (i >= j)) {
                plotcorrInternal()
            }
            else if (type == "upper" && (i <= j)) {
                plotcorrInternal()
            }
        }
    }
    invisible()
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param fit 
##' @param ... 
##' @return 
##' @author Casper
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
