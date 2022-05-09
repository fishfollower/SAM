##' Plot sam object 
##' @method plot sam
##' @param  x fitted object as returned from the \code{\link{sam.fit}} function.
##' @param  ... extra arguments (not possible to use add=TRUE --- instead collect to a list of fits using e.g the c(...), and then plot that collected object).
##' @importFrom graphics par
##' @details gives a 3 plot overview plot og ssb, fbar, and recruits. These plots are available individually via the functions \code{\link{ssbplot}}, \code{fbarplot}, and \code{\link{recplot}}.  
##' @export
plot.sam<-function(x, ...){
  dots <- list(...) 
  if("add" %in% names(dots) && dots$add==TRUE){stop("Not possible to use add=TRUE here --- please collect to a list of fits using c(...)")}
  op<-par(mfrow=c(3,1))
  ssbplot(x,...)
  fbarplot(x,...)
  recplot(x,...)
  par(op)
}

##' Plot samforecast object 
##' @method plot samforecast
##' @param  x fitted object as returned from the \code{\link{sam.fit}} function
##' @param  ... extra arguments
##' @importFrom graphics par
##' @details ...
##' @export
plot.samforecast<-function(x, ...){
  op<-par(mfrow=c(3,1))
  ssbplot(x,...)
  fbarplot(x, drop=0,...)
  recplot(x,...)
  par(op)
}

##' Plot hcr object 
##' @method plot hcr
##' @param  x output from the \code{\link{hcr}} function
##' @param  ... extra arguments
##' @importFrom graphics par
##' @details ...
##' @export
plot.hcr<-function(x, ...){
  op<-par(mfrow=c(3,1))
  ssbplot(x,...)
  fbarplot(x, drop=0,...)
  recplot(x,...)
  par(op)
}

##' Collect sam objects 
##' @method c sam
##' @param  ... one or more sam fits (as returned from the \code{\link{sam.fit}} function) to be combined 
##' @importFrom graphics par
##' @details ...
##' @export
c.sam<-function(...){
  ret<-list(...)
  class(ret)<-"samset"
  ret
}

##' Plot sam object 
##' @method plot samset
##' @param  x fitted object as returned from the \code{\link{sam.fit}} function.
##' @param  ... extra arguments
##' @importFrom graphics par
##' @details ...
##' @export
plot.samset<-function(x, ...){
  op<-par(mfrow=c(3,1))
  ssbplot(x,...)
  fbarplot(x,...)
  recplot(x,...)
  par(op)
}

##' Compute process residuals (single joint sample) 
##' @param fit the fitted object as returned from the \code{\link{sam.fit}} function
##' @param map map from original fit 
##' @param ... extra arguments (not currently used)
##' @return an object of class \code{samres}
##' @details Single joint sample residuals of log(F) and log(N)
##' @importFrom TMB sdreport
##' @export
procres <- function(fit, map = fit$obj$env$map, ...){
  pp<-fit$pl
  attr(pp,"what") <- NULL
  pp$missing <- NULL
  fakefile <- file()
  sink(fakefile)
  saveConf(fit$conf, file="")
  sink()
  conf <- loadConf(fit$data, fakefile, patch=TRUE)
  close(fakefile)
  fit.co<-sam.fit(fit$data, conf, pp, run=FALSE, map = map)
  fit.co$obj$env$data$resFlag<-1
  fit.co$obj$retape()
  sdrep <- sdreport(fit.co$obj,fit$opt$par)  
  ages <- as.integer(colnames(fit.co$data$natMor))

  if (exists(".Random.seed")){
    oldseed <- get(".Random.seed", .GlobalEnv)
    oldRNGkind <- RNGkind()
  }
  set.seed(123456)
  idx <- which(names(sdrep$value)=="resN")
  resN <- rmvnorm(1,mu=sdrep$value[idx], Sigma=sdrep$cov[idx,idx])
  resN <- matrix(resN, nrow=nrow(fit.co$pl$logN))
  resN <- data.frame(year=fit.co$data$years[as.vector(col(resN))],
                     fleet=1,
                     age=ages[as.vector(row(resN))],
                     residual=as.vector(resN))
  idx <- which(names(sdrep$value)=="resF")
  resF <- rmvnorm(1,mu=sdrep$value[idx], Sigma=sdrep$cov[idx,idx])
  resF <- matrix(resF, nrow=nrow(fit.co$pl$logF))

  fleets <- which(fit.co$data$fleetTypes==0)
  nfleets <- length(fleets)
  f<-function(i){
      if(i==0){
          return(resN)
      }else{
        iF<-fit.co$conf$keyLogFsta[fleets[i],]
        sub<-resF[iF[iF>=0]+1,]  
        iF.std <- iF
        iF.std[iF.std>=0]<-iF.std[iF.std>=0]-min(iF.std[iF.std>=0])+1
        return(data.frame(year=fit.co$data$years[as.vector(col(sub))],
               fleet=1+i,
               age=ages[iF.std[iF.std>=0]][as.vector(row(sub))],
               residual=as.vector(sub)))
      }
  }
  resNF<-lapply(0:nfleets, f)
  ret <- do.call(rbind, resNF)
  attr(ret, "fleetNames") <- c("Joint sample residuals log(N)", attr(fit.co$data, "fleetNames")[fleets])
  class(ret) <- "samres"
  if (exists("oldseed")){
    do.call("RNGkind",as.list(oldRNGkind))
    assign(".Random.seed", oldseed, .GlobalEnv)
  }
  return(ret)
}

##' Plot sam residuals 
##' @method plot samres 
##' @param x an object of type 'samres' as returned from functions \code{\link{residuals.sam}} or \code{\link{procres}}.
##' @param type either "bubble" (default) or "summary"
##' @param ... extra arguments
##' @details In the "bubble" type red indicate negative residuals and blue positive. The area of the circles scales with the absolute size of the residuals.  
##' @importFrom graphics abline
##' @importFrom stats acf na.pass qqnorm
##' @export
##' @examples
##' \dontrun{
##' data(nscodData)
##' data(nscodConf)
##' data(nscodParameters)
##' fit <- sam.fit(nscodData, nscodConf, nscodParameters)
##' par(ask=FALSE)
##' plot(residuals(fit))
##' }
plot.samres<-function(x, type="bubble",...){
  if(type=="bubble"){
    add_legend <- function(x, ...) {
      opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
                  mar=c(0, 0, 0, 0), new=TRUE)
      on.exit(par(opar))
      plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
      zscale <- pretty(x$residual,min.n=4)
      uu<-par("usr")
      yy<-rep(uu[3]+.03*(uu[4]-uu[3]), length(zscale))
      xx<-seq(uu[1]+.10*(uu[2]-uu[1]),uu[1]+.4*(uu[2]-uu[1]), length=length(zscale))
      text(xx,yy,labels=zscale)
      colb <- ifelse(zscale<0, rgb(1, 0, 0, alpha=.5), rgb(0, 0, 1, alpha=.5))
      bs<-1
      if("bubblescale"%in%names(list(...))) bs <- list(...)$bubblescale
      points(xx,yy,cex=sqrt(abs(zscale))/max(sqrt(abs(zscale)), na.rm=TRUE)*5*bs, pch=19, col=colb)
    }
    neg.age <- (x$age < -1.0e-6)
    x$age[neg.age] <- mean(x$age[!neg.age],na.rm=TRUE)
    plotby(x$year, x$age, x$residual, by=attr(x,"fleetNames")[x$fleet], xlab="Year", ylab="Age", ...)
    add_legend(x, ...)
  } else if(type=="summary"){
      tmp <- xtabs( x$residual ~ x$year + x$age + x$fleet)
      nfleets <- dim(tmp)[3]
      op<-par(mfrow=rev(n2mfrow(3*nfleets)),mar=c(4,5,1,1))
      on.exit(par(op))
      fnames <- attr(x,"fleetNames")
      for(f in 1:nfleets){
          sel<-which(x$fleet==f)
          plot(factor(x$age[sel]),x$residual[sel])
          mtext(fnames[f],side=2,line=3) 
          tmp2 <- t(tmp[,,f])
          tmp2 <- rbind(tmp2, matrix(NA,nrow=nrow(tmp2),ncol=ncol(tmp2)))
          acf(as.vector(tmp2),na.action=na.pass,lag.max=dim(tmp)[2])
          qqnorm(x$residual[sel])
          abline(0,1)
      }
  } else stop(paste0("unknown type: ",type))
}

##' Print sam object 
##' @method print sam 
##' @param  x the fitted object as returned from the \code{\link{sam.fit}} function
##' @param  ... extra arguments
##' @details prints the log-likelihood and the main convergence criteria
##' @export
print.sam<-function(x, ...){
  cat("SAM model: log likelihood is", logLik.sam(x,...),"Convergence", ifelse(0==x$opt$convergence, "OK\n", "failed\n"))
}

##' Print samres object 
##' @method print samres 
##' @param  x a sam residual object as returned from either \code{\link{residuals.sam}} or \code{\link{procres}}
##' @param  ... extra arguments
##' @details prints the residuals as a data.frame 
##' @export
print.samres<-function(x, ...){
  class(x)<-NULL
  print(as.data.frame(x))
}

##' Print hcr object 
##' @method print hcr 
##' @param  x a sam hcr object as returned by \code{\link{hcr}}
##' @param  ... extra arguments
##' @details prints the HCR forecast
##' @export
print.hcr<-function(x, ...){
  print(x$forecast)
}


##' Print referencepoint object 
##' @method print sam_referencepoints 
##' @param  x a sam referencepoint object as returned by \code{\link{referencepoints}}
##' @param  ... extra arguments
##' @details prints the F reference point table
##' @export
print.sam_referencepoints<-function(x, tables = c("F","Biomass","Yield"),digits = 4, format = "f"){
    if(is.na(any(match(tables, names(x$tables)))))
        stop("Table does not exist")
    tbntxt <- ifelse(tables == "F", "Fishing mortality rate", paste("Equilibrium",tables))
    for(tt in seq_along(tables)){
        t0 <- x$tables[[tables[tt]]]
        t1 <- matrix(formatC(t0,digits = digits, format = format),nrow(t0), ncol(t0))
        rownames(t1) <- rownames(t0)
        colnames(t1) <- colnames(t0)   
        tab <- as.table(t1)
        dn <- attr(tab,"dimnames")
        names(dn) <- c(tbntxt[tt],"")
        attr(tab,"dimnames") <- dn
        print(tab, justify = "right", right = TRUE)
    }
}


##' Log likelihood of sam object 
##' @method logLik sam 
##' @param  object sam fitted object as returned from the \code{\link{sam.fit}} function
##' @param  ... extra arguments
##' @details log likelihood of a sam model run 
##' @export
logLik.sam<-function(object, ...){
  ret<- -object$opt$objective
  attr(ret,"df")<-length(object$opt$par)
  class(ret)<-"logLik"
  ret
}

##' Extract fixed coefficients of sam object 
##' @method coef sam 
##' @param  object sam fitted object as returned from the \code{\link{sam.fit}} function
##' @param  ... extra arguments
##' @details fixed coefficients of sam object
##' @importFrom stats coef
##' @export
coef.sam <- function(object, ...){
  ret <- object$sdrep$par.fixed
  attr(ret,"cov") <- object$sdrep$cov.fixed
  attr(ret,"sd") <- sqrt(diag(object$sdrep$cov.fixed))
  class(ret)<-"samcoef"
  ret
}

##' Print samcoef object 
##' @method print samcoef 
##' @param  x samcoef object as returned from the \code{\link{coef.sam}} function
##' @param  ... extra arguments
##' @export
print.samcoef<-function(x, ...){
  y<-as.vector(x)
  names(y)<-names(x)
  print(y)
}


##' Extract number of observations from sam object 
##' @method nobs sam 
##' @param object sam fitted object as returned from the \code{\link{sam.fit}} function
##' @param ... extra arguments
##' @importFrom stats nobs
##' @details ...
##' @export
nobs.sam<-function(object, ...){
  as.integer(object$data$nobs)
}

##' Extract residuals from sam object 
##' @method residuals sam 
##' @param object sam fitted object as returned from the \code{\link{sam.fit}} function
##' @param discrete logical if model contain discrete observations  
##' @param ... extra arguments for TMB's oneStepPredict
##' @details one-observation-ahead quantile residuals are calculated
##' @importFrom stats residuals
##' @importFrom TMB oneStepPredict
##' @details ...
##' @export
residuals.sam<-function(object, discrete=FALSE, ...){
  pp<-object$pl
  attr(pp,"what") <- NULL
  pp$missing <- NULL
  fakefile <- file()
  sink(fakefile)
  saveConf(object$conf, file="")
  sink()
  conf <- loadConf(object$data, fakefile, patch=TRUE)
  close(fakefile)
  object.co<-sam.fit(object$data, conf, pp, run=FALSE, map = object$obj$env$map)
  cat("One-observation-ahead residuals. Total number of observations: ", nobs(object), "\n")  
  res <- oneStepPredict(object.co$obj, observation.name="logobs", data.term.indicator="keep", discrete=discrete,...)
  cat("One-observation-ahead residuals. Done\n")  
  ret <- cbind(object.co$data$aux, res)
  attr(ret,"fleetNames") <- attr(object.co$data, "fleetNames")
  class(ret)<-"samres"
  ret
}

##' Summary of sam object 
##' @method summary sam 
##' @param object sam fitted object as returned from the \code{\link{sam.fit}} function
##' @param ... extra arguments 
##' @details summary table containing recruits, SSB, and Fbar
##' @export
summary.sam<-function(object, ...){
  dots<-list(...)  
  ret <- cbind(round(rectable(object,...)), round(ssbtable(object)), round(fbartable(object),3))
  add <- 0
  if(!is.null(dots$lagR)){
    if(dots$lagR==TRUE){
      add <- 1
    }
  }
  colnames(ret)[1] <- paste("R(age ", object$conf$minAge+add, ")", sep="") 
  colnames(ret)[4] <- "SSB"
  colnames(ret)[7] <- paste("Fbar(",object$conf$fbarRange[1], "-", object$conf$fbarRange[2], ")", sep="")
  ret
}

##' Simulate from a sam object 
##' @method simulate sam 
##' @param object sam fitted object as returned from the \code{\link{sam.fit}} function
##' @param nsim number of response lists to simulate. Defaults to 1.
##' @param seed random number seed
##' @param full.data logical, should each inner list contain a full list of data. Defaults to TRUE  
##' @param keep.process Keep logN and logF processes when full.data = TRUE?
##' @param ... extra arguments
##' @importFrom stats simulate
##' @details simulates data sets from the model fitted and conditioned on the random effects estimated 
##' @return returns a list of lists. The outer list has length \code{nsim}. Each inner list contains simulated values of \code{logF}, \code{logN}, and \code{obs} with dimensions equal to those parameters.
##' @export
simulate.sam<-function(object, nsim=1, seed=NULL, full.data=TRUE, keep.process = FALSE, ...){
    if(!is.null(seed)) set.seed(seed)
    pl <- object$pl
    map <- object$obj$env$map
    with.map <- intersect(names(pl),names(map))
    applyMap <- function(par.name){
        tapply(pl[[par.name]],map[[par.name]],mean)
    }
    pl[with.map] <- sapply(with.map, applyMap, simplify = FALSE)
    est <- unlist(pl)
    if(full.data){
        ret <- replicate(nsim, {
            sval <- object$obj$simulate(est)
            ret <- c(object$data[names(object$data)!="logobs"],#all the old data
                     sval["logobs"])#simulated observations
            if(keep.process){
                ret$logN <- sval$logN
                ret$logF <- sval$logF
            }
            ret
        },
        simplify=FALSE)
        ret<-lapply(ret, function(x){attr(x,"fleetNames") <- attr(object$data,"fleetNames");x})
    }else{
  	ret <- replicate(nsim, object$obj$simulate(est), simplify=FALSE)
    }
    ret
}


##' Plot sam object 
##' @method plot samypr
##' @param  x ...
##' @param  ... extra arguments 
##' @importFrom graphics par title
##' @details ...
##' @export
plot.samypr<-function(x, ...){
  par(mar=c(5.1,4.1,4.1,5.1))
  plot(x$fbar, x$yield, type='l', xlab=x$fbarlab, ylab='Yield per recruit', ...)
  lines(c(x$fmax,x$fmax), c(par('usr')[1],x$yield[x$fmaxIdx]), lwd=3, col='red')
  lines(c(x$f01,x$f01), c(par('usr')[1],x$yield[x$f01Idx]), lwd=3, col='blue')  
  ssbscale <- max(x$yield)/max(x$ssb)

  lines(x$fbar, ssbscale*x$ssb, lty='dotted')
  ssbtick <- pretty(x$ssb)
  ssbat <- ssbtick*ssbscale
  axis(4,at=ssbat, labels=ssbtick)
  mtext('SSB per recruit', side=4, line=2)

  lines(c(x$f35,x$f35), c(par('usr')[1],x$ssb[x$f35Idx]*ssbscale), lwd=3, col='green')

  title(eval(substitute(expression(F[max]==fmax~ ~ ~ ~ ~F[0.10]==f01~ ~ ~ ~ ~F[0.35*SPR]==f35), 
                        list(fmax=round(x$fmax,2), f01=round(x$f01,2), f35=round(x$f35,2)))))
}


##' Print samypr object 
##' @method print samypr 
##' @param  x an object as returned from the ypr function
##' @param  ... extra arguments
##' @details ...
##' @export
print.samypr <- function(x, ...){
  idx <- c(x$fmaxIdx, x$f01Idx, x$f35Idx)
  ret <- cbind(x$fbar[idx],x$ssb[idx],x$yield[idx])
  rownames(ret) <- c("Fmax", "F01", "F35")
  colnames(ret) <- c("Fbar", "SSB", "Yield")
  print(ret)
}

##' Print samforecast object 
##' @method print samforecast 
##' @param  x an object as returned from the forecast function
##' @param  ... extra arguments
##' @details ...
##' @export
print.samforecast<-function(x, ...){
  print(attr(x,"tab"))
}

##' Print samset object 
##' @method print samset 
##' @param  x a list of sam models
##' @param  ... extra arguments
##' @details ...
##' @importFrom stats logLik
##' @export 
print.samset<-function(x,...){
  if("jitflag"%in%names(attributes(x))){  
    fit<-attr(x,"fit")
    maxabsdiff <- apply(abs(do.call(cbind, lapply(x, function(f)unlist(f$pl)-unlist(fit$pl)))),1,max)
    maxlist <- relist(maxabsdiff, fit$pl)
    ret <- as.data.frame(unlist(lapply(maxlist,function(x)if(length(x)>0)max(x) else NULL)))
    fbar <- max(unlist(lapply(x, function(f)abs(fbartable(f)[,1]-fbartable(fit)[,1]))))
    ssb <- max(unlist(lapply(x, function(f)abs(ssbtable(f)[,1]-ssbtable(fit)[,1]))))
    rec <- max(unlist(lapply(x, function(f)abs(rectable(f)[,1]-rectable(fit)[,1]))))
    catch <- max(unlist(lapply(x, function(f)abs(catchtable(f)[,1]-catchtable(fit)[,1]))))
    logLik <- max(abs(unlist(lapply(x, logLik))-logLik(fit)))
    ret <- rbind(ret, ssb=ssb,  fbar=fbar, rec=rec, catch=catch, logLik=logLik)
    names(ret) <- "max(|delta|)" 
    print(ret)
  }else{
    print.default(x,...)
  }
}
