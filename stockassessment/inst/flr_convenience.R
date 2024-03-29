as.matrix.FLQuant <- function(v){
    require(FLCore)
    d <- dim(v@.Data)
    dn <- dimnames(v@.Data)
    xx <- t(matrix(v@.Data, d[1],d[2]))
    colnames(xx) <- dn[[1]]
    rownames(xx) <- dn[[2]]
    xx
}

as.FLStock.sam <- function(fit, unit.w = "kg", name = "", desc = "", predicted = FALSE){
    require(FLCore)
    toFLQ <- function(x, unit = "NA"){
        FLQuant(t(x),dimnames=list(age=as.numeric(colnames(x)),
                                   year=as.numeric(rownames(x)),
                                   unit="unique",
                                   season="all",
                                   area="unique",
                                   iter=1),
                units = unit)
    }
    na2zero <- function(x){
        x[is.na(x)] <- 0
        x
    }
    resize <- function(x, ages, years, replicate = FALSE){
        x1 <- matrix(NA,length(years),length(ages))
        colnames(x1) <- ages
        rownames(x1) <- years
        x1[match(rownames(x),rownames(x1)),match(colnames(x),colnames(x1))] <- x
        if(replicate){
            ## Years
            ii <- which(!match(rownames(x1),rownames(x), FALSE))
            if(length(ii) > 0){
                isAfter <- as.numeric(rownames(x1)[ii]) > max(as.numeric(rownames(x)))
                x1[ii,] <- x[c(1,nrow(x))[as.numeric(isAfter)+1],]
            }
            ## Ages
            ii <- which(!match(colnames(x1),colnames(x), FALSE))
            if(length(ii) > 0){
                isAfter <- as.numeric(colnames(x1)[ii]) > max(as.numeric(colnames(x)))
                x1[,ii] <- x[,c(1,ncol(x))[as.numeric(isAfter)+1]]
            }
        }
        x1
    }
    ## catch.n
    ages <- fit$conf$minAge:fit$conf$maxAge
    years <- fit$data$years
    CN <- toFLQ(resize(na2zero(getFleet(fit,1, predicted)),ages,years))
    LF <- toFLQ(resize(fit$data$landFrac,ages,years,TRUE))
    ## mat
    MA <- toFLQ(resize(fit$data$propMat,ages,years,TRUE))
    ## stock.wt
    SW <- toFLQ(resize(fit$data$stockMeanWeight,ages,years,TRUE), unit.w)
    ## stock.n
    SN <- toFLQ(resize(ntable(fit),ages,years))
    ## catch.wt
    CW <- toFLQ(resize(fit$data$catchMeanWeight,ages,years,TRUE), unit.w)
    ## discards.wt
    DW <- toFLQ(resize(fit$data$disMeanWeight,ages,years,TRUE), unit.w)
    ## discards.n
    DN <- CN * (1 - LF)
    ## landings.wt
    LW <- toFLQ(resize(fit$data$landMeanWeight,ages,years,TRUE), unit.w)
    ## landings.n
    LN <- CN * LF
    ## m
    M <- toFLQ(resize(fit$data$natMor,ages,years,TRUE))
    ## harvest.spwn
    PF <- toFLQ(resize(fit$data$propF,ages,years,TRUE))
    ## m.spwn
    PM <- toFLQ(resize(fit$data$propM,ages,years,TRUE))
    ## harvest
    F <- toFLQ(resize(faytable(fit),ages,years,FALSE), "f")
    FLStock(catch = FLCore::apply(CN * CW,2,sum), # Total catch weight (‘FLQuant’)
            catch.n = CN, # Catch numbers (‘FLQuant’)
            catch.wt = CW, # Mean catch weights (‘FLQuant’)
            discards = FLCore::apply(DW * DN,2,sum), # Total discards weight (‘FLQuant’)
            discards.n = DN, # Discard numbers (‘FLQuant’)
            discards.wt = DW, # Mean discard weights (‘FLQuant’)
            landings = FLCore::apply(LW * LN,2,sum), # Total landings weight (‘FLQuant’)
            landings.n = LN, # Landing numbers (‘FLQuant’)
            landings.wt = LW, # Landing weights (‘FLQuant’)
            stock = FLCore::apply(SW * SN,2,sum), # Total stock weight (‘FLQuant’)
            stock.n = SN, # Stock numbers (‘FLQuant’)
            stock.wt = SW, # Mean stock weights (‘FLQuant’)
            m = M, # Natural mortality (‘FLQuant’)
            mat = MA, # Proportion mature (‘FLQuant’)
            harvest = F, # Harvest rate or fishing mortality. The units of this slot should be set to 'hr' or 'f' accordingly (‘FLQuant’)
            harvest.spwn = PF, # Proportion of harvest/fishing mortality before spawning (‘FLQuant’)
            m.spwn = PM, # Proportion of natural mortality before spawning (‘FLQuant’)
            name = name, # Name of the stock (‘character’)
            desc = desc, # Description of the stock (‘character’)
            range = c(min = fit$conf$minAge,
                      max = fit$conf$maxAge,
                      plusgroup = ifelse(fit$conf$maxAgePlusGroup[1],fit$conf$maxAge,NA),
                      minyear = min(fit$data$years),
                      maxyear = max(fit$data$years),
                      minfbar = fit$conf$fbarRange[2],
                      maxfbar = fit$conf$fbarRange[2]
                      ) # Named numeric vector containing the quant and year ranges, the plusgroup and the quant range that the average fishing mortality should be calculated over (‘numeric’)
            )
}

as.sam.fit <- function(x,...){ UseMethod("as.sam.fit")    }

getSamFitData <- function(x, makePlusGroup = TRUE){
    CN <- as.matrix(catch.n(x))
    CN[CN==0] <- NA
    MA <- as.matrix(mat(x))
    SW <- as.matrix(stock.wt(x))
    CW <- as.matrix(catch.wt(x))
    DW <- as.matrix(discards.wt(x))
    LW <- as.matrix(landings.wt(x))
    M <- as.matrix(m(x))
    PF <- as.matrix(harvest.spwn(x))
    PM <- as.matrix(m.spwn(x))
    LF <- as.matrix(landings.n(x) / catch.n(x))
    LF[is.nan(LF)] <- 1
    ii <- which(apply(CN,2,function(x)mean(is.na(x))) > 0.10)          
    if(makePlusGroup && length(ii) > 0){
        if(is.numeric(makePlusGroup)){
            AgeIPlus <- which(colnames(CN == makePlusGroup))
        }else{
            AgeIPlus <- min(ii)
        }
        if(AgeIPlus > 1){
            reduceMat <- function(x,Ai,FUN,...){
                i1 <- seq(1,Ai-1,1)
                i2 <- seq(Ai,ncol(x),1)
                inam <- seq(1,Ai,1)
                x2 <- apply(x[,i2, drop = FALSE], 1, function(x)FUN(na.omit(x),...))
                y <- cbind(x[,i1,drop=FALSE], x2)
                colnames(y) <- colnames(x)[inam]
                y
            }
            CN <- reduceMat(CN,AgeIPlus,sum)
            CN[CN==0] <- NA
            MA <- reduceMat(MA,AgeIPlus,mean)
            SW <- reduceMat(SW,AgeIPlus,mean)
            CW <- reduceMat(CW,AgeIPlus,mean)
            DW <- reduceMat(DW,AgeIPlus,mean)
            LW <- reduceMat(LW,AgeIPlus,mean)
            M <- reduceMat(M,AgeIPlus,mean)
            PF <- reduceMat(PF,AgeIPlus,mean)
            PM <- reduceMat(PM,AgeIPlus,mean)
            LF <- reduceMat(LF,AgeIPlus,mean)
        }else{
            warning("Calculate plus group is the first age, skipping!")
        }
    }else if(makePlusGroup && length(ii) == 0){
        message("Constructed plus group not needed.")
    }
    dat <- setup.sam.data(fleets = NULL,
                          surveys = NULL,
                          residual.fleet = CN, 
                          prop.mature = MA,
                          stock.mean.weight = SW,
                          catch.mean.weight = CW, 
                          dis.mean.weight = DW, 
                          land.mean.weight = LW,
                          natural.mortality = M, 
                          prop.f = PF, 
                          prop.m = PM,
                          land.frac = LF,
                          recapture = NULL)
    dat
}
as.sam.fit.FLStock <- function(x, sr, rec_pars, makePlusGroup = TRUE){
    if(missing(sr))
        sr <- 0
    dat <- getSamFitData(x, makePlusGroup)
    cnf <- defcon(dat)
    cnf$keyLogFsta[1,] <- seq_len(ncol(cnf$keyLogFsta))-1
    f1 <- max(min(dat$minAgePerFleet), x@range["minfbar"])
    f2 <- min(max(dat$maxAgePerFleet), x@range["maxfbar"])
    cnf$fbarRange <- c(f1,f2)
    cnf$stockRecruitmentModelCode <- sr
    pl <- defpar(dat,cnf)
    if(!missing(rec_pars))
        if(length(rec_pars) == length(pl$rec_pars)){
            pl$rec_pars <- rec_pars
        }else{
            warning("rec_pars wrong length, ignored!")
        }
    sam.fit(dat,cnf,pl,run=FALSE)
}

as.sam.tmbargs.FLStock <- function(x, sr, rec_pars, makePlusGroup = TRUE){
    sim.condRE <- TRUE
    penalizeSpline <- FALSE
    if(missing(sr))
        sr <- 0
    data <- getSamFitData(x, makePlusGroup)
    conf <- defcon(data)
    conf$keyLogFsta[1,] <- seq_len(ncol(conf$keyLogFsta))-1
    f1 <- max(min(data$minAgePerFleet), x@range["minfbar"])
    f2 <- min(max(data$maxAgePerFleet), x@range["maxfbar"])
    conf$fbarRange <- c(f1,f2)
    conf$stockRecruitmentModelCode <- sr
    parameters <- defpar(data,conf)
    if(!missing(rec_pars))
        if(length(rec_pars) == length(parameters$rec_pars)){
            parameters$rec_pars <- rec_pars
        }else{
            warning("rec_pars wrong length, ignored!")
        }
     if(length(conf$maxAgePlusGroup)==1){
    tmp <- conf$maxAgePlusGroup    
    conf$maxAgePlusGroup <- defcon(data)$maxAgePlusGroup
    conf$maxAgePlusGroup[1] <- tmp
  }
  definit <- defpar(data, conf)
  if(!identical(parameters,relist(unlist(parameters), skeleton=definit))){
    warning("Initial values are not consistent, so running with default init values from defpar()")
    parameters<-definit
  }
  data<-clean.void.catches(data,conf)
  
  confTmp = defcon(data)
  for(i in 1:length(confTmp)){
    if(!names(confTmp)[i] %in% names(conf)){
      conf[[length(conf)+1]] = confTmp[[i]]
      names(conf)[length(conf)] = names(confTmp)[i]
    }
  }
  
  tmball <- c(data, conf, list(simFlag=rep(as.integer(sim.condRE),length = 2)))
  if(is.null(tmball$resFlag)){tmball$resFlag <- 0}  
  nmissing <- sum(is.na(data$logobs))
  parameters$missing <- numeric(nmissing)
  ran <- c("logN", "logF", "missing", "logSW", "logCW", "logitMO", "logNM")
  if(penalizeSpline)
      ran <- c(ran, "rec_pars")
  
  args <- c(list(data = tmball,
                 parameters = parameters,
                 random = ran,
                 DLL = "stockassessment"))

  mapRP <- list(logFScaleMSY = factor(NA),
                implicitFunctionDelta = factor(NA),
                logScaleFmsy = factor(NA),
                logScaleFmypyl = factor(NA),
                logScaleFmdy = factor(NA),
                logScaleFmax = factor(NA),
                logScaleFxdYPR = factor(rep(NA,length(args$parameters$logScaleFxdYPR))),
                logScaleFxB0 = factor(rep(NA,length(args$parameters$logScaleFxB0))),
                logScaleFcrash = factor(NA),
                logScaleFext = factor(NA),
                logScaleFxPercent = factor(rep(NA,length(args$parameters$logScaleFxPercent))),
                logScaleFlim = factor(NA),
                splinePenalty = factor(ifelse(penalizeSpline,1,NA))
                )
  if(is.null(args$map) || !is.list(args$map) || length(args$map) == 0){
      args$map <- mapRP
  }else{
      args$map <- c(args$map, mapRP)
  }
  
  if(!is.null(conf$hockeyStickCurve))
      if(is.null(args$map$rec_pars) &
         !is.na(conf$hockeyStickCurve) &
         conf$stockRecruitmentModelCode == 63)
          args$map$rec_pars = factor(c(1,2,NA))
    args
}
