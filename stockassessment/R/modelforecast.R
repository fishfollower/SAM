##' Parallel replicate for modelforecast
##'
##' @param n number of replicates
##' @param expr expression
##' @param simplify simplify passes to sapply
##' @param ncores number of cores
##' @return output
##' @importFrom parallel makeCluster clusterSetRNGStream parSapply stopCluster
.SAM_replicate <- function(n, expr, simplify = "array", ncores = 1, env = parent.frame(n+1)){
    if(ncores > 1){
        cl <- parallel::makeCluster(ncores)
        on.exit(parallel::stopCluster(cl))
        parallel::clusterSetRNGStream(cl)
        eval(expression(2+2),env)
        parallel::clusterExport(cl, "env", environment())
        v <- parallel::parSapply(cl, integer(n), eval(substitute(function(...) expr),env), 
                                 simplify = simplify)        
    }else{
        v <- sapply(integer(n), eval(substitute(function(...) expr),env), 
                    simplify = simplify)
    }
    v

}

.forecastDefault <- function(){
    ## list(specification = "DEFAULT_NOT_TO_BE_USED",
    ##      Amin = NA_real_,
    ##      Amax = NA_real_,
    ##      fleet = NA_real_,
    ##      relative = 0,
    ##      cstr = as.numeric(99),
    ##      target = NA_real_)
    list()
}

.parseRel <- function(type,s){
    if(s == ""){                        #Absolute value
        v <- -3
        attr(v,"lag") <- 1
        return(v)
    }
    if(s == "*"){                      #Relative to last year
        v <- -2
        attr(v,"lag") <- 1
        return(v)
    }
    ## Check relative to same type of value
    if(grepl("(\\*)(+|-)([[:digit:]]+)",s)){
        v <- -2
        sign <- 1
        if(gsub("(\\*)(+|-)([[:digit:]]+)","\\2",s)=="+")
            sign <- -1
        attr(v,"lag") <- sign * as.numeric(gsub("(\\*)(+|-)([[:digit:]]+)","\\3",s))
    }
    if(as.character(type) != gsub("(\\*)(F|C|SSB|TSB|L)(.*)","\\2",s))
        stop(sprintf("%s values can only be relative to %s values",as.character(type),as.character(type)))
    ## Find relative type
    if(grepl("^\\*(F|C|SSB|TSB|L)$",s)){ #Relative to total
        v <- -1
        attr(v,"lag") <- 1
        return(-1)
    }else if(grepl("^\\*(F|C|SSB|TSB|L)(\\[[[:digit:]]+\\])$",s)){
        v <- as.numeric(gsub("(^\\*(F|C|SSB|TSB|L))(\\[)([[:digit:]]+)(\\]$)","\\4",s))-1
        attr(v,"lag") <- 1
        return(v)
    }
    stop("Wrong specification of relative value. Relative values can only specify a fleet")
}

.parseForecast <- function(s, FbarRange, fleetTypes, ageRange, useNonLinearityCorrection){
    nCF <- sum(fleetTypes == 0)
    ## Split constraints
    sL <- strsplit(s,"[[:space:]]*&+[[:space:]]*")
    forecastEnum <- c("F","C","SSB","TSB","L","KeepRelF","HCR")
    getAgeFleetSettings <- function(spec, type){
        if(spec == ""){                 # Not specified
            Amin <- as.numeric(ifelse(type=="F",FbarRange[1],ageRange[1]))
            Amax <- as.numeric(ifelse(type=="F",FbarRange[2],ageRange[2]))
            fleet <- -1
        }else if(grepl("\\[[[:digit:]]+,[[:digit:]]+-[[:digit:]]+\\]",spec)){ # Fully specified
            Amin <- as.numeric(gsub("(\\[)([[:digit:]]+)(,)([[:digit:]]+)(-)([[:digit:]]+)(\\])","\\4",spec))
            Amax <- as.numeric(gsub("(\\[)([[:digit:]]+)(,)([[:digit:]]+)(-)([[:digit:]]+)(\\])","\\6",spec))
            fleet <- as.numeric(gsub("(\\[)([[:digit:]]+)(,)([[:digit:]]+)(-)([[:digit:]]+)(\\])","\\2",spec))-1
        }else if(grepl("\\[[[:digit:]]+\\]",spec)){ # Fleet specified
            Amin <- as.numeric(ifelse(type=="F",FbarRange[1],ageRange[1]))
            Amax <- as.numeric(ifelse(type=="F",FbarRange[2],ageRange[2]))
            fleet <- as.numeric(gsub("(\\[)([[:digit:]]+)(\\])","\\2",spec))-1
        }else if(grepl("\\[[[:digit:]]+-[[:digit:]]+\\]",spec)){ # Age range specified
            Amin <- as.numeric(gsub("(\\[)([[:digit:]]+)(-)([[:digit:]]+)(\\])","\\2",spec))
            Amax <- as.numeric(gsub("(\\[)([[:digit:]]+)(-)([[:digit:]]+)(\\])","\\4",spec))
            fleet <- -1
        }else{ #Wrongly specified
            stop("Wrong specification of forecast")
        }
        list(Amin = Amin, Amax = Amax, fleet = fleet)
    }
    parseHCR <- function(ss){
        ## HCR[settings]=value(F|C|L)[settings]~value(SSB|TSB)[settings]lag
        ## Match 1: HCR
        ## Match 2: HCR settings
        ## Match 3: =
        ## Match 4: Target value 
        ## Match 5: Target type (F, C, L)
        ## Match 6: Target settings (age range, fleet)
        ## Match 7: ~
        ## Match 8: Trigger value
        ## Match 9: Trigger type (SSB, TSB)
        ## Match 10: Trigger settings (age range)
        ## Match 11: Trigger lag
        frmt <- "^(HCR)(\\[.+\\])?(=)([[:digit:]]*\\.?[[:digit:]]+)(F|C|L)?(\\[.+\\])?(~)([[:digit:]]*\\.?[[:digit:]]+)(SSB|TSB|B0)?(\\[.+\\])?(-[[:digit:]]*)?$"
        if(!grepl(frmt,ss))
            stop("Wrong specification of HCR")
        type <- factor(gsub(frmt,"\\1",ss),forecastEnum)
        rel <- -3
        ## Get target info
        target <- log(as.numeric(gsub(frmt,"\\4",ss)))
        if(is.na(target))
            stop("wrong specification of forecast target")
        targetType <- gsub(frmt,"\\5",ss)
        if(targetType == "")
            targetType <- "F"
        targetTypeI <- match(targetType, c("F","C","L"))-1
        if(is.na(targetTypeI))
            stop("wrong specification of forecast target type")
        targetSpec <- gsub(frmt,"\\6",ss)
        targetSettings <- getAgeFleetSettings(targetSpec,targetType)
        ## Get trigger info (Need to trick gsub)     
        frmtTrig <- "^([[:digit:]]*\\.?[[:digit:]]+)(SSB|TSB|B0)?(\\[.+\\])?(-[[:digit:]]*)?$"
        ssTrig <- gsub(".+~","",ss)
        if(!grepl(frmtTrig,ssTrig))
            stop("Wrong specification of HCR trigger")     
        triggerValue <- as.numeric(gsub(frmtTrig,"\\1",ssTrig))
        triggerType <- gsub(frmtTrig,"\\2",ssTrig)
        if(triggerType == "")
            triggerType <- "SSB"
        triggerTypeI <- match(triggerType, c("SSB","TSB","B0"))-1
        triggerSpec <- gsub(frmtTrig,"\\3",ssTrig)
        triggerSettings <- getAgeFleetSettings(triggerSpec,triggerType)
        triggerLag <- gsub(frmtTrig,"\\4",ssTrig)
        if(triggerLag == ""){
            triggerLag <- 0
        }else{
            triggerLag <- -as.numeric(triggerLag)
        }
        ## Get other HCR info
        hcrSpec <- lapply(strsplit(gsub("(\\[|\\])","",gsub(frmt,"\\2",ss)),",")[[1]],function(x) strsplit(x,"=")[[1]])
        hcrVal <- lapply(hcrSpec, function(x) as.numeric(x[2]))
        names(hcrVal) <- sapply(hcrSpec, function(x) x[1])
        hcrConf <- list(FO = 1e-9, FC = 1e-9, BO = 0, BC = 0)
        nms <- intersect(names(hcrConf),names(hcrVal))
        hcrConf[nms] <- hcrVal[nms]
        list(specification = ss,
             Amin = targetSettings$Amin,
             Amax = targetSettings$Amax,
             fleet = targetSettings$fleet,
             relative = rel,
             cstr = as.numeric(type)-1,
             target = target,
             settings = c(biomassLag = triggerLag,
                          biomassType=triggerTypeI,
                          targetType = targetTypeI,
                          triggerAmin = triggerSettings$Amin,
                          triggerAmax = triggerSettings$Amax,
                          Forigin = hcrConf$FO,
                          Fcap = hcrConf$FC,
                          Borigin = hcrConf$BO,
                          Bcap = hcrConf$BC,
                          Btrigger = triggerValue),
             useNonLinearityCorrection = as.logical(useNonLinearityCorrection)
             )
    }
    parseOne <- function(ss){
        if(grepl("^HCR",ss))
            return(parseHCR(ss))
        frmt <- "^(F|C|SSB|TSB|L)(\\[.+\\])?(=)([[:digit:]]*\\.?[[:digit:]]*)(\\*.*)?$"
        if(!grepl(frmt,ss))
            stop("Wrong specification of forecast")      
        type <- factor(gsub(frmt,"\\1",ss),forecastEnum)
        spec <- gsub(frmt,"\\2",ss)
        target <- log(as.numeric(gsub(frmt,"\\4",ss)))
        if(is.na(target))
            stop("wrong specification of forecast target")
        rel <- .parseRel(type,gsub(frmt,"\\5",ss))
        lag <- attr(rel,"lag")
        if(is.null(lag))
            lag <- 1
        if(is.na(rel) || !is.numeric(rel) || rel < -3){
            stop("wrong specification of relative value")
        }
        fa <- getAgeFleetSettings(spec, type)
        Amin <- fa$Amin
        Amax <- fa$Amax
        fleet <- fa$fleet
        if((type == "SSB" || type == "TSB") && fleet >= 0){
            warning("For SSB and TSB constraints, fleet should not be set")
            fleet <- -1
        }
        if(Amin < ageRange[1] || Amax > ageRange[2])
            stop("Specified ages not in age range of model")
        if(fleet >= nCF)
            stop("fleet number is too high")
        if(rel >= nCF)
            stop("relative fleet number is too high")
        if(rel == fleet)
            stop("A fleet cannot be relative to itself")
        v <- list(specification = ss,
                  Amin = Amin,
                  Amax = Amax,
                  fleet = fleet,
                  relative = rel,
                  cstr = as.numeric(type)-1,
                  target = target,
                  settings = c(compareLag = lag),
                  useNonLinearityCorrection = as.logical(useNonLinearityCorrection))
        v
    }
    cstr <- lapply(sL,function(x)lapply(x,parseOne))
    ## Add constraints to fully specify
    fleetsSpec <- lapply(cstr,function(xx) factor(sapply(xx,function(yy)yy$fleet),c(-1,seq_len(nCF)-1)))
    if(any(unlist(sapply(fleetsSpec,is.na))))
        stop("Fleet numbers should be from 1 to the number of residual catch fleets")
    ft <- lapply(fleetsSpec,table)
    if(any(sapply(ft, function(x) any(x>1) || sum(x) > nCF)))
        stop("Each fleet/total can only have one constraint and there can not be more constraints than fleets")
    makeRelCon <- function(f1,f2){
        list(specification = sprintf("AUTOGENERATED~F%d/F%d",f2,f1),
             Amin = FbarRange[1],
             Amax = FbarRange[2],
             fleet = as.numeric(f1)-1,
             relative = as.numeric(f2)-1,
             cstr = as.numeric(factor("KeepRelF",forecastEnum))-1,
             target = 0,
             settings = c(compareLag = 1),
             useNonLinearityCorrection = as.logical(useNonLinearityCorrection))  
    }
    cToAdd <- lapply(ft,function(x){
        ii <- unname(which(x[-1]==0))
        if(length(ii) > 1){
            return(lapply(ii[-1],makeRelCon,f1=ii[1]))
        }
        list()
    })
    lapply(seq_along(cstr), function(i){
        c(cstr[[i]], cToAdd[[i]])
    })
}




##' Model based forecast function
##' @rdname modelforecast 
##' @param fit an assessment object of type sam, as returned from the function sam.fit
##' @param ... other variables used by the methods
##' @export
modelforecast <- function(fit, ...){
    UseMethod("modelforecast")
}

##' Model based forecast function
##' @param fit SAM model fit
##' @param constraints a character vector of forecast constraint specifications
##' @param fscale a vector of f-scales. See details.  
##' @param catchval a vector of target catches. See details "old specification".
##' @param fval a vector of target f values. See details "old specification".
##' @param nextssb a vector target SSB values the following year. See details "old specification".
##' @param landval a vector of target catches. See details "old specification".   
##' @param nosim number of simulations. If 0, the Laplace approximation is used for forecasting.
##' @param year.base starting year default last year in assessment. Currently it is only supported to use last assessment year or the year before  
##' @param ave.years vector of years to average for weights, maturity, M and such  
##' @param overwriteBioModel Overwrite GMRF models with ave.years?
##' @param rec.years vector of years to use to resample recruitment from. If the vector is empty, the stock recruitment model is used.
##' @param label optional label to appear in short table
##' @param overwriteSelYears if a vector of years is specified, then the average selectivity of those years is used (not recommended)
##' @param deterministicF option to set F variance to (almost) zero (not recommended)
##' @param processNoiseF option to turn off process noise in F
##' @param fixedFdeviation Use a fixed F deviation from target?
##' @param useFHessian Use the covariance of F estimates instead of the estimated process covariance for forecasting?
##' @param resampleFirst Resample base year when nosim > 0?
##' @param customSel supply a custom selection vector that will then be used as fixed selection in all years after the final assessment year (not recommended)
##' @param lagR if the second youngest age should be reported as recruits
##' @param splitLD if TRUE the result is split in landing and discards
##' @param addTSB if TRUE the total stock biomass (TSB) is added
##' @param biasCorrect Do bias correction of reported variables. Can be turned off to reduce running time (not recommended).
##' @param returnAllYears If TRUE, all years are bias corrected. Otherwise, only forecast years are corrected.
##' @param returnObj Only return TMB object?
##' @param progress Show progress bar for simulations?
##' @param estimate the summary function used (typically mean or median) for simulations
##' @param silent Passed to MakeADFun. Should the TMB object be silent?
##' @param newton_config Configuration for newton optimizer to find F values. See ?TMB::newton for details. Use NULL for TMB defaults.
##' @param custom_pl Parameter list. By default, the parameter list from fit is used.
##' @param useNonLinearityCorrection Should a non linearity correction be added to transformation of logF? See Details - Non-linearity correction.
##' @details
##' Function to forecast the model under specified catch constraints. In the forecast, catch constraints are used to set the mean of the \eqn{log(F)} process for each simulation. Therefore, catch constraints are not matched exactly in individual simulations. Likewise, the summary of a specific set of simulations will not match exactly due to random variability.
##' By default, recruitment is forecasted using the estimated recruitment model. If a vector of recruitment years is given, recruitment is forecasted using a log-normal distribution with the same mean and variance as the recruitment in the years given. This is different from the forecast function, which samples from the recruitment estimates.
##' Catch scenarios are specified by a vector of target constraints. The first value determines F in the year after the base year.
##'
##' @section Forecast constraints:
##'
##' \subsection{F based constraints:}{
##' Forecasts for F values are specified by the format "F[f,a0-a1]=x" where f is the residual catch fleet and a0-a1 is an age range. For example, "F[2,2-4]=0.3" specifies that the average F for the second fleet over ages 2-4 should be 0.3.
##' If an "*" is added to the target value, the target will be relative to the year before. For example, "F[2,2-4]=0.9*" specifies that the average F for the second fleet over ages 2-4 should be 90% of the year before. Further, the target for a fleet can be relative to the total by adding "*F" or to another fleet by adding "*F[f]" where f is the fleet number. The same age range will always be used.
##' If the fleet is omitted (e.g., F[2-4]), the target is for the total F.
##' If the age range is omitted (e.g., F[2]), the fbar range of the model is used.
##' Likewise, both fleet and age range can be omited (e.g., F=0.3) to specify a value for total F with the range used in the model.
##'
##' For example:
##' \describe{
##' \item{"F=0.2"}{Will set the median average total fishing mortality rate to 0.2}
##' \item{"F[1]=0.2"}{Will set the median average fishing mortality rate of the first fleet to 0.2}
##' \item{"F[2-4]=0.2"}{Will set the median average total fishing mortality rate over ages 2 to 4 to 0.2}
##' \item{"F[3,2-4]=0.2"}{Will set the median average fishing mortality rate over ages 2 to 4 for the third fleet to 0.2}
##' }
##' }
##'
##' \subsection{Catch/Landing based constraints:}{
##' Forecasts for catch and landing values are specified by the format "C[f,a0-a1]=x" for catch and "L[f,a0-a1]" for landings. If the age range is omitted, all modelled ages are used. Otherwise, the format is similar to F based scenarios.
##' If an "*" is added to the target value, the target will be relative to the year before.
##' Further, the catch target for a fleet can be relative to the total by adding "*C" or to another fleet by adding "*C[f]" where f is the fleet number. The same age range will always be used. Likewise, relative landing targets can be specified using "*", "*L", or "*L[f]" for targets relative to last year, the total, or fleet f, respectively.
##'
##' For example:
##' \describe{
##' \item{"C=100000"}{Will scale F such that the total predicted catch is 100000}
##' \item{"C[1]=100000"}{Will scale F such that the predicted catch of the first fleet is 100000}
##' \item{"C[2-4]=100000"}{Will scale F such that the total predicted catch for ages 2 to 4 is 100000}
##' \item{"C[3,2-4]=100000"}{Will scale F such that the predicted catch for ages 2 to 4 in the third fleet is 100000}
##' \item{"L=100000"}{Will scale F such that the total predicted landing is 100000}
##' \item{"L[1]=100000"}{Will scale F such that the predicted landing of the first fleet is 100000}
##' \item{"L[2-4]=100000"}{Will scale F such that the total predicted landing for ages 2 to 4 is 100000}
##' \item{"L[3,2-4]=100000"}{Will scale F such that the predicted landing for ages 2 to 4 in the third fleet is 100000}
##' }
##'}
##' 
##' \subsection{Next year's SSB/TSB based constraints:}{
##' Forecasts for spawning stock biomass (SSB) and total stock biomass (TSB) values are specified by the format "SSB[a0-a1]=x" for SSB and "TSB[a0-a1]" for TSB. For setting F in year y, the relevant biomass for year y+1 is predicted for the constraint. If spawning is not at the beginning of the year, F is assumed to be the same for year y and y+1 in the prediction.
##' The format is similar to catch/landing based scenarios. However, fleets have no effect. If an age range is omitted, the full age range of the model is used.
##' If an "*" is added to the target value, the target will be relative to the year before. That is, when setting F in year y, the predicted biomass in year y+1 will be relative to the biomass in year y-1.
##' Note that since SSB and TSB used for catch constraints are predicted, the input constraint will differ from the output SSB and TSB estimates due to process variability.
##'
##' For example:
##' \describe{
##' \item{SSB=200000}{Will scale F such that the predicted SSB at the beginning of the next year is 200000}
##' \item{SSB[3-9]=200000}{Will scale F such that the predicted SSB for ages 3 to 9 at the beginning of the next year is 200000}
##' \item{TSB=200000}{Will scale F such that the predicted TSB at the beginning of the next year is 200000}
##' \item{TSB[3-9]=200000}{Will scale F such that the predicted TSB for ages 3 to 9 at the beginning of the next year is 200000}
##'}
##' }
##' 
##' \subsection{Harvest control rule based constraints:}{
##' Harvest control rules can be specified for forecasts using the format "HCR=x~y" where x is the target and y is the biomass trigger (see ?hcr for full details on the form of the harvest control rule). Further, the target can be specified as an F target ("HCR=xF~y"), catch target ("HCR=xC~y"), or landing target ("HCR=xL~y"). Likewise the trigger can either be for SSB ("HCR=x~ySSB") or TSB ("HCR=x~yTSB"). Age ranges can be set for both triggers and targets and a fleet can be set for the target. The notation and defaults are similar to the F based and SSB/TSB based constraints, respectively.
##'When setting F in year y, the projected biomass in year y is used by default. To use the (at this time known) biomass in a previous year, a time lag can be specified. To specify a time lag of, e.g., 1 year for SSB the format is "HCR=x~ySSB-1".
##' Finally, the origin and cap for the HCR can be set using "HCR[FO=a,FC=b,BO=d,BC=e]=x~y", where FO is the F (or catch or landing) value at origin, BO is the biomass at origin, FC is the F (or catch or landing) value when the HCR is capped and BC is the biomass at which the HCR is capped. See ?hcr for further details on the shape of the HCR.
##' For a HCR similar to the ICES advice rule, the specification is on the form "HCR[BC=Blim] = fmsy~MSYBtrigger". Note that, unlike an ICES advice rule, the HCR does not do a forecast to determine if fishing can continue below Blim.
##'
##' For example:
##' \describe{
##' \item{HCR=0.9~100000}{Will apply a harvest control rule with an F target of 0.9 and a biomass trigger of 100000 on SSB}
##' \item{HCR=10000C~100000}{Will apply a harvest control rule with a catch target of 10000 and a biomass trigger of 100000 on SSB}
##' \item{HCR=0.9~100000SSB}{Will apply a harvest control rule with an F target of 0.9 and a biomass trigger of 100000 on SSB}
##' \item{HCR=0.9F[1,2-4]~100000SSB}{Will apply a harvest control rule with an F target on the first fleet ages 2-4 of 0.9 and a biomass trigger of 100000 on SSB}
##' \item{HCR=0.9~100000TSB[0-4]}{Will apply a harvest control rule with an F target of 0.9 and a biomass trigger of 100000 on TSB for ages 0 to 4}
##' \item{HCR[FC=1e-9,BC=20000]=0.9~100000}{Will apply a harvest control rule with an F target of 0.9 and a biomass trigger of 100000 on SSB where biomass values below 20000 will give an F of 1e-9}
##' \item{HCR[FO=0,BO=30000]=0.9~100000}{Will apply a harvest control rule with an F target of 0.9 and a biomass trigger of 100000 on SSB where the slope on which F is reduced goes to zero F at a biomass of 30000}
##' }
##'}
##' 
##' \subsection{Combining constraints:}{
##' Constraints for different fleets can be combined by "&".
##' For example, "F[2-4]=0.5 & C[2]=10000" specifies that total Fbar over ages 2-4 should be 0.5 while the catch for the second residual catch fleet should be 10,000t.
##' The constraints cannot affect within-fleet selectivity. Therefore, a fleet can at most have one constraint per year, and the total number of constraints cannot exceed the number of catch fleets. That is, if a constraint is given for the sum of fleets, there must be at least one fleet without any constraints.
##' For fleets where no constraints are given, a constraint is set to keep their relative Fs constant.
##' }
##'
##' \subsection{Values relative to previous year:}{
##' Catch constraints specified as specific values are inherently different from catch constraints specified as relative values, even if they lead to the same F. Catch constraints specified as relative values will propagate the uncertainty in, e.g, F from previous years whereas constraints specified as specific values will not. This is different from the \link{forecast} function where, for example, a forecast using fval is the same as a forecast using fscale, if they lead to the same F. 
##'}
##'
##'
##' ##' \subsection{Process variability:}{
##' In the forecast, constraints are used to set the predicted F value in year y based on information available until year y-1. Therefore, constraints using predicted values for year y, such as catch, will not be matched exactly by the realized catch due to process variability in F, N, biological processes and catch itself.
##'}
##'
##' @section Non-linearity correction:
##' In the model forecasts, constraints are calculated to set the mean of the log(F) process, corresponding to the median F-at-ages. Typically, the constraints are non-linear functions of log(F)-at-age. Therefore, when stochasticity is added to log(F) (i.e., deterministicF=FALSE), target values will correspond to a transformation of the median, and not the median of the transformation. For example, a target for the average fishing mortality (Fbar) will correspond to the average of the median F at age, which will be different from the median Fbar.
##'
##' The "useNonLinearityCorrection" argument can be used to shift the target from a function of the mean log(F) (median F) towards the log-mean of the function of log(F), which is approximately the median of the function of log(F). 
##' 
##' @section Old specification:
##' It is also possible to specify forecast constraints in a way similar to the \link{forecast} function. 
##' There are four ways to specify a scenario. If e.g. four F values are specified (e.g. fval=c(.1,.2,.3,4)), then the first value is used in the year after the last assessment year (base.year + 1), and the three following in the three following years. Alternatively F's can be specified by a scale, or a target catch. Only one option can be used per year. So for instance to set a catch in the first year and an F-scale in the following one would write catchval=c(10000,NA,NA,NA), fscale=c(NA,1,1,1). If only NA's are specified in a year, the F model is used for forecasting. The length of the vector specifies how many years forward the scenarios run. Unlike the forecast function, no value should be given for the base year.
##' Internally, the old specification is translated such that "fval=x" becomes "F=x", "fscale=x" becomes "F=x*", "catchval=x" becomes "C=x", "nextssb=x" becomes "SSB=x", and "landval=x" becomes "L=x".
##'
##' @section Forecasts using Laplace approximation or simulations:
##' Forecasts can be made using either a Laplace approximation projection (nosim=0) or simulations (nosim > 0). When using the Laplace approximation, the most likely projected trajectory of the processes along with a confidence interval is returned. In contrast, simulation based forecasts will return individual simulated trajectories and summarize using the function given as the estimate argument along with an interval covering 95% of the simulations.
##'
##' @section Warnings:
##' Long term forecasts with random walk recruitment can lead to unstable behaviour and difficulties finding suitable F values for the constraints.
##' If no suitable F value can be found, an error message will be shown, and F values will be NA or NaN. Likewise, forecasts leading to high F values in some years (or large changes from one year to another) may cause problems for the optimization as they will be used as starting values for the next years.
##' Since the model works on log space, all target values should be strictly positive. Values too close to zero may cause problems.
##'
##'
##' 
##' @return an object of type samforecast
##' @seealso forecast
##' @importFrom methods formalArgs
##' @importFrom stats var
##' @importFrom utils tail txtProgressBar setTxtProgressBar
##' @rdname modelforecast
##' @method modelforecast sam
##' @export
modelforecast.sam <- function(fit,
                              constraints = NULL,
                              fscale = NULL,
                              catchval = NULL,
                              fval = NULL,
                              nextssb = NULL,
                              landval = NULL,
                              nosim = 0,
                              year.base = max(fit$data$years),
                              ave.years = max(fit$data$years)+(-9:0),
                              overwriteBioModel = FALSE,
                              rec.years = c(), #max(fit$data$years)+(-9:0),
                              label = NULL,
                              overwriteSelYears = NULL,
                              deterministicF = FALSE,
                              processNoiseF = FALSE,
                              fixedFdeviation = FALSE,                              
                              useFHessian = FALSE,
                              resampleFirst = !is.null(nosim) && nosim > 0,
                              fixFirstN = TRUE,
                              customSel = NULL,
                              lagR = FALSE,
                              splitLD = FALSE,
                              addTSB = FALSE,
                              biasCorrect = FALSE,
                              returnAllYears = FALSE,
                              returnObj = FALSE,
                              progress = TRUE,
                              estimate = median,
                              silent = TRUE,
                              newton_config = NULL,
                              custom_pl = NULL,
                              useNonLinearityCorrection = (nosim > 0 && !deterministicF),
                              ncores = 1,
                              ...
                              ){
    ## Check for hcr, findMSY, hcrConf, hcrCurrentSSB
    dots <- list(...)
    findMSY  <-  NULL
    if(!is.na(match("findMSY",names(dots))))
        findMSY <- dots[[match("findMSY",names(dots))]]
    hcr  <-  NULL
    if(!is.na(match("hcr",names(dots))))
        hcr <- dots[[match("hcr",names(dots))]]
    hcrConf = numeric(0)
    if(!is.na(match("hcrConf",names(dots))))
        hcrConf <- dots[[match("hcrConf",names(dots))]]
    hcrCurrentSSB = 0
    if(!is.na(match("hcrCurrentSSB",names(dots))))
        hcrCurrentSSB <- dots[[match("hcrCurrentSSB",names(dots))]]
    nCatchAverageYears  <-  1
    if(!is.na(match("nCatchAverageYears",names(dots))))
        nCatchAverageYears <- dots[[match("nCatchAverageYears",names(dots))]]
    fastFixedF <- FALSE
    if(!is.na(match("fastFixedF",names(dots))))
        fastFixedF <- dots[[match("fastFixedF",names(dots))]]
    

    if(!is.null(nosim) && nosim > 0){ 
        estimateLabel <- paste(deparse(substitute(estimate), 500L), collapse = " ")
    }else{
        estimateLabel <- "mostLikelyTrajectory"
    }
 
    
    if(progress && !returnObj && !is.null(nosim) && nosim > 0){
        pb <- utils::txtProgressBar(min = 0, max = nosim+3, style = 3)
        incpb <- function() utils::setTxtProgressBar(pb, pb$getVal()+1)
    }else{
        incpb <- function(){ return(invisible(NULL)) }
    }

    ## Handle year.base < max(fit$data$years)
    if(year.base > max(fit$data$years)){
        stop("")
    }else if(year.base < max(fit$data$years)){
        ## warning("year.base is ignored for now")
    }

    ## if(length(ave.years) == 0){
    useModelBio <- !overwriteBioModel
    ##     ave.years = max(fit$data$years)+(-9:0)
    ## }else{
    ##     useModelBio <- FALSE
    ## }

    ## Checks
    if(deterministicF && length(fscale) > 0 && any(!is.na(fscale)) && is.null(customSel))
        warning("Forecasted F values may depend on the last estimated F vector and can therefore have uncertainty. Provide a custom selectivity to circumvent this.")
    
    
    ## Get number of forecast years
    lengthVec <- c(length(constraints),
                   length(fscale),
                   length(catchval),
                   length(fval),length(nextssb),length(landval),
                   length(hcr))
    if(any(lengthVec > 0 & lengthVec < max(lengthVec)))
        stop("All target vectors must have the same length")    
    nYears <- max(length(constraints),length(fscale),length(catchval),length(fval),length(nextssb),length(landval),length(findMSY), length(hcr))

    ## Convert input to an F model code and a target value
    if(is.null(constraints))
        constraints <- rep(NA_character_, nYears)    

    if(is.null(fscale))
        fscale <- rep(NA_real_, nYears)
    if(any(!is.na(constraints) & !is.na(fscale)))
        warning("fscale specified for years with other constraints. Using previously defined constraints.")
    constraints[is.na(constraints) & !is.na(fscale)] <- sprintf("F=%f*",fscale[is.na(constraints) & !is.na(fscale)])

    if(is.null(catchval))
        catchval <- rep(NA_real_, nYears)
    if(any(!is.na(constraints) & !is.na(catchval)))
        warning("catchval specified for years with other constraints. Using previously defined constraints.")
    constraints[is.na(constraints) & !is.na(catchval)] <- sprintf("C=%f",catchval[is.na(constraints) & !is.na(catchval)])

    if(is.null(fval))
        fval <- rep(NA_real_, nYears)
    if(any(!is.na(constraints) & !is.na(fval)))
        warning("fval specified for years with other constraints. Using previously defined constraints.")
    constraints[is.na(constraints) & !is.na(fval)] <- sprintf("F=%f",fval[is.na(constraints) & !is.na(fval)])
    
    if(is.null(nextssb))
        nextssb <- rep(NA_real_, nYears)
    if(any(!is.na(constraints) & !is.na(nextssb)))
        warning("nextssb specified for years with other constraints. Using previously defined constraints.")
constraints[is.na(constraints) & !is.na(nextssb)] <- sprintf("SSB=%f",nextssb[is.na(constraints) & !is.na(nextssb)])

    if(is.null(landval))
        landval <- rep(NA_real_, nYears)
    if(any(!is.na(constraints) & !is.na(landval)))
        warning("landval specified for years with other constraints. Using previously defined constraints.")
    constraints[is.na(constraints) & !is.na(landval)] <- sprintf("L=%f",landval[is.na(constraints) & !is.na(landval)])

    if(is.null(findMSY))
        findMSY <- rep(NA_real_, nYears)
    if(any(!is.na(constraints) & !is.na(findMSY)))
        warning("findMSY specified for years with other constraints. findMSY will overwrite other values.")

    if(is.null(hcr))
        hcr <- rep(NA_real_, nYears)
    if(any(!is.na(constraints) & !is.na(hcr)))
        warning("hcr specified for years with other constraints. hcr will overwrite other values.")
          
    ## tab <- rbind(fscale,fval,catchval,nextssb,landval, findMSY, hcr)
    ## FModel <- apply(tab,2, function(x){
    ##     y <- which(!is.na(x))
    ##     switch(as.character(length(y)),
    ##            "0"=0,
    ##            "1"=y,
    ##            stop("At most one target can be specified per year."))
    ## })
    ## target <- apply(tab,2, function(x){
    ##     y <- which(!is.na(x))
    ##     switch(as.character(length(y)),
    ##            "0"=NA_real_,
    ##            "1"=x[y],
    ##            stop("At most one target can be specified per year."))
    ## })
    FModel <- rep(0,nYears)
    FModel[!is.na(constraints)] <- 1
    if(fastFixedF)
        FModel[!is.na(constraints)] <- 2
    FModel[!is.na(findMSY)] <- 3
    FModel[!is.na(hcr)] <- 4
    
    nYears <- length(constraints)
    cstr <- replicate(nYears, .forecastDefault(), simplify = FALSE)
    cstr[!is.na(constraints)] <- .parseForecast(constraints[!is.na(constraints)], fit$conf$fbarRange, fit$data$fleetTypes, c(fit$conf$minAge,fit$conf$maxAge), useNonLinearityCorrection)

    ## Use custom selectivity?
    if(is.null(customSel)){
        customSel <- numeric(0)
        if(!is.null(overwriteSelYears)){
            fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)  
            Ftab <- faytable(fit)
            customSel <- colMeans(Ftab[as.integer(rownames(Ftab)) %in% overwriteSelYears, , drop=FALSE])
            customSel <- customSel/mean(customSel[fromto[1]:fromto[2]])
        }
    }else{
        if(!is.null(overwriteSelYears))
            warning("overwriteSelYears is ignored when customSel is given.")
        fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)  
        customSel <- customSel/mean(customSel[fromto[1]:fromto[2]])
    }

    ## Get recruitment model
    if(length(rec.years) == 0){
        recModel <- rep(0,nYears)
        logRecruitmentMedian <- NA_real_
        logRecruitmentVar <- NA_real_
    }else{
        rectab <- rectable(fit)
        recpool <- rectab[rownames(rectab)%in%rec.years,1]
        recModel <- rep(1,nYears)
        logRecruitmentMedian <- log(median(recpool))
        logRecruitmentVar <- stats::var(log(recpool))
    }

    ## Get F process time scale model
    ## By default, scale as random walk
    fsdTimeScaleModel <- rep(0,nYears)
    if(deterministicF){ ## 'Zero' variance of F process
        fsdTimeScaleModel <- rep(2,nYears)
    }else if(fixedFdeviation){
        fsdTimeScaleModel <- rep(3,nYears)
    }else if(!processNoiseF){ ## Constant variance of F process
        fsdTimeScaleModel <- rep(1,nYears)
    }
    ## When F model is used, fsdTimeScaleModel should be 1
    fsdTimeScaleModel[FModel == 0] <- 1


    ## Convert average years to indices
    ave.yearsIn <- ave.years
    ave.years <- match(ave.years, fit$data$years) - 1
    if(any(is.na(ave.years)))
        stop("ave.years has years without data.")

    ## Find base year number
    preYears <- match(year.base, fit$data$years)
    postYears <- nYears - (fit$data$noYears - preYears)

    ## Prepare forecast
    obj0 <- fit$obj
    invisible(obj0$fn(fit$opt$par))
    invisible(obj0$gr(fit$opt$par))

    args <- as.list(obj0$env)[methods::formalArgs(TMB::MakeADFun)[methods::formalArgs(TMB::MakeADFun) != "..."]]
    args$silent <- silent
    if(is.null(custom_pl)){
        pl <- fit$pl
    }else{
        pl <- custom_pl
    }
    pl$logF <- cbind(pl$logF,matrix(pl$logF[,ncol(pl$logF)],nrow(pl$logF),postYears))
    pl$logN <- cbind(pl$logN,matrix(pl$logN[,ncol(pl$logN)],nrow(pl$logN),postYears))
    if(useModelBio){
        splitArray <- function(a){
            nr <- dim(a)[1]; nc <- dim(a)[2]; na <- dim(a)[3]
            lapply(split(a,rep(seq_len(na),each=nr*nc)), matrix, nrow = nr, ncol = nc)
        }
        extendBio <- function(x) rbind(x,matrix(x[nrow(x)],postYears,ncol(x)))
        if(nrow(pl$logitMO) > 0){            
            pl$logitMO <- extendBio(pl$logitMO)
        }else{
            #warning(sprintf("No MO random effects. Using data average over %s.",paste(ave.yearsIn,collapse=", ")))
        }
        if(nrow(pl$logNM) > 0){            
            pl$logNM <- extendBio(pl$logNM)
        }else{
            #warning(sprintf("No NM random effects. Using data average over %s.",paste(ave.yearsIn,collapse=", ")))
        }
        if(nrow(pl$logSW) > 0){            
            pl$logSW <- extendBio(pl$logSW)
        }else{
            #warning(sprintf("No SW random effects. Using data average over %s.",paste(ave.yearsIn,collapse=", ")))
        }
        if(dim(pl$logCW)[1] > 0){            
            pl$logCW <- simplify2array(lapply(splitArray(pl$logCW), extendBio))
        }else{
            #warning(sprintf("No CW random effects. Using data average over %s.",paste(ave.yearsIn,collapse=", ")))
        }
    }
    d0 <- dim(pl$logitFseason)
    lfsOld <- pl$logitFseason
    pl$logitFseason <- array(0, c(d0[1],postYears+d0[2], d0[3]))
    pl$logitFseason[,1:d0[2],] <- lfsOld
    args$parameters <- pl
    args$random <- unique(names(obj0$env$par[obj0$env$random]))
    args$data$reportingLevel <- 0

    if(useFHessian){
        if(year.base==max(fit$data$years)){
            est <- fit$sdrep$estY
            FEstCov <- fit$sdrep$covY[grepl("LogF$",names(est)),grepl("LogF$",names(est))]
        }else if(year.base==(max(fit$data$years)-1)){
            est <- fit$sdrep$estYm1
            FEstCov <- fit$sdrep$covYm1[grepl("LogF$",names(est)),grepl("LogF$",names(est))]
        }else{
            stop("year.base not implemented yet more than one year before end of assessment.")           
        }
    }else{
        FEstCov <- matrix(0,0,0)
    }
    
    args$data$forecast <- list(nYears = as.numeric(nYears),
                               preYears = as.numeric(preYears),
                               nCatchAverageYears = as.numeric(nCatchAverageYears),
                               aveYears = as.numeric(ave.years),
                               forecastYear = as.numeric(c(rep(0,preYears),seq(1,nYears,length=nYears))),
                               FModel = as.numeric(FModel),
                               ##target = as.numeric(target),
                               constraints = cstr,                               
                               cfg = newton_config,
                               selectivity = as.numeric(customSel),
                               recModel = as.numeric(recModel),
                               logRecruitmentMedian = as.numeric(logRecruitmentMedian),
                               logRecruitmentVar = as.numeric(logRecruitmentVar),
                               fsdTimeScaleModel = as.numeric(fsdTimeScaleModel),
                               simFlag = c(0,0),
                               hcrConf = hcrConf,
                               hcrCurrentSSB = hcrCurrentSSB,
                               Fdeviation = rnorm(nrow(pl$logF)),
                               FdeviationCov = diag(1,nrow(pl$logF),nrow(pl$logF)),
                               FEstCov = FEstCov,
                               fixFirstN = fixFirstN)

    if(any(!is.na(findMSY))){
        args$map$logFScaleMSY <- NULL
        return(args)
    }
    ## Create forecast object
    if(!is.null(nosim) && nosim > 0){
        args$type <- "Fun"
        args$integrate <- NULL
        args$random <- NULL
    }

    ## Done with initial work chunk
    incpb()

    obj <- do.call(TMB::MakeADFun, args)
    ## Done with initial MakeADFun
    incpb()

    if(as.integer(returnObj)==1)
        return(obj)

    fleetHasF <- apply(fit$conf$keyLogFsta>-1,1,any)
    if(!is.null(nosim) && nosim > 0){
        if(year.base==max(fit$data$years)){
            est <- fit$sdrep$estY
            cov <- fit$sdrep$covY
        }else if(year.base==(max(fit$data$years)-1)){
            est <- fit$sdrep$estYm1
            cov <- fit$sdrep$covYm1
        }else{
            stop("year.base not implemented yet more than one year before end of assessment.")           
        }
        names(est) <- gsub("(^.*[lL]ast)(Log[NF]$)","\\2",names(est))
        i0 <- which(fit$data$year == year.base)
        plMap <- pl
        map <- fit$obj$env$map
        with.map <- intersect(names(plMap), names(map))
        applyMap <- function(par.name) {
            tapply(plMap[[par.name]], map[[par.name]], mean)
        }
        plMap[with.map] <- sapply(with.map, applyMap, simplify = FALSE)
        sniii <- 1
        doSim <- function(re_constraint = NULL, re_pl = NULL){
            obj2 <- obj
            if(!is.null(re_constraint)){
                ## Check length of constraints?
                cstr <- replicate(nYears, .forecastDefault(), simplify = FALSE)
                cstr[!is.na(re_constraint)] <- .parseForecast(re_constraint[!is.na(re_constraint)], fit$conf$fbarRange, fit$data$fleetTypes, c(fit$conf$minAge,fit$conf$maxAge), useNonLinearityCorrection)                
                obj2$env$data$forecast$constraints <- cstr
            }
            sim0 <- 0*est 
            if(resampleFirst){
                sim0 <- rmvnorm(1, mu=0*est, Sigma=cov)
            }
            dList0 <- split(as.vector(sim0), names(est))
            estList0 <- split(as.vector(sim0+est), names(est))
            if(!is.null(re_pl)){
                plMap2 <- re_pl
                map <- fit$obj$env$map
                with.map <- intersect(names(plMap2), names(map))
                applyMap <- function(par.name) {
                    tapply(plMap2[[par.name]], map[[par.name]], mean)
                }
                plMap2[with.map] <- sapply(with.map, applyMap, simplify = FALSE)
                p <- unlist(plMap2)
                names(p) <- rep(names(plMap2), times = sapply(plMap2,length))
            }else{
                p <- unlist(plMap)
                names(p) <- rep(names(plMap), times = sapply(plMap,length))
            }            
            ## Only works when year.base is last assessment year
            indxN <- matrix(which(names(p) %in% "logN"),nrow=length(estList0$LogN))[,i0]
            indxF <- matrix(which(names(p) %in% "logF"),nrow=length(estList0$LogF))[,i0]
            p[indxN] <- estList0$LogN
            p[indxF] <- estList0$LogF
            obj2$env$data$forecast$Fdeviation[] <- dList0$LogF
            obj2$env$data$forecast$FdeviationCov <- cov[names(est) %in% "LogF",names(est) %in% "LogF"]
            v <- obj2$simulate(par = p)
            sniii <<- sniii+1
            incpb()
            return(v)
        }
        if(as.integer(returnObj)==2)
            return(doSim)
        simvals <- .SAM_replicate(nosim, doSim(), simplify = FALSE, ncores = ncores, env = environment(doSim))
        simlist <- vector("list",length(FModel) + 1)
        for(i in 0:(length(FModel))){
            y<-year.base+i
            ii <- i0 + i
            simlist[[i+1]] <- list(sim = do.call("rbind",lapply(simvals,function(x) c(x$logN[,ii], x$logF[,ii]))),
                                   fbar = sapply(simvals,function(x) exp(x$logfbar[ii])),
                                   catch = sapply(simvals,function(x) exp(x$logCatch[ii])),
                                   ssb = sapply(simvals,function(x) exp(x$logssb[ii])),
                                   rec = sapply(simvals,function(x) exp(x$logN[1,ii])),
                                   cwF = rep(NA, nosim),
                                   catchatage = do.call("cbind",lapply(simvals,function(x) exp(x$logCatchAge[,ii]))),
                                   catchbyfleet = do.call("cbind",lapply(simvals,function(x) t(exp(x$logCatchByFleet[ii,fleetHasF,drop=FALSE])))),
                                   fbarbyfleet = do.call("cbind",lapply(simvals,function(x) exp(x$logFbarByFleet[fleetHasF,ii,drop=FALSE]))),
                                   land = sapply(simvals,function(x) exp(x$logLand[ii])),
                                   fbarL = sapply(simvals,function(x) exp(x$logfbarL[ii])),
                                   tsb = sapply(simvals,function(x) exp(x$logtsb[ii])),
                                   logEmpiricalSPR  = sapply(simvals,function(x) (x$logEmpiricalSPR[ii])),
                                   logEmpiricalYPR  = sapply(simvals,function(x) (x$logEmpiricalYPR[ii])),
                                   logEmpiricalYPR_L  = sapply(simvals,function(x) (x$logEmpiricalYPR_L[ii])),
                                   logEmpiricalYPR_D  = sapply(simvals,function(x) (x$logEmpiricalYPR_D[ii])),
                                   year=y)
            rownames(simlist[[i+1]]$catchatage) <- seq(fit$conf$minAge,fit$conf$maxAge,1)
        }
        attr(simlist, "fit")<-fit
        ## Similar to stockassessment::forecast
        collect <- function(x){
            est <- estimate(x)
            quan <- unname(quantile(x, c(.025,.975), na.rm = TRUE))
            v <- c(estimate=est, low=quan[1], high=quan[2])
            names(v)  <- c(estimateLabel, "low","high")
            v
        }
        fbar <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$fbar))),3)
        fbarL <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$fbarL))),3)  
        rec <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$rec))))
        ssb <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$ssb))))
        tsb <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$tsb))))
        catch <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$catch))))
        land <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$land))))  
        caytable<-round(do.call(rbind, lapply(simlist, function(xx)apply(xx$catchatage,1,collect))))
        cbftable<-round(do.call(rbind, lapply(simlist, function(xx)apply(xx$catchbyfleet,1,collect))))
        fbftable<-round(do.call(rbind, lapply(simlist, function(xx)apply(xx$fbarbyfleet,1,collect))),3)
        colnames(fbftable) <- attr(fit$data,"fleetNames")[fleetHasF]
        colnames(cbftable) <- attr(fit$data,"fleetNames")[fleetHasF]

        tab <- cbind(fbar,rec,ssb,catch)
        if(splitLD){
            tab<-cbind(tab,fbarL,fbar-fbarL,land,catch-land)
        }
        if(addTSB){
            tab<-cbind(tab,tsb)
        }
        ## if(!missing(customWeights)) tab <- cbind(tab,cwF=round(do.call(rbind, lapply(simlist, function(xx)collect(xx$cwF))),3))
        rownames(tab) <- unlist(lapply(simlist, function(xx)xx$year))
        nam <- c(estimateLabel,"low","high")
        basename<-c("fbar:","rec:","ssb:","catch:")
        if(splitLD){
            basename<-c(basename,"fbarL:","fbarD:","Land:","Discard:")    
        }
        if(addTSB){
            basename<-c(basename,"tsb:")    
        }
        ## if(!missing(customWeights)){
        ##     basename<-c(basename,"cwF:")    
        ## }
        colnames(tab)<-paste0(rep(basename, each=length(nam)), nam)
        
        attr(simlist, "tab")<-tab
        shorttab<-t(tab[,grep(estimateLabel,colnames(tab))])
        rownames(shorttab)<-sub(paste0(":",estimateLabel),"",paste0(label,if(!is.null(label))":",rownames(shorttab)))
        attr(simlist, "shorttab")<-shorttab
        attr(simlist, "label") <- label
        attr(simlist, "caytable")<-caytable        
        attr(simlist, "catchby")<-cbftable
        attr(simlist, "fbarby")<-fbftable
        attr(simlist, "nosim") <- nosim
        class(simlist) <- "samforecast"
        attr(simlist,"estimateLabel") <- estimateLabel
        attr(simlist,"useNonLinearityCorrection") <- useNonLinearityCorrection
        ## Done with reporting
        incpb()
        if(progress)
            close(pb)
        return(simlist)
    }else{
        invisible(obj$fn(fit$opt$par))
        ## Get results
        sdr <- TMB::sdreport(obj, fit$opt$par, svd_solve(fit$sdrep$cov.fixed),
                             bias.correct= biasCorrect,
                             skip.delta.method = biasCorrect,
                             bias.correct.control = list(sd = TRUE,
                                                         split = lapply(obj$env$ADreportIndex()[c("logfbar","logssb","logR","logCatch","logtsb","logLagR","logLand","logDis","loglandfbar","logdisfbar")[c(TRUE,TRUE,TRUE,TRUE,addTSB,lagR,splitLD,splitLD,splitLD,splitLD)]], utils::tail, n = nYears + 1 + (fit$data$noYears-1) * as.numeric(returnAllYears))
                                                         )
                             )
        ssdr <- summary(sdr)
        invisible(obj$fn(fit$opt$par))       
        rp <- obj$report(obj$env$last.par)

        simlist <- list()
        for(i in 0:(length(FModel))){
            y<-year.base+i 
            simlist[[i+1]] <- list(sim=NA, fbar=NA, catch=NA, ssb=NA, rec=NA,
                                   cwF=NA, catchatage=NA, catchbyfleet=NA, fbarbyfleet=NA, land=NA, fbarL=NA, tsb=NA, logEmpiricalSPR=NA, logEmpiricalYPR=NA, year=y)
        }

        attr(simlist, "fit")<-fit

        toCI <- function(x, trans = exp){
            trans(x %*% matrix(c(1,0,1,-2,1,2), nrow = 2, ncol =3))
        }

        indx <- 1:2
        if(biasCorrect)
            indx <- 3:4
        
        fbar <- toCI(ssdr[rownames(ssdr) %in% "logfbar",indx])
        rec <- toCI(ssdr[rownames(ssdr) %in% "logR",indx])
        if(lagR)
            rec <- toCI(ssdr[rownames(ssdr) %in% "logLagR",indx])
        ssb <- toCI(ssdr[rownames(ssdr) %in% "logssb",indx])
        catch <- toCI(ssdr[rownames(ssdr) %in% "logCatch",indx])
        cayvec <- toCI(ssdr[rownames(ssdr) %in% "logCatchAge",indx])
        cbfvec <- toCI(ssdr[rownames(ssdr) %in% "logCatchByFleet",indx])
        fbfvec <- toCI(ssdr[rownames(ssdr) %in% "logFbarByFleet",indx])
       
        tab <- cbind(fbar,rec,ssb,catch)
        if(splitLD){
            land <- toCI(ssdr[rownames(ssdr) %in% "logLand",indx])
            fbarL <- toCI(ssdr[rownames(ssdr) %in% "loglandfbar",indx])
            dis <- toCI(ssdr[rownames(ssdr) %in% "logDis",indx])
            fbarD <- toCI(ssdr[rownames(ssdr) %in% "logdisfbar",indx])
            tab<-cbind(tab,fbarL,fbarD,land,dis)
        }
        if(addTSB){
            tsb <- toCI(ssdr[rownames(ssdr) %in% "logtsb",indx])
            tab<-cbind(tab,tsb)
        }
       
        ## Handle cwF !!
        
        ## Row and column names on table
        futureYears <- unlist(lapply(simlist, function(xx)xx$year))
        rownames(tab) <- sort(unique(c(fit$data$years, futureYears)))
        fullTable <- tab
        tab <- tab[as.numeric(rownames(tab)) %in% futureYears, , drop = FALSE]
        
        nam <- c(estimateLabel,"low","high")
        basename<-c("fbar:","rec:","ssb:","catch:")
        if(splitLD){
            basename<-c(basename,"fbarL:","fbarD:","Land:","Discard:")    
        }
        if(addTSB){
            basename<-c(basename,"tsb:")    
        }
        ## if(!missing(customWeights)){
        ##     basename<-c(basename,"cwF:")    
        ## }

        colnames(tab)<-paste0(rep(basename, each=length(nam)), nam)
        if(returnAllYears || !biasCorrect)
            colnames(fullTable) <- colnames(tab)
        attr(simlist, "tab")<-tab
        shorttab<-t(tab[,grep(estimateLabel,colnames(tab))])
        rownames(shorttab)<-sub(paste0(":",estimateLabel),"",paste0(label,if(!is.null(label))":",rownames(shorttab)))

        ## Make caycbftable and fbftable
#### cay
        ii <- matrix(seq_along(rp$logCatchAge), nrow(rp$logCatchAge), ncol(rp$logCatchAge))[,as.numeric(rownames(fullTable)) %in% futureYears,drop=FALSE]
        caytable <- do.call("rbind",lapply(seq_len(ncol(ii)), function(yy){
            v <- t(cayvec[ii[,yy],,drop=FALSE])
        }))
        #### cbf
        ii <- matrix(seq_along(rp$logCatchByFleet), nrow(rp$logCatchByFleet), ncol(rp$logCatchByFleet))[as.numeric(rownames(fullTable)) %in% futureYears,fleetHasF,drop=FALSE]
        cbftable <- do.call("rbind",lapply(seq_len(nrow(ii)), function(yy){
            v <- t(cbfvec[ii[yy,],,drop=FALSE])
        }))
        #### fay
        ii <- matrix(seq_along(rp$logFbarByFleet), nrow(rp$logFbarByFleet), ncol(rp$logFbarByFleet))[fleetHasF,as.numeric(rownames(fullTable)) %in% futureYears,drop=FALSE]
        fbftable <- do.call("rbind",lapply(seq_len(ncol(ii)), function(yy){
            v <- t(fbfvec[ii[,yy],,drop=FALSE])
        }))
        #### dimnames
        rownames(caytable) <- rep(c(estimateLabel,"low","high"), length.out = nrow(caytable))
        rownames(cbftable) <- rownames(fbftable) <- rep(c(estimateLabel,"low","high"), length.out = nrow(cbftable))
        colnames(cbftable) <- colnames(fbftable) <- attr(fit$data,"fleetNames")[fleetHasF]
        colnames(caytable) <- seq(fit$conf$minAge,fit$conf$maxAge,1)

        ##

        
        attr(simlist, "shorttab")<-shorttab
        attr(simlist, "label") <- label
        attr(simlist, "caytable")<-caytable        
        attr(simlist, "catchby")<-cbftable
        attr(simlist, "fbarby")<-fbftable
        attr(simlist, "fulltab") <- fullTable
        attr(simlist,"nosim") <- 0
        class(simlist) <- "samforecast"
        attr(simlist,"estimateLabel") <- estimateLabel
        attr(simlist,"useNonLinearityCorrection") <- useNonLinearityCorrection
        return(simlist)    
    }
}
