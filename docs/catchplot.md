

<!-- toc -->

October 16, 2017

# DESCRIPTION

```
Package: stockassessment
Title: State-Space Assessment Model
Version: 0.5.2
Date: 2017-10-11
Authors@R: c(person("Anders","Nielsen",role=c("aut","cre"), email="an@aqua.dtu.dk"),
	    person("Casper","Berg",role="aut"),
            person("Kasper","Kristensen",role="aut"), 
            person("Mollie", "Brooks", role="aut"),
	    person(c("Christoffer","Moesgaard"), "Albertsen", role="aut"))
Description: Fitting SAM...
License: GPL-2
Imports: TMB, ellipse, MASS
LinkingTo: TMB, RcppEigen
Suggests: testthat
URL: https://github.com/fishfollower/SAM
LazyData: TRUE
BugReports: https://github.com/fishfollower/SAM/issues
RoxygenNote: 6.0.1```


# `addforecast`: SAM add forecasts

## Description


 SAM add forecasts


## Usage

```r
addforecast(fit, what, dotcol = "black", dotpch = 19, dotcex = 1.5,
  intervalcol = gray(0.5, alpha = 0.5))
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     the object returned from sam.fit
```what```     |     what to plot
```dotcol```     |     color for dot
```dotpch```     |     pch for dot
```dotcex```     |     cex for dot
```intervalcol```     |     color for interval

## Details


 internal plotting fun


# `catchplot`: SAM catch plot

## Description


 SAM catch plot


## Usage

```r
catchplot(fit, obs.show = TRUE, drop = NULL, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     the object returned from sam.fit
```obs.show```     |     if observations are to be shown also
```drop```     |     number of years to be left unplotted at the end. Default (NULL) is to not show years at the end with no catch information
```...```     |     extra arguments transferred to plot including the following: list()  `add` logical, plotting is to be added on existing plot list()  `ci` logical, confidence intervals should be plotted list()  `cicol` color to plot the confidence polygon

## Details


 Plot of estimated (and optionally observed) total catch in weight


# `catchtable`: Catch table

## Description


 Catch table


## Usage

```r
catchtable(fit, obs.show = FALSE)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     ...
```obs.show```     |     logical add a column with catch sum of product rowsums(C*W)

## Details


 ...


# `clean.void.catches`: remove void catches

## Description


 remove void catches


## Usage

```r
clean.void.catches(dat, conf)
```


## Arguments

Argument      |Description
------------- |----------------
```dat```     |     data for the sam model as returned from the setup.sam.data function
```conf```     |     model configuration which can be set up using the [`defcon`](defcon.html) function and then modified

## Value


 an updated dataset without the catches where F is fixed to zero


# `coef.sam`: Extract fixed coefficients of sam object

## Description


 Extract fixed coefficients of sam object


## Usage

```r
list(list("coef"), list("sam"))(object, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```object```     |     sam fitted object (result from sam.fit)
```...```     |     extra arguments

## Details


 ...


# `corplotcommon`: Common function for plotting correlation matrices.

## Description


 Common function for plotting correlation matrices.


## Usage

```r
corplotcommon(x, fn, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |     a list of correlation matrices
```fn```     |     a vector of fleet names
```...```     |     extra arguments to plotcorr

# `corplot`: Plots between-age correlations by fleet, either estimated or empirical using residuals.

## Description


 Plots between-age correlations by fleet, either estimated or empirical using residuals.


## Usage

```r
corplot(x, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |     Either a sam fit as returned by sam.fit OR the object returned from residuals.sam
```...```     |     extra arguments to plot

# `c.sam`: Collect sam objects

## Description


 Collect sam objects


## Usage

```r
list(list("c"), list("sam"))(...)
```


## Arguments

Argument      |Description
------------- |----------------
```...```     |     sam fits to be combined

## Details


 ...


# `defcon`: Setup basic minimal configuration for sam assessment

## Description


 Setup basic minimal configuration for sam assessment


## Usage

```r
defcon(dat)
```


## Arguments

Argument      |Description
------------- |----------------
```dat```     |     sam data object

## Details


 ...


## Value


 a list containing
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  


# `defpar`: Setup minimal initial parameters

## Description


 Setup minimal initial parameters


## Usage

```r
defpar(dat, conf)
```


## Arguments

Argument      |Description
------------- |----------------
```dat```     |     sam data object
```conf```     |     configuration list

## Details


 ...


## Value


 a list containing the following
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  
 

*  


# `empirobscorrplot`: Plots the residual between-age correlation matrices by fleet.

## Description


 Plots the residual between-age correlation matrices by fleet.


## Usage

```r
empirobscorrplot(res, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```res```     |     the object returned from residuals.sam
```...```     |     extra arguments to plot

# `faytable`: F-at-age table

## Description


 F-at-age table


## Usage

```r
faytable(fit)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     ...

## Details


 ...


# `fbarplot`: SAM Fbar plot

## Description


 SAM Fbar plot


## Usage

```r
fbarplot(fit, partial = (class(fit) == "sam"), drop = NULL,
  pcol = "lightblue", page = NULL, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     the object returned from sam.fit
```partial```     |     true if included partial F's are to be plotted
```drop```     |     number of years to be left unplotted at the end. Default (NULL) is to not show years at the end with no catch information
```pcol```     |     color of partial lines
```page```     |     partial ages to plot
```...```     |     extra arguments transferred to plot including the following: list()  `add` logical, plotting is to be added on existing plot list()  `ci` logical, confidence intervals should be plotted list()  `cicol` color to plot the confidence polygon

## Details


 Plot the defined fbar.


# `fbartable`: Fbar table

## Description


 Fbar table


## Usage

```r
fbartable(fit)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     ...

## Details


 ...


# `fitplot`: Plots fit to data

## Description


 Plots fit to data


## Usage

```r
fitplot(fit, log = TRUE, fleets = unique(fit$data$aux[, "fleet"]), ...)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     the object returned from sam.fit
```log```     |     should the plot be against log-obs
```fleets```     |     an integer vector of fleets to plot. Default is all of them
```...```     |     extra arguments to plot

# `forecast`: forecast function to do shortterm

## Description


 forecast function to do shortterm


## Usage

```r
forecast(fit, fscale = NULL, catchval = NULL, fval = NULL, nosim = 1000,
  year.base = max(fit$data$years), ave.years = max(fit$data$years) + (-4:0),
  rec.years = max(fit$data$years) + (-9:0), label = NULL,
  overwriteSelYears = NULL, deterministic = FALSE)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     an assessment object of type sam, as returned from the function sam.fit
```fscale```     |     a vector of f-scales. See details.
```catchval```     |     a vector of target catches. See details.
```fval```     |     a vector of target f values. See details.
```nosim```     |     number of simulations default is 1000
```year.base```     |     starting year default last year in assessment. Currently it is only supported to use last assessment year or the year before
```ave.years```     |     vector of years to average for weights, maturity, M and such
```rec.years```     |     vector of years to use to resample recruitment from
```label```     |     optional label to appear in short table
```overwriteSelYears```     |     if a vector of years is specified, then the average selectivity of those years is used (not recommended)
```deterministic```     |     option to turn all process noise off (not recommended, as it will likely cause bias)

## Details


 There are three ways to specify a scenario. If e.g. four F values are specified (e.g. fval=c(.1,.2,.3,4)), then the first value is used in the last assessment year (base.year), and the three following in the three following years. Alternatively F's can be specified by a scale, or a target catch. Only one option can be used per year. So for instance to set a catch in the first year and an F-scale in the following one would write catchval=c(10000,NA,NA,NA), fscale=c(NA,1,1,1). The length of the vector specifies how many years forward the scenarios run.


## Value


 an object of type samforecast


# `getLowerBounds`: Bounds

## Description


 Bounds


## Usage

```r
getLowerBounds(parameters)
```


## Arguments

Argument      |Description
------------- |----------------
```parameters```     |     initial values for the model in a format similar to what is returned from the defpar function

## Value


 a named list


# `getUpperBounds`: Bounds

## Description


 Bounds


## Usage

```r
getUpperBounds(parameters)
```


## Arguments

Argument      |Description
------------- |----------------
```parameters```     |     initial values for the model in a format similar to what is returned from the defpar function

## Value


 a named list


# `is.whole.positive.number`: Function to test if x is ...

## Description


 Function to test if x is ...


## Usage

```r
is.whole.positive.number(x, tol = .Machine$double.eps^0.5)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |     number
```tol```     |     precision

## Details


 ...


# `jit`: Jitter runs

## Description


 Jitter runs


## Usage

```r
jit(fit, nojit = 10, par = defpar(fit$data, fit$conf), sd = 0.25,
  ncores = detectCores())
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     a fitted model object as returned from sam.fit
```nojit```     |     a list of vectors. Each element in the list specifies a run where the fleets mentioned are omitted
```par```     |     initial values to jitter around. The defaule ones are returned from the defpar function
```sd```     |     the standard deviation used to jitter the initial values (most parameters are on a log scale, so similar to cv)
```ncores```     |     the number of cores to attemp to use

## Details


 ...


## Value


 A "samset" object, which is basically a list of sam fits


# `leaveout`: leaveout run

## Description


 leaveout run


## Usage

```r
leaveout(fit, fleet = as.list(2:fit$data$noFleets), ncores = detectCores(),
  ...)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     a fitted model object as returned from sam.fit
```fleet```     |     a list of vectors. Each element in the list specifies a run where the fleets mentioned are omitted
```ncores```     |     the number of cores to attemp to use
```...```     |     extra arguments to [`sam.fit`](sam.fit.html)

## Details


 ...


# `loadConf`: Loads a model configuration from a file

## Description


 Loads a model configuration from a file


## Usage

```r
loadConf(dat, file)
```


## Arguments

Argument      |Description
------------- |----------------
```dat```     |     sam data list as returned from the function setup.sam.data
```file```     |     the file to read the configuration from

## Details


 function useful loading a model configuration. Such a configuration can be saved via the saveConf function


# `logLik.sam`: Log likelihood of sam object

## Description


 Log likelihood of sam object


## Usage

```r
list(list("logLik"), list("sam"))(object, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```object```     |     sam fitted object (result from sam.fit)
```...```     |     extra arguments

## Details


 ...


# `modelDescription`: Description of model

## Description


 Description of model


## Usage

```r
modelDescription(fit, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     returned object from sam.fit
```...```     |     Additional parameters to be passed to ...

## Details


 ...


# `modeltable`: model table

## Description


 model table


## Usage

```r
modeltable(fits)
```


## Arguments

Argument      |Description
------------- |----------------
```fits```     |     A sam fit as returned from the sam.fit function, or a collection c(fit1, fit2, ...) of such fits

## Details


 ...


# `modelVersionInfo`: Description of model

## Description


 Description of model


## Usage

```r
modelVersionInfo(fit, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     returned object from sam.fit
```...```     |     Additional parameters to be passed to ...

## Details


 Writes a string to install the version of the package which was used to run the model.


# `nobs.sam`: Extract number of observations from sam object

## Description


 Extract number of observations from sam object


## Usage

```r
list(list("nobs"), list("sam"))(object, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```object```     |     sam fitted object (result from sam.fit)
```...```     |     extra arguments

## Details


 ...


# `nscodConf`: 
 nscodConf


## Description


 nscodConf


## Format


 The format is:
 $ minAge
 $ maxAge
 $ maxAgePlusGroup
 $ keyLogFsta
 $ corFlag
 $ keyLogFpar
 $ keyQpow
 $ keyVarF
 $ keyVarLogN
 $ keyVarObs
 $ stockRecruitmentModelCode
 $ noScaledYears
 $ keyScaledYears
 $ keyParScaledYA
 $ fbarRange


## Usage

```r
data("nscodConf")```


## Details


 ...


## References


 ...


## Examples

```r 
 data(nscodConf)
 ## maybe str(nscodConf) ; plot(nscodConf) ...
 ``` 

# `nscodData`: 
 nscodData


## Description


 nscodData


## Format


 The format is:
 $ noFleets
 $ fleetTypes
 $ sampleTimes
 $ noYears
 $ years
 $ nobs
 $ idx1
 $ idx2
 $ aux
 $ logobs
 $ propMat
 $ stockMeanWeight
 $ catchMeanWeight
 $ natMor
 $ landFrac
 $ disMeanWeight
 $ landMeanWeight
 $ propF
 $ propM


## Usage

```r
data("nscodData")```


## Details


 ...


## References


 ...


## Examples

```r 
 data(nscodData)
 ## maybe str(nscodData) ; plot(nscodData) ...
 ``` 

# `nscodParameters`: 
 nscodParameters


## Description


 nscodParameters


## Format


 The format is:
 List of 14
 $ logFpar
 $ logQpow
 $ logSdLogFsta
 $ logSdLogN
 $ logSdLogObs
 $ rec_loga
 $ rec_logb
 $ itrans_rho
 $ logScale
 $ logScaleSSB
 $ logPowSSB
 $ logSdSSB
 $ logF
 $ logN


## Usage

```r
data("nscodParameters")```


## Details


 ...


## References


 ...


## Examples

```r 
 data(nscodParameters)
 ## maybe str(nscodParameters) ; plot(nscodParameters) ...
 ``` 

# `ntable`: N table

## Description


 N table


## Usage

```r
ntable(fit)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     ...

## Details


 ...


# `obscorrplot`: Plots the estimated correlation matrices by fleet.

## Description


 Plots the estimated correlation matrices by fleet.


## Usage

```r
obscorrplot(fit, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     the object returned from sam.fit
```...```     |     extra arguments to plot

# `obscov`: Extract observation covariance matrices from a SAM fit

## Description


 Extract observation covariance matrices from a SAM fit


## Usage

```r
obscov(fit, corr = FALSE)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     the object returned from sam.fit
```corr```     |     if TRUE return correlation matrices rather than covariances

## Value


 a list of matrices


# `parplot`: SAM parameter plot

## Description


 SAM parameter plot


## Usage

```r
parplot(fit, cor.report.limit = 0.95, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     the object returned from sam.fit
```cor.report.limit```     |     correlations with absolute value > this number is reported in the plot
```...```     |     extra arguments transferred to plot

## Details


 Plot of all estimated model parameters (fixed effects). Shown with confidence interval.


# `partable`: parameter table

## Description


 parameter table


## Usage

```r
partable(fit)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     ...

## Details


 ...


# `plotby`: Plot by one or two

## Description


 Plot by one or two


## Usage

```r
plotby(x = NULL, y = NULL, z = NULL, x.line = NULL, y.line = NULL,
  z.line = NULL, by = NULL, bubblescale = 1, x.common = !is.null(x),
  y.common = !is.null(y), z.common = !is.null(z), xlab = NULL,
  ylab = NULL, xlim = NULL, ylim = NULL, zmax = NULL, axes = TRUE,
  ...)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |     numeric vector of points to be plotted
```y```     |     numeric vector of points to be plotted
```z```     |     numeric vector of points to be plotted
```x.line```     |     numeric vector of points of line to be added
```y.line```     |     numeric vector of points of line to be added
```z.line```     |     numeric vector of points of line to be added
```by```     |     vector or two column matrix to create sub sets from
```bubblescale```     |     scaling of bubble size
```x.common```     |     logical: use same x-axis for all plots
```y.common```     |     logical: use same y-axis for all plots
```z.common```     |     logical: use same z-axis for all plots
```xlab```     |     normal graphical parameter
```ylab```     |     normal graphical parameter
```xlim```     |     normal graphical parameter
```ylim```     |     normal graphical parameter
```zmax```     |     internally used to scale bubbles similarly
```axes```     |     normal graphical parameter
```...```     |     additional graphical parameters

## Details


 Function used for splitting plots e.g. used to plot residuals


## Examples

```r 
 exdat<-expand.grid(age=1:5, year=1950:2016, fleet=1:3)
 exdat$perfectres<-rnorm(nrow(exdat))
 attach(exdat)
 par(ask=FALSE)
 plotby(year,age,perfectres, by=fleet)
 detach(exdat)
 ``` 

# `plot.samforecast`: Plot samforecast object

## Description


 Plot samforecast object


## Usage

```r
list(list("plot"), list("samforecast"))(x, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |     ...
```...```     |     extra arguments

## Details


 ...


# `plot.sam`: Plot sam object

## Description


 Plot sam object


## Usage

```r
list(list("plot"), list("sam"))(x, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |     ...
```...```     |     extra arguments (not possible to use add=TRUE --- please collect to a list of fits using e.g the c(...), and then plot that collected object)

## Details


 ...


# `plot.samres`: Plot sam residuals

## Description


 Plot sam residuals


## Usage

```r
list(list("plot"), list("samres"))(x, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |     ...
```...```     |     extra arguments

## Details


 ...


## Examples

```r 
 list("\n", "data(nscodData)\n", "data(nscodConf)\n", "data(nscodParameters)\n", "fit <- sam.fit(nscodData, nscodConf, nscodParameters)\n", "par(ask=FALSE)\n", "plot(residuals(fit))\n") 
 ``` 

# `plot.samset`: Plot sam object

## Description


 Plot sam object


## Usage

```r
list(list("plot"), list("samset"))(x, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |     ...
```...```     |     extra arguments

## Details


 ...


# `plot.samypr`: Plot sam object

## Description


 Plot sam object


## Usage

```r
list(list("plot"), list("samypr"))(x, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |     ...
```...```     |     extra arguments

## Details


 ...


# `print.samcoef`: Print samcoef object

## Description


 Print samcoef object


## Usage

```r
list(list("print"), list("samcoef"))(x, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |     ...
```...```     |     extra arguments

## Details


 ...


# `print.samforecast`: Print samforecast object

## Description


 Print samforecast object


## Usage

```r
list(list("print"), list("samforecast"))(x, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |     an object as returned from the forecast function
```...```     |     extra arguments

## Details


 ...


# `print.sam`: Print sam object

## Description


 Print sam object


## Usage

```r
list(list("print"), list("sam"))(x, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |     ...
```...```     |     extra arguments

## Details


 ...


# `print.samres`: Print samres object

## Description


 Print samres object


## Usage

```r
list(list("print"), list("samres"))(x, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |     ...
```...```     |     extra arguments

## Details


 ...


# `print.samset`: Print samset object

## Description


 Print samset object


## Usage

```r
list(list("print"), list("samset"))(x, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |     a list of sam models
```...```     |     extra arguments

## Details


 ...


# `print.samypr`: Print samypr object

## Description


 Print samypr object


## Usage

```r
list(list("print"), list("samypr"))(x, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |     an object as returned from the ypr function
```...```     |     extra arguments

## Details


 ...


# `procres`: Compute process residuals (single joint sample)

## Description


 Compute process residuals (single joint sample)


## Usage

```r
procres(fit, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     the fitted object as returned from the sam.fit function
```...```     |     extra arguments (not currently used)

## Details


 ...


## Value


 an object of class `samres` 


# `read.ices`: Function to read ICES/CEFAS data files and validate if input makes sense

## Description


 Function to read ICES/CEFAS data files and validate if input makes sense


## Usage

```r
read.ices(filen)
```


## Arguments

Argument      |Description
------------- |----------------
```filen```     |     The filename

## Details


 First two lines are ignored and can be used for comments.
 Can read formats 1 full, 2 row, 3 scalar, and 5 column
 
 Tests:
 Formatcode is valid, years and ages are pos. integers
 minimum <= maximum for years and ages
 number of rows and coulmns match year and age ranges
 data contains only numbers.
 
 Returns: A validated data matrix.


# `read.surveys`: Function to read ices survey format

## Description


 Function to read ices survey format


## Usage

```r
read.surveys(filen)
```


## Arguments

Argument      |Description
------------- |----------------
```filen```     |     the file

## Details


 ...


# `read.table.nowarn`: Function to supress incomplete final line warning

## Description


 Function to supress incomplete final line warning


## Usage

```r
read.table.nowarn(...)
```


## Arguments

Argument      |Description
------------- |----------------
```...```     |     arguments

## Details


 ...


# `recplot`: SAM Recruits plot

## Description


 SAM Recruits plot


## Usage

```r
recplot(fit, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     the object returned from sam.fit
```...```     |     extra arguments transferred to plot including the following: list()  `add` logical, plotting is to be added on existing plot list()  `ci` logical, confidence intervals should be plotted list()  `cicol` color to plot the confidence polygon

## Details


 Plot of numbers of recruits (youngest age class)


# `rectable`: Recruit table

## Description


 Recruit table


## Usage

```r
rectable(fit)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     ...

## Details


 ...


# `reduce`: reduce helper function to reduce data

## Description


 reduce helper function to reduce data


## Usage

```r
reduce(data, year = NULL, fleet = NULL, age = NULL, conf = NULL)
```


## Arguments

Argument      |Description
------------- |----------------
```data```     |     a data object as returned by the function setup.sam.data
```year```     |     a vector of years to be excluded.
```fleet```     |     a vector of fleets to be excluded.
```age```     |     a vector of ages fleets to be excluded.
```conf```     |     an optional corresponding configuration to be modified along with the data change. Modified is returned as attribute "conf"

## Details


 When more than one vector is supplied they need to be of same length, as only the pairs are excluded


# `residuals.sam`: Extract residuals from sam object

## Description


 Extract residuals from sam object


## Usage

```r
list(list("residuals"), list("sam"))(object, discrete = FALSE, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```object```     |     sam fitted object (result from sam.fit)
```discrete```     |     logical if model contain discrete observations
```...```     |     extra arguments for TMB's oneStepPredict

## Details


 ...


# `retro`: retro run

## Description


 retro run


## Usage

```r
retro(fit, year = NULL, ncores = detectCores(), ...)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     a fitted model object as returned from sam.fit
```year```     |     either 1) a single integer n in which case runs where all fleets are reduced by 1, 2, ..., n are returned, 2) a vector of years in which case runs where years from and later are excluded for all fleets, and 3 a matrix of years were each column is a fleet and each column corresponds to a run where the years and later are excluded.
```ncores```     |     the number of cores to attemp to use
```...```     |     extra arguments to [`sam.fit`](sam.fit.html)

## Details


 ...


# `runwithout`: runwithout helper function

## Description


 runwithout helper function


## Usage

```r
runwithout(fit, year = NULL, fleet = NULL, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     a fitted model object as returned from sam.fit
```year```     |     a vector of years to be excluded.  When both fleet and year are supplied they need to be of same length, as only the pairs are excluded
```fleet```     |     a vector of fleets to be excluded.  When both fleet and year are supplied they need to be of same length, as only the pairs are excluded
```...```     |     extra arguments to sam.fit

## Details


 ...


# `sam.fit`: Fit SAM model

## Description


 Fit SAM model


## Usage

```r
sam.fit(data, conf, parameters, newtonsteps = 3, rm.unidentified = FALSE,
  run = TRUE, lower = getLowerBounds(parameters),
  upper = getUpperBounds(parameters), sim.condRE = TRUE, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```data```     |     data for the sam model as returned from the setup.sam.data function
```conf```     |     model configuration which can be set up using the [`defcon`](defcon.html) function and then modified.
```parameters```     |     initial values which can be set up using the [`defpar`](defpar.html) function and then modified.
```newtonsteps```     |     optional extra true newton steps
```rm.unidentified```     |     option to eliminate unidentified model parameters based on gradient in initial value (somewhat experimental)
```run```     |     if FALSE return AD object without running the optimization
```lower```     |     named list with lower bounds for optimization (only met before extra newton steps)
```upper```     |     named list with upper bounds for optimization (only met before extra newton steps)
```sim.condRE```     |     logical with default `TRUE` . Simulated observations will be conditional on estimated values of F and N, rather than also simulating F and N forward from their initial values.
```...```     |     extra arguments to MakeADFun

## Details


 ...


## Value


 an object of class `sam` 


## Examples

```r 
 data(nscodData)
 data(nscodConf)
 data(nscodParameters)
 fit <- sam.fit(nscodData, nscodConf, nscodParameters)
 ``` 

# `saveConf`: Saves a model configuration list to a file

## Description


 Saves a model configuration list to a file


## Usage

```r
saveConf(x, file = "", overwrite = FALSE)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |     sam configuration list as returned from defcon or loadConf
```file```     |     the file to save the configuration to
```overwrite```     |     logical if an existing file should be overwritten (FALSE by default)

## Details


 function useful for saving a model configuration. A saved configuration can be read back in via the loadConf function


# `setSeq`: small helper function

## Description


 small helper function


## Usage

```r
setSeq(min, max)
```


## Arguments

Argument      |Description
------------- |----------------
```min```     |     from
```max```     |     to

## Details


 internal


# `setS`: small helper function

## Description


 small helper function


## Usage

```r
setS(x)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |     vector if indices

## Details


 internal


# `setup.sam.data`: Combine the data sources to SAM readable object

## Description


 Combine the data sources to SAM readable object


## Usage

```r
setup.sam.data(fleets = NULL, surveys = NULL, residual.fleet = NULL,
  prop.mature = NULL, stock.mean.weight = NULL, catch.mean.weight = NULL,
  dis.mean.weight = NULL, land.mean.weight = NULL,
  natural.mortality = NULL, prop.f = NULL, prop.m = NULL,
  land.frac = NULL, recapture = NULL)
```


## Arguments

Argument      |Description
------------- |----------------
```fleets```     |     comm fleets vith effort (currently unimplemented)
```surveys```     |     surveys
```residual.fleet```     |     total catch minus commercial
```prop.mature```     |     pm
```stock.mean.weight```     |     sw
```catch.mean.weight```     |     cw
```dis.mean.weight```     |     dw
```land.mean.weight```     |     lw
```natural.mortality```     |     nm
```prop.f```     |     ...
```prop.m```     |     ...
```land.frac```     |     ...
```recapture```     |     ...

## Details


 ...


# `simulate.sam`: Simulate from a sam object

## Description


 Simulate from a sam object


## Usage

```r
list(list("simulate"), list("sam"))(object, nsim = 1, seed = NULL, full.data = TRUE,
  ...)
```


## Arguments

Argument      |Description
------------- |----------------
```object```     |     sam fitted object (result from sam.fit)
```nsim```     |     number of response lists to simulate. Defaults to 1.
```seed```     |     random number seed
```full.data```     |     logical, should each inner list contain a full list of data. Defaults to TRUE
```...```     |     extra arguments

## Details


 ...


## Value


 returns a list of lists. The outer list has length `nsim` . Each inner list contains simulated values of `logF` , `logN` , and `obs` with dimensions equal to those parameters.


# `srplot`: Plots the stock recruitment

## Description


 Plots the stock recruitment


## Usage

```r
srplot(fit, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     the object returned from sam.fit
```...```     |     extra arguments to plot

# `ssbplot`: SAM SSB plot

## Description


 SAM SSB plot


## Usage

```r
ssbplot(fit, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     the object returned from sam.fit
```...```     |     extra arguments transferred to plot including the following: list()  `add` logical, plotting is to be added on existing plot list()  `ci` logical, confidence intervals should be plotted list()  `cicol` color to plot the confidence polygon

## Details


 Plot of spawning stock biomass


# `ssbtable`: SSB table

## Description


 SSB table


## Usage

```r
ssbtable(fit)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     ...

## Details


 ...


# `summary.sam`: Summary of sam object

## Description


 Summary of sam object


## Usage

```r
list(list("summary"), list("sam"))(object, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```object```     |     sam fitted object (result from sam.fit)
```...```     |     extra arguments

## Details


 ...


# `tsbplot`: SAM TSB plot

## Description


 SAM TSB plot


## Usage

```r
tsbplot(fit, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     the object returned from sam.fit
```...```     |     extra arguments transferred to plot including the following: list()  `add` logical, plotting is to be added on existing plot list()  `ci` logical, confidence intervals should be plotted list()  `cicol` color to plot the confidence polygon

## Details


 Plot of total stock biomass


# `tsbtable`: TSB table

## Description


 TSB table


## Usage

```r
tsbtable(fit)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     ...

## Details


 ...


# `ypr`: Yield per recruit calculation

## Description


 Yield per recruit calculation


## Usage

```r
ypr(fit, Flimit = 2, Fdelta = 0.01, aveYears = min(15,
  length(fit$data$years)), ageLimit = 100)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     the object returned from sam.fit
```Flimit```     |     Upper limit for Fbar
```Fdelta```     |     increments on the Fbar axis
```aveYears```     |     Number of years back to use when calculating averages (selection, weights, ...)
```ageLimit```     |     Oldest age used (should be high)

