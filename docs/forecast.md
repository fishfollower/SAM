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


