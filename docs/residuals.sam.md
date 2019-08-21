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
```object```     |     sam fitted object as returned from the [`sam.fit`](sam.fit.html) function
```discrete```     |     logical if model contain discrete observations
```...```     |     extra arguments for TMB's oneStepPredict

## Details


 one-observation-ahead quantile residuals are calculated
 
 ...


