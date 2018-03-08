# `fitplot`: Plots fit to data

## Description


 Plots fit to data


## Usage

```r
fitplot(fit, log = TRUE, ...)
list(list("fitplot"), list("sam"))(fit, log = TRUE, fleets = unique(fit$data$aux[,
  "fleet"]), ...)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     the object returned from sam.fit
```log```     |     should the plot be against log-obs
```...```     |     extra arguments to plot
```fleets```     |     an integer vector of fleets to plot. Default is all of them

