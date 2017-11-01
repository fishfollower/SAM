# `plotit`: Plot helper

## Description


 Plot helper


## Usage

```r
plotit(fit, what, x = fit$data$years, ylab = what, xlab = "Years",
  ex = numeric(0), trans = function(x) x, add = FALSE, ci = TRUE,
  cicol = gray(0.5, alpha = 0.5), addCI = if (class(fit) == "samset") {    
  rep(FALSE, length(fit)) } else {     NA }, drop = 0,
  unnamed.basename = "current", xlim = NULL, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     the fitted object from sam.fit of a set of such fits c(fit1,fit2)
```what```     |     quoted name of object to extract
```x```     |     x-values
```ylab```     |     label on y-axis
```xlab```     |     label on x-axis
```ex```     |     extra y's to make room for
```trans```     |     function to transform values by
```add```     |     logical, plotting is to be added on existing plot
```ci```     |     logical, confidence intervals should be plotted
```cicol```     |     color to plot the confidence polygon
```addCI```     |     A logical vector indicating if confidence intervals should be plotted for the added fits.
```drop```     |     number of years to be left unplotted at the end.
```unnamed.basename```     |     the name to assign an unnamed basefit
```xlim```     |     ...
```...```     |     extra arguments transferred to plot

## Details


 The basic plotting used bu many of the plotting functions (e.g. ssbplot, fbarplot ...)


