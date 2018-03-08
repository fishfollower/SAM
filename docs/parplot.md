# `parplot`: SAM parameter plot

## Description


 SAM parameter plot


## Usage

```r
parplot(fit, cor.report.limit = 0.95, ...)
list(list("parplot"), list("sam"))(fit, cor.report.limit = 0.95, ...)
list(list("parplot"), list("samset"))(fit, cor.report.limit = 0.95, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     the object returned from sam.fit
```cor.report.limit```     |     correlations with absolute value > this number is reported in the plot
```...```     |     extra arguments transferred to plot

## Details


 Plot of all estimated model parameters (fixed effects). Shown with confidence interval.


