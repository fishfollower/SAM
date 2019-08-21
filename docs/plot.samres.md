# `plot.samres`: Plot sam residuals

## Description


 Plot sam residuals


## Usage

```r
list(list("plot"), list("samres"))(x, type = "bubble", ...)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |     an object of type 'samres' as returned from functions [`residuals.sam`](residuals.sam.html) or [`procres`](procres.html) .
```type```     |     either "bubble" (default) or "summary"
```...```     |     extra arguments

## Details


 In the "bubble" type red indicate negative residuals and blue positive. The area of the circles scales with the absolute size of the residuals.


## Examples

```r 
 list("\n", "data(nscodData)\n", "data(nscodConf)\n", "data(nscodParameters)\n", "fit <- sam.fit(nscodData, nscodConf, nscodParameters)\n", "par(ask=FALSE)\n", "plot(residuals(fit))\n") 
 ``` 

