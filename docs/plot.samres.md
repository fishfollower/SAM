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

