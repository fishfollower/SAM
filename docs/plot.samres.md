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
```x```     |     an object of type 'samres' as returned from residuals or procres.
```type```     |     either "bubble" (default) or "summary"
```...```     |     extra arguments

## Details


 ...


## Examples

```r 
 list("\n", "data(nscodData)\n", "data(nscodConf)\n", "data(nscodParameters)\n", "fit <- sam.fit(nscodData, nscodConf, nscodParameters)\n", "par(ask=FALSE)\n", "plot(residuals(fit))\n") 
 ``` 

