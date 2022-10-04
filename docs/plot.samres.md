# `plot.samres`

Plot sam residuals


## Description

Plot sam residuals


## Usage

```r
list(list("plot"), list("samres"))(x, type = "bubble", ...)
```


## Arguments

Argument      |Description
------------- |----------------
`x`     |     an object of type 'samres' as returned from functions [`residuals.sam`](#residuals.sam) or [`procres`](#procres) .
`type`     |     either "bubble" (default) or "summary"
`...`     |     extra arguments


## Details

In the "bubble" type red indicate negative residuals and blue positive. The area of the circles scales with the absolute size of the residuals.


## Examples

```r
data(nscodData)
data(nscodConf)
data(nscodParameters)
fit <- sam.fit(nscodData, nscodConf, nscodParameters)
par(ask=FALSE)
plot(residuals(fit))
```


