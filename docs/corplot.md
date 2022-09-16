# `corplot`

Plots between-age correlations by fleet, either estimated or empirical using residuals.


## Description

Plots between-age correlations by fleet, either estimated or empirical using residuals.


## Usage

```r
corplot(x, ...)
list(list("corplot"), list("sam"))(x, ...)
list(list("corplot"), list("samres"))(x, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`x`     |     Either a sam fit as returned by sam.fit OR the object returned from residuals.sam
`...`     |     extra arguments to plot


