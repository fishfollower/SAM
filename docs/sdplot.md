# `sdplot`

Plots the sd of the log observations as estimated in SAM in increasing order


## Description

Plots the sd of the log observations as estimated in SAM in increasing order


## Usage

```r
sdplot(fit, barcol = NULL, marg = NULL, ylim = NULL, ...)
list(list("sdplot"), list("sam"))(fit, barcol = NULL, marg = NULL, ylim = NULL, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     the object returned from sam.fit
`barcol`     |     color for each fleet and age
`marg`     |     margin for plot (mar in par())
`ylim`     |     bounds for y-axis
`...`     |     extra arguments to plot


