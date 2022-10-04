# `yprplot`

SAM YPR plot


## Description

SAM YPR plot


## Usage

```r
yprplot(fit, ...)
list(list("yprplot"), list("default"))(fit, ...)
list(list("yprplot"), list("samforecast"))(fit, ...)
list(list("yprplot"), list("hcr"))(fit, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     the object returned from sam.fit
`...`     |     extra arguments transferred to plot including the following: list()  `add` logical, plotting is to be added on existing plot list()  `ci` logical, confidence intervals should be plotted list()  `cicol` color to plot the confidence polygon


## Details

Plot of deterministic equilibrium yield per recruit assuming biological parameters and selectivity for that year remains unchanged in the future.


