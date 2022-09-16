# `rmaxplot`

SAM rmax plot


## Description

SAM rmax plot


## Usage

```r
rmaxplot(fit, ...)
list(list("rmaxplot"), list("default"))(fit, ...)
list(list("rmaxplot"), list("samforecast"))(fit, ...)
list(list("rmaxplot"), list("hcr"))(fit, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     the object returned from sam.fit
`...`     |     extra arguments transferred to plot including the following: list()  `add` logical, plotting is to be added on existing plot list()  `ci` logical, confidence intervals should be plotted list()  `cicol` color to plot the confidence polygon


## Details

Plot of life expectancy


