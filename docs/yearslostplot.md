# `yearslostplot`

SAM years lost to fishing plot


## Description

SAM years lost to fishing plot


## Usage

```r
yearslostplot(fit, cause, ...)
list(list("yearslostplot"), list("default"))(fit, cause = c("Fishing", "Other", "LifeExpectancy"), ...)
list(list("yearslostplot"), list("samforecast"))(fit, cause = c("Fishing", "Other", "LifeExpectancy"), ...)
list(list("yearslostplot"), list("hcr"))(fit, cause = c("Fishing", "Other", "LifeExpectancy"), ...)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     the object returned from sam.fit
`cause`     |     Fisning, Other, or LifeExpectancy
`...`     |     extra arguments transferred to plot including the following: list()  `add` logical, plotting is to be added on existing plot list()  `ci` logical, confidence intervals should be plotted list()  `cicol` color to plot the confidence polygon


## Details

Plot of years lost to fishing


