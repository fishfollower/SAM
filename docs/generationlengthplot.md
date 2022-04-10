# `generationlengthplot`

SAM generation length plot


## Description

SAM generation length plot


## Usage

```r
generationlengthplot(fit, ...)
list(list("generationlengthplot"), list("default"))(fit, ...)
list(list("generationlengthplot"), list("samforecast"))(fit, ...)
list(list("generationlengthplot"), list("hcr"))(fit, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     the object returned from sam.fit
`...`     |     extra arguments transferred to plot including the following: list()  `add` logical, plotting is to be added on existing plot list()  `ci` logical, confidence intervals should be plotted list()  `cicol` color to plot the confidence polygon


## Details

Plot of life expectancy


