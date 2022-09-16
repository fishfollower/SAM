# `b0plot`

SAM equilibrium biomass in the absence of fishing plot


## Description

SAM equilibrium biomass in the absence of fishing plot


## Usage

```r
b0plot(fit, ...)
list(list("b0plot"), list("default"))(fit, ...)
list(list("b0plot"), list("samforecast"))(fit, ...)
list(list("b0plot"), list("hcr"))(fit, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     the object returned from sam.fit
`...`     |     extra arguments transferred to plot including the following: list()  `add` logical, plotting is to be added on existing plot list()  `ci` logical, confidence intervals should be plotted list()  `cicol` color to plot the confidence polygon


## Details

Plot of deterministic equilibrium biomass in the absence of fishing assuming biological parameters and selectivity for that year remains unchanged in the future.


