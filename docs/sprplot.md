# `sprplot`

SAM SPR plot


## Description

SAM SPR plot


## Usage

```r
sprplot(fit, ...)
list(list("sprplot"), list("default"))(fit, ...)
list(list("sprplot"), list("samforecast"))(fit, ...)
list(list("sprplot"), list("hcr"))(fit, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     the object returned from sam.fit
`...`     |     extra arguments transferred to plot including the following: list()  `add` logical, plotting is to be added on existing plot list()  `ci` logical, confidence intervals should be plotted list()  `cicol` color to plot the confidence polygon


## Details

Plot of deterministic equilibrium spawners per recruit assuming biological parameters and selectivity for that year remains unchanged in the future.


