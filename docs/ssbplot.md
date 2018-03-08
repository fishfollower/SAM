# `ssbplot`: SAM SSB plot

## Description


 SAM SSB plot


## Usage

```r
ssbplot(fit, ...)
list(list("ssbplot"), list("default"))(fit, ...)
list(list("ssbplot"), list("samforecast"))(fit, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     the object returned from sam.fit
```...```     |     extra arguments transferred to plot including the following: list()  `add` logical, plotting is to be added on existing plot list()  `ci` logical, confidence intervals should be plotted list()  `cicol` color to plot the confidence polygon

## Details


 Plot of spawning stock biomass


