# `runwithout`

runwithout helper function


## Description

runwithout helper function


## Usage

```r
runwithout(fit, year, fleet, ...)
list(list("runwithout"), list("sam"))(fit, year = NULL, fleet = NULL, map = fit$obj$env$map, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     a fitted model object as returned from sam.fit
`year`     |     a vector of years to be excluded.  When both fleet and year are supplied they need to be of same length, as only the pairs are excluded
`fleet`     |     a vector of fleets to be excluded.  When both fleet and year are supplied they need to be of same length, as only the pairs are excluded
`...`     |     extra arguments to sam.fit
`map`     |     map to use


## Details

...


