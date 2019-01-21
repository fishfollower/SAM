# `retro`: retro run

## Description


 retro run


## Usage

```r
retro(fit, year = NULL, ncores = detectCores(), ...)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     a fitted model object as returned from sam.fit
```year```     |     either 1) a single integer n in which case runs where all fleets are reduced by 1, 2, ..., n are returned, 2) a vector of years in which case runs where years from and later are excluded for all fleets, and 3 a matrix of years were each column is a fleet and each column corresponds to a run where the years and later are excluded.
```ncores```     |     the number of cores to attempt to use
```...```     |     extra arguments to [`sam.fit`](sam.fit.html)

## Details


 ...


