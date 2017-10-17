# `jit`: Jitter runs

## Description


 Jitter runs


## Usage

```r
jit(fit, nojit = 10, par = defpar(fit$data, fit$conf), sd = 0.25,
  ncores = detectCores())
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     a fitted model object as returned from sam.fit
```nojit```     |     a list of vectors. Each element in the list specifies a run where the fleets mentioned are omitted
```par```     |     initial values to jitter around. The defaule ones are returned from the defpar function
```sd```     |     the standard deviation used to jitter the initial values (most parameters are on a log scale, so similar to cv)
```ncores```     |     the number of cores to attemp to use

## Details


 ...


## Value


 A "samset" object, which is basically a list of sam fits


