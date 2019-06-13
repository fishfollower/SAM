# `defpar`: Setup initial values for all model parameters and random effects.

## Description


 Setup initial values for all model parameters and random effects.


## Usage

```r
defpar(dat, conf)
```


## Arguments

Argument      |Description
------------- |----------------
```dat```     |     sam data object as returned from the function `setup.sam.data`
```conf```     |     sam configuration list, which could be read from a configuration file via the `loadConf` function. A default/dummy configuration can be generated via the `defcon` function.

## Details


 The model parameters and random effects are not initialized in any clever way - most are simply set to zero. If convergence problems occour different initial values can be tested, but it is more likely a problem with the model configuration.


## Value


 a list containing initial values for all model parameters and random effects in the model.


