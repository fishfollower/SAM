# `simulate.sam`: Simulate from a sam object

## Description


 Simulate from a sam object


## Usage

```r
list(list("simulate"), list("sam"))(object, nsim = 1, seed = NULL,
  full.data = TRUE, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```object```     |     sam fitted object (result from sam.fit)
```nsim```     |     number of response lists to simulate. Defaults to 1.
```seed```     |     random number seed
```full.data```     |     logical, should each inner list contain a full list of data. Defaults to TRUE
```...```     |     extra arguments

## Details


 ...


## Value


 returns a list of lists. The outer list has length `nsim` . Each inner list contains simulated values of `logF` , `logN` , and `obs` with dimensions equal to those parameters.


