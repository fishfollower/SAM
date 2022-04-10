# `ICESvalues`

Function to estimate ICES values


## Description

The function calculates MSYBtrigger and Bpa from estimated reference points through a simulation forecast. If a sam object is used, other reference points (e.g., Fmsy) are estimated first.


## Usage

```r
ICESvalues(x, nosim, nyears, ntail, ...)
list(list("ICESvalues"), list("sam"))(x, nosim, nyears, ...)
list(list("ICESvalues"), list("sam_referencepoints"))(
  x,
  nosim = 1000,
  nyears = 100,
  ntail = 10,
  quantile_buffer_Fmsy = 0.5,
  calculate_Fp05 = TRUE,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`x`     |     object to calculate for
`nosim`     |     Number of simulations
`nyears`     |     Number of years to forecast before equilibrium
`ntail`     |     Number of years to use for calculations
`...`     |     Other parameters passed to [referencepoints](#referencepoints)


## Value

ICES values


## Author

Christoffer Moesgaard Albertsen


