# `ypr`

Yield per recruit calculation


## Description

Yield per recruit calculation


## Usage

```r
ypr(
  fit,
  Flimit = 2,
  Fdelta = 0.01,
  aveYears = min(15, length(fit$data$years)),
  ageLimit = 100,
  sprProp = 0.35,
  ...
)
list(list("ypr"), list("sam"))(
  fit,
  Flimit = 2,
  Fdelta = 0.01,
  aveYears = min(15, length(fit$data$years)),
  ageLimit = 100,
  sprProp = 0.35,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     the object returned from sam.fit
`Flimit`     |     Upper limit for Fbar
`Fdelta`     |     increments on the Fbar axis
`aveYears`     |     Number of years back to use when calculating averages (selection, weights, ...)
`ageLimit`     |     Oldest age used (should be high)
`sprProp`     |     Proportion of SPR at F=0, for example 0.35 if F0.35SPR
`...`     |     extra arguments not currently used


