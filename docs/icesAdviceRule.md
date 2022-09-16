# `icesAdviceRule`

Forecast with an ICES advice rule


## Description

Forecast with an ICES advice rule


## Usage

```r
icesAdviceRule(
  x,
  Fmsy,
  MSYBtrigger,
  Blim,
  nosim = 10000,
  ave.years = max(x$data$years) + (-4:0),
  rec.years = numeric(0),
  preForecast = list(),
  currentSSB = FALSE,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`x`     |     Fitted assessment model
`Fmsy`     |     ICES Fmsy which is used as target F
`MSYBtrigger`     |     ICES MSYBtrigger below which F is reduced
`Blim`     |     ICES Blim below which F is set to zero.
`nosim`     |     Number of simulations to do. If NULL a model forecast based on the Laplace approximation is used
`ave.years`     |     vector of years to average for weights, maturity, M and such
`rec.years`     |     vector of years to use to resample recruitment from. If an empty vector is given, recruitment is based on the fitted model.
`preForecast`     |     list of forecast parameters (i.e., fval, fscale, catchval, landval, or nextssb) to use before the HCR
`currentSSB`     |     if TRUE, SSB at the begining of the control rule year is used. If FALSE, SSB at the begining of the previous year is used.
`...`     |     Other arguments passes to hcr


## Value

hcr object


## Seealso

[hcr](#hcr)


## References

ICES (2021) Advice on fishing opportunities. DOI: 10.17895/ices.advice.7720


