# `referencepoints`

Estimate reference points


## Description

Estimate reference points


## Usage

```r
referencepoints(
  fit,
  nYears,
  Fsequence,
  aveYears,
  selYears,
  SPRpercent,
  catchType,
  MSYreduction,
  newtonSteps = 3,
  optN = 100,
  jacobianHScale = 0.5,
  ...
)
list(list("referencepoints"), list("sam"))(
  fit,
  nYears = 100,
  Fsequence = seq(0, 4, len = 200),
  aveYears = max(fit$data$years) + (-9:0),
  selYears = max(fit$data$years),
  SPRpercent = c(0.35),
  dYPRpercent = c(0.1),
  B0percent = c(0.2),
  catchType = "catch",
  MSYreduction = c(0.05),
  newtonSteps = 3,
  optN = 20,
  jacobianHScale = 0.5,
  fixRP = c(),
  RecCorrection = "median",
  biasCorrect = FALSE,
  nlminb.control = list(eval.max = 1000, iter.max = 1000),
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     an object to calculate reference points for
`nYears`     |     Number of years to use in per-recruit calculations
`Fsequence`     |     Sequence of F values used to report per-recruit and equilibrium values
`aveYears`     |     Vector of year indices used to calculate average natural mortality, weights, etc. Following ICES guidelines, the default is the last 10 years (starting at 0)
`selYears`     |     Vector of year indices used to calculate selectivity (starting at 0)
`SPRpercent`     |     Vector of x values for F[x * 100%] reference points. Default is 0.35.
`catchType`     |     Catch type used: (total) catch, landings, discard.
`MSYreduction`     |     Vector of proportions for MSY ranges. Default is 0.05 giving an MSY range corresponding to no more than a 5% yield reduction.
`newtonSteps`     |     Number of additional Newton steps at the end of the reference point optimization.
`optN`     |     N used for numerical optimizers to find equilibrium biomass
`jacobianHScale`     |     Scale step size in jacobian calculation
`...`     |     not used
`dYPRpercent`     |     Defunct
`B0percent`     |     Defunct
`fixRP`     |     Defunct
`RecCorrection`     |     Defunct
`biasCorrect`     |     Defunct
`nlminb.control`     |     Defunct


## Value

a sam_referencepoints fit


## Seealso

[forecastMSY](#forecastmsy)


## References

Albertsen, C. M. and Trijoulet, V. (2020) Model-based estimates of reference points in an age-based state-space stock assessment model. Fisheries Research, 230, 105618. doi: 10.1016/j.fishres.2020.105618


