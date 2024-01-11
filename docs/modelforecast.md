# `modelforecast`

Model based forecast function


## Description

Model based forecast function
 
 Model based forecast function


## Usage

```r
modelforecast(fit, ...)
list(list("modelforecast"), list("sam"))(
  fit,
  constraints = NULL,
  fscale = NULL,
  catchval = NULL,
  fval = NULL,
  nextssb = NULL,
  landval = NULL,
  nosim = 0,
  year.base = max(fit$data$years),
  ave.years = max(fit$data$years) + (-9:0),
  overwriteBioModel = FALSE,
  rec.years = c(),
  label = NULL,
  overwriteSelYears = NULL,
  deterministicF = FALSE,
  processNoiseF = FALSE,
  fixedFdeviation = FALSE,
  useFHessian = FALSE,
  resampleFirst = !is.null(nosim) && nosim > 0,
  useModelLastN = TRUE,
  customSel = NULL,
  lagR = FALSE,
  splitLD = FALSE,
  addTSB = FALSE,
  biasCorrect = FALSE,
  returnAllYears = FALSE,
  returnObj = FALSE,
  progress = TRUE,
  estimate = median,
  silent = TRUE,
  newton_config = NULL,
  custom_pl = NULL,
  useNonLinearityCorrection = (nosim > 0 && !deterministicF),
  ncores = 1,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     SAM model fit
`...`     |     other variables used by the methods
`constraints`     |     a character vector of forecast constraint specifications
`fscale`     |     a vector of f-scales. See details.
`catchval`     |     a vector of target catches. See details "old specification".
`fval`     |     a vector of target f values. See details "old specification".
`nextssb`     |     a vector target SSB values the following year. See details "old specification".
`landval`     |     a vector of target catches. See details "old specification".
`nosim`     |     number of simulations. If 0, the Laplace approximation is used for forecasting.
`year.base`     |     starting year default last year in assessment. Currently it is only supported to use last assessment year or the year before
`ave.years`     |     vector of years to average for weights, maturity, M and such
`overwriteBioModel`     |     Overwrite GMRF models with ave.years?
`rec.years`     |     vector of years to use to resample recruitment from. If the vector is empty, the stock recruitment model is used.
`label`     |     optional label to appear in short table
`overwriteSelYears`     |     if a vector of years is specified, then the average selectivity of those years is used (not recommended)
`deterministicF`     |     option to set F variance to (almost) zero (not recommended)
`processNoiseF`     |     option to turn off process noise in F
`fixedFdeviation`     |     Use a fixed F deviation from target?
`useFHessian`     |     Use the covariance of F estimates instead of the estimated process covariance for forecasting?
`resampleFirst`     |     Resample base year when nosim > 0?
`customSel`     |     supply a custom selection vector that will then be used as fixed selection in all years after the final assessment year (not recommended)
`lagR`     |     if the second youngest age should be reported as recruits
`splitLD`     |     if TRUE the result is split in landing and discards
`addTSB`     |     if TRUE the total stock biomass (TSB) is added
`biasCorrect`     |     Do bias correction of reported variables. Can be turned off to reduce running time (not recommended).
`returnAllYears`     |     If TRUE, all years are bias corrected. Otherwise, only forecast years are corrected.
`returnObj`     |     Only return TMB object?
`progress`     |     Show progress bar for simulations?
`estimate`     |     the summary function used (typically mean or median) for simulations
`silent`     |     Passed to MakeADFun. Should the TMB object be silent?
`newton_config`     |     Configuration for newton optimizer to find F values. See ?TMB::newton for details. Use NULL for TMB defaults.
`custom_pl`     |     Parameter list. By default, the parameter list from fit is used.
`useNonLinearityCorrection`     |     Should a non linearity correction be added to transformation of logF? See Details - Non-linearity correction.


## Details

Function to forecast the model under specified catch constraints. In the forecast, catch constraints are used to set the mean of the $log(F)$ process for each simulation. Therefore, catch constraints are not matched exactly in individual simulations. Likewise, the summary of a specific set of simulations will not match exactly due to random variability.
 By default, recruitment is forecasted using the estimated recruitment model. If a vector of recruitment years is given, recruitment is forecasted using a log-normal distribution with the same mean and variance as the recruitment in the years given. This is different from the forecast function, which samples from the recruitment estimates.
 Catch scenarios are specified by a vector of target constraints. The first value determines F in the year after the base year.


## Value

an object of type samforecast


## Seealso

forecast


