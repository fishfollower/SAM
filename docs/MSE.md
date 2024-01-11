# `MSE`

Management strategy evaluation using SAM models


## Description

Management strategy evaluation using SAM models


## Usage

```r
MSE(
  OM,
  EM,
  nYears,
  AdviceForecastSettings,
  AdviceYears = 1,
  AdviceLag = 0,
  initialAdvice = NA,
  implementationError = function(x) x,
  knotRange = 3,
  intermediateFleets = numeric(0),
  OMselectivityFixed = FALSE,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`OM`     |     sam.fit that will work as operating model
`EM`     |     sam.fit that will work as estimation model
`nYears`     |     Number of years to run simulation
`AdviceYears`     |     Number of years advice given at a time. How advice is given is determined by forecastSettings
`AdviceLag`     |     Lag between assessment and advice
`initialAdvice`     |     Advice in the first AdviceLag years
`implementationError`     |     Function to add implementation error (i.e, transform advice to target catch)
`knotRange`     |     Range of spline knot values to try
`intermediateFleets`     |     Fleets that are available in the (first) intermediate year
`...`     |     arguments passed on to addSimulatedYears
`forecastSettings`     |     Settings to do forecast that determines advice


## Value

a list with MSE result


