# `addRecruitmentCurve`

Add stock-recruitment curve to srplot


## Description

Add stock-recruitment curve to srplot


## Usage

```r
addRecruitmentCurve(
  fit,
  CI = TRUE,
  col = rgb(0.6, 0, 0),
  cicol = rgb(0.6, 0, 0, 0.3),
  plot = TRUE,
  PI = FALSE,
  picol = rgb(0.6, 0, 0),
  pilty = 2,
  ...
)
list(list("addRecruitmentCurve"), list("sam"))(
  fit,
  CI = TRUE,
  col = rgb(0.6, 0, 0),
  cicol = rgb(0.6, 0, 0, 0.3),
  plot = TRUE,
  PI = FALSE,
  picol = rgb(0.6, 0, 0),
  pilty = 2,
  year = NA_real_,
  lastR = NA_real_,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     Object to show SR-curve for
`CI`     |     Add confidence intervals?
`col`     |     Color of fitted line
`cicol`     |     Color of confidence intervals
`plot`     |     Add the curve to a plot?
`PI`     |     Add prediction intervals?
`picol`     |     Color of prediction interval line
`pilty`     |     Line type of prediction interval line
`...`     |     not used


## Seealso

srplot


## Author

Christoffer Moesgaard Albertsen


