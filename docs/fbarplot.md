# `fbarplot`

SAM Fbar plot


## Description

SAM Fbar plot


## Usage

```r
fbarplot(fit, ...)
list(list("fbarplot"), list("sam"))(
  fit,
  partial = TRUE,
  drop = NULL,
  pcol = "lightblue",
  page = NULL,
  plot = TRUE,
  effectiveF = any(!fit$conf$seasonTimes %in% c(0, 1)),
  ...
)
list(list("fbarplot"), list("samset"))(
  fit,
  partial = FALSE,
  drop = NULL,
  pcol = "lightblue",
  page = NULL,
  ...
)
list(list("fbarplot"), list("samforecast"))(
  fit,
  partial = FALSE,
  drop = NULL,
  pcol = "lightblue",
  page = NULL,
  ...
)
list(list("fbarplot"), list("hcr"))(
  fit,
  partial = FALSE,
  drop = NULL,
  pcol = "lightblue",
  page = NULL,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     the object returned from sam.fit
`...`     |     extra arguments transferred to plot including the following: list()  `add` logical, plotting is to be added on existing plot list()  `ci` logical, confidence intervals should be plotted list()  `cicol` color to plot the confidence polygon
`partial`     |     true if included partial F's are to be plotted
`drop`     |     number of years to be left unplotted at the end. Default (NULL) is to not show years at the end with no catch information
`pcol`     |     color of partial lines
`page`     |     partial ages to plot
`plot`     |     true if fbar should be plotted
`effectiveF`     |     If TRUE, effective full year F based on catch and survival is plotted. If FALSE, full year F based on survival is plotted.


## Details

Plot the defined fbar.


