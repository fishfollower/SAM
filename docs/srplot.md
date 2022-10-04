# `srplot`

Plots the stock recruitment


## Description

Plots the stock recruitment


## Usage

```r
srplot(fit, ...)
list(list("srplot"), list("sam"))(
  fit,
  textcol = "red",
  years = TRUE,
  linetype = "l",
  linecol = "black",
  polycol = do.call("rgb", c(as.list(col2rgb("black")[, 1]), list(alpha = 0.1))),
  polyborder = do.call("rgb", c(as.list(col2rgb("black")[, 1]), list(alpha = 0.3))),
  polylty = 3,
  polylwd = 1,
  xlim,
  ylim,
  add = FALSE,
  CIlevel = 0.95,
  addCurve = TRUE,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     the object returned from sam.fit
`...`     |     extra arguments to plot
`textcol`     |     color of years on plot
`years`     |     the plotting symbols are the years
`linetype`     |     type for the plot (default line)
`linecol`     |     color of lines between points
`polycol`     |     Inner color of error ellipses
`polyborder`     |     Border color of error ellipses
`polylty`     |     Border line type of error ellipses
`polylwd`     |     Border line width of error ellipses
`xlim`     |     bounds for x-axis
`ylim`     |     bounds for y-axis
`add`     |     false if a new plot should be created
`CIlevel`     |     Confidence level for error ellipses on stock-recruitment pairs
`addCurve`     |     Call addRecruitmentCurve?


