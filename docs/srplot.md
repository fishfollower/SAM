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
  xlim,
  ylim,
  add = FALSE,
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
`xlim`     |     bounds for x-axis
`ylim`     |     bounds for y-axis
`add`     |     false if a new plot should be created


