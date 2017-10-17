# `plotby`: Plot by one or two

## Description


 Plot by one or two


## Usage

```r
plotby(x = NULL, y = NULL, z = NULL, x.line = NULL, y.line = NULL,
  z.line = NULL, by = NULL, bubblescale = 1, x.common = !is.null(x),
  y.common = !is.null(y), z.common = !is.null(z), xlab = NULL,
  ylab = NULL, xlim = NULL, ylim = NULL, zmax = NULL, axes = TRUE,
  ...)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |     numeric vector of points to be plotted
```y```     |     numeric vector of points to be plotted
```z```     |     numeric vector of points to be plotted
```x.line```     |     numeric vector of points of line to be added
```y.line```     |     numeric vector of points of line to be added
```z.line```     |     numeric vector of points of line to be added
```by```     |     vector or two column matrix to create sub sets from
```bubblescale```     |     scaling of bubble size
```x.common```     |     logical: use same x-axis for all plots
```y.common```     |     logical: use same y-axis for all plots
```z.common```     |     logical: use same z-axis for all plots
```xlab```     |     normal graphical parameter
```ylab```     |     normal graphical parameter
```xlim```     |     normal graphical parameter
```ylim```     |     normal graphical parameter
```zmax```     |     internally used to scale bubbles similarly
```axes```     |     normal graphical parameter
```...```     |     additional graphical parameters

## Details


 Function used for splitting plots e.g. used to plot residuals


## Examples

```r 
 exdat<-expand.grid(age=1:5, year=1950:2016, fleet=1:3)
 exdat$perfectres<-rnorm(nrow(exdat))
 attach(exdat)
 par(ask=FALSE)
 plotby(year,age,perfectres, by=fleet)
 detach(exdat)
 ``` 

