# `componentplot`

Area plot of spawning components


## Description

Area plot of spawning components


## Usage

```r
componentplot(fit, ...)
list(list("componentplot"), list("sam"))(
  fit,
  onlyComponentYears = FALSE,
  ylab = "Composition",
  colSet = c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100",
    "#CC6677", "#882255", "#AA4499"),
  legend.pos = "bottom",
  bg = "white",
  ncol = length(cf),
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     sam fit
`...`     |     passed to legend
`onlyComponentYears`     |     If true, x axis is limited to the range with spawning component data. Otherwise, the model years are used.
`ylab`     |     Label for y axis
`colSet`     |     Colors
`legend.pos`     |     Legend position. See ?legend
`bg`     |     Background of legend. See ?legend
`ncol`     |     Number of columns in legend. See ?legend


## Value

Nothing


## Author

Christoffer Moesgaard Albertsen


