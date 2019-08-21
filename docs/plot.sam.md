# `plot.sam`: Plot sam object

## Description


 Plot sam object


## Usage

```r
list(list("plot"), list("sam"))(x, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```x```     |     fitted object as returned from the [`sam.fit`](sam.fit.html) function.
```...```     |     extra arguments (not possible to use add=TRUE --- instead collect to a list of fits using e.g the c(...), and then plot that collected object).

## Details


 gives a 3 plot overview plot og ssb, fbar, and recruits. These plots are available individually via the functions [`ssbplot`](ssbplot.html) , `fbarplot` , and [`recplot`](recplot.html) .


