# `fbarplot`: SAM Fbar plot

## Description


 SAM Fbar plot


## Usage

```r
fbarplot(fit, partial = (class(fit) == "sam"), drop = NULL,
  pcol = "lightblue", page = NULL, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     the object returned from sam.fit
```partial```     |     true if included partial F's are to be plotted
```drop```     |     number of years to be left unplotted at the end. Default (NULL) is to not show years at the end with no catch information
```pcol```     |     color of partial lines
```page```     |     partial ages to plot
```...```     |     extra arguments transferred to plot including the following: list()  `add` logical, plotting is to be added on existing plot list()  `ci` logical, confidence intervals should be plotted list()  `cicol` color to plot the confidence polygon

## Details


 Plot the defined fbar.


