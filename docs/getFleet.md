# `getFleet`

Extract a fleet observed or predicted value from a fitted object


## Description

Extract a fleet observed or predicted value from a fitted object


## Usage

```r
getFleet(fit, fleet, pred = "FALSE")
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     A fitted object as returned from sam.fit
`fleet`     |     The fleet number
`pred`     |     Should it be predicted value, default is observed


## Details

Extract for example the observed or predicted catch at age of fleet "fleet"


## Value

A matrix of observed or predicted values for fleet "fleet"


