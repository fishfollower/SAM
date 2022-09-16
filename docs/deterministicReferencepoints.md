# `deterministicReferencepoints`

Function to calculate reference points for the embedded deterministic model of a SAM fit


## Description

The function estimates reference points based on deterministic per-recruit calculations with no process variance.
 The following reference points are implemented:
 list("\n", "   ", list(list("F=x"), list("F fixed to x, e.g., ", list("\"F=0.3\""))), "\n", "   ", list(list("StatusQuo"), list("F in the last year of the assessment")), "\n", "   ", list(list("StatusQuo-y"), list("F in the y years before the last in the assessment, e.g., ", list("\"StatusQuo-1\""))), "\n", "   ", list(list("MSY"), list("F that maximizes yield")), "\n", "   ", list(list("0.xMSY"), list("Fs that gives 0.x*100% of MSY, e.g., ", list("\"0.95MSY\""))), "\n", "   ", list(list("Max"), 
    list("F that maximizes yield per recruit")), "\n", "   ", list(list("0.xdYPR"), list("F such that the derivative of yield per recruit is 0.x times the derivative at F=0, e.g., ", list("\"0.1dYPR\""))), "\n", "   ", list(list("0.xSPR"), list("F such that spawners per recruit is 0.x times spawners per recruit at F=0, e.g., ", list("\"0.35SPR\""))), "\n", "   ", list(list("0.xB0"), list("F such that biomass is 0.x times the biomass at F=0, e.g., ", list("\"0.2B0\""))), "\n")


## Usage

```r
deterministicReferencepoints(fit, referencepoints, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     A fitted SAM model
`referencepoints`     |     list of reference points to calculate (See details)
`...`     |     other arguments not used


## Value

List of estimated reference points


## Examples

```r
deterministicReferencepoints(fit, c("MSY","0.95MSY","Max","0.35SPR","0.1dYPR","StatusQuo-3"))
```


