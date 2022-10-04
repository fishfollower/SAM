# `recplot`

SAM Recruits plot


## Description

SAM Recruits plot


## Usage

```r
recplot(fit, lagR = FALSE, ...)
list(list("recplot"), list("sam"))(fit, lagR = FALSE, ...)
list(list("recplot"), list("samset"))(fit, lagR = FALSE, ...)
list(list("recplot"), list("samforecast"))(fit, lagR = FALSE, ...)
list(list("recplot"), list("hcr"))(fit, lagR = FALSE, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     the object returned from sam.fit
`lagR`     |     use the age after the youngest as R
`...`     |     extra arguments transferred to plot including the following: list()  `add` logical, plotting is to be added on existing plot list()  `ci` logical, confidence intervals should be plotted list()  `cicol` color to plot the confidence polygon


## Details

Plot of numbers of recruits (youngest age class)


