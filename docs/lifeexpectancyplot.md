# `lifeexpectancyplot`

SAM life expectancy plot


## Description

SAM life expectancy plot


## Usage

```r
lifeexpectancyplot(fit, atRecruit = TRUE, ...)
list(list("lifeexpectancyplot"), list("default"))(fit, atRecruit = TRUE, ylimAdd = fit$conf$maxAge, ...)
list(list("lifeexpectancyplot"), list("samforecast"))(fit, atRecruit = TRUE, ylimAdd = fit$conf$maxAge, ...)
list(list("lifeexpectancyplot"), list("hcr"))(fit, atRecruit = TRUE, ylimAdd = fit$conf$maxAge, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     the object returned from sam.fit
`atRecruit`     |     If true, show life expectancy given survival until minAge, otherwise show life expectancy at birth
`...`     |     extra arguments transferred to plot including the following: list()  `add` logical, plotting is to be added on existing plot list()  `ci` logical, confidence intervals should be plotted list()  `cicol` color to plot the confidence polygon
`ylimAdd`     |     values to add when calculating ylim for the plot


## Details

Plot of life expectancy


