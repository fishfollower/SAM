# `catchplot`

SAM catch plot


## Description

SAM catch plot


## Usage

```r
catchplot(fit, obs.show = TRUE, drop = NULL, ...)
list(list("catchplot"), list("sam"))(fit, obs.show = TRUE, drop = NULL, ...)
list(list("catchplot"), list("samset"))(fit, obs.show = TRUE, drop = NULL, ...)
list(list("catchplot"), list("samforecast"))(fit, obs.show = TRUE, drop = NULL, ...)
list(list("catchplot"), list("hcr"))(fit, obs.show = TRUE, drop = NULL, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     the object returned from sam.fit
`obs.show`     |     if observations are to be shown also
`drop`     |     number of years to be left unplotted at the end. Default (NULL) is to not show years at the end with no catch information
`...`     |     extra arguments transferred to plot including the following: list()  `add` logical, plotting is to be added on existing plot list()  `ci` logical, confidence intervals should be plotted list()  `cicol` color to plot the confidence polygon


## Details

Plot of estimated (and optionally observed) total catch in weight


