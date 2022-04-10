# `obscov`

Extract observation covariance matrices from a SAM fit


## Description

Extract observation covariance matrices from a SAM fit


## Usage

```r
obscov(fit, corr = FALSE, ...)
list(list("obscov"), list("sam"))(fit, corr = FALSE, ...)
list(list("obscov"), list("samset"))(fit, corr = FALSE, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     the object returned from sam.fit
`corr`     |     if TRUE return correlation matrices rather than covariances
`...`     |     extra arguments not currently used


## Value

a list of matrices


