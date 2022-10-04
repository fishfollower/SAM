# `bc`

Spline basis for use with formula interface


## Description

Spline basis for use with formula interface


## Usage

```r
bc(x, df = 3L, knots = NULL, Boundary.knots = range(x), intercept = FALSE)
```


## Arguments

Argument      |Description
------------- |----------------
`x`     |     Points to evaluate the basis in
`df`     |     Degrees of freedom
`knots`     |     Internal knots. If NULL, they are selected from quantiles of x.
`Boundary.knots`     |     Boundary knots. Defaults to range of x
`intercept`     |     Include an intercept in basis?


## Value

A spline basis


