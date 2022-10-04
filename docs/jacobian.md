# `jacobian`

Calculate jacobian of a function


## Description

Calculate jacobian of a function


## Usage

```r
jacobian(
  func,
  x,
  h = abs(1e-04 * x) + 1e-04 * (abs(x) < sqrt(.Machine$double.eps/7e-07)),
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`func`     |     function
`x`     |     parameter values
`h`     |     step size
`...`     |     passed to func


## Value

jacobian matrix


