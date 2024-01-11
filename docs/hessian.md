# `hessian`

Calculate hessian of a function


## Description

Calculate hessian of a function


## Usage

```r
hessian(
  func,
  x,
  h = abs(1e-04 * x) + 1e-04 * (abs(x) < sqrt(.Machine$double.eps/7e-07)),
  columns = seq_along(x),
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


## Note

Could be made more accurate in some cases


