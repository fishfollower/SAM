# `leaveout`

leaveout run


## Description

leaveout run


## Usage

```r
leaveout(
  fit,
  fleet = as.list(2:fit$data$noFleets),
  ncores = detectCores(),
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     a fitted model object as returned from sam.fit
`fleet`     |     a list of vectors. Each element in the list specifies a run where the fleets mentioned are omitted
`ncores`     |     the number of cores to attemp to use
`...`     |     extra arguments to [`sam.fit`](#sam.fit)


## Details

...


