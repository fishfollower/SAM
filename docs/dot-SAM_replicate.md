# `.SAM_replicate`

Parallel replicate for modelforecast


## Description

Parallel replicate for modelforecast


## Usage

```r
.SAM_replicate(
  n,
  expr,
  simplify = "array",
  ncores = 1,
  env = parent.frame(n + 1),
  par_precall = NULL,
  type = ifelse(.Platform$OS.type == "unix", "mclapply", "PSOCK")
)
```


## Arguments

Argument      |Description
------------- |----------------
`n`     |     number of replicates
`expr`     |     expression
`simplify`     |     simplify passes to sapply
`ncores`     |     number of cores


## Value

output


