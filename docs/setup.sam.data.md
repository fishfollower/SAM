# `setup.sam.data`

Combine the data sources to SAM readable object


## Description

Combine the data sources to SAM readable object


## Usage

```r
setup.sam.data(
  fleets = NULL,
  surveys = NULL,
  residual.fleets = NULL,
  prop.mature = NULL,
  stock.mean.weight = NULL,
  catch.mean.weight = NULL,
  dis.mean.weight = NULL,
  land.mean.weight = NULL,
  natural.mortality = NULL,
  prop.f = NULL,
  prop.m = NULL,
  land.frac = NULL,
  recapture = NULL,
  sum.residual.fleets = NULL,
  aux.fleets = NULL,
  keep.all.ages = FALSE,
  average.sampleTimes.survey = TRUE,
  fleetnames.remove.space = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`fleets`     |     comm fleets vith effort (currently unimplemented)
`surveys`     |     surveys
`residual.fleets`     |     fleet, or list of fleets without effort information
`prop.mature`     |     pm
`stock.mean.weight`     |     sw
`catch.mean.weight`     |     cw
`dis.mean.weight`     |     dw
`land.mean.weight`     |     lw
`natural.mortality`     |     nm
`prop.f`     |     ...
`prop.m`     |     ...
`land.frac`     |     ...
`recapture`     |     ...
`sum.residual.fleets`     |     ...
`aux.fleets`     |     ...
`keep.all.ages`     |     ...
`average.sampleTimes.survey`     |     Should sample times for surveys be averaged?
`fleetnames.remove.space`     |     Should white space in fleet names be removed?


## Details

...


