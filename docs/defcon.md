# `defcon`

Setup basic minimal configuration for sam assessment


## Description

Setup basic minimal configuration for sam assessment


## Usage

```r
defcon(dat, level = 1)
```


## Arguments

Argument      |Description
------------- |----------------
`dat`     |     sam data object
`level`     |     1 or 2 (1 most basic configuration, 2 configuration with AR correlation structure on surveys)


## Details

The configuration returned by defcon is intended as a help to set up a syntactically correct configuration for the sam model. The dimensions are set from the data (years, age-classes, and fleet types available). The configuration is intended to be fairly simplistic in the hope that the model configured will at least converge (not guaranteed). Most importantly: No model validation has been performed, so it should not be assumed that the returned model configuration will result in a sensible assessment of the stock. The actual model configuration is the responsibility of the user.


## Value

a list containing the elements needed to configure a sam model (e.g. minAge, maxAge, maxAgePlusGroup, keyLogFsta, ...).


