# `sam.fit`: Fit SAM model

## Description


 Fit SAM model


## Usage

```r
sam.fit(data, conf, parameters, newtonsteps = 3, rm.unidentified = FALSE,
  run = TRUE, lower = getLowerBounds(parameters),
  upper = getUpperBounds(parameters), sim.condRE = TRUE, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```data```     |     data for the sam model as returned from the setup.sam.data function
```conf```     |     model configuration which can be set up using the [`defcon`](defcon.html) function and then modified.
```parameters```     |     initial values which can be set up using the [`defpar`](defpar.html) function and then modified.
```newtonsteps```     |     optional extra true newton steps
```rm.unidentified```     |     option to eliminate unidentified model parameters based on gradient in initial value (somewhat experimental)
```run```     |     if FALSE return AD object without running the optimization
```lower```     |     named list with lower bounds for optimization (only met before extra newton steps)
```upper```     |     named list with upper bounds for optimization (only met before extra newton steps)
```sim.condRE```     |     logical with default `TRUE` . Simulated observations will be conditional on estimated values of F and N, rather than also simulating F and N forward from their initial values.
```...```     |     extra arguments to MakeADFun

## Details


 ...


## Value


 an object of class `sam` 


## Examples

```r 
 data(nscodData)
 data(nscodConf)
 data(nscodParameters)
 fit <- sam.fit(nscodData, nscodConf, nscodParameters)
 ``` 

