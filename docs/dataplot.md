# `dataplot`: SAM Data plot

## Description


 SAM Data plot


## Usage

```r
dataplot(fit, col = NULL, fleet_type = NULL, fleet_names = NULL)
list(list("dataplot"), list("sam"))(fit, col = NULL, fleet_type = NULL,
  fleet_names = NULL)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     the object returned from sam.fit
```col```     |     color to use for each fleet, default is two sequential colors list()
```fleet_type```     |     character vector giving the type of data per fleet. The default uses fit$data$fleetTypes as follows: list()  `fit$data$fleetTypes==0` "Catch at age" list()  `fit$data$fleetTypes==1` "Catch at age with effort" list()  `fit$data$fleetTypes==2 or 6` "Index at age" list()  `fit$data$fleetTypes==3` "Biomass or catch index" list()  `fit$data$fleetTypes==5` "Tagging data" list()  `fit$data$fleetTypes==7` "Sum of fleets"
```fleet_names```     |     character vector giving fleet names. The default is given by attr(fit$data,"fleetNames")

## Details


 Plot data available for the stock


