# `reduce`: reduce helper function to reduce data

## Description


 reduce helper function to reduce data


## Usage

```r
reduce(data, year = NULL, fleet = NULL, age = NULL, conf = NULL)
```


## Arguments

Argument      |Description
------------- |----------------
```data```     |     a data object as returned by the function setup.sam.data
```year```     |     a vector of years to be excluded.
```fleet```     |     a vector of fleets to be excluded.
```age```     |     a vector of ages fleets to be excluded.
```conf```     |     an optional corresponding configuration to be modified along with the data change. Modified is returned as attribute "conf"

## Details


 When more than one vector is supplied they need to be of same length, as only the pairs are excluded


