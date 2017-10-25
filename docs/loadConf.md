# `loadConf`: Loads a model configuration from a file

## Description


 Loads a model configuration from a file


## Usage

```r
loadConf(dat, file, patch = FALSE)
```


## Arguments

Argument      |Description
------------- |----------------
```dat```     |     sam data list as returned from the function setup.sam.data
```file```     |     the file to read the configuration from
```patch```     |     logical if TRUE missing entries will be automatically filled with default values

## Details


 function useful loading a model configuration. Such a configuration can be saved via the saveConf function


