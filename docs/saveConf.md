# `saveConf`

Saves a model configuration list to a file


## Description

Saves a model configuration list to a file


## Usage

```r
saveConf(x, file = "", overwrite = FALSE)
```


## Arguments

Argument      |Description
------------- |----------------
`x`     |     sam configuration list as returned from defcon or loadConf
`file`     |     the file to save the configuration to
`overwrite`     |     logical if an existing file should be overwritten (FALSE by default)


## Details

function useful for saving a model configuration. A saved configuration can be read back in via the loadConf function


