# `write.ices`

Write ICES/CEFAS data file from matrix


## Description

Write ICES/CEFAS data file from matrix


## Usage

```r
write.ices(x, fileout, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`x`     |     a matrix where rownames are taken as years and colnames are taken as ages
`fileout`     |     file name or connection
`...`     |     Arguments to be passed to write


## Details

Takes the data and writes them in the ICES/CEFAS format. It is assumed that rows represent consecutive years and cols consecutive ages


