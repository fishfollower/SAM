# `read.ices`

Function to read ICES/CEFAS data files and validate if input makes sense


## Description

Function to read ICES/CEFAS data files and validate if input makes sense


## Usage

```r
read.ices(filen)
```


## Arguments

Argument      |Description
------------- |----------------
`filen`     |     The filename


## Details

First two lines are ignored and can be used for comments.
 Can read formats 1 full, 2 row, 3 scalar, and 5 column
 
 Tests:
 Formatcode is valid, years and ages are pos. integers
 minimum <= maximum for years and ages
 number of rows and coulmns match year and age ranges
 data contains only numbers.
 
 Returns: A validated data matrix.


