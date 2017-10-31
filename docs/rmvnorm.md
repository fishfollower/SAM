# `rmvnorm`: rmvnorm helper function to draw multivariate normal samples

## Description


 rmvnorm helper function to draw multivariate normal samples


## Usage

```r
rmvnorm(n = 1, mu, Sigma)
```


## Arguments

Argument      |Description
------------- |----------------
```n```     |     the number of samples.
```mu```     |     the mean vector.
```Sigma```     |     a positive-definite symmetric matrix specifying the covariance matrix.

## Details


 Generates samples via the Cholesky decomposition, which is less platform dependent than eigenvalue decomposition.


## Value


 If n = 1 a vector of the same length as mu, otherwise an n by length(mu) matrix with one sample in each row.


