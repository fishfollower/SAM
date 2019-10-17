library(stockassessment)

tol <- sqrt(.Machine$double.eps)

set.seed(298374)

jacobian <- function(func, x, ...){
         r <- .Call("jacobian",
                   function(x)func(x,...),
                   x,
                   globalenv(),
                   30L,
                   0.1 * 10^floor(log10(abs(x))) + 1e-4,##abs(1e-4 * x) + 1e-4 * (abs(x) < 1e-8),
                   1e-12)
        do.call("cbind",r[-1])
}
grad <- function(func, x, ...){
         r <- .Call("jacobian",
                   function(x)func(x,...),
                   x,
                   globalenv(),
                   30L,
                   0.1 * 10^floor(log10(abs(x))) + 1e-4,##abs(1e-4 * x) + 1e-4 * (abs(x) < 1e-8),
                   1e-12)
         v <- do.call("cbind",r[-1])
         if(nrow(v) == 1)
             return(as.vector(v))
         return(diag(v))         
}


gr.sym <- function(x, fn){
    eval(D(parse(text=as.character(functionBody(fn))[2]),"x"))
}

check <- function(x, fn, gr = grad){
    r <- (gr(fn,x) / gr.sym(x,fn) - 1)
    b <- max(abs(r)) < tol
    attr(b,"error") <- max(abs(r))
    b
}

## Test
fn0 <- function(x){
    sin(x)
}
x0 <- runif(1000, -pi, pi)
v <- check(x0,fn0)
cat(isTRUE(v),"\n",
    sep = "",
    file = "res.out")

## Test
fn1 <- function(x){
    cos(x)
}
x1 <- runif(1000, -pi, pi)
v <- check(x1,fn1)
cat(isTRUE(v),"\n",
    sep = "",
    file = "res.out",
    append = TRUE)

## Test
fn2 <- function(x){
    exp(x)
}
x2 <- runif(1000, -pi, pi)
v <- check(x2,fn2)
cat(isTRUE(v),"\n",
    sep = "",
    file = "res.out",
    append = TRUE)


## Test
fn3 <- function(x){
    exp(x) / sqrt(cos(x)^3 + sin(x)^3)
}
x3 <- runif(1000, 0, 2)
v <- check(x3,fn3)
cat(isTRUE(v),"\n",
    sep = "",
    file = "res.out",
    append = TRUE)


## Test
fn4 <- function(x){
    1/x
}
x4 <- runif(1000, -5, 5)
v <- check(x4,fn4)
cat(isTRUE(v),"\n",
    sep = "",
    file = "res.out",
    append = TRUE)



## Test
fn5 <- function(x){
    1/x^2
}
x5 <- runif(1000, -5, 5)
v <- check(x5,fn5)
cat(isTRUE(v),"\n",
    sep = "",
    file = "res.out",
    append = TRUE)


## Test
fn6 <- function(x){
    exp(x) / (sin(x) - x^2)
}
x6 <- runif(1000, -5, 5)
v <- check(x6,fn6)
cat(isTRUE(v),"\n",
    sep = "",
    file = "res.out",
    append = TRUE)



## Test
fn7 <- function(x){
    asin(x)
}
x7 <- runif(1000, 0, 1)
suppressWarnings(v <- check(x7,fn7)) ## default step size produces NaN for x close to 1
cat(isTRUE(v),"\n",
    sep = "",
    file = "res.out",
    append = TRUE)



## Test
fn8 <- function(x){
    (x-100)^2 + 1.e-06 * (x - 300)^3
}
x8 <- runif(1000, -5, 5)
v <- check(x8,fn8)
cat(isTRUE(v),"\n",
    sep = "",
    file = "res.out",
    append = TRUE)

## Test
fn9 <- function(x){
    digamma(x) * x
}
x9 <- runif(1000, 0, 5)
v <- check(x9,fn9)
cat(isTRUE(v),"\n",
    sep = "",
    file = "res.out",
    append = TRUE)

## Test
fn10 <- function(x){
    log1p(dnorm(x))
}
x10 <- runif(1000, -5, 5)
v <- check(x10,fn10)
cat(isTRUE(v),"\n",
    sep = "",
    file = "res.out",
    append = TRUE)



### Compare to other methods
if(FALSE){

    gr1 <- function(fn, x){
        e <- new.env()
        e$x <- x
        diag(attr(stats::numericDeriv(quote(fn(x)),"x",e),"gradient"))
    }

    gr2 <- function(fn, x){
        numDeriv::grad(fn,x)
    }

    compare <- function(fn,x){
        v <- c("stats" = attr(check(x,fn,gr1),"error"),
               "numDeriv" = attr(check(x,fn,gr2),"error"),
               "SAM" = attr(check(x,fn,grad),"error"))
        cat("Not worse than numDeriv:",v[3] <= v[2],"\n")
        cat("Not much worse than numDeriv:",v[3] < 1.05 * v[2],"\n")
        cat("Not worse than stats:",v[3] <= v[1],"\n")
        sort(v)
    }
          

compare(fn0,x0)
compare(fn1,x1)
compare(fn2,x2)
compare(fn3,x3)
compare(fn4,x4)
compare(fn5,x5)
compare(fn6,x6)
compare(fn7,x7)
compare(fn8,x8)
compare(fn9,x9)
compare(fn10,x10)


}
