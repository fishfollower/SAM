
##' Calculate jacobian of a function
##'
##' @param func function
##' @param x parameter values
##' @param h step size
##' @param ... passed to func
##' @return jacobian matrix
jacobian <- function(func, x,
                     h = abs(1e-04 * x) + 1e-04 * (abs(x) < sqrt(.Machine$double.eps/7e-07)),
                     subset = seq_along(x),
                     ## 0.1 * 10^floor(log10(abs(x))) + 1e-4,
                     ...){
    subset <- as.integer(subset)
    if(any(subset < 1 | subset > length(x)))
        stop("Issue with subset")
        r <- .Call(C_jacobian,
                   function(x)func(x,...),
                   x,
                   globalenv(),
                   30L,
                   h,##abs(1e-4 * x) + 1e-4 * (abs(x) < 1e-8),
                   1e-12,
                   as.integer(subset-1))
        do.call("cbind",r[-1])
}
##' Calculate gradient of a function
##'
##' @param func function
##' @param x parameter values
##' @param h step size
##' @param ... passed to func
##' @return gradient vector
grad <- function(func, x,
                 h = abs(1e-04 * x) + 1e-04 * (abs(x) < sqrt(.Machine$double.eps/7e-07)),
                 ##0.1 * 10^floor(log10(abs(x))) + 1e-4,
                 ...){
         r <- .Call(C_jacobian,
                   function(x)func(x,...),
                   x,
                   globalenv(),
                   30L,
                   h,##abs(1e-4 * x) + 1e-4 * (abs(x) < 1e-8),
                   1e-12)
         v <- do.call("cbind",r[-1])
         if(nrow(v) == 1)
             return(as.vector(v))
         return(diag(v))         
}

svd_solve <- function(x){
    ss <- svd(x)
    ss$v %*% diag(1/ss$d, length(ss$d), length(ss$d)) %*% t(ss$u)
}


##' Calculate jacobian of a function
##'
##' @param func function
##' @param x parameter values
##' @param h step size
##' @param ... passed to func
##' @return jacobian matrix
##' @note Could be made more accurate in some cases
hessian <- function(func, x,
                    h = abs(1e-04 * x) + 1e-04 * (abs(x) < sqrt(.Machine$double.eps/7e-07)),
                    columns = seq_along(x),
                    ...){
    if(any(columns < 1 | columns > length(x)))
        stop("Issue with columns")
    r <- .Call(C_hessian,
               function(x)func(x,...),
               x,
               globalenv(),
               h,
               as.integer(columns-1))
    r
}

hessian_gr <- function(gr, x,
                       h = abs(1e-04 * x) + 1e-04 * (abs(x) < sqrt(.Machine$double.eps/7e-07)),
                       rows = seq_along(x),
                       ...){
    t(jacobian(function(x)gr(x, ...), x, h=h, subset=rows))
}
