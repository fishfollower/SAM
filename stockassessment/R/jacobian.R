
##' Calculate jacobian of a function
##'
##' @param func function
##' @param x parameter values
##' @param h step size
##' @param ... passed to func
##' @return jacobian matrix
jacobian <- function(func, x,
                     h = abs(1e-04 * x) + 1e-04 * (abs(x) < sqrt(.Machine$double.eps/7e-07)),
                     maxit = 30L,
                     tol = 1e-12,
                     subset = seq_along(x),
                     ## 0.1 * 10^floor(log10(abs(x))) + 1e-4,
                     ...){
    subset <- as.integer(subset)
    if(any(subset < 1 | subset > length(x)))
        stop("Issue with subset")
        r <- .Call(C_jacobian,
                   function(x)func(x,...),
                   as.numeric(x),
                   globalenv(),
                   as.integer(maxit),
                   h,##abs(1e-4 * x) + 1e-4 * (abs(x) < 1e-8),
                   tol,
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
                   as.numeric(x),
                   globalenv(),
                   30L,
                   h,##abs(1e-4 * x) + 1e-4 * (abs(x) < 1e-8),
                   1e-12,
                   as.integer(seq_along(x)-1))
         v <- do.call("cbind",r[-1])
         if(nrow(v) == 1)
             return(as.vector(v))
         return(diag(v))         
}

svd_solve <- function(x){
    ss <- svd(x)
    ss$v %*% diag(1/ss$d, length(ss$d), length(ss$d)) %*% t(ss$u)
}


##' Calculate hessian of a function
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
               as.numeric(x),
               globalenv(),
               h,
               as.integer(columns-1))
    r
}

hessian_gr <- function(gr, x,
                       h = abs(1e-04 * x) + 1e-04 * (abs(x) < sqrt(.Machine$double.eps/7e-07)),
                       maxit = 4L,
                       tol = 1e-4,
                       subset = seq_along(x),
                       method = c("romberg","central","forward"),
                       symmetrize = TRUE,
                       ...){
    method <- match.arg(method)

    if(any(subset < 1 | subset > length(x)))
        stop("Issue with subset")
  
    if(method == "romberg"){
        r <- (t(jacobian(function(x)gr(x, ...), x, h=h,maxit=maxit,tol=tol, subset=subset)))
    }else if(method == "central"){
        r <- .Call(C_hessian_gr_central,
                     function(x)gr(x,...),
                     as.numeric(x),
                     globalenv(),
                     h,
                     as.integer(subset-1))
    }else if(method == "forward"){
        r <- .Call(C_hessian_gr_forward,
                     function(x)gr(x,...),
                     as.numeric(x),
                     globalenv(),
                     h,
                     as.integer(subset-1))
    }else{
        stop("Wrong method")
    }
    add <- function(a,b){
        ifelse(is.na(a),b,a)+ifelse(is.na(b),a,b)        
    }
    if(symmetrize)
        r[] <- 0.5 * add(r[],t(r)[])
    r
}
