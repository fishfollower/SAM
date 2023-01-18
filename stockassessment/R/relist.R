#' @export
#' @importFrom utils relist
#' @method relist array
relist.array <- function (flesh, skeleton = attr(flesh, "skeleton")) 
{
    d <- dim(skeleton)
    if(sum(d)==0)
        return(array(0,dim=d, dimnames = dimnames(skeleton)))
    if (is.numeric(skeleton[1])) 
        return(array(flesh, dim=d, dimnames = dimnames(skeleton)))
    ## stop("Only implemented for numeric arrays")
    result <- skeleton
    ind <- 1L
    for (i in seq_along(skeleton)) {
        size <- length(unlist(skeleton[[i]]))
        result[[i]] <- relist(flesh[seq.int(ind, length.out = size)], 
            skeleton[[i]])
        ind <- ind + size
    }
    result
}
