##' Description of model  
##' @param fit returned object from sam.fit
##' @param ... Additional parameters to be passed to ...
##' @importFrom utils packageVersion
##' @details ...
##' @export
modelDescription <-function (fit,...){
   ret<-list()
   ret$modelVersion<-paste('The model is a state-space stock assessment (SAM from the package "stockassessment" version ',
                           packageVersion("stockassessment"), ".", sep="")
   ret$modelIntuition<-'The model works by assuming that stock-sizes at age (N) and fishing mortalities (F) are unobserved ...' 
   cat(ret$modelVersion,ret$modelIntuition)
}
