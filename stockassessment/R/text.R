##' Description of model  
##' @param fit returned object from sam.fit
##' @param ... Additional parameters to be passed to ...
##' @importFrom utils packageVersion
##' @details ...
##' @export
modelDescription <-function (fit,...){
  "+" = function(x,y) {
      if(is.character(x) || is.character(y)) {
        return(paste(x , y, sep=""))
      } else {
        .Primitive("+")(x,y)
      }
  } 
  ret<-list()
  ret$modelVersion <- 'The model is a state-space stock assessment (SAM from the package "stockassessment" version ' +
                      packageVersion("stockassessment") + ").\n\n"
  ret$modelIntro <- 'The model works by assuming that stock-sizes at age (N) and fishing mortalities at age (F) ' +
                      'are unobserved processes.'
  ret$ages <-  'The first age group is age ' + fit$conf$minAge + ' and the last age group is age ' + fit$conf$maxAge +
                      ifelse(fit$conf$maxAgePlusGroup==1,'+. ','. ') 
  ret$data <- 'The data period covers ' + fit$data$noYears + ' years (from ' + min(fit$data$years) + ' to ' +
              max(fit$data$years) + '). ' 
  cat(ret$modelVersion+ret$modelIntro+ret$ages+ret$data)
}
