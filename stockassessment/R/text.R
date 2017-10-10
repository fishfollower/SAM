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
                      'are unobserved processes. '
  ret$ages <- 'The first age group is age ' + fit$conf$minAge + ' and the last age group is age ' + fit$conf$maxAge +
               ifelse(fit$conf$maxAgePlusGroup==1,'+. ','. ') 
  ret$data <- 'The data period covers ' + fit$data$noYears + ' years (from ' + min(fit$data$years) + ' to ' +
               max(fit$data$years) + '). '+'The data contains '+fit$data$noFleets+' fleets. '
  #ret$Fmod <-
  #ret$Nmod <-
  cat(ret$modelVersion+ret$modelIntro+ret$ages+ret$data, ...)
}

##' Description of model  
##' @param fit returned object from sam.fit
##' @param ... Additional parameters to be passed to ...
##' @details Writes a string to install the version of the package which was used to run the model.
##' @export
modelVersionInfo <-function (fit,...){
  ret<-c(
    '# The fit was run with a specific version of stockassessment package.',
    '# If in the mean time version on your system has been updated',
    '# you can revert back to the version used by inserting this:',
    '',
    paste0('devtools::install_github("fishfollower/SAM/stockassessment@',attr(fit,"RemoteSha"),'")'),
    '',
    '# right before the stockassessment package is loaded'
  )
  writeLines(ret, ...)
}
