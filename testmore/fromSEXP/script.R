library(stockassessment)

src <- '
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <TMB.hpp>
#include <SAM.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
return 0;
}

extern "C" {

SEXP testData(SEXP x) {
   dataSet<double> d1(x);
   dataSet<double> d2 = d1;
   return Rf_mkString("OK");
}

SEXP testConf(SEXP x) {
   confSet c1(x);
   confSet c2 = c1;
   return Rf_mkString("OK");
}

R_CallMethodDef callMethods[] = {
 {"testData", (DL_FUNC) &testData, 1},
 {"testConf", (DL_FUNC) &testConf, 1},
 {NULL, NULL, 0}
};

  void R_init_%s(DllInfo *info)
  {
    /* Register the .C and .Call routines.
       No .Fortran() or .External() routines,
       so pass those arrays as NULL.
    */
    R_registerRoutines(info,
		       NULL, callMethods,
		       NULL, NULL);
    R_useDynamicSymbols(info, (Rboolean)FALSE);
  }


}
'

f <- tempfile(fileext = ".cpp")
if(file.exists(TMB::dynlib(sub("\\.cpp$","",f))))
    file.remove(TMB::dynlib(sub("\\.cpp$","",f)))
            
cat(sprintf(src,gsub("((^.+/)|(\\.cpp$))","",f)),file = f)

foutp <- tempfile(fileext = ".txt")
zz <- file(foutp, open = "wt")
sink(zz, type = "output")
sink(zz, type ="message")
TMB::compile(f,
             CXXFLAGS = paste0("-std=c++14 -Wno-ignored-attributes -I", system.file("include", package = "stockassessment")),
             framework="TMBad",
             libinit = FALSE)
sink(type ="message")
sink(type = "output")
close(zz)

if(file.exists(TMB::dynlib(sub("\\.cpp$","",f)))){

dyn.load(TMB::dynlib(sub("\\.cpp$","",f)))


data(nscodData)
data(nscodConf)
par<-defpar(nscodData,nscodConf)

fit<-sam.fit(nscodData,nscodConf,par,run = FALSE)



cat(tryCatch(.Call("testData", x=fit$obj$env$data), error = function(e){paste("NOT_OK:",e)}),
    "\n",
    sep = "",
    file = "res.out")
cat(tryCatch(.Call("testConf", x=fit$obj$env$data), error = function(e){paste("NOT_OK:",e)}),
    "\n",
    file = "res.out",
    sep = "",
    append = TRUE)


}else{
    cat("Compilation failed",
        "\n",
        sep = "",
        file = "res.out")
}
