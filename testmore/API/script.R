library(stockassessment)

src <- '
#include <SAM_API.hpp>

extern "C" {
SEXP addOne(SEXP x){
  double xr = Rf_asReal(x);
  return Rf_ScalarReal(xr + 1.0);
}
}

  static const
  R_CallMethodDef callMethods[] = {
    SAM_CALLDEFS,
    {"addOne", (DL_FUNC)&addOne, 1},
    {NULL,NULL,0}
  };

  void R_init_testAPI(DllInfo *info)
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
'


f <- file.path(tempdir(), "testAPI.cpp")
cat(src,file = f)

foutp <- tempfile(fileext = ".txt")
zz <- file(foutp, open = "wt")
sink(zz, type = "output")
sink(zz, type ="message")
TMB::compile(f,
             CXXFLAGS=sprintf("-I%s",system.file("include",package="stockassessment")),
             libinit = FALSE)
sink(type ="message")
sink(type = "output")
close(zz)


dyn.load(TMB::dynlib(file.path(tempdir(),"testAPI")))

cat(is.loaded("addOne","testAPI"), "\n", file = "res.out")
cat(is.loaded("sam_hcr","testAPI"), "\n", file = "res.out", append = TRUE)
cat(is.loaded("sam_jacobian","testAPI"), "\n", file = "res.out", append = TRUE)
cat(is.loaded("sam_perRecruit","testAPI"), "\n", file = "res.out", append = TRUE)
cat(is.loaded("sam_stockRecruitmentModel","testAPI"), "\n", file = "res.out", append = TRUE)
cat(is.loaded("sam_Se_sbh","testAPI"), "\n", file = "res.out", append = TRUE)
cat(is.loaded("sam_Se_sl","testAPI"), "\n", file = "res.out", append = TRUE)
cat(is.loaded("sam_bcspline","testAPI"), "\n", file = "res.out", append = TRUE)
cat(is.loaded("sam_ibcspline","testAPI"), "\n", file = "res.out", append = TRUE)
cat(is.loaded("sam_ibcdspline","testAPI"), "\n", file = "res.out", append = TRUE)
cat(is.loaded("sam_ibcispline","testAPI"), "\n", file = "res.out", append = TRUE)


cat(isTRUE(all.equal(.Call("addOne",1), 2)), "\n", file = "res.out", append = TRUE)
cat(isTRUE(all.equal(.Call("addOne",10), 11)), "\n", file = "res.out", append = TRUE)
cat(isTRUE(all.equal(do.call("rbind",.Call("sam_jacobian",function(x) x^2, 4, globalenv(), 20, 0.1, 1e-6)[-1])[1,1],8)), "\n", file = "res.out", append = TRUE)

