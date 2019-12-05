#define WITH_LIBTMB
#include <TMB.hpp>

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


extern "C" {
#include <R_ext/Rdynload.h>

  SEXP perRecruitR(SEXP Fbar, SEXP dat, SEXP conf, SEXP pl, SEXP sel, SEXP aveYears, SEXP nYears);
  SEXP stockRecruitmentModelR(SEXP ssb, SEXP rec_pars, SEXP code);
  SEXP jacobian(SEXP fn, SEXP par, SEXP rho, SEXP maxit, SEXP h, SEXP tolerance);
  SEXP Se_sbhR(SEXP lambda, SEXP a, SEXP b, SEXP g);
  SEXP Se_slR(SEXP lambda, SEXP a, SEXP b, SEXP g);
  
#define CALLDEF(name,n) {#name, (DL_FUNC) &name, n}
  
  static const
  R_CallMethodDef callMethods[] = {

    // TMB
    #ifdef TMB_CALLDEFS
    TMB_CALLDEFS,
    #else
    CALLDEF(MakeADFunObject, 4),
    CALLDEF(InfoADFunObject, 1),
    CALLDEF(EvalADFunObject, 3),
    CALLDEF(MakeDoubleFunObject, 3),
    CALLDEF(EvalDoubleFunObject, 3),
    CALLDEF(getParameterOrder, 3),
    CALLDEF(MakeADGradObject, 3),
    CALLDEF(MakeADHessObject2, 4),
    CALLDEF(usingAtomics, 0),
    CALLDEF(TMBconfig, 2),
    #endif
    
    CALLDEF(perRecruitR,7),
    CALLDEF(stockRecruitmentModelR,3),
    CALLDEF(jacobian,6),
    CALLDEF(Se_sbhR,4),
    CALLDEF(Se_slR,4),

    {NULL,NULL,0}
  };

  void R_init_stockassessment(DllInfo *info)
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
