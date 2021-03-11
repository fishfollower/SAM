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
  SEXP hcrR(SEXP ssb, SEXP hcrConf);
  SEXP jacobian(SEXP fn, SEXP par, SEXP rho, SEXP maxit, SEXP h, SEXP tolerance);
  SEXP Se_sbhR(SEXP lambda, SEXP a, SEXP b, SEXP g);
  SEXP Se_slR(SEXP lambda, SEXP a, SEXP b, SEXP g);
  SEXP bcsplineR(SEXP x, SEXP knots, SEXP pars);
  SEXP ibcsplineR(SEXP x, SEXP knots, SEXP pars);
  SEXP ibcdsplineR(SEXP x, SEXP knots, SEXP pars);
  SEXP ibcisplineR(SEXP x, SEXP knots, SEXP pars);


  
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
    CALLDEF(hcrR,2),
    CALLDEF(jacobian,6),
    CALLDEF(Se_sbhR,4),
    CALLDEF(Se_slR,4),
    CALLDEF(bcsplineR,3),
    CALLDEF(ibcsplineR,3),
    CALLDEF(ibcdsplineR,3),
    CALLDEF(ibcisplineR,3),
    {NULL,NULL,0}
  };

#define CALLABLE(name) R_RegisterCCallable("stockassessment", #name, (DL_FUNC) &name)
  
  void R_init_stockassessment(DllInfo *info)
  {
    /* Register the .C and .Call routines.
       No .Fortran() or .External() routines,
       so pass those arrays as NULL.
    */
    R_registerRoutines(info,
		       NULL, callMethods,
		       NULL, NULL);

    CALLABLE(hcrR);
    CALLABLE(jacobian);
    CALLABLE(perRecruitR);
    CALLABLE(stockRecruitmentModelR);
    CALLABLE(Se_sbhR);
    CALLABLE(Se_slR);
    CALLABLE(bcsplineR);
    CALLABLE(ibcsplineR);
    CALLABLE(ibcdsplineR);
    CALLABLE(ibcisplineR);

    R_useDynamicSymbols(info, (Rboolean)FALSE);
  }


}
