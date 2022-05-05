#define WITH_LIBTMB
#include <TMB.hpp>

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


extern "C" {
#include <R_ext/Rdynload.h>

  SEXP perRecruitR(SEXP logFbar, SEXP tmbdat, SEXP pl, SEXP sel, SEXP aveYears, SEXP nYears, SEXP CT);
  SEXP perRecruitSR(SEXP logFbar, SEXP dat, SEXP conf, SEXP pl, SEXP sel, SEXP aveYears, SEXP nYears, SEXP CT, SEXP logNinit);
  SEXP stockRecruitmentModelR(SEXP ssb, SEXP rec_pars, SEXP code, SEXP constRecBreaks, SEXP year, SEXP lastR);
  SEXP logSRR(SEXP logssb, SEXP rec_pars, SEXP code, SEXP constRecBreaks, SEXP year, SEXP lastR);
  SEXP hcrR(SEXP ssb, SEXP hcrConf);
  SEXP jacobian(SEXP fn, SEXP par, SEXP rho, SEXP maxit, SEXP h, SEXP tolerance);
  SEXP bcsplineR(SEXP x, SEXP knots, SEXP pars);
  SEXP ibcsplineR(SEXP x, SEXP knots, SEXP pars);
  SEXP ibcdsplineR(SEXP x, SEXP knots, SEXP pars);
  SEXP ibcisplineR(SEXP x, SEXP knots, SEXP pars);
  SEXP iibcisplineR(SEXP x, SEXP knots, SEXP pars);
  SEXP splinebasis_bcR(SEXP x, SEXP knots);
  SEXP splinebasis_ibcR(SEXP x, SEXP knots);
  SEXP splinebasis_iibcR(SEXP x, SEXP knots);
    
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
    CALLDEF(perRecruitSR,8),
    CALLDEF(stockRecruitmentModelR,6),
    CALLDEF(logSRR,6),
    CALLDEF(hcrR,2),
    CALLDEF(jacobian,6),
    CALLDEF(bcsplineR,3),
    CALLDEF(ibcsplineR,3),
    CALLDEF(ibcdsplineR,3),
    CALLDEF(ibcisplineR,3),
    CALLDEF(iibcisplineR,3),
    CALLDEF(splinebasis_bcR,2),
    CALLDEF(splinebasis_ibcR,2),
    CALLDEF(splinebasis_iibcR,2),
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
    CALLABLE(perRecruitSR);
    CALLABLE(stockRecruitmentModelR);
    CALLABLE(logSRR);
    CALLABLE(bcsplineR);
    CALLABLE(ibcsplineR);
    CALLABLE(ibcdsplineR);
    CALLABLE(ibcisplineR);
    CALLABLE(iibcisplineR);
    CALLABLE(splinebasis_bcR);
    CALLABLE(splinebasis_ibcR);
    CALLABLE(splinebasis_iibcR);
 
    R_useDynamicSymbols(info, (Rboolean)FALSE);
  }


}
