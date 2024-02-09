// #define WTIH_LIBTMB
// #include "TMB.h"

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


extern "C" {
#include <R_ext/Rdynload.h>

  SEXP MakeADFunObject(SEXP data, SEXP parameters, SEXP report, SEXP control);
  SEXP FreeADFunObject(SEXP f);
  SEXP InfoADFunObject(SEXP f);
  SEXP tmbad_print(SEXP f, SEXP control);
  SEXP EvalADFunObject(SEXP f, SEXP theta, SEXP control);
  SEXP TransformADFunObject(SEXP f, SEXP control);
  SEXP MakeDoubleFunObject(SEXP data, SEXP parameters, SEXP report, SEXP control);
  SEXP EvalDoubleFunObject(SEXP f, SEXP theta, SEXP control);
  SEXP getParameterOrder(SEXP data, SEXP parameters, SEXP report, SEXP control);
  SEXP MakeADGradObject(SEXP data, SEXP parameters, SEXP report, SEXP control);
  SEXP MakeADHessObject2(SEXP data, SEXP parameters, SEXP report, SEXP control);
  SEXP usingAtomics();
  SEXP getFramework();
  SEXP getSetGlobalPtr(SEXP ptr);
  SEXP TMBconfig(SEXP envir, SEXP cmd);  
  
  SEXP perRecruitR(SEXP logFbar, SEXP tmbdat, SEXP pl, SEXP sel, SEXP aveYears, SEXP nYears, SEXP CT);
  SEXP perRecruitSR(SEXP logFbar, SEXP dat, SEXP conf, SEXP pl, SEXP sel, SEXP aveYears, SEXP nYears, SEXP CT, SEXP logNinit);
  SEXP perRecruitSR_Calc(SEXP logFbar, SEXP dat, SEXP conf, SEXP pl, SEXP sel, SEXP aveYears, SEXP nYears, SEXP CT, SEXP logNinit, SEXP DT);
  SEXP stockRecruitmentModelR(SEXP ssb, SEXP rec_pars, SEXP code, SEXP constRecBreaks, SEXP year, SEXP lastR);
  SEXP logSRR(SEXP logssb, SEXP rec_pars, SEXP code, SEXP constRecBreaks, SEXP year, SEXP lastR);
  SEXP hcrR(SEXP ssb, SEXP hcrConf);
  SEXP jacobian(SEXP fn, SEXP par, SEXP rho, SEXP maxit, SEXP h, SEXP tolerance, SEXP subset);
  SEXP hessian(SEXP fn, SEXP par, SEXP rho, SEXP h, SEXP columns);
  SEXP hessian_gr_central(SEXP fn, SEXP par, SEXP rho, SEXP h, SEXP subset);
  SEXP hessian_gr_forward(SEXP fn, SEXP par, SEXP rho, SEXP h, SEXP subset);
  SEXP bcsplineR(SEXP x, SEXP knots, SEXP pars);
  SEXP ibcsplineR(SEXP x, SEXP knots, SEXP pars);
  SEXP ibcdsplineR(SEXP x, SEXP knots, SEXP pars);
  SEXP ibcisplineR(SEXP x, SEXP knots, SEXP pars);
  SEXP iibcisplineR(SEXP x, SEXP knots, SEXP pars);
  SEXP splinebasis_bcR(SEXP x, SEXP knots);
  SEXP splinebasis_ibcR(SEXP x, SEXP knots);
  SEXP splinebasis_iibcR(SEXP x, SEXP knots);
  SEXP recruitmentProperties(SEXP tmbdat, SEXP pl);
  SEXP roots2coefficients(SEXP roots);
  SEXP roots2ARpar(SEXP roots);
  SEXP logitroots2ARpar(SEXP x);
  
#define CALLDEF(name,n) {#name, (DL_FUNC) &name, n}
  
  static const
  R_CallMethodDef callMethods[] = {

				   // TMB
				   {"MakeADFunObject",     (DL_FUNC) &MakeADFunObject,     4},   
				   {"FreeADFunObject",     (DL_FUNC) &FreeADFunObject,     1},   
				   {"InfoADFunObject",     (DL_FUNC) &InfoADFunObject,     1},   
				   {"tmbad_print",         (DL_FUNC) &tmbad_print,         2},   
				   {"EvalADFunObject",     (DL_FUNC) &EvalADFunObject,     3},   
				   {"TransformADFunObject",(DL_FUNC) &TransformADFunObject,2},   
				   {"MakeDoubleFunObject", (DL_FUNC) &MakeDoubleFunObject, 4},   
				   {"EvalDoubleFunObject", (DL_FUNC) &EvalDoubleFunObject, 3},   
				   {"getParameterOrder",   (DL_FUNC) &getParameterOrder,   4},   
				   {"MakeADGradObject",    (DL_FUNC) &MakeADGradObject,    4},   
				   {"MakeADHessObject2",   (DL_FUNC) &MakeADHessObject2,   4},   
				   {"usingAtomics",        (DL_FUNC) &usingAtomics,        0},   
				   {"getFramework",        (DL_FUNC) &getFramework,        0},
				   {"getSetGlobalPtr",     (DL_FUNC) &getFramework,        1},
				   {"TMBconfig",           (DL_FUNC) &TMBconfig,           2},
      
    CALLDEF(perRecruitR,7),
    CALLDEF(perRecruitSR,8),
    CALLDEF(perRecruitSR_Calc,9),
    CALLDEF(stockRecruitmentModelR,6),
    CALLDEF(logSRR,6),
    CALLDEF(hcrR,2),
    CALLDEF(jacobian,7),
    CALLDEF(hessian,5),
    CALLDEF(hessian_gr_central,5),
    CALLDEF(hessian_gr_forward,5),
    CALLDEF(bcsplineR,3),
    CALLDEF(ibcsplineR,3),
    CALLDEF(ibcdsplineR,3),
    CALLDEF(ibcisplineR,3),
    CALLDEF(iibcisplineR,3),
    CALLDEF(splinebasis_bcR,2),
    CALLDEF(splinebasis_ibcR,2),
    CALLDEF(splinebasis_iibcR,2),
    CALLDEF(recruitmentProperties,2),
				   CALLDEF(roots2coefficients,1),
				   CALLDEF(roots2ARpar,1),
				   CALLDEF(logitroots2ARpar,1),
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
    CALLABLE(hessian);
    CALLABLE(hessian_gr_central);
    CALLABLE(hessian_gr_forward);
    CALLABLE(perRecruitR);
    CALLABLE(perRecruitSR);
    CALLABLE(perRecruitSR_Calc);
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
    CALLABLE(recruitmentProperties);
    CALLABLE(roots2coefficients);
    CALLABLE(roots2ARpar);
    CALLABLE(logitroots2ARpar);
 
    R_useDynamicSymbols(info, (Rboolean)FALSE);
    R_forceSymbols(info, (Rboolean)FALSE);
  }


}
