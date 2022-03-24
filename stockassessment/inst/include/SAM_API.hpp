#ifndef SAM_R_API
#define SAM_R_API

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define SAM_CALLDEFS							  \
  {"sam_hcr",                   (DL_FUNC) &sam_hcr,                   2}, \
  {"sam_jacobian",              (DL_FUNC) &sam_jacobian,              6}, \
  {"sam_perRecruit",            (DL_FUNC) &sam_perRecruitR,            7}, \
  {"sam_stochPerRecruit",       (DL_FUNC) &sam_stockPerRecruitR,       8}, \
  {"sam_stockRecruitmentModel", (DL_FUNC) &sam_stockRecruitmentModel, 3}, \
  {"sam_Se_sbh",                (DL_FUNC) &sam_Se_sbh,                3}, \
  {"sam_Se_sl",                 (DL_FUNC) &sam_Se_sl,                 3}, \
  {"sam_bcspline",              (DL_FUNC) &sam_bcspline,              3}, \
  {"sam_ibcspline",             (DL_FUNC) &sam_ibcspline,             3}, \
  {"sam_ibcdspline",            (DL_FUNC) &sam_ibcdspline,            3}, \
  {"sam_ibcispline",            (DL_FUNC) &sam_ibcispline,            3} 

extern "C" {


  SEXP sam_hcr(SEXP ssb, SEXP hcrConf){
    static SEXP(*fp)(SEXP, SEXP) = (SEXP(*)(SEXP, SEXP)) R_GetCCallable("stockassessment", "hcrR");
    if (fp==NULL){
      Rf_error("hcr not found");
      return R_NilValue;
    }
    return fp(ssb, hcrConf);
  }
  
  SEXP sam_jacobian(SEXP fn, SEXP par, SEXP rho, SEXP maxit, SEXP h, SEXP tolerance){
    static SEXP(*fp)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("stockassessment", "jacobian");
    if (fp==NULL){
      Rf_error("jacobian not found");
      return R_NilValue;
    }
    return fp(fn, par, rho, maxit, h, tolerance);
  }

  SEXP sam_perRecruitR(SEXP logFbar, SEXP tmbdat, SEXP pl, SEXP sel, SEXP aveYears, SEXP nYears, SEXP CT){
    static SEXP(*fp)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("stockassessment", "perRecruitR");
    if (fp==NULL){
      Rf_error("perRecruit not found");
      return R_NilValue;
    }
    return fp(logFbar, tmbdat, pl, sel, aveYears, nYears, CT);
  }


  SEXP sam_stochPerRecruitR(SEXP logFbar, SEXP tmbdat, SEXP pl, SEXP sel, SEXP aveYears, SEXP nYears, SEXP CT, SEXP logNinit){
    static SEXP(*fp)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("stockassessment", "stochPerRecruitR");
    if (fp==NULL){
      Rf_error("stochPerRecruit not found");
      return R_NilValue;
    }
    return fp(logFbar, tmbdat, pl, sel, aveYears, nYears, CT, logNinit);
  }
  
  SEXP sam_stockRecruitmentModel(SEXP ssb, SEXP rec_pars, SEXP code){
    static SEXP(*fp)(SEXP, SEXP, SEXP) = (SEXP(*)(SEXP, SEXP, SEXP)) R_GetCCallable("stockassessment", "stockRecruitmentModelR");
    if (fp==NULL){
      Rf_error("stockRecruitmentModel not found");
      return R_NilValue;
    }
    return fp(ssb, rec_pars, code);
  }

  SEXP sam_Se_sbh(SEXP lambda, SEXP a, SEXP b, SEXP g){
    static SEXP(*fp)(SEXP, SEXP, SEXP, SEXP) = (SEXP(*)(SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("stockassessment", "Se_sbhR");
    if (fp==NULL){
      Rf_error("Se_sbh not found");
      return R_NilValue;
    }
    return fp(lambda, a, b, g);
  }

  SEXP sam_Se_sl(SEXP lambda, SEXP a, SEXP b, SEXP g){
    static SEXP(*fp)(SEXP, SEXP, SEXP, SEXP) = (SEXP(*)(SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("stockassessment", "Se_slR");
    if (fp==NULL){
      Rf_error("Se_sl not found");
      return R_NilValue;
    }
    return fp(lambda, a, b, g);
  }

  SEXP sam_bcspline(SEXP x, SEXP knots, SEXP pars){
    static SEXP(*fp)(SEXP, SEXP, SEXP) = (SEXP(*)(SEXP, SEXP, SEXP)) R_GetCCallable("stockassessment", "bcsplineR");
    if (fp==NULL){
      Rf_error("bcspline not found");
      return R_NilValue;
    }
    return fp(x, knots, pars);
  }

  SEXP sam_ibcspline(SEXP x, SEXP knots, SEXP pars){
    static SEXP(*fp)(SEXP, SEXP, SEXP) = (SEXP(*)(SEXP, SEXP, SEXP)) R_GetCCallable("stockassessment", "ibcsplineR");
    if (fp==NULL){
      Rf_error("ibcspline not found");
      return R_NilValue;
    }
    return fp(x, knots, pars);
  }

  SEXP sam_ibcdspline(SEXP x, SEXP knots, SEXP pars){
    static SEXP(*fp)(SEXP, SEXP, SEXP) = (SEXP(*)(SEXP, SEXP, SEXP)) R_GetCCallable("stockassessment", "ibcdsplineR");
    if (fp==NULL){
      Rf_error("ibcdspline not found");
      return R_NilValue;
    }
    return fp(x, knots, pars);
  }

  SEXP sam_ibcispline(SEXP x, SEXP knots, SEXP pars){
    static SEXP(*fp)(SEXP, SEXP, SEXP) = (SEXP(*)(SEXP, SEXP, SEXP)) R_GetCCallable("stockassessment", "ibcisplineR");
    if (fp==NULL){
      Rf_error("ibcispline not found");
      return R_NilValue;
    }
    return fp(x, knots, pars);
  }
  
}

#endif
