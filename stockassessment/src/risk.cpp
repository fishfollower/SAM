#define WITH_SAM_LIB
#include "SAM.h"

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


extern "C" {

  
  SEXP logRiskHazard(SEXP x, SEXP m, SEXP logk, SEXP loga, SEXP logb, SEXP model){
    if(!(Rf_isReal(x) && Rf_length(x) > 0))
      Rf_error("x must be a non-zero length numeric vector");
    if(!(Rf_isReal(m) && Rf_length(m) > 0))
      Rf_error("m must be a non-zero length numeric vector");
    if(!(Rf_isReal(logk) && Rf_length(logk) > 0))
      Rf_error("logk must be a non-zero length numeric vector");
    if(!(Rf_isReal(loga) && Rf_length(loga) > 0))
      Rf_error("loga must be a non-zero length numeric vector");
    if(!(Rf_isReal(logb) && Rf_length(logb) > 0))
      Rf_error("logb must be a non-zero length numeric vector");
    if(!(Rf_isInteger(model) && Rf_length(model) > 0))
      Rf_error("model must be a non-zero length integer vector");

    int n1 = Rf_length(x);
    int n2 = Rf_length(m);
    int n3 = Rf_length(logk);
    int n4 = Rf_length(loga);
    int n5 = Rf_length(logb);
    int n6 = Rf_length(model);
    int n = std::max({n1,n2,n3,n4,n5,n6});
    SEXP res = PROTECT(Rf_allocVector(REALSXP, n));
    double* p_r = REAL(res);
    double* p_x = REAL(x);
    double* p_m = REAL(m);
    double* p_logk = REAL(logk);
    double* p_loga = REAL(loga);
    double* p_logb = REAL(logb);
    int* p_mod = INTEGER(model);
    for(int i = 0; i < n; ++i){
      p_r[i] = logRiskHazard<double>(p_x[i % n1],
				     p_m[i % n2],
				     p_logk[i % n3],
				     p_loga[i % n4],
				     p_logb[i % n5],
				     p_mod[i % n6]);
    }
    UNPROTECT(1);
    return res;  
  }

}
