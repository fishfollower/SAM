#define WITH_SAM_LIB
#include "SAM.h"

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


extern "C" {
  
  SEXP roots2coefficients(SEXP roots){
    if(!(Rf_isReal(roots) && Rf_length(roots) > 0))
      Rf_error("par must be a non-zero length numeric vector");

    int n = Rf_length(roots);
    SEXP coef = PROTECT(Rf_allocVector(REALSXP, n+1));
    double* rR = REAL(roots);
    double* rC = REAL(coef);
    rC[n] = 1.0;

    for(int i = 1; i <= n; ++i){
      for(int j = n - i - 1; j < n; ++j){
	if(j >= 0)
	  rC[j] += (-1) * rR[i-1] * rC[j+1];
      }
    }
    UNPROTECT(1);
    return coef;  
  }

  SEXP roots2ARpar(SEXP roots){
    SEXP coef = PROTECT(roots2coefficients(roots));
    SEXP phi = PROTECT(Rf_allocVector(REALSXP, Rf_length(roots)));
    int n = Rf_length(roots);
    for(int i = 0; i < n; ++i)
      REAL(phi)[i] = -REAL(coef)[(n-1)-i];
    UNPROTECT(2);
    return phi;
  }

  SEXP logitroots2ARpar(SEXP x){
    return asSEXP(logitroots2ARpar<double>(asVector<double>(x)));
  }

}
