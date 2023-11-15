#define WITH_LIBTMB
#include "TMB.h"

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <algorithm>

typedef struct r_function {

  SEXP fcall;
  SEXP env;
  int nIn;
  int nOut;

  vector<double> operator()(vector<double> x){
    return this->operator()(asSEXP(x));
  }

    vector<double> operator()(SEXP x){
    SEXP ans;
    SETCADR(fcall, x);    
    PROTECT(ans = Rf_duplicate(Rf_eval(fcall, env)));
    UNPROTECT(1);
    return asVector<double>(ans);
  }

} r_function, *Rfunction;


extern "C" {

  /*
    Implementation of:
    C.J.F. Ridders,
    Accurate computation of F′(x) and F′(x) F″(x),
    Advances in Engineering Software (1978),
    Volume 4, Issue 2,
    1982,
    Pages 75-76,
    ISSN 0141-1195,
    https://doi.org/10.1016/S0141-1195(82)80057-0.
   */

  SEXP jacobian(SEXP fn, SEXP par, SEXP rho, SEXP maxit, SEXP h, SEXP tolerance, SEXP subset){
    SEXP ans;


    // Check input
    //// fn must be a function
    if(!Rf_isFunction(fn))
      Rf_error("fn must be a function");
    ///// par is a numeric vector
    if(!(Rf_isReal(par) && Rf_length(par) > 0))
      Rf_error("par must be a non-zero length numeric vector");
    ///// rho is an environment
    if(!Rf_isEnvironment(rho))
      Rf_error("rho must be an environment");
    ///// maxit is an integer
    if(!(Rf_isInteger(maxit) || Rf_isReal(maxit)) && Rf_length(maxit) == 1)
      Rf_error("maxit must be a single integer");
    ///// h is a vector of the same length as par
    if(!(Rf_isReal(h) && Rf_length(h) == Rf_length(par)))
      Rf_error("h must be a numeric vector with the same length as par");
    //// tolerance is a scalar
    if(!(Rf_isReal(tolerance) && Rf_length(tolerance) == 1))
      Rf_error("tolerance must be a numeric scalar");

    
    Rfunction RF = (Rfunction) R_alloc(1, sizeof(r_function));
    RF->fcall = PROTECT(Rf_lang2(fn, R_NilValue));
    RF->env = rho;
    RF->nIn = Rf_length(par);

    vector<double> funres = RF->operator()(par);
    // Check if values are finite
    // for(int i = 0; i < funres.size(); ++i)
    //   if(!R_finite(funres(i)))
    if(!funres.unaryExpr(&R_finite).all())
      Rf_error("Function evaluation returned non-finite values");
    
    RF->nOut = funres.size();
    
    ans = PROTECT(Rf_allocVector(VECSXP, RF->nIn+1));
    SET_VECTOR_ELT(ans, 0, asSEXP(funres));

    vector<double> dpar = asVector<double>(par);

    vector<double> h0 = asVector<double>(h);
    
    int ntab = Rf_asInteger(maxit);
    double tol = Rf_asReal(tolerance);
    // Loop over parameters
    for(int p0 = 0; p0 < Rf_length(subset); ++p0){
      int p = INTEGER(subset)[p0];
      vector<double> dir(RF->nIn);
      dir.setZero();
      dir(p) = 1.0;

      double err = R_PosInf;
      double hh = 2.0 * h0(p);
      
      vector<double> v1(RF->nIn);
      vector<double> v2(RF->nIn);

      do {
	hh /= 2.0;
	v1 = RF->operator()(dpar + dir * hh);
	v2 = RF->operator()(dpar - dir * hh);
      }while(!v1.unaryExpr(&R_finite).all() || !v2.unaryExpr(&R_finite).all());

      array<double> tab(RF->nOut, ntab, ntab);
      tab.setZero();

      tab.col(0).col(0) = (v1 - v2) / (2.0 * hh);

      for(int i = 1; i < ntab; ++i){
	do {
	  hh /= 2.0;
	  v1 = RF->operator()(dpar + dir * hh);
	  v2 = RF->operator()(dpar - dir * hh);
	}while(!v1.unaryExpr(&R_finite).all() || !v2.unaryExpr(&R_finite).all());
	
	tab.col(i).col(0) = (v1 - v2) / (2.0 * hh);
	double f = 1.0;
	for(int m = 1; m <= i; ++m){
	  f *= 4.0;
	  tab.col(i).col(m) = ((vector<double>)tab.col(i).col(m-1) * f - (vector<double>)tab.col(i-1).col(m-1)) / (f - 1.0);
	  double errtmp = std::max((tab.col(i).col(m) - tab.col(i).col(m-1)).abs().maxCoeff(),
			 (tab.col(i).col(m) - tab.col(i-1).col(i-1)).abs().maxCoeff());
	  if(errtmp < err){
	    SET_VECTOR_ELT(ans, p+1, asSEXP((vector<double>)tab.col(i).col(m)));
	    err = errtmp;
	    if(err < tol)
	      goto endloop;
	  }
	}
	if((tab.col(i).col(i) - tab.col(i-1).col(i-1)).abs().maxCoeff() > 2.0 * err)
	  break;
      }
    endloop:
      Rf_setAttrib(VECTOR_ELT(ans,p+1), Rf_install("error"), asSEXP(err));
    }    
    
    UNPROTECT(2);
    return ans;
  
  }


  SEXP hessian(SEXP fn, SEXP par, SEXP rho, SEXP h, SEXP columns){
    SEXP ans;


    // Check input
    //// fn must be a function
    if(!Rf_isFunction(fn))
      Rf_error("fn must be a function");
    ///// par is a numeric vector
    if(!(Rf_isReal(par) && Rf_length(par) > 0))
      Rf_error("par must be a non-zero length numeric vector");
    ///// rho is an environment
    if(!Rf_isEnvironment(rho))
      Rf_error("rho must be an environment");
    ///// h is a vector of the same length as par
    if(!(Rf_isReal(h) && Rf_length(h) == Rf_length(par)))
      Rf_error("h must be a numeric vector with the same length as par");
   
    
    Rfunction RF = (Rfunction) R_alloc(1, sizeof(r_function));
    RF->fcall = PROTECT(Rf_lang2(fn, R_NilValue));
    RF->env = rho;
    RF->nIn = Rf_length(par);

    vector<double> funres = RF->operator()(par);
    // Check if values are finite
    // for(int i = 0; i < funres.size(); ++i)
    //   if(!R_finite(funres(i)))
    if(!funres.unaryExpr(&R_finite).all())
      Rf_error("Function evaluation returned non-finite values");
    
    RF->nOut = funres.size();
    if(RF->nOut != 1)
      Rf_error("Hessian expects a scalar valued function");

    int npar = RF->nIn;
    ans = PROTECT(Rf_allocMatrix(REALSXP, npar, npar));
    double* pa = REAL(ans);

    if(Rf_length(columns) < npar){
      for(int i = 0; i < npar; ++i)
	for(int j = 0; j < npar; ++j)
	  pa[i * npar + j] = R_NaReal;
    }
    
    vector<double> dpar = asVector<double>(par);

    vector<double> h0 = asVector<double>(h);
    
    double f0 = funres(0);
    // Diagonal
    for(int i = 0; i < npar; ++i){
      vector<double> dir(npar);
      dir.setZero();
      dir(i) = h0(i);
      double f1 = RF->operator()(dpar + dir)(0);
      double f2 = RF->operator()(dpar - dir)(0);
      pa[i * npar + i] = (f1 + f2 - 2.0 * f0) / ( h0(i) * h0(i) );
    }
    for(int k = 0; k < Rf_length(columns); ++k){ // Columns
      int j = INTEGER(columns)[k];
      for(int i = k + 1; i < npar; ++i){ // Rows, lower tri
	vector<double> dir(npar);
	dir.setZero();
	dir(i) = h0(i);
	dir(j) = h0(j);
	double f1 = RF->operator()(dpar + dir)(0);
	double f2 = RF->operator()(dpar - dir)(0);
	double hdi = pa[i * npar + i] * h0(i) * h0(i);
	double hdj = pa[j * npar + j] * h0(j) * h0(j);
	pa[j * npar + i] = (f1 + f2 - 2.0 * f0 - hdi - hdj) / ( 2.0 * h0(i) * h0(j) );
	pa[i * npar + j] = pa[j * npar + i];
      }
    }
    UNPROTECT(2);
    return ans;  
  }
}
