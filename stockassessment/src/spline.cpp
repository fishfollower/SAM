#define WITH_LIBTMB
#include <TMB.hpp>
#include "../inst/include/SAM.hpp"

  
extern "C" {


  SEXP bcsplineR(SEXP x, SEXP knots, SEXP pars){
    double r = bcspline(Rf_asReal(x),
			asVector<double>(knots),
			asVector<double>(pars));
    return asSEXP(r);
  }
  SEXP ibcsplineR(SEXP x, SEXP knots, SEXP pars){
    double r = ibcspline(Rf_asReal(x),
			 asVector<double>(knots),
			 asVector<double>(pars));
    return asSEXP(r);
  }
  SEXP ibcdsplineR(SEXP x, SEXP knots, SEXP pars){
    double r = ibcdspline(Rf_asReal(x),
			  asVector<double>(knots),
			  asVector<double>(pars));
    return asSEXP(r);
  }
  SEXP ibcisplineR(SEXP x, SEXP knots, SEXP pars){
    double r = ibcdspline(Rf_asReal(x),
			  asVector<double>(knots),
			  asVector<double>(pars));
    return asSEXP(r);
  }

  SEXP iibcisplineR(SEXP x, SEXP knots, SEXP pars){
    double r = ibcdspline(Rf_asReal(x),
			  asVector<double>(knots),
			  asVector<double>(pars));
    return asSEXP(r);
  }

}
