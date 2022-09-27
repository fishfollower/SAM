// #define WITH_LIBTMB
// #include <TMB.hpp>
// #include "../inst/include/SAM.hpp"
// #include "TMB.h"

#define WITH_SAM_LIB
#include "SAM.h"



matrix<double> splinebasis_bc(vector<double> x,
			      vector<double> knots){
   matrix<double> sg = spline_helper::getSigAndGam(knots);
   matrix<double> res(x.size(), knots.size());
   res.setZero();
   for(int k = 0; k < knots.size(); ++k){
     for(int i = 0; i < x.size(); ++i){
       double tmp = spline_helper::dkwnorm((double)x(i), (double)knots(k), (double)sg(k,0), (double)sg(k,1), false);
       res(i,k) = tmp * sg(k,0);
     }
   }
   return res;
}


matrix<double> splinebasis_ibc(vector<double> x,
			      vector<double> knots){
   matrix<double> sg = spline_helper::getSigAndGam(knots);
   matrix<double> res(x.size(), knots.size());
   res.setZero();
   for(int k = 0; k < knots.size(); ++k){
     double v0 = spline_helper::pkwnorm((double)knots(0), (double)knots(k), (double)sg(k,0), (double)sg(k,1));
     for(int i = 0; i < x.size(); ++i){
       double tmp = spline_helper::pkwnorm((double)x(i), (double)knots(k), (double)sg(k,0), (double)sg(k,1));
       res(i,k) = tmp - v0;
     }
   }
   return res;
}

matrix<double> splinebasis_iibc(vector<double> x,
			      vector<double> knots){
   matrix<double> sg = spline_helper::getSigAndGam(knots);
   matrix<double> res(x.size(), knots.size());
   res.setZero();
   for(int k = 0; k < knots.size(); ++k){
     double v0 = spline_helper::pkwnorm((double)knots(0), (double)knots(k), (double)sg(k,0), (double)sg(k,1));
     for(int i = 0; i < x.size(); ++i){
       double tmp = spline_helper::ipkwnorm((double)x(i), (double)knots(k), (double)sg(k,0), (double)sg(k,1), 0.0);
       res(i,k) = tmp - v0 * x(i);
     }
   }
   return res;
}

  

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
    double r = ibcispline(Rf_asReal(x),
			  asVector<double>(knots),
			  asVector<double>(pars));
    return asSEXP(r);
  }

  SEXP iibcisplineR(SEXP x, SEXP knots, SEXP pars){
    double r = iibcispline(Rf_asReal(x),
			  asVector<double>(knots),
			  asVector<double>(pars));
    return asSEXP(r);
  }


  SEXP splinebasis_bcR(SEXP x, SEXP knots){
    return asSEXP(splinebasis_bc(asVector<double>(x),
				 asVector<double>(knots)));
  }

  SEXP splinebasis_ibcR(SEXP x, SEXP knots){
    return asSEXP(splinebasis_ibc(asVector<double>(x),
				  asVector<double>(knots)));
  }

  SEXP splinebasis_iibcR(SEXP x, SEXP knots){
    return asSEXP(splinebasis_iibc(asVector<double>(x),
				   asVector<double>(knots)));
  }


}
