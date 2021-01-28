#define WITH_LIBTMB
#include <TMB.hpp>
#include "../inst/include/SAM.hpp"


extern "C" {
  
  SEXP hcrR(SEXP ssb, SEXP hcrConf){
    vector<double> s = asVector<double>(ssb);
    vector<double> hc = asVector<double>(hcrConf);
    vector<double> r(s.size());
    r.setZero();
    for(int i = 0; i < s.size(); ++i)
      r(i) = hcr(s(i), hc);
    return asSEXP(exp(r));
  }


  SEXP perRecruitR(SEXP Fbar, SEXP dat, SEXP conf, SEXP pl, SEXP sel, SEXP aveYears, SEXP nYears){
    dataSet<double> d0(dat);
    confSet c0(conf);
    paraSet<double> p0(pl);
    vector<double> s0 = asVector<double>(sel);
    vector<double> ls0(s0.size());
    for(int i = 0; i < ls0.size(); ++i)
      ls0 = log(s0);
    vector<int> a0 = asVector<int>(aveYears);
    double Fbar0 = Rf_asReal(Fbar);
    int nY0 = Rf_asInteger(nYears);
 
    PERREC_t<double> y = perRecruit<double, double>(Fbar0, d0, c0, p0, s0, a0, nY0);
    const char *resNms[] = {"logF", "logYPR", "logSPR", "logSe", "logRe", "logYe", ""}; // Must end with ""
    SEXP res;
    PROTECT(res = Rf_mkNamed(VECSXP, resNms));
    SET_VECTOR_ELT(res, 0, asSEXP(y.logFbar));
    SET_VECTOR_ELT(res, 1, asSEXP(y.logYPR));
    SET_VECTOR_ELT(res, 2, asSEXP(y.logSPR));
    SET_VECTOR_ELT(res, 3, asSEXP(y.logSe));
    SET_VECTOR_ELT(res, 4, asSEXP(y.logRe));
    SET_VECTOR_ELT(res, 5, asSEXP(y.logYe));

    UNPROTECT(1);    
    return res;

  }


  SEXP stockRecruitmentModelR(SEXP ssb, SEXP rec_pars, SEXP code){
    double b = Rf_asReal(ssb);
    vector<double> rp = asVector<double>(rec_pars);
    int srmc = Rf_asInteger(code);
	
    double v = exp(functionalStockRecruitment(b, rp, srmc));

#ifdef CPPAD_FRAMEWORK
    vector<AD<double> > rp2(rp.size() + 1);
    for(int i = 0; i < rp.size(); ++i)
      rp2(i) = rp(i);
    rp2(rp.size()) = b;
    CppAD::Independent(rp2);
    // vector<AD<double> > x( 1 );
    // x[0] = b;
    // CppAD::Independent(x);
    vector<AD<double> > y( 1 );
    y[0] = exp(functionalStockRecruitment(rp2(rp.size()), (vector<AD<double> >)rp2.head(rp.size()), srmc));
    CppAD::ADFun<double> F(rp2, y);
    vector<double> x_eval( rp.size() + 1 );
    for(int i = 0; i < rp.size(); ++i)
      x_eval(i) = rp(i);
    x_eval[rp.size()] = b;
    vector<double> r = F.Jacobian(x_eval);
#endif
#ifdef TMBAD_FRAMEWORK
   
    F_dFunctionalSR2 Fd = {rp.size(),srmc};
    vector<double> x(rp.size() + 1);
    for(int i = 0; i < rp.size(); ++i)
      x(i) = rp(i);
    x(rp.size()) = b;
    TMBad::ADFun<> G(TMBad::StdWrap<F_dFunctionalSR2,vector<TMBad::ad_aug> >(Fd), x);
    // TMBad::ADFun<> G(Fd,x);
    G = G.JacFun();
    vector<double> r = G(x);
#endif

    const char *resNms[] = {"Recruits", "Gradient", ""}; // Must end with ""
    SEXP res;
    PROTECT(res = Rf_mkNamed(VECSXP, resNms));
    SET_VECTOR_ELT(res, 0, asSEXP(v));
    SET_VECTOR_ELT(res, 1, asSEXP(r));
 
    UNPROTECT(1);    
    return res;
      
  }

  SEXP Se_sbhR(SEXP lambda, SEXP a, SEXP b, SEXP g){
    double r = Se_sbh(Rf_asReal(lambda), Rf_asReal(a), Rf_asReal(b), Rf_asReal(g));
    return asSEXP(r);
  }

  SEXP Se_slR(SEXP lambda, SEXP a, SEXP b, SEXP g){
    double r = Se_sl(Rf_asReal(lambda), Rf_asReal(a), Rf_asReal(b), Rf_asReal(g));
    return asSEXP(r);
  }

}
