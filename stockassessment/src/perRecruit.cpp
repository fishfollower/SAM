#define WITH_SAM_LIB
#include "SAM.h"

// #include "../inst/include/SAM.hpp"


// #include <R.h>
// #include <Rmath.h>
// #include <Rinternals.h>
// #include <R_ext/Rdynload.h>



PERREC_t<double> perRecruit_S(double logFbar, dataSet<double>& dat, confSet& conf, paraSet<double>& par, vector<double>& logSel, vector<int> aveYears, int nYears, int CT, vector<double> logNinit){
  if(nYears < 0)
    Rf_error("nYears must be non-negative.");
  if(aveYears.size() == 0)
    Rf_error("aveYears must be given.");
  if(logSel.size() != conf.keyLogFsta.maxCoeff()+1)
    Rf_error("Wrong size of selectivity vector");
  if(logNinit.size() != conf.maxAge - conf.minAge + 1)
    Rf_error("Wrong size of initial N");

  // Extend arrays
  dataSet<double> newDat = dat;
  int nMYears = dat.noYears;
  // propMat
  bool guessNY = nYears < 1;
  if(guessNY){		// Calculate nYears from M
    array<double> pm = newDat.propMat;
    extendArray(pm, nMYears, 1, aveYears, false);
    double Mlt = pm(pm.dim[0]-1,pm.dim[1]-1); // Last year last age
    double Flt = exp(logFbar + logSel(logSel.size()-1)); // Last age
    double Zlt = Mlt + Flt;
    nYears = std::max(pm.dim[1] - (-60) / Zlt, (double)pm.dim[1] * 20.0);       
  }
  extendArray(newDat.propMat, nMYears, nYears, aveYears, false);
  // stockMeanWeight
  extendArray(newDat.stockMeanWeight, nMYears, nYears, aveYears, false);
  // catchMeanWeight
  extendArray(newDat.catchMeanWeight, nMYears, nYears, aveYears, false);
  // natMor
  extendArray(newDat.natMor, nMYears, nYears, aveYears, false);
  // landFrac
  extendArray(newDat.landFrac, nMYears, nYears, aveYears, false);
  // disMeanWeight
  extendArray(newDat.disMeanWeight, nMYears, nYears, aveYears, false);
  // landMeanWeight
  extendArray(newDat.landMeanWeight, nMYears, nYears, aveYears, false);
  // propF
  extendArray(newDat.propF, nMYears, nYears, aveYears, false);
  // propM
  extendArray(newDat.propM, nMYears, nYears, aveYears, false);
  newDat.noYears = nYears;

  Recruitment<double> recruit = makeRecruitmentFunction(conf, par);

  // Make logF array
  array<double> logF(logSel.size(), nYears);
  logF.setZero();
  for(int i = 0; i < nYears; ++i)
    logF.col(i) = logSel + logFbar;

  // Make logN array - start with one recruit
  int nAge = conf.maxAge - conf.minAge + 1;
  array<double> logN(nAge, nYears);
  logN.setConstant(R_NegInf);
  logN(0,0) = 0.0;
  

  matrix<double> nvar = get_nvar(newDat, conf, par, logN, logF);
  vector<double> fracMixN(conf.fracMixN.size());
  for(int i=0; i<conf.fracMixN.size(); ++i){fracMixN(i)=conf.fracMixN(i);}
  MVMIX_t<double> neg_log_densityN(nvar,fracMixN);
  MortalitySet<double> mort(newDat, conf, par, logF);
  for(int i = 1; i < nYears; ++i){
    vector<double> predN = predNFun(newDat, conf, par, logN, logF, recruit, mort, i);
    logN.col(i) = predN + neg_log_densityN.simulate();
    // Remove recruitment
    logN(0,i) = R_NegInf;    
  }
  
  vector<double> cat(nYears);
  cat.setZero();
  typename referencepointSet<double>::CatchType catchType = static_cast<typename referencepointSet<double>::CatchType>(CT);
  switch(catchType){
  case referencepointSet<double>::totalCatch:
    cat = catchFun(newDat, conf, logN, logF,mort);
    break;
  case referencepointSet<double>::landings:
    cat = landFun(newDat, conf, logN, logF,mort);
    break;
  case referencepointSet<double>::discard:
    cat = disFun(newDat, conf, logN, logF,mort);
    break;
  default:
    Rf_error("Unknown reference point catch type.");
    break;
  }
  double logYPR = log(sum(cat) + SAM_Zero);//

  double logYLTF = log(yearsLostFishing_i(newDat, conf, logF, newDat.natMor.dim(0)-1, conf.minAge, conf.maxAge) + SAM_Zero);
  double logLifeExpectancy = log(temporaryLifeExpectancy_i(newDat, conf, logF, newDat.natMor.dim(0)-1, conf.minAge, 10 * conf.maxAge) + (double)conf.minAge + SAM_Zero);
  vector<double> ssb = ssbFun(newDat, conf, logN, logF, mort);
  double logSPR = log(sum(ssb) + SAM_Zero); //log(sum(ssb)); log(sum(ssb) + (T)exp(-12.0));


  ////////////////////////////////////////////////////////////////////////////////
  // Survival calculations                                                      //
  ////////////////////////////////////////////////////////////////////////////////
  double discYPR = 0.0;
  // Custom recursive calculation inspired by survival.hpp
  double yl_logp = 0.0;
  double yl_q = 0.0;
  double yl = 0.0;
 
  for(int aa = 0; aa < cat.size(); ++aa){
    int j = std::min(std::max(aa-conf.minAge,0), newDat.natMor.cols()-1);		// Cohort age index
    int i = std::min(aa, newDat.natMor.rows()-1);
    double M = newDat.natMor(i,j);
    double F = 0.0;    
    for(int f = 0; f < conf.keyLogFsta.dim(0); ++f)
      if(conf.keyLogFsta(f,j)>(-1) && aa >= conf.minAge)
	F += exp(logF(conf.keyLogFsta(f,j),i));
    double Z = M + F;
    yl += yl_q + exp(yl_logp) * F / Z * (1.0 - 1.0 / Z * (1.0 - exp(-Z)));    
    discYPR += cat(aa) * exp(-yl);
    // discYe += catYe(aa) * exp(-yl);
    yl_q += exp(yl_logp) * F / Z * (1.0 - exp(-Z));
    yl_logp += -Z;
  }
  double logDiscYPR = log(discYPR);
  ////////////////////////////////////////////////////////////////////////////////
  // End survival calculations                                                  //
  ////////////////////////////////////////////////////////////////////////////////


  
  if(conf.stockRecruitmentModelCode == 0){//  ||
    // conf.stockRecruitmentModelCode == 3){
    PERREC_t<double> res = {logFbar, // logFbar
			  logYPR,	   // logYPR
			  logSPR,	   // logSPR
			  R_NaReal,	   // logSe
			  R_NaReal,		 // logRe
			  R_NaReal,// logYe
			  R_NaReal,// dSR0
			  logLifeExpectancy, // logLifeExpectancy
			  logYLTF, // logYearsLost
			  logDiscYPR, // logDiscYPR
			  R_NaReal	 // logDiscYe
    }; 
    return res;
  }

 

  // Simulate biomass
  //  int nYearsEqMax = (guessNY) ? 1000 : nYears;
  array<double> logNeq(nAge, nYears);
  logNeq.setConstant(R_NegInf);
  logNeq.col(0) = logNinit;
 

  // MortalitySet<double> mort(newDat, conf, par, logF);
  if(guessNY){
    // Improve Initial value by starting stochastic simulation at deterministic equilibrium, but avoid a crashed equilibrium.
    vector<double> logNinit2 = logNinit;
    array<double> logFInit(logF.rows(), 2);
    for(int i = 0; i < logFInit.cols(); ++i)
      logF.col(i) = logSel + logFbar;
    // array<double> noF(logF.rows(),1);
    // noF.setConstant(R_NegInf);
    double imp = 100.0;
    while(imp > 1e-6){
      array<double> tmp(nAge,2);
      tmp.setZero();
      tmp.col(0) = logNinit2;
      // Try unfished equilibrium? Or adaptively reduce F if crashing?
      vector<double> lnitmp = predNFun(newDat, conf, par, tmp, logFInit,recruit,mort, 1);
      if(conf.minAge == 0){
	tmp.col(1) = lnitmp;
	lnitmp(0) = predNFun(newDat, conf, par, tmp, logFInit,recruit,mort, 1)(0);
      }
      if(lnitmp(0) < logNinit(0) + log(0.01)){ // Probably heading for crash
	imp = 0.0;
      }else{
	imp = (logNinit2 - lnitmp).abs().maxCoeff();
      }
      logNinit2 = lnitmp;    
    }
    logNeq.col(0) = logNinit2;
    // for(int i = 0; i < logNinit2.size(); ++i) Rprintf("%d: %.03f\n",i,logNinit2(i));
    // Rprintf("\n");
  }
   
  matrix<double> nvareq = get_nvar(newDat, conf, par, logNeq, logF);
  MVMIX_t<double> neg_log_densityNeq(nvareq,fracMixN);
  
  // vector<double> runMean = exp(logNeq0.col(0));
  // runMean.setZero();
  // int nYearsEq = 1;
  for(int i = 1; i < nYears; ++i){
    //Harvest control rule?
    vector<double> predN = predNFun(newDat, conf, par, logNeq, logF,recruit,mort, i);
    vector<double> noiseN = neg_log_densityNeq.simulate();
    logNeq.col(i) = predN + noiseN;
    if(conf.minAge == 0){
      // In this case, predicted recruitment was wrong since it depends on SSB the same year
      // Predict recruitment again with updated ssb (assuming maturity at age 0 is 0):
      double predRec = predNFun(newDat, conf, par, logNeq, logF,recruit,mort, i)(0);
      // Overwrite recruitment, but keep simulated noise to retain correlation:
      logNeq(0,i) = predRec + noiseN(0);      
    }
  }

  vector<double> catYe(nYears);
  catYe.setZero();
  switch(catchType){
  case referencepointSet<double>::totalCatch:
    catYe = catchFun(newDat, conf, logNeq, logF,mort);
    break;
  case referencepointSet<double>::landings:
    catYe = landFun(newDat, conf, logNeq, logF,mort);
    break;
  case referencepointSet<double>::discard:
    catYe = disFun(newDat, conf, logNeq, logF,mort);
    break;
  default:
    Rf_error("Unknown reference point catch type.");
    break;
  }

  vector<double> ssbeq = ssbFun(newDat, conf, logNeq, logF,mort);
  double logSe = log(ssbeq(nYears-1));
  double logYe = log(catYe(nYears-1));
  double logRe = logNeq(0,nYears-1);
  double logDiscYe = logDiscYPR + logRe;

  double dsr0 = recruit.dSR(-30.0);
  
  // Return
  PERREC_t<double> res = {logFbar, // logFbar
			logYPR,	// logYPR
			logSPR,	// logSPR
			logSe,	// logSe
			logRe,	// logRe
			logYe,	// logYe
			dsr0,
			logLifeExpectancy, // logLifeExpectancy
			logYLTF, // logYearsLost
			logDiscYPR,
			logDiscYe,
  };	// DSR0



  return res;
}



struct F_dFunctionalSR2 {
  int srmc;
  double year;
  double lastR;
  vector<double> crb;

  template <template<class> class V, class T>
  T operator()(const V<T> &x){
    confSet conf;
    conf.stockRecruitmentModelCode = srmc;
    conf.constRecBreaks = crb;
    vector<T> rec_pars(x.size()-1);
    for(int i = 0; i < rec_pars.size(); ++i)
      rec_pars(i) = x(i);
    T logssb = x(x.size()-1);
    paraSet<T> par;
    par.rec_pars = rec_pars;
    Recruitment<T> rec = makeRecruitmentFunction(conf, par);
    return rec(logssb, T(year), T(lastR));
  
  }
};

int rec_hasEquilibrium(Recruitment<double>& rec){
    double lse = rec.logSe(0.0);
    return R_finite(lse);
};

int rec_isCompensatory(Recruitment<double>& rec){
  double h = 0.1;
  double x = -5.0;
  bool isComp = !rec.isAutoregressive() && !rec.isTimevarying();
  while(isComp && x < 20){
    if(rec.dSR(x) >= exp(rec(x,R_NaReal,R_NaReal)-x))
      isComp = false;
    x += h;
  }
  return isComp;
};

int rec_hasOvercompensation(Recruitment<double>& rec){
  double h = 0.01;
  double x = -5.0;
  bool isOComp = !rec.isAutoregressive() && !rec.isTimevarying();
  while(!isOComp && x < 20){
    if(rec.dSR(x) < 0)
      isOComp = true;
    x += h;
  }
  return isOComp;
};

int rec_hasMaxAtFiniteS(Recruitment<double>& rec){
  return R_finite(rec.logSAtMaxR());
};

int rec_hasFiniteMax(Recruitment<double>& rec){
  return R_finite(rec.logMaxR());
};

int rec_hasFiniteMaxGradient(Recruitment<double>& rec){
  double g = rec.maxGradient();
  return R_finite(g) && g < 1.0e8;
};

int rec_isAutoregressive(Recruitment<double>& rec){
  return rec.isAutoregressive();
}

int rec_isTimevarying(Recruitment<double>& rec){
  return rec.isTimevarying();
}

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


  SEXP perRecruitR(SEXP logFbar, SEXP tmbdat, SEXP pl, SEXP sel, SEXP aveYears, SEXP nYears, SEXP CT){
    dataSet<double> d0(tmbdat);
    confSet c0(tmbdat);
    paraSet<double> p0(pl);
    vector<double> s0 = asVector<double>(sel);
    vector<double> ls0(s0.size());
    for(int i = 0; i < ls0.size(); ++i)
      ls0 = log(s0);
    vector<int> a0 = asVector<int>(aveYears);
    double logFbar0 = Rf_asReal(logFbar);
    int nY0 = Rf_asInteger(nYears);
    // int RC0 = Rf_asInteger(RC);
    int CT0 = Rf_asInteger(CT);
    PERREC_t<double> y = perRecruit_D<double>(logFbar0, d0, c0, p0, ls0, a0, nY0, CT0);
    const char *resNms[] = {"logF", "logYPR", "logSPR", "logSe", "logRe", "logYe", "dSR0", "logLifeExpectancy", "logYearsLost","logDiscYe","logDiscYPR", ""}; // Must end with ""
    SEXP res;
    PROTECT(res = Rf_mkNamed(VECSXP, resNms));
    SET_VECTOR_ELT(res, 0, asSEXP(y.logFbar));
    SET_VECTOR_ELT(res, 1, asSEXP(y.logYPR));
    SET_VECTOR_ELT(res, 2, asSEXP(y.logSPR));
    SET_VECTOR_ELT(res, 3, asSEXP(y.logSe));
    SET_VECTOR_ELT(res, 4, asSEXP(y.logRe));
    SET_VECTOR_ELT(res, 5, asSEXP(y.logYe));
    SET_VECTOR_ELT(res, 6, asSEXP(y.dSR0));
    SET_VECTOR_ELT(res, 7, asSEXP(y.logLifeExpectancy));
    SET_VECTOR_ELT(res, 8, asSEXP(y.logYearsLost));
    SET_VECTOR_ELT(res, 9, asSEXP(y.logDiscYe));
    SET_VECTOR_ELT(res, 10, asSEXP(y.logDiscYPR));

    UNPROTECT(1);    
    return res;

  }

  SEXP perRecruitSR(SEXP logFbar, SEXP tmbdat, SEXP pl, SEXP sel, SEXP aveYears, SEXP nYears, SEXP CT, SEXP logNinit){
    dataSet<double> d0(tmbdat);
    confSet c0(tmbdat);
    paraSet<double> p0(pl);
    vector<double> s0 = asVector<double>(sel);
    vector<double> ls0(s0.size());
    for(int i = 0; i < ls0.size(); ++i)
      ls0 = log(s0);
    vector<int> a0 = asVector<int>(aveYears);
    double logFbar0 = Rf_asReal(logFbar);
    int nY0 = Rf_asInteger(nYears);
    // int RC0 = Rf_asInteger(RC);
    int CT0 = Rf_asInteger(CT);
    vector<double> logNinit0 = asVector<double>(logNinit);
    GetRNGstate();
    PERREC_t<double> y = perRecruit_S(logFbar0, d0, c0, p0, ls0, a0, nY0, CT0, logNinit0);
    PutRNGstate();
    const char *resNms[] = {"logF", "logYPR", "logSPR", "logSe", "logRe", "logYe", "dSR0", "logLifeExpectancy", "logYearsLost","logDiscYe","logDiscYPR", ""}; // Must end with ""
    SEXP res;
    PROTECT(res = Rf_mkNamed(VECSXP, resNms));
    SET_VECTOR_ELT(res, 0, asSEXP(y.logFbar));
    SET_VECTOR_ELT(res, 1, asSEXP(y.logYPR));
    SET_VECTOR_ELT(res, 2, asSEXP(y.logSPR));
    SET_VECTOR_ELT(res, 3, asSEXP(y.logSe));
    SET_VECTOR_ELT(res, 4, asSEXP(y.logRe));
    SET_VECTOR_ELT(res, 5, asSEXP(y.logYe));
    SET_VECTOR_ELT(res, 6, asSEXP(y.dSR0));
    SET_VECTOR_ELT(res, 7, asSEXP(y.logLifeExpectancy));
    SET_VECTOR_ELT(res, 8, asSEXP(y.logYearsLost));
    SET_VECTOR_ELT(res, 9, asSEXP(y.logDiscYe));
    SET_VECTOR_ELT(res, 10, asSEXP(y.logDiscYPR));
    
    UNPROTECT(1);    
    return res;

  }



  
  SEXP logSRR(SEXP logssb, SEXP rec_pars, SEXP code, SEXP constRecBreaks, SEXP year, SEXP lastR){
    // Fake paraSet and confSet
    vector<double> rp = asVector<double>(rec_pars);
    paraSet<double> par;
    par.rec_pars = rp;
    int srmc = Rf_asInteger(code);
    vector<double> crb = asVector<double>(constRecBreaks);
    confSet conf;
    conf.stockRecruitmentModelCode = srmc;
    conf.constRecBreaks = crb;
    // Make recruitment
    Recruitment<double> rec = makeRecruitmentFunction(conf, par);
    // Calculate
    int n = Rf_length(logssb);
    SEXP v = PROTECT(Rf_allocVector(REALSXP, n));
    double* LS = REAL(logssb);
    double* LR = REAL(v);
    double* y = REAL(year);
    double* lastRec = REAL(lastR);
    for(int i = 0; i < n; ++i)
      LR[i] = rec(LS[i], y[i], lastRec[i]);
    UNPROTECT(1);
    return v;
  }
 

  SEXP stockRecruitmentModelR(SEXP logssb, SEXP rec_pars, SEXP code, SEXP constRecBreaks, SEXP year, SEXP lastR){
    double b = Rf_asReal(logssb);
    double y = Rf_asReal(year);
    double lr = Rf_asReal(lastR);
     // Fake paraSet and confSet
    vector<double> rp = asVector<double>(rec_pars);
    paraSet<double> par;
    par.rec_pars = rp;
    int srmc = Rf_asInteger(code);
    vector<double> crb = asVector<double>(constRecBreaks);
    confSet conf;
    conf.stockRecruitmentModelCode = srmc;
    conf.constRecBreaks = crb;
    // Make recruitment
    Recruitment<double> rec = makeRecruitmentFunction(conf, par);
    // Calculate	
    double v = rec(b, y, lr);
    F_dFunctionalSR2 Fd = {srmc,y,lr,crb};
    vector<double> x(rp.size() + 1);
    x.setZero();
    for(int i = 0; i < rp.size(); ++i)
      x(i) = rp(i);
    x(x.size()-1) = b;
    TMBad::ADFun<> G(TMBad::StdWrap<F_dFunctionalSR2,vector<TMBad::ad_aug> >(Fd), x);
    G = G.JacFun();
    vector<double> r = G(x);
 
    // vector<double> r = x;
    const char *resNms[] = {"logRecruits", "Gradient","dRdS", ""}; // Must end with ""
    SEXP res;
    PROTECT(res = Rf_mkNamed(VECSXP, resNms));
    SET_VECTOR_ELT(res, 0, asSEXP(v));
    SET_VECTOR_ELT(res, 1, asSEXP(r));
    SET_VECTOR_ELT(res, 2, asSEXP(rec.dSR(b)));
    UNPROTECT(1);    
    return res;
      
  }
  
  SEXP recruitmentProperties(SEXP tmbdat, SEXP pl){
    confSet c0(tmbdat);
    paraSet<double> p0(pl);
    Recruitment<double> rec = makeRecruitmentFunction(c0,p0);
    const char *resNms[] = {"name","hasEquilibrium", "isCompensatory", "hasMaxAtFiniteS", "isAutoregressive","isTimevarying","hasOvercompensation","hasFiniteMax","hasFiniteMaxGradient", "logSAtMaxR", "logMaxR", "maxGradient", ""}; // Must end with ""
    SEXP res;
    PROTECT(res = Rf_mkNamed(VECSXP, resNms));
    SET_VECTOR_ELT(res, 0, Rf_mkString(rec.name));
    SET_VECTOR_ELT(res, 1, Rf_ScalarLogical(rec_hasEquilibrium(rec)));
    SET_VECTOR_ELT(res, 2, Rf_ScalarLogical(rec_isCompensatory(rec)));
    SET_VECTOR_ELT(res, 3, Rf_ScalarLogical(rec_hasMaxAtFiniteS(rec)));
    SET_VECTOR_ELT(res, 4, Rf_ScalarLogical(rec_isAutoregressive(rec)));
    SET_VECTOR_ELT(res, 5, Rf_ScalarLogical(rec_isTimevarying(rec)));
    SET_VECTOR_ELT(res, 6, Rf_ScalarLogical(rec_hasOvercompensation(rec))); // Must be same as hasMaxAtFiniteS
    SET_VECTOR_ELT(res, 7, Rf_ScalarLogical(rec_hasFiniteMax(rec)));
    SET_VECTOR_ELT(res, 8, Rf_ScalarLogical(rec_hasFiniteMaxGradient(rec)));
    SET_VECTOR_ELT(res, 9, Rf_ScalarReal(rec.logSAtMaxR()));
    SET_VECTOR_ELT(res, 10, Rf_ScalarReal(rec.logMaxR()));
    SET_VECTOR_ELT(res, 11, Rf_ScalarReal(rec.maxGradient()));
    UNPROTECT(1);
    return res;
  }

}
