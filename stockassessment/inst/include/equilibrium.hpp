#pragma once
#ifndef SAM_EQUILIBRIUM_HPP
#define SAM_EQUILIBRIUM_HPP


template<class Type>
struct PERREC_t {
  Type logFbar;
  Type logYPR;
  Type logSPR;
  Type logSe;
  Type logRe;
  Type logYe;
  Type dSR0;
  Type logLifeExpectancy;
  Type logYearsLost;
  Type logDiscYPR;
  Type logDiscYe;  
};


///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Deterministic /////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

template<class Type>
PERREC_t<Type> perRecruit_D(const Type& logFbar, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, vector<Type>& logSel, vector<int>& aveYears, int nYears = 300, int CT = 0){
  if(nYears <= 0)
    Rf_error("nYears must be greater than 0.");
  if(aveYears.size() == 0)
    Rf_error("aveYears must be given.");
  if(logSel.size() != conf.keyLogFsta.maxCoeff()+1)
    Rf_error("Wrong size of selectivity vector");

  // Prepare data
  dataSet<Type> newDat = dat;
  int nMYears = dat.noYears;
  // propMat
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

  Recruitment<Type> recruit = makeRecruitmentFunction(conf, par);
  Recruitment<Type> rec0 = Recruitment<Type>("zero",new Rec_None<Type>());

  // Make logF array
  array<Type> logF(logSel.size(), nYears);
  logF.setZero();
  for(int i = 0; i < nYears; ++i)
    logF.col(i) = logSel + logFbar;

  // Make logN array - start with one recruit
  int nAge = conf.maxAge - conf.minAge + 1;
  array<Type> logN(nAge, nYears);
  logN.setConstant(R_NegInf);
  logN(0,0) = 0.0;

  MortalitySet<Type> mort(newDat, conf, par, logF);
  
  // Run loop over years
  for(int i = 1; i < nYears; ++i){
    logN.col(i) = predNFun(newDat, conf, par, logN, logF, rec0, mort, i);
    //logN(0,i) = R_NegInf;
  }
 
  // Calculate yield
  vector<Type> cat(nYears);
  cat.setZero();
  typename referencepointSet<Type>::CatchType catchType = static_cast<typename referencepointSet<Type>::CatchType>(CT);
  switch(catchType){
  case referencepointSet<Type>::totalCatch:
    cat = catchFun(newDat, conf, logN, logF, mort);
    break;
  case referencepointSet<Type>::landings:
    cat = landFun(newDat, conf, logN, logF, mort);
    break;
  case referencepointSet<Type>::discard:
    cat = disFun(newDat, conf, logN, logF, mort);
    break;
  default:
    Rf_error("Unknown reference point catch type.");
    break;
  }
  Type logYPR = log(sum(cat) + SAM_Zero);//

  Type logYLTF = log(yearsLostFishing_i(newDat, conf, logF, newDat.natMor.dim(0)-1, conf.minAge, conf.maxAge) + SAM_Zero);
  Type logLifeExpectancy = log(temporaryLifeExpectancy_i(newDat, conf, logF, newDat.natMor.dim(0)-1, conf.minAge, 10 * conf.maxAge) + (Type)conf.minAge + SAM_Zero);

  // Calculate spawners
  vector<Type> ssb = ssbFun(newDat, conf, logN, logF, mort);
  Type logSPR = log(sum(ssb) + SAM_Zero);

  ////////////////////////////////////////////////////////////////////////////////
  // Survival calculations                                                      //
  ////////////////////////////////////////////////////////////////////////////////
  Type discYPR = 0.0;
  // T discYe = 0.0;
  // Custom recursive calculation inspired by survival.hpp
  Type yl_logp = 0.0;
  Type yl_q = 0.0;
  Type yl = 0.0;
  for(int aa = 0; aa < cat.size(); ++aa){
    int j = std::min(std::max(aa-conf.minAge,0), newDat.natMor.cols()-1);		// Cohort age index
    int i = std::min(aa, newDat.natMor.rows()-1);
    Type M = newDat.natMor(i,j);
    Type F = 0.0;
    for(int f = 0; f < conf.keyLogFsta.dim(0); ++f)
      if(conf.keyLogFsta(f,j)>(-1) && aa >= conf.minAge)
	F += exp(logF(conf.keyLogFsta(f,j),i));
   Type Z = M + F;
    yl += yl_q + exp(yl_logp) * F / Z * (1.0 - 1.0 / Z * (1.0 - exp(-Z)));
    discYPR += cat(aa) * exp(-yl);
    // discYe += catYe(aa) * exp(-yl);
    yl_q += exp(yl_logp) * F / Z * (1.0 - exp(-Z));
    yl_logp += -Z;
  }
  Type logDiscYPR = log(discYPR + SAM_Zero);
  // T logDiscYe = log(discYe + (T)exp(SAM_NegInf));
  ////////////////////////////////////////////////////////////////////////////////
  // End survival calculations                                                  //
  ////////////////////////////////////////////////////////////////////////////////

  Type logSe = recruit.logSe(logSPR);
  Type dSR0 = recruit.dSR((Type)SAM_NegInf);
  Type logRe = logSe - logSPR;
  Type logYe = logSe - logSPR + logYPR;
  Type logDiscYe = logDiscYPR + logRe;
 
   // Return
  PERREC_t<Type> res = {logFbar, // logFbar
		     logYPR,	// logYPR
		     logSPR,	// logSPR
		     logSe,	// logSe
		     logRe,	// logRe
		     logYe,	// logYe
		     dSR0,      // DSR0
		     logLifeExpectancy, // logLifeExpectancy
		     logYLTF, // logYearsLost
		     logDiscYPR,
		     logDiscYe,
  };

  return res;
   
}


////////// Convenience functions to get Yield per recruit and derivative //////////

template<class Type>
struct Funct_YPR {
  dataSet<Type> dat;
  confSet conf;
  paraSet<Type> par;
  referencepointSet<Type> rp;
  
  Funct_YPR() = default;
  
  Funct_YPR(const dataSet<Type>& dat_,
  	    const confSet& conf_,
  	    const paraSet<Type>& par_,
  	    const referencepointSet<Type>& rp_) : 
    dat(dat_),
    conf(conf_),
    par(par_),
    rp(rp_){};

  // template<class T>
  // T operator()(vector<T> logFbar){
  Type operator()(vector<Type> logFbar){
    // dataSet<T> d2(dat);
    // paraSet<T> p2(par);
    vector<Type> ls = rp.getLogSelectivity();
    // vector<T> ls(ls0);
    // PERREC_t<T> r = perRecruit_D(logFbar(0), d2, conf, p2, ls, rp.aveYears, rp.nYears, rp.catchType);
    PERREC_t<Type> r = perRecruit_D(logFbar(0), dat, conf, par, ls, rp.aveYears, rp.nYears, rp.catchType);
    return exp(r.logYPR);
  }  
};

template<class Type>
Type dYPR(Type logFbar, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, referencepointSet<Type> rp){
  Funct_YPR<Type> f(dat,conf,par,rp);
  vector<Type> u(1); u(0) = logFbar;
  //vector<Type> g = autodiff::gradient(f, u);
  // autodiff::gradient gives memory not mapped error
  // Using numeric gradient instead
  Type h = 0.001;
  Type v = -f((vector<Type>)(u + 2.0 * h)) + 8.0 * f((vector<Type>)(u + h)) - 8.0 * f((vector<Type>)(u - h)) + f((vector<Type>)(u - 2.0 * h));
  Type g = v / (12.0 * h);
  // Return diff(YPR(f))|_f=exp(logFbar)
  return g / exp(logFbar);
}

// For calculations about future
template<class Type>
Type yieldPerRecruit_i(const Type& logFbar, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, referencepointSet<Type>& rp, bool give_log = false){
  vector<Type> ls = rp.getLogSelectivity();
  PERREC_t<Type> r =  perRecruit_D(logFbar, dat, conf, par, ls, rp.aveYears, rp.nYears, rp.catchType);
  if(give_log)
    return r.logYPR;
  return exp(r.logYPR);
}

// For calculations about assessment period for one year
template<class Type>
Type yieldPerRecruit_i(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, int i, int CT, int nYears, bool give_log = false){
  // Make referencepointSet
  referencepointSet<Type> rp(nYears, CT, i, logF, conf);
  //rp.setLogSelectivity(logF,conf);
  Type logFbar = rp.logFbar(logF, conf);  
  return yieldPerRecruit_i(logFbar, dat, conf, par, rp, give_log);
}

template<class Type>
vector<Type> yieldPerRecruit(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, bool give_log = false){
  vector<Type> r(dat.catchMeanWeight.dim(0));
  r.setZero();
  for(int i = 0; i < r.size(); ++i){
    r(i) = yieldPerRecruit_i(dat, conf, par, logF, i, 0, 300, give_log);
  }
  return r;
}

////////// Convenience functions to get Spawners per recruit and derivative //////////



// For calculations about future
template<class Type>
Type spawnersPerRecruit_i(const Type& logFbar, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, referencepointSet<Type>& rp, bool give_log = false){
  vector<Type> ls = rp.getLogSelectivity();
  PERREC_t<Type> r =  perRecruit_D(logFbar, dat, conf, par, ls, rp.aveYears, rp.nYears, rp.catchType);
  if(give_log)
    return r.logSPR;
  return exp(r.logSPR);
}

// For calculations about assessment period for one year
template<class Type>
Type spawnersPerRecruit_i(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, int i, int CT, int nYears, bool give_log = false){
  // Make referencepointSet
  referencepointSet<Type> rp(nYears, CT, i, logF, conf);
  Type logFbar = rp.logFbar(logF, conf);
  return spawnersPerRecruit_i(logFbar, dat, conf, par, rp, give_log);
}

template<class Type>
vector<Type> spawnersPerRecruit(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, bool give_log = false){
  vector<Type> r(dat.catchMeanWeight.dim(0));
  r.setZero();
  for(int i = 0; i < r.size(); ++i)
    r(i) = spawnersPerRecruit_i(dat, conf, par, logF, i, 0, 300, give_log);
  return r;
}

// For calculations about assessment period for one year
template<class Type>
Type SPR0_i(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, int i, int CT, int nYears, bool give_log = false){
  // Make referencepointSet
  referencepointSet<Type> rp(nYears, CT, i, logF, conf);
  Type logFbar = R_NegInf; //rp.logFbar(logF, conf);
  return spawnersPerRecruit_i(logFbar, dat, conf, par, rp, give_log);
}


////////// Convenience functions to get equilibrium biomass //////////



// For calculations about future
template<class Type>
Type equilibriumBiomass_i(const Type& logFbar, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, referencepointSet<Type>& rp, bool give_log = false){
  vector<Type> ls = rp.getLogSelectivity();
  PERREC_t<Type> r =  perRecruit_D(logFbar, dat, conf, par, ls, rp.aveYears, rp.nYears, rp.catchType);
  if(give_log)
    return r.logSe;
  return exp(r.logSe);
}

// For calculations about assessment period for one year
template<class Type>
Type equilibriumBiomass_i(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, int i, int CT, int nYears, bool give_log = false){
  // Make referencepointSet
  referencepointSet<Type> rp(nYears, CT, i, logF, conf);
  Type logFbar = rp.logFbar(logF, conf);
  return equilibriumBiomass_i(logFbar, dat, conf, par, rp, give_log);
}

template<class Type>
vector<Type> equilibriumBiomass(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, bool give_log = false){
  vector<Type> r(dat.catchMeanWeight.dim(0));
  r.setZero();
  for(int i = 0; i < r.size(); ++i)
    r(i) = equilibriumBiomass_i(dat, conf, par, logF, i, 0, 300, give_log);
  return r;
}


////////// Convenience functions to get equilibrium biomass //////////



// For calculations about future
template<class Type>
Type B0_i(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, referencepointSet<Type>& rp, bool give_log = false){
  vector<Type> ls = rp.getLogSelectivity();
  PERREC_t<Type> r =  perRecruit_D(Type(R_NegInf), dat, conf, par, ls, rp.aveYears, rp.nYears, rp.catchType);
  if(give_log)
    return r.logSe;
  return exp(r.logSe);
}

// For calculations about assessment period for one year
template<class Type>
Type B0_i(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, int i, int CT, int nYears, bool give_log = false){
  // Make referencepointSet
  referencepointSet<Type> rp(nYears, CT, i, logF, conf);
  return B0_i(dat, conf, par, rp, give_log);
}

template<class Type>
vector<Type> B0(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, bool give_log = false){
  vector<Type> r(dat.catchMeanWeight.dim(0));
  r.setZero();
  for(int i = 0; i < r.size(); ++i)
    r(i) = B0_i(dat, conf, par, logF, i, 0, 300, give_log);
  return r;
}

  


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Stochastic //////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////



template<class Type>
PERREC_t<Type> perRecruit_S(Type logFbar, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, vector<Type>& logSel, vector<int> aveYears, int nYears, int CT, vector<Type> logNinit){
  if(nYears < 0)
    Rf_error("nYears must be non-negative.");
  if(aveYears.size() == 0)
    Rf_error("aveYears must be given.");
  if(logSel.size() != conf.keyLogFsta.maxCoeff()+1)
    Rf_error("Wrong size of selectivity vector");
  if(logNinit.size() != conf.maxAge - conf.minAge + 1)
    Rf_error("Wrong size of initial N");

  // Extend arrays
  dataSet<Type> newDat = dat;
  int nMYears = dat.noYears;
  // propMat
  bool guessNY = nYears < 1;
  if(guessNY){		// Calculate nYears from M
    array<Type> pm = newDat.propMat;
    extendArray(pm, nMYears, 1, aveYears, false);
    Type Mlt = pm(pm.dim[0]-1,pm.dim[1]-1); // Last year last age
    Type Flt = exp(logFbar + logSel(logSel.size()-1)); // Last age
    Type Zlt = Mlt + Flt;
    nYears = std::max(pm.dim[1] - (-60) / Zlt, (Type)pm.dim[1] * 20.0);       
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

  Recruitment<Type> recruit = makeRecruitmentFunction(conf, par);

  // Make logF array
  array<Type> logF(logSel.size(), nYears);
  logF.setZero();
  for(int i = 0; i < nYears; ++i)
    logF.col(i) = logSel + logFbar;

  // Make logN array - start with one recruit
  int nAge = conf.maxAge - conf.minAge + 1;
  array<Type> logN(nAge, nYears);
  logN.setConstant(R_NegInf);
  logN(0,0) = 0.0;
  

  matrix<Type> nvar = get_nvar(newDat, conf, par, logN, logF);
  vector<Type> fracMixN(conf.fracMixN.size());
  for(int i=0; i<conf.fracMixN.size(); ++i){fracMixN(i)=conf.fracMixN(i);}
  MVMIX_t<Type> neg_log_densityN(nvar,fracMixN);
  MortalitySet<Type> mort(newDat, conf, par, logF);
  for(int i = 1; i < nYears; ++i){
    vector<Type> predN = predNFun(newDat, conf, par, logN, logF, recruit, mort, i);
    logN.col(i) = predN + neg_log_densityN.simulate();
    // Remove recruitment
    logN(0,i) = R_NegInf;    
  }
  
  vector<Type> cat(nYears);
  cat.setZero();
  typename referencepointSet<Type>::CatchType catchType = static_cast<typename referencepointSet<Type>::CatchType>(CT);
  switch(catchType){
  case referencepointSet<Type>::totalCatch:
    cat = catchFun(newDat, conf, logN, logF,mort);
    break;
  case referencepointSet<Type>::landings:
    cat = landFun(newDat, conf, logN, logF,mort);
    break;
  case referencepointSet<Type>::discard:
    cat = disFun(newDat, conf, logN, logF,mort);
    break;
  default:
    Rf_error("Unknown reference point catch type.");
    break;
  }
  Type logYPR = log(sum(cat) + SAM_Zero);//

  Type logYLTF = log(yearsLostFishing_i(newDat, conf, logF, newDat.natMor.dim(0)-1, conf.minAge, conf.maxAge) + SAM_Zero);
  Type logLifeExpectancy = log(temporaryLifeExpectancy_i(newDat, conf, logF, newDat.natMor.dim(0)-1, conf.minAge, 10 * conf.maxAge) + (Type)conf.minAge + SAM_Zero);
  vector<Type> ssb = ssbFun(newDat, conf, logN, logF, mort);
  Type logSPR = log(sum(ssb) + SAM_Zero); //log(sum(ssb)); log(sum(ssb) + (T)exp(-12.0));


  ////////////////////////////////////////////////////////////////////////////////
  // Survival calculations                                                      //
  ////////////////////////////////////////////////////////////////////////////////
  Type discYPR = 0.0;
  // Custom recursive calculation inspired by survival.hpp
  Type yl_logp = 0.0;
  Type yl_q = 0.0;
  Type yl = 0.0;
 
  for(int aa = 0; aa < cat.size(); ++aa){
    int j = std::min(std::max(aa-conf.minAge,0), newDat.natMor.cols()-1);		// Cohort age index
    int i = std::min(aa, newDat.natMor.rows()-1);
    Type M = newDat.natMor(i,j);
    Type F = 0.0;    
    for(int f = 0; f < conf.keyLogFsta.dim(0); ++f)
      if(conf.keyLogFsta(f,j)>(-1) && aa >= conf.minAge)
	F += exp(logF(conf.keyLogFsta(f,j),i));
    Type Z = M + F;
    yl += yl_q + exp(yl_logp) * F / Z * (1.0 - 1.0 / Z * (1.0 - exp(-Z)));    
    discYPR += cat(aa) * exp(-yl);
    // discYe += catYe(aa) * exp(-yl);
    yl_q += exp(yl_logp) * F / Z * (1.0 - exp(-Z));
    yl_logp += -Z;
  }
  Type logDiscYPR = log(discYPR);
  ////////////////////////////////////////////////////////////////////////////////
  // End survival calculations                                                  //
  ////////////////////////////////////////////////////////////////////////////////


  
  if(conf.stockRecruitmentModelCode == 0){//  ||
    // conf.stockRecruitmentModelCode == 3){
    PERREC_t<Type> res = {logFbar, // logFbar
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
  array<Type> logNeq(nAge, nYears);
  logNeq.setConstant(R_NegInf);
  logNeq.col(0) = logNinit;
 
  // array<Type> logFeq0(logSel.size(), nYearsEqMax);
  // logFeq0.setZero();
  // for(int i = 0; i < nYearsEqMax; ++i)
  //   logFeq0.col(i) = logSel + logFbar;

  // MortalitySet<Type> mort(newDat, conf, par, logF);
  if(guessNY){
    // Improve Initial value by starting stochastic simulation at deterministic equilibrium, but avoid a crashed equilibrium.
    vector<Type> logNinit2 = logNinit;
    array<Type> logFInit(logF.rows(), 2);
    for(int i = 0; i < logFInit.cols(); ++i)
      logF.col(i) = logSel + logFbar;
    // array<Type> noF(logF.rows(),1);
    // noF.setConstant(R_NegInf);
    Type imp = 100.0;
    while(imp > 1e-6){
      array<Type> tmp(nAge,2);
      tmp.setZero();
      tmp.col(0) = logNinit2;
      // Try unfished equilibrium? Or adaptively reduce F if crashing?
      vector<Type> lnitmp = predNFun(newDat, conf, par, tmp, logFInit,recruit,mort, 1);
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
   
  matrix<Type> nvareq = get_nvar(newDat, conf, par, logNeq, logF);
  MVMIX_t<Type> neg_log_densityNeq(nvareq,fracMixN);
  
  // vector<Type> runMean = exp(logNeq0.col(0));
  // runMean.setZero();
  // int nYearsEq = 1;
  for(int i = 1; i < nYears; ++i){
    //Harvest control rule?
    vector<Type> predN = predNFun(newDat, conf, par, logNeq, logF,recruit,mort, i);
    vector<Type> noiseN = neg_log_densityNeq.simulate();
    logNeq.col(i) = predN + noiseN;
    if(conf.minAge == 0){
      // In this case, predicted recruitment was wrong since it depends on SSB the same year
      // Predict recruitment again with updated ssb (assuming maturity at age 0 is 0):
      Type predRec = predNFun(newDat, conf, par, logNeq, logF,recruit,mort, i)(0);
      // Overwrite recruitment, but keep simulated noise to retain correlation:
      logNeq(0,i) = predRec + noiseN(0);      
    }
  }

  vector<Type> catYe(nYears);
  catYe.setZero();
  switch(catchType){
  case referencepointSet<Type>::totalCatch:
    catYe = catchFun(newDat, conf, logNeq, logF,mort);
    break;
  case referencepointSet<Type>::landings:
    catYe = landFun(newDat, conf, logNeq, logF,mort);
    break;
  case referencepointSet<Type>::discard:
    catYe = disFun(newDat, conf, logNeq, logF,mort);
    break;
  default:
    Rf_error("Unknown reference point catch type.");
    break;
  }

  vector<Type> ssbeq = ssbFun(newDat, conf, logNeq, logF,mort);
  Type logSe = log(ssbeq(nYears-1));
  Type logYe = log(catYe(nYears-1));
  Type logRe = logNeq(0,nYears-1);
  Type logDiscYe = logDiscYPR + logRe;


  // Type dsr0 = 10000.0;
  // if(conf.stockRecruitmentModelCode != 0 &&
  //    conf.stockRecruitmentModelCode != 3 &&
  //    conf.stockRecruitmentModelCode != 62 &&
  //    conf.stockRecruitmentModelCode != 65 &&
  //    (conf.stockRecruitmentModelCode != 68) && // || newPar.rec_pars[2] < 0) &&
  //    (conf.stockRecruitmentModelCode != 69) &&
  //    conf.stockRecruitmentModelCode != 90 &&
  //    conf.stockRecruitmentModelCode != 91 &&
  //    conf.stockRecruitmentModelCode != 92){ // || newPar.rec_pars[2] < 0 ))
  //   dsr0 = dFunctionalSR(Type(SAM_Zero), par.rec_pars, conf.stockRecruitmentModelCode);
  // }
  // if(conf.stockRecruitmentModelCode == 68 || conf.stockRecruitmentModelCode == 69){
  //   dsr0 = CppAD::CondExpLt(par.rec_pars[2],
  // 			    (Type)0.0,
  // 			    dFunctionalSR(Type(SAM_Zero),
  // 					  par.rec_pars,
  // 					  conf.stockRecruitmentModelCode),
  // 			    dsr0);
  // }

  // if(conf.stockRecruitmentModelCode == 90 ||
  //    conf.stockRecruitmentModelCode == 91 ||
  //    conf.stockRecruitmentModelCode == 92){
  //   // dSplineSR uses logssb
  //   dsr0 = dSplineSR((Type)SAM_Zero,
  // 		     (vector<Type>)conf.constRecBreaks.template cast<Type>(),
  // 		     par.rec_pars,
  // 		     conf.stockRecruitmentModelCode);
  // }
  Type dsr0 = recruit.dSR(-30.0);
  
  // Return
  PERREC_t<Type> res = {logFbar, // logFbar
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


  // PERREC_t<Type> res = {R_NaReal, // logFbar
  // 		       R_NaReal,	   // logYPR
  // 		       R_NaReal,	   // logSPR
  // 		       R_NaReal,	   // logSe
  // 		       R_NaReal,		 // logRe
  // 		       R_NaReal,// logYe
  // 		       R_NaReal,// dSR0
  // 		       R_NaReal, // logLifeExpectancy
  // 		       R_NaReal, // logYearsLost
  // 		       R_NaReal, // logDiscYPR
  // 		       R_NaReal	 // logDiscYe
  // }; 

  return res;
}






// Spawners per recruit

// Yield per recruit

// Spawning biomass

// Recruitment

// Yield

// B0


#endif
