SAM_DEPENDS(define)
SAM_DEPENDS(incidence)
SAM_DEPENDS(recruitment)
SAM_DEPENDS(forecast)
SAM_DEPENDS(refpointset)
SAM_DEPENDS(derived)
SAM_DEPENDS(predn)
SAM_DEPENDS(survival)

HEADER(
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
});

SAM_SPECIALIZATION(struct PERREC_t<double>);
SAM_SPECIALIZATION(struct PERREC_t<TMBad::ad_aug>);

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Deterministic /////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

template<class Type>
PERREC_t<Type> perRecruit_D(const Type& logFbar, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, vector<Type>& logSel, vector<int>& aveYears, int nYears DEFARG(= 300), int CT DEFARG(= 0))SOURCE({
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
  extendArray(newDat.propMat, nMYears, nYears, aveYears, par.meanLogitMO, conf.keyMatureMean, 1, false);
   // stockMeanWeight
  extendArray(newDat.stockMeanWeight, nMYears, nYears, aveYears, par.meanLogSW, conf.keyStockWeightMean, 0, false);
  // catchMeanWeight
  extendArray(newDat.catchMeanWeight, nMYears, nYears, aveYears, par.meanLogCW, conf.keyCatchWeightMean, 0, false);
  // natMor
  extendArray(newDat.natMor, nMYears, nYears, aveYears, par.meanLogNM, conf.keyMortalityMean, 0, false);
  // landFrac (No biopar process)
  extendArray(newDat.landFrac, nMYears, nYears, aveYears, false);
  // disMeanWeight (No biopar process)
  extendArray(newDat.disMeanWeight, nMYears, nYears, aveYears, false);
  // landMeanWeight (No biopar process)
  extendArray(newDat.landMeanWeight, nMYears, nYears, aveYears, false);
  // propF (No biopar process)
  extendArray(newDat.propF, nMYears, nYears, aveYears, false);
  // propM (No biopar process)
  extendArray(newDat.propM, nMYears, nYears, aveYears, false);
  newDat.noYears = nYears;

  Recruitment<Type> recruit = makeRecruitmentFunction(conf, par);
  confSet conf2(conf); conf2.stockRecruitmentModelCode = -1;
  Recruitment<Type> rec0 = makeRecruitmentFunction(conf2, par);

  // Make logF array
  array<Type> logF(logSel.size(), nYears);
  logF.setZero();
  for(int i = 0; i < nYears; ++i)
    logF.col(i) = logSel + logFbar;

  // Make logitF season array
  array<Type> logitFseason(par.seasonMu.rows(), par.seasonMu.cols(),nYears);
  logitFseason.setZero();
  for(int i = 0; i < nYears; ++i)
    for(int j = 0; j < par.seasonMu.cols(); ++j)
      for(int k = 0; k < par.seasonMu.rows(); ++k)
	logitFseason(k,j,i) = par.seasonMu(k,j);
      

  // Make logN array - start with one recruit
  int nAge = conf.maxAge - conf.minAge + 1;
  array<Type> logN(nAge, nYears);
  logN.setConstant(R_NegInf);
  logN(0,0) = 0.0;

  MortalitySet<Type> mort(newDat, conf, par, logF, logitFseason);
  
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
   
  });

SAM_SPECIALIZATION(PERREC_t<double> perRecruit_D(const double&, dataSet<double>&, confSet&, paraSet<double>&, vector<double>&, vector<int>&, int, int));
SAM_SPECIALIZATION(PERREC_t<TMBad::ad_aug> perRecruit_D(const TMBad::ad_aug&, dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, vector<TMBad::ad_aug>&, vector<int>&, int, int));

////////// Convenience functions to get Yield per recruit and derivative //////////

#ifndef WITH_SAM_LIB

namespace equilibrium_fun {
template<class Type>
struct Funct_YPR {
  dataSet<Type> dat;
  confSet conf;
  paraSet<Type> par;
  referencepointSet<Type> rp;
  
  Funct_YPR() : dat(), conf(), par(), rp() {};
  
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
}
#endif

template<class Type>
Type dYPR(Type logFbar, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, referencepointSet<Type>& rp)SOURCE({
    equilibrium_fun::Funct_YPR<Type> f(dat,conf,par,rp);
  vector<Type> u(1); u(0) = logFbar;
  //vector<Type> g = autodiff::gradient(f, u);
  // autodiff::gradient gives memory not mapped error
  // Using numeric gradient instead
  Type h = 0.001;
  Type v = -f((vector<Type>)(u + 2.0 * h)) + 8.0 * f((vector<Type>)(u + h)) - 8.0 * f((vector<Type>)(u - h)) + f((vector<Type>)(u - 2.0 * h));
  Type g = v / (12.0 * h);
  // Return diff(YPR(f))|_f=exp(logFbar)
  return g / exp(logFbar);
  });

SAM_SPECIALIZATION(double dYPR(double, dataSet<double>&, confSet&, paraSet<double>&, referencepointSet<double>&));
SAM_SPECIALIZATION(TMBad::ad_aug dYPR(TMBad::ad_aug, dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, referencepointSet<TMBad::ad_aug>&));

// #ifndef WITH_SAM_LIB
// namespace equilibrium_fun {
// For calculations about future
template<class Type>
Type yieldPerRecruit_i(const Type& logFbar, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, referencepointSet<Type>& rp, bool give_log DEFARG(= false))SOURCE({
  vector<Type> ls = rp.getLogSelectivity();
  PERREC_t<Type> r =  perRecruit_D(logFbar, dat, conf, par, ls, rp.aveYears, rp.nYears, rp.catchType);
  if(give_log)
    return r.logYPR;
  return exp(r.logYPR);
  });

SAM_SPECIALIZATION(double yieldPerRecruit_i(const double&, dataSet<double>&, confSet&, paraSet<double>&, referencepointSet<double>&, bool));
SAM_SPECIALIZATION(TMBad::ad_aug yieldPerRecruit_i(const TMBad::ad_aug&, dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, referencepointSet<TMBad::ad_aug>&, bool));


// For calculations about assessment period for one year
template<class Type>
Type yieldPerRecruit_i(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, int i, int CT, int nYears, bool give_log DEFARG(= false))SOURCE({
  // Make referencepointSet
  referencepointSet<Type> rp(nYears, CT, i, logF, conf);
  //rp.setLogSelectivity(logF,conf);
  Type logFbar = rp.logFbar(logF, conf);  
  return yieldPerRecruit_i(logFbar, dat, conf, par, rp, give_log);
  })

  
SAM_SPECIALIZATION(double yieldPerRecruit_i(dataSet<double>&, confSet&, paraSet<double>&, array<double>&, int, int, int, bool));
SAM_SPECIALIZATION(TMBad::ad_aug yieldPerRecruit_i(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&, int, int, int, bool));

// }
// #endif

template<class Type>
vector<Type> yieldPerRecruit(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, bool give_log DEFARG(= false))SOURCE({
  vector<Type> r(dat.catchMeanWeight.dim(0));
  r.setZero();
  for(int i = 0; i < r.size(); ++i){
    r(i) = yieldPerRecruit_i(dat, conf, par, logF, i, 0, 300, give_log);
  }
  return r;
  })


  SAM_SPECIALIZATION(vector<double> yieldPerRecruit(dataSet<double>&, confSet&, paraSet<double>&, array<double>&, bool));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> yieldPerRecruit(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&, bool));
  

////////// Convenience functions to get Spawners per recruit and derivative //////////


// #ifndef WITH_SAM_LIB
// namespace equilibrium_fun {
// For calculations about future
template<class Type>
Type spawnersPerRecruit_i(const Type& logFbar, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, referencepointSet<Type>& rp, bool give_log DEFARG(= false))SOURCE({
  vector<Type> ls = rp.getLogSelectivity();
  PERREC_t<Type> r =  perRecruit_D(logFbar, dat, conf, par, ls, rp.aveYears, rp.nYears, rp.catchType);
  if(give_log)
    return r.logSPR;
  return exp(r.logSPR);
  });

SAM_SPECIALIZATION(double spawnersPerRecruit_i(const double&, dataSet<double>&, confSet&, paraSet<double>&, referencepointSet<double>&, bool));
SAM_SPECIALIZATION(TMBad::ad_aug spawnersPerRecruit_i(const TMBad::ad_aug&, dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, referencepointSet<TMBad::ad_aug>&, bool));

// For calculations about assessment period for one year
template<class Type>
Type spawnersPerRecruit_i(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, int i, int CT, int nYears, bool give_log DEFARG(= false))SOURCE({
  // Make referencepointSet
  referencepointSet<Type> rp(nYears, CT, i, logF, conf);
  Type logFbar = rp.logFbar(logF, conf);
  return spawnersPerRecruit_i(logFbar, dat, conf, par, rp, give_log);
  });

SAM_SPECIALIZATION(double spawnersPerRecruit_i(dataSet<double>&, confSet&, paraSet<double>&, array<double>&, int, int, int, bool));
SAM_SPECIALIZATION(TMBad::ad_aug spawnersPerRecruit_i(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&, int, int, int, bool));

// }
// #endif

template<class Type>
vector<Type> spawnersPerRecruit(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, bool give_log DEFARG(= false))SOURCE({
  vector<Type> r(dat.catchMeanWeight.dim(0));
  r.setZero();
  for(int i = 0; i < r.size(); ++i)
    r(i) = spawnersPerRecruit_i(dat, conf, par, logF, i, 0, 300, give_log);
  return r;
  });

SAM_SPECIALIZATION(vector<double> spawnersPerRecruit(dataSet<double>&, confSet&, paraSet<double>&, array<double>&, bool));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> spawnersPerRecruit(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&, bool));

// For calculations about assessment period for one year
template<class Type>
Type SPR0_i(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, int i, int CT, int nYears, bool give_log DEFARG(= false))SOURCE({
  // Make referencepointSet
  referencepointSet<Type> rp(nYears, CT, i, logF, conf);
  Type logFbar = R_NegInf; //rp.logFbar(logF, conf);
  return spawnersPerRecruit_i(logFbar, dat, conf, par, rp, give_log);
  })

SAM_SPECIALIZATION(double SPR0_i(dataSet<double>&, confSet&, paraSet<double>&, array<double>&,  int, int, int, bool));
SAM_SPECIALIZATION(TMBad::ad_aug SPR0_i(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&, int, int, int, bool));

////////// Convenience functions to get equilibrium biomass //////////


// #ifndef WITH_SAM_LIB
// namespace equilibrium_fun {
// For calculations about future
template<class Type>
Type equilibriumBiomass_i(const Type& logFbar, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, referencepointSet<Type>& rp, bool give_log DEFARG(= false))SOURCE({
  vector<Type> ls = rp.getLogSelectivity();
  PERREC_t<Type> r =  perRecruit_D(logFbar, dat, conf, par, ls, rp.aveYears, rp.nYears, rp.catchType);
  if(give_log)
    return r.logSe;
  return exp(r.logSe);
  })

SAM_SPECIALIZATION(double equilibriumBiomass_i(const double&, dataSet<double>&, confSet&, paraSet<double>&, referencepointSet<double>&, bool));
SAM_SPECIALIZATION(TMBad::ad_aug equilibriumBiomass_i(const TMBad::ad_aug&, dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, referencepointSet<TMBad::ad_aug>&, bool));


// For calculations about assessment period for one year
template<class Type>
Type equilibriumBiomass_i(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, int i, int CT, int nYears, bool give_log DEFARG(= false))SOURCE({
  // Make referencepointSet
  referencepointSet<Type> rp(nYears, CT, i, logF, conf);
  Type logFbar = rp.logFbar(logF, conf);
  return equilibriumBiomass_i(logFbar, dat, conf, par, rp, give_log);
  })


SAM_SPECIALIZATION(double equilibriumBiomass_i(dataSet<double>&, confSet&, paraSet<double>&, array<double>&, int, int, int, bool));
SAM_SPECIALIZATION(TMBad::ad_aug equilibriumBiomass_i(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&, int, int, int, bool));

// }
// #endif

template<class Type>
vector<Type> equilibriumBiomass(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, bool give_log DEFARG(= false))SOURCE({
  vector<Type> r(dat.catchMeanWeight.dim(0));
  r.setZero();
  for(int i = 0; i < r.size(); ++i)
    r(i) = equilibriumBiomass_i(dat, conf, par, logF, i, 0, 300, give_log);
  return r;
  })


SAM_SPECIALIZATION(vector<double> equilibriumBiomass(dataSet<double>&, confSet&, paraSet<double>&, array<double>&, bool));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> equilibriumBiomass(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&, bool));


////////// Convenience functions to get equilibrium biomass //////////


// #ifndef WITH_SAM_LIB
// namespace equilibrium_fun {
// For calculations about future
template<class Type>
Type B0_i(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, referencepointSet<Type>& rp, bool give_log DEFARG(= false))SOURCE({
  vector<Type> ls = rp.getLogSelectivity();
  PERREC_t<Type> r =  perRecruit_D(Type(R_NegInf), dat, conf, par, ls, rp.aveYears, rp.nYears, rp.catchType);
  if(give_log)
    return r.logSe;
  return exp(r.logSe);
  })


  SAM_SPECIALIZATION(double B0_i(dataSet<double>&, confSet&, paraSet<double>&, referencepointSet<double>&, bool));
  SAM_SPECIALIZATION(TMBad::ad_aug B0_i(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, referencepointSet<TMBad::ad_aug>&, bool));


// For calculations about assessment period for one year
template<class Type>
Type B0_i(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, int i, int CT, int nYears, bool give_log DEFARG(= false))SOURCE({
  // Make referencepointSet
  referencepointSet<Type> rp(nYears, CT, i, logF, conf);
  return B0_i(dat, conf, par, rp, give_log);
  })


SAM_SPECIALIZATION(double B0_i(dataSet<double>&, confSet&, paraSet<double>&, array<double>&, int, int, int, bool));
  SAM_SPECIALIZATION(TMBad::ad_aug B0_i(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&, int, int, int, bool));

// }
// #endif

template<class Type>
vector<Type> B0(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, bool give_log DEFARG(= false))SOURCE({
  vector<Type> r(dat.catchMeanWeight.dim(0));
  r.setZero();
  for(int i = 0; i < r.size(); ++i)
    r(i) = B0_i(dat, conf, par, logF, i, 0, 300, give_log);
  return r;
  })

  SAM_SPECIALIZATION(vector<double> B0(dataSet<double>&, confSet&, paraSet<double>&, array<double>&, bool));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> B0(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&, bool));


