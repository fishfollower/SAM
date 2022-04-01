#pragma once
#ifndef SAM_HCR_HPP
#define SAM_HCR_HPP

template <class Type>
Type hcr_min(Type a, Type b){
  return 0.5 * (a + b - sqrt(1e-4 + (a-b) * (a-b)));
}
template <class Type>
Type hcr_max(Type a, Type b){
  return 0.5 * (a + b + sqrt(1e-4 + (a-b) * (a-b)));
}

template <class Type>
Type hcr(Type ssb, vector<Type> hcrConf){
  Type Ftarget = hcrConf(0);
  Type Flim = hcrConf(1);
  Type Flow = hcrConf(2);
  Type Blim = hcrConf(3);
  Type Blow = hcrConf(4);
  Type Btrigger = hcrConf(5);

  Type newF = CppAD::CondExpLt(ssb,
			       Blow,
			       Flow,
			       hcr_min(Ftarget, hcr_max(Flim, Flim + (ssb - Blim) * (Ftarget - Flim) / (Btrigger - Blim))));
  return log(hcr_max(newF, (Type)exp(-10)));	  
}

template <class Type>
void forecastSimulation(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, forecastSet<Type>& forecast, array<Type>& logN, array<Type>& logF, Recruitment<Type> recruit, objective_function<Type> *of){
  // Only for forecast simulation
  if(forecast.nYears == 0 || !(isDouble<Type>::value) || !(of->do_simulate))
    return;

  // General setup
  // int stateDimF=logF.dim[0];
  // int timeSteps=logF.dim[1];
  // int stateDimN=conf.keyLogFsta.dim[1];

  // Setup for F
  matrix<Type> fvar = get_fvar(dat, conf, par, logF);
  MVMIX_t<Type> neg_log_densityF(fvar,Type(conf.fracMixF));

  // Setup for N
  matrix<Type> nvar = get_nvar(dat, conf, par, logN, logF);
  vector<Type> fracMixN(conf.fracMixN.size());
  for(int i=0; i<conf.fracMixN.size(); ++i){fracMixN(i)=conf.fracMixN(i);}
  MVMIX_t<Type> neg_log_densityN(nvar,fracMixN);


  int nYears = forecast.nYears;
  for(int i = 0; i < nYears; ++i){
    int indx = forecast.forecastYear.size() - nYears + i;
    // Update forecast
    forecast.updateForecast(i, logF, logN, dat, conf, par, recruit);
    // Simulate F
    // int forecastIndex = CppAD::Integer(forecast.forecastYear(i))-1;
    if(forecast.simFlag(0) == 0){
      Type timeScale = forecast.forecastCalculatedLogSdCorrection(i);
      logF.col(indx) = (vector<Type>)forecast.forecastCalculatedMedian.col(i) + neg_log_densityF.simulate() * timeScale;
    }
    // Simulate N
    if(forecast.simFlag(1) == 0){
      vector<Type> predN = predNFun(dat,conf,par,logN,logF,recruit,indx);
      vector<Type> Nscale(logN.rows());
      Nscale.setConstant((Type)1.0);
      if(forecast.recModel(CppAD::Integer(forecast.forecastYear(indx))-1) != forecast.asRecModel){
	Nscale(0) = sqrt(forecast.logRecruitmentVar) / sqrt(nvar(0,0));
	predN(0) = forecast.logRecruitmentMedian;
      }
      // logN.col(indx) = predN + neg_log_densityN.simulate();// * Nscale;
      vector<Type> noiseN = neg_log_densityN.simulate();
      logN.col(indx) = predN + noiseN;
      if(conf.minAge == 0 &&
	 forecast.recModel(CppAD::Integer(forecast.forecastYear(indx))-1) == forecast.asRecModel)
	logN(0,indx) = predNFun(dat,conf,par,logN,logF,recruit,i)(0) + noiseN(0);
    }
  }
  return;
}


#endif
