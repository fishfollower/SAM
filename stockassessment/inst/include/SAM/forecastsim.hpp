SAM_DEPENDS(define)
SAM_DEPENDS(recruitment)
SAM_DEPENDS(incidence)
SAM_DEPENDS(mvmix)
SAM_DEPENDS(derived)
SAM_DEPENDS(f)
SAM_DEPENDS(predn)
SAM_DEPENDS(n)


template <class Type>
void forecastSimulation(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, forecastSet<Type>& forecast, array<Type>& logN, array<Type>& logF, Recruitment<Type>& recruit, MortalitySet<Type>& mort, objective_function<Type> *of)SOURCE({
  // Only for forecast simulation
  if(forecast.nYears == 0 || !(isDouble<Type>::value) || !(of->do_simulate))
    return;

  // MortalitySet<Type> mort2(mort);
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
    forecast.updateForecast(i, logF, logN, dat, conf, par, recruit, mort);
    // Simulate F
    // int forecastIndex = CppAD::Integer(forecast.forecastYear(i))-1;
    if(forecast.simFlag(0) == 0){
      Type timeScale = forecast.forecastCalculatedLogSdCorrection(i);
      logF.col(indx) = (vector<Type>)forecast.forecastCalculatedMedian.col(i) + neg_log_densityF.simulate() * timeScale;
      mort.updateYear(dat,conf,par,logF,indx);
    }
    // Simulate N
    if(forecast.simFlag(1) == 0){
      vector<Type> predN = predNFun(dat,conf,par,logN,logF,recruit,mort,indx);
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
	logN(0,indx) = predNFun(dat,conf,par,logN,logF,recruit,mort,indx)(0) + noiseN(0);
    }
  }
  return;
})

SAM_SPECIALIZATION(void forecastSimulation(dataSet<double>&, confSet&, paraSet<double>&, forecastSet<double>&, array<double>&, array<double>&, Recruitment<double>&, MortalitySet<double>&, objective_function<double>*));
SAM_SPECIALIZATION(void forecastSimulation(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, forecastSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, Recruitment<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>&, objective_function<TMBad::ad_aug>*));
