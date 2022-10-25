SAM_DEPENDS(define)
SAM_DEPENDS(incidence)
SAM_DEPENDS(mvmix)
SAM_DEPENDS(predn)
SAM_DEPENDS(forecast)


template <class Type>
matrix<Type> get_nvar(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, array<Type> &logF)SOURCE({
  int stateDimN=logN.dim[0];
  // int timeSteps=logN.dim[1];
  matrix<Type> nvar(stateDimN,stateDimN);
  vector<Type> varLogN=exp(par.logSdLogN*Type(2.0));
  for(int i=0; i<stateDimN; ++i){
    for(int j=0; j<stateDimN; ++j){
      if(i!=j){nvar(i,j)=0.0;}else{nvar(i,j)=varLogN(conf.keyVarLogN(i));}
    }
  }
  return nvar;
  });

SAM_SPECIALIZATION(matrix<double> get_nvar(dataSet<double>&, confSet&, paraSet<double>&, array<double>&, array<double>&));
SAM_SPECIALIZATION(matrix<TMBad::ad_aug> get_nvar(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&));
  

template <class Type>
Type nllN(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, forecastSet<Type>& forecast, array<Type> &logN, array<Type> &logF, Recruitment<Type> &recruit, MortalitySet<Type>& mort, data_indicator<vector<Type>,Type> &keep, objective_function<Type> *of)SOURCE({
  Type nll=0;
  int stateDimN=logN.dim[0];
  int timeSteps=logN.dim[1];
  array<Type> resN(stateDimN,timeSteps-1);
  // matrix<Type> nvar(stateDimN,stateDimN);
  // vector<Type> varLogN=exp(par.logSdLogN*Type(2.0));
  // for(int i=0; i<stateDimN; ++i){
  //   for(int j=0; j<stateDimN; ++j){
  //     if(i!=j){nvar(i,j)=0.0;}else{nvar(i,j)=varLogN(conf.keyVarLogN(i));}
  //   }
  // }
  matrix<Type> nvar = get_nvar(dat, conf, par, logN, logF);
  vector<Type> fracMixN(conf.fracMixN.size());
  for(int i=0; i<conf.fracMixN.size(); ++i){fracMixN(i)=conf.fracMixN(i);}
  MVMIX_t<Type> neg_log_densityN(nvar,fracMixN);
  Eigen::LLT< Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> > lltCovN(nvar);
  matrix<Type> LN = lltCovN.matrixL();
  matrix<Type> LinvN = LN.inverse();

  matrix<Type> pn(stateDimN,timeSteps);
  pn.setZero();

  if(conf.initState){
    resN.col(0) = LinvN*(vector<Type>(logN.col(0)-par.initN));    
    nll+= neg_log_densityN(logN.col(0)-par.initN) ;//density::MVNORM(diagonalMatrix(Type(0.1),stateDimN))(logN.col(0)-par.initN); //neg_log_densityN(logN.col(0)-par.initN); // N-Process likelihood 
    SIMULATE_F(of){
      logN.col(0) = par.initN + neg_log_densityN.simulate(); //SCALE(N01(par.initN),1.0); // + neg_log_densityN.simulate();
    }
  }
  
  for(int i = 1; i < timeSteps; ++i){
    vector<Type> predN = predNFun(dat,conf,par,logN,logF,recruit,mort, i);
    pn.col(i) = predN;
    if(forecast.nYears > 0 &&
       forecast.forecastYear(i) > 0 &&
       forecast.recModel(CppAD::Integer(forecast.forecastYear(i))-1) != forecast.asRecModel){
      // Forecast
      vector<Type> Nscale(logN.rows());
      Nscale.setZero();
      Nscale += 1.0;
      Nscale(0) = sqrt(forecast.logRecruitmentVar) / sqrt(nvar(0,0));
      vector<Type> predNTmp = predN;
      predNTmp(0) = forecast.logRecruitmentMedian;
      // MVMIX_t<Type> nllTmp(nvar,Type(conf.fracMixN));
      nll+=neg_log_densityN((logN.col(i)-predNTmp) / Nscale) + (log(Nscale)).sum();
      SIMULATE_F(of){
    	if(forecast.simFlag(1) == 0){
    	  logN.col(i) = predNTmp + neg_log_densityN.simulate() * Nscale;
    	}
      }
    }else{
      resN.col(i-1) = LinvN*(vector<Type>(logN.col(i)-predN));    
      nll+=neg_log_densityN(logN.col(i)-predN); // N-Process likelihood
      SIMULATE_F(of){
    	if(forecast.nYears > 0 &&
    	   forecast.forecastYear(i) > 0){
    	  // In forecast
    	  if(forecast.simFlag(1)==0){
    	    vector<Type> noiseN = neg_log_densityN.simulate();
    	    logN.col(i) = predN + noiseN;
	    if(conf.minAge == 0){
	      logN(0,i) = predNFun(dat,conf,par,logN,logF,recruit,mort,i)(0) + noiseN(0);
	    }
    	  }
    	}else{
    	  if(conf.simFlag(1)==0){
	    vector<Type> noiseN = neg_log_densityN.simulate();
    	    logN.col(i) = predN + noiseN;
	    // Handle recruitment if minAge == 0, assuming propMat(-,0)=0
	    if(conf.minAge == 0){
	      logN(0,i) = predNFun(dat,conf,par,logN,logF,recruit,mort,i)(0) + noiseN(0);
	    }
    	  }
    	}
      }
    }
  }
  if(conf.resFlag==1){
    ADREPORT_F(resN,of);
  }
  if(CppAD::Variable(keep.sum()) && conf.initState == 0){ // add wide prior for first state, but _only_ when computing ooa residuals
    Type huge = 10;
    for (int i = 0; i < stateDimN; i++) nll -= dnorm(logN(i, 0), Type(0), huge, true);  
  } 
  REPORT_F(nvar,of);
  REPORT_F(pn,of);
  return nll;
  });

SAM_SPECIALIZATION(double nllN(dataSet<double>&, confSet&, paraSet<double>&, forecastSet<double>&, array<double>&, array<double>&, Recruitment<double>&, MortalitySet<double>&, data_indicator<vector<double>,double>&, objective_function<double>*));
SAM_SPECIALIZATION(TMBad::ad_aug nllN(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, forecastSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, Recruitment<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>&, data_indicator<vector<TMBad::ad_aug>,TMBad::ad_aug>&, objective_function<TMBad::ad_aug>*));
