#ifndef SAM_N_HPP
#define SAM_N_HPP

template <class Type>
matrix<Type> get_nvar(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, array<Type> &logF){
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
}

template <class Type>
Type nllN(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, array<Type> &logF, data_indicator<vector<Type>,Type> &keep, objective_function<Type> *of){ 
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
  MVMIX_t<Type> neg_log_densityN(nvar,Type(conf.fracMixN));
  Eigen::LLT< Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> > lltCovN(nvar);
  matrix<Type> LN = lltCovN.matrixL();
  matrix<Type> LinvN = LN.inverse();

  for(int i = 1; i < timeSteps; ++i){ 
    vector<Type> predN = predNFun(dat,conf,par,logN,logF,i);
    if(dat.forecast.nYears > 0 &&
       dat.forecast.recModel(CppAD::Integer(dat.forecast.forecastYear(i))-1) != dat.forecast.asRecModel &&
       dat.forecast.forecastYear(i) > 0){
      // Forecast
      vector<Type> Nscale(logN.rows());
      Nscale.setZero();
      Nscale += 1.0;
      Nscale(0) = sqrt(dat.forecast.logRecruitmentVar) / sqrt(nvar(0,0));
      vector<Type> predNTmp = predN;
      predNTmp(0) = dat.forecast.logRecruitmentMedian;
      // MVMIX_t<Type> nllTmp(nvar,Type(conf.fracMixN));
      nll+=neg_log_densityN((logN.col(i)-predNTmp) / Nscale) + (log(Nscale)).sum();
      SIMULATE_F(of){
    	if(dat.forecast.simFlag(1) == 0){
    	  logN.col(i) = predNTmp + neg_log_densityN.simulate() * Nscale;
    	}
      }
    }else{
      resN.col(i-1) = LinvN*(vector<Type>(logN.col(i)-predN));    
      nll+=neg_log_densityN(logN.col(i)-predN); // N-Process likelihood 
      SIMULATE_F(of){
    	if(dat.forecast.nYears > 0 &&
    	   dat.forecast.forecastYear(i) > 0){
    	  // In forecast
    	  if(dat.forecast.simFlag(1)==0){
    	    logN.col(i) = predN + neg_log_densityN.simulate();
    	  }
    	}else{
    	  if(conf.simFlag(1)==0){
    	    logN.col(i) = predN + neg_log_densityN.simulate();
    	  }
    	}
      }
    }
  }
  if(conf.resFlag==1){
    ADREPORT_F(resN,of);
  }
  if(CppAD::Variable(keep.sum())){ // add wide prior for first state, but _only_ when computing ooa residuals
    Type huge = 10;
    for (int i = 0; i < stateDimN; i++) nll -= dnorm(logN(i, 0), Type(0), huge, true);  
  } 
  return nll;
}

#endif
