SAM_DEPENDS(convenience)
SAM_DEPENDS(define)
SAM_DEPENDS(forecast)


template <class Type>
Type nllSeason(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, forecastSet<Type>& forecast, array<Type> &logitFseason, data_indicator<vector<Type>,Type> &keep, objective_function<Type> *of)SOURCE({

    int nSeasonPar = logitFseason.dim[0];
    int timeSteps = logitFseason.dim[1];
   int nProcesses = logitFseason.dim[2];
 
    Type nll = 0.0;
    
    for(int p = 0; p < nProcesses; ++p){
      for(int i = 1; i < timeSteps; ++i){
      	for(int s = 0; s < nSeasonPar; ++s){
	  Type b = toInterval((Type)par.seasonLogitRho(s,p), Type(0.0), Type(1.0), Type(1.0));
	  Type mu = par.seasonMu(s,p);
	  Type pred = mu + b * (logitFseason(s,i-1,p) - mu);

	  nll -= dnorm(logitFseason(s, i, p), pred, exp(par.seasonLogSd(s,p)), true);
	  if(!(forecast.nYears > 0 && forecast.forecastYear(i) > 0)){
	  // if(forecast.nYears == 0){
	    SIMULATE_F(of){
	      if(conf.simFlag(0)==0){
	  	// Do pre-forecast simulation here
	  	logitFseason(s,i,p) = rnorm(pred, exp(par.seasonLogSd(s,p)));
	      }
	    }
	  // }else if(forecast.nYears > 0 && forecast.forecastYear(i) > 0){
	  //   SIMULATE_F(of){
	  //     if(conf.simFlag(0)==0){
	  // 	// Do pre-forecast simulation here
	  // 	logitFseason(s,i,p) = rnorm(1, pred, exp(par.seasonLogSd(s,p)))(0);
	  //     }
	  //   }
	  }
	}
      }
    }
    return nll;
  }
  )



SAM_SPECIALIZATION(double nllSeason(dataSet<double>&, confSet&, paraSet<double>&, forecastSet<double>&, array<double>&, data_indicator<vector<double>,double>&, objective_function<double>*));
SAM_SPECIALIZATION(TMBad::ad_aug nllSeason(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, forecastSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&, data_indicator<vector<TMBad::ad_aug>,TMBad::ad_aug>&, objective_function<TMBad::ad_aug>*));
