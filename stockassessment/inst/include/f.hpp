#pragma once
#ifndef SAM_F_HPP
#define SAM_F_HPP

template <class Type>
Type trans(Type x){
  return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);
}

template <class Type>
Type jacobiUVtrans( array<Type> logF){
  //  int nr=logF.rows(); //Cange to .rows eller .cols
  //  int nc=logF.cols();
  int nr=logF.cols(); //Cange to .rows eller .cols
  int nc=logF.rows();
  matrix<Type> A(nc,nc);
  for(int i=0; i<nc; ++i){
    for(int j=0; j<nc; ++j){
      A(i,j) = -1;
    }
  }
  for(int i=0; i<nc; ++i){
    A(0,i)=1;
  }
  for(int i=1; i<nc; ++i){
    A(i,i-1)=nc-1;
  }
  A/=nc;
  
  return nr*log(CppAD::abs(A.determinant()));
}

template<class Type>
matrix<Type> get_fvar(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logF){
  using CppAD::abs;
  int stateDimF=logF.dim[0];
  int timeSteps=logF.dim[1];
  int stateDimN=conf.keyLogFsta.dim[1];
  vector<Type> sdLogFsta=exp(par.logSdLogFsta);
  array<Type> resF(stateDimF,timeSteps-1);
  matrix<Type> fvar(stateDimF,stateDimF);
  matrix<Type> fcor(stateDimF,stateDimF);
  vector<Type> fsd(stateDimF);  

  if(conf.corFlag==0){
    fcor.setZero();
  }

  for(int i=0; i<stateDimF; ++i){
    fcor(i,i)=1.0;
  }

  if(conf.corFlag==1){
    for(int i=0; i<stateDimF; ++i){
      for(int j=0; j<i; ++j){
        fcor(i,j)=trans(par.itrans_rho(0));
        fcor(j,i)=fcor(i,j);
      }
    } 
  }

  if(conf.corFlag==2){
    for(int i=0; i<stateDimF; ++i){
      for(int j=0; j<i; ++j){
        fcor(i,j)=pow(trans(par.itrans_rho(0)),abs(Type(i-j)));
        fcor(j,i)=fcor(i,j);
      }
    } 
  }

  int i,j;
  for(i=0; i<stateDimF; ++i){
    for(j=0; j<stateDimN; ++j){
      if(conf.keyLogFsta(0,j)==i)break;
    }
    fsd(i)=sdLogFsta(conf.keyVarF(0,j));
  }
 
  for(i=0; i<stateDimF; ++i){
    for(j=0; j<stateDimF; ++j){
      fvar(i,j)=fsd(i)*fsd(j)*fcor(i,j);
    }
  }
  return fvar;
}

template <class Type>
Type nllF(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, forecastSet<Type>& forecast, array<Type> &logF, data_indicator<vector<Type>,Type> &keep, objective_function<Type> *of){
  Type nll=0; 
  int stateDimF=logF.dim[0];
  int timeSteps=logF.dim[1];
  // int stateDimN=conf.keyLogFsta.dim[1];
  vector<Type> sdLogFsta=exp(par.logSdLogFsta);
  array<Type> resF(stateDimF,timeSteps-1);

  
  if(conf.corFlag==3){
    return(nllFseparable(dat, conf, par, forecast, logF, keep ,of));
  }
 
  //density::MVNORM_t<Type> neg_log_densityF(fvar);
  matrix<Type> fvar = get_fvar(dat, conf, par, logF);
  MVMIX_t<Type> neg_log_densityF(fvar,Type(conf.fracMixF));
  Eigen::LLT< Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> > lltCovF(fvar);
  matrix<Type> LF = lltCovF.matrixL();
  matrix<Type> LinvF = LF.inverse();

  for(int i=1;i<timeSteps;i++){
    resF.col(i-1) = LinvF*(vector<Type>(logF.col(i)-logF.col(i-1)));

    if(forecast.nYears > 0 && forecast.forecastYear(i) > 0){
      // Forecast
      int forecastIndex = CppAD::Integer(forecast.forecastYear(i))-1;
      Type timeScale = forecast.forecastCalculatedLogSdCorrection(forecastIndex);

      nll += neg_log_densityF((logF.col(i) - (vector<Type>)forecast.forecastCalculatedMedian.col(forecastIndex)) / timeScale) + log(timeScale) * Type(stateDimF);

      SIMULATE_F(of){
    	if(forecast.simFlag(0) == 0){
    	  logF.col(i) = (vector<Type>)forecast.forecastCalculatedMedian.col(forecastIndex) + neg_log_densityF.simulate() * timeScale;
    	}
      }
    }else{
      nll+=neg_log_densityF(logF.col(i)-logF.col(i-1)); // F-Process likelihood
      SIMULATE_F(of){
	if(conf.simFlag(0)==0){
	  logF.col(i)=logF.col(i-1)+neg_log_densityF.simulate();
	}
      }
    }
  }

  if(CppAD::Variable(keep.sum())){ // add wide prior for first state, but _only_ when computing ooa residuals
    Type huge = 10;
    for (int i = 0; i < stateDimF; i++) nll -= dnorm(logF(i, 0), Type(0), huge, true);  
  } 

  if(conf.resFlag==1){
    ADREPORT_F(resF,of);
  }
  REPORT_F(fvar,of);
  return nll;
}



template <class Type>
Type nllFseparable(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, forecastSet<Type>& forecast, array<Type> &logF, data_indicator<vector<Type>,Type> &keep, objective_function<Type> *of){
  
  int stateDimF=logF.dim[0];
  int timeSteps=logF.dim[1];

  matrix<Type> SigmaU(stateDimF-1,stateDimF-1);
  SigmaU.setZero();
  vector<Type> sdU(stateDimF-1);  
  vector<Type> sdV(1);  
  for(int i=0; i<sdU.size(); ++i){
    sdU(i) = exp(par.sepFlogSd(0));
  }
  sdV(0) = exp(par.sepFlogSd(1));
  
  Type rhoU = trans(par.sepFlogitRho(0));
  Type rhoV = trans(par.sepFlogitRho(1));
    
  Type nll=0; 

  matrix<Type> logU(timeSteps,stateDimF-1);
  logU.setZero();
  vector<Type> logV(timeSteps);
  logV.setZero();
  for(int i=0; i<timeSteps; ++i){
    if(forecast.nYears > 0 && forecast.forecastYear(i) > 0){
      Rf_warning("Forecast with separable F is experimental");
      int forecastIndex = CppAD::Integer(forecast.forecastYear(i))-1;
      vector<Type> logFtmp = (vector<Type>)forecast.forecastCalculatedMedian.col(forecastIndex);
      logV(i)=(logFtmp).mean();
      for(int j=0; j<stateDimF-1; ++j){
	logU(i,j)=logFtmp(j)-logV(i);
      }
    }else{
      logV(i)=(logF).col(i).mean();
      // }
      // for(int i=0; i<timeSteps; ++i){
      for(int j=0; j<stateDimF-1; ++j){
	logU(i,j)=logF(j,i)-logV(i);
      }      
      logV(i)=(logF).col(i).mean();
    }
  }

  SigmaU.diagonal() = sdU*sdU;
  
  density::MVNORM_t<Type> nldens(SigmaU);
  for(int y=1; y<timeSteps; ++y){
    vector<Type> diff=vector<Type>(logU.row(y))-rhoU*vector<Type>(logU.row(y-1))- par.sepFalpha.segment(0,par.sepFalpha.size()-1);
    nll += nldens(diff);

    SIMULATE_F(of){
      if(conf.simFlag(0)==0){
        vector<Type> uu = nldens.simulate();
        Type sumUZero = 0;
        for(int j=0; j<stateDimF-1; ++j){
	  logU(y,j)=rhoU*logU(y-1,j) +uu(j)+ par.sepFalpha(j);
	  logF(j,y) = logU(y,j);
	  sumUZero += logU(y,j);
        }
        logF(stateDimF-1,y) = -sumUZero;
      }
    }
  }
  for(int y=1; y<timeSteps; ++y){
    nll += -dnorm(logV(y),rhoV* logV(y-1) - par.sepFalpha(par.sepFalpha.size()-1) ,sdV(0),true);
    SIMULATE_F(of){
      if(conf.simFlag(0)==0){
        logV(y)=rhoV*logV(y-1)+ rnorm( Type(0) , sdV(0))+ par.sepFalpha(par.sepFalpha.size()-1); 
        for(int j=0; j<stateDimF; ++j){
          logF(j,y) =  logF(j,y)+ logV(y) ;
        }
      }
    }
  }
  nll += -jacobiUVtrans(logF);
  
  return nll;
}


#endif
