SAM_DEPENDS(convenience)
SAM_DEPENDS(define)
SAM_DEPENDS(forecast)
SAM_DEPENDS(mvmix)

#ifndef WITH_SAM_LIB
namespace f_fun {

  template <class Type>
  Type jacobiUVtrans( array<Type>& logF) SOURCE({
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
  
    return nr*log(fabs(A.determinant()));
    });
  
  SAM_SPECIALIZATION(double jacobiUVtrans(array<double>&));
  SAM_SPECIALIZATION(TMBad::ad_aug jacobiUVtrans(array<TMBad::ad_aug>&));


  
  template <class Type>
  Type nllFseparable(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, forecastSet<Type>& forecast, array<Type> &logF, data_indicator<vector<Type>,Type> &keep, objective_function<Type> *of)SOURCE({
  
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
  
    Type rhoU = toInterval((Type)par.sepFlogitRho(0),Type(-1.0),Type(1.0),Type(2.0));
    Type rhoV = toInterval((Type)par.sepFlogitRho(1),Type(-1.0),Type(1.0),Type(2.0));
    
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
    })

  SAM_SPECIALIZATION(double nllFseparable(dataSet<double>&, confSet&, paraSet<double>&, forecastSet<double>&, array<double>&, data_indicator<vector<double>,double>&, objective_function<double>*));
  SAM_SPECIALIZATION(TMBad::ad_aug nllFseparable(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, forecastSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&, data_indicator<vector<TMBad::ad_aug>,TMBad::ad_aug>&, objective_function<TMBad::ad_aug>*));



  
} // End of namespace F functionns
#endif

  
template<class Type>
matrix<Type> get_fvar(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logF)SOURCE({
  // using CppAD::abs;
  int stateDimF=logF.dim[0];
  int timeSteps=logF.dim[1];
  int noFleets=conf.keyLogFsta.dim[0];
  int stateDimN=conf.keyLogFsta.dim[1];
  vector<Type> sdLogFsta=exp(par.logSdLogFsta);
  array<Type> resF(stateDimF,timeSteps-1);
  matrix<Type> fvar(stateDimF,stateDimF);
  matrix<Type> fcor(stateDimF,stateDimF);
  vector<Type> fsd(stateDimF);  
  vector<int> statesFleets(stateDimF);
  
  //Fill statesFleets:
  for(int f=0; f<noFleets;f++){
    for(int i=0;i<stateDimF;i++){
      for(int j=0;j<stateDimN;j++){
	if(conf.keyLogFsta(f,j)==i){
	  statesFleets(i)=f;
	  break;
	}  
      }
    }
  }
  
  fcor.setZero();
  // for(int i=0; i<stateDimF; ++i){
  //   fcor(i,i)=1.0;
  // }
  int count=0; //if corFlag varies between 0-2, itrans_rho is shorter than comm fleet length
  for(int f=0;f<noFleets;f++){
    bool nxtPar = false;
    for(int i=0; i<stateDimF; ++i){
      fcor(i,i)=1.0;
      for(int j=0; j<i; ++j){
        if(statesFleets(i)==f && statesFleets(j)==f){
	  switch(conf.corFlag(f)){
	  case 0:		// Independent
	    fcor(j,i) = 0.0;
	    break;
	  case 10:		// Independent
	    fcor(j,i) = 0.0;
	    break;
	  case 1:		// Compound symmetry
	    fcor(j,i)=toInterval((Type)par.itrans_rho(count),Type(-1.0),Type(1.0),Type(2.0));
	    nxtPar = true;
	    break;
	  case 11:		// Compound symmetry
	    fcor(j,i)=toInterval((Type)par.itrans_rho(count),Type(-1.0),Type(1.0),Type(2.0));
	    nxtPar = true;
	    break;
	  case 2:		// AR(1) structure
	    fcor(j,i)=pow(toInterval((Type)par.itrans_rho(count),Type(-1.0),Type(1.0),Type(2.0)),abs(i-j));
	    nxtPar = true;
	    break;
	  case 12:		// AR(1) structure
	    fcor(j,i)=pow(toInterval((Type)par.itrans_rho(count),Type(-1.0),Type(1.0),Type(2.0)),abs(i-j));
	    nxtPar = true;
	    break;
	    // case 3: separable structure
	  case 4:		// (almost) Perfect correlation
            fcor(j,i)= 0.99;
	    break;
	  default:
	    Rf_error("F correlation not implemented");
	    break;
	  }
	  fcor(i,j) = fcor(j,i);
        }
      }
    }
    if(nxtPar)
      count++;
  }

  for(int i=0; i<stateDimF; ++i){
    bool stop = false;
    int ff = 0;
    int j = 0;
    for(ff=0; ff<noFleets; ff++){
      for(j=0; j<stateDimN; j++){
        if(conf.keyLogFsta(ff,j)==i){
          stop=true;
          break;
        } 
      }
      if(stop)break;
    }
    fsd(i)=sdLogFsta(conf.keyVarF(ff,j));
  }
 
  for(int i=0; i<stateDimF; ++i){
    for(int j=0; j<stateDimF; ++j){
      fvar(i,j)=fsd(i)*fsd(j)*fcor(i,j);
    }
  }
  return fvar;
  })

SAM_SPECIALIZATION(matrix<double> get_fvar(dataSet<double>&, confSet&, paraSet<double>&, array<double>&));
SAM_SPECIALIZATION(matrix<TMBad::ad_aug> get_fvar(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&));



template<class Type>
vector<Type> get_fphi(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logF)SOURCE({
  // using CppAD::abs;
  int stateDimF=logF.dim[0];
  int timeSteps=logF.dim[1];
  int noFleets=conf.keyLogFsta.dim[0];
  int stateDimN=conf.keyLogFsta.dim[1];
  vector<int> statesFleets(stateDimF);
  
  //Fill statesFleets:
  for(int f=0; f<noFleets;f++){
    for(int i=0;i<stateDimF;i++){
      for(int j=0;j<stateDimN;j++){
	if(conf.keyLogFsta(f,j)==i){
	  statesFleets(i)=f;
	  break;
	}  
      }
    }
  }
  
  vector<Type> phi(stateDimF);
  phi = Type(1.0);
  int whichPar=-1;
  int lastSuchFleet = -1;
  for(int i=0; i<stateDimF; ++i){
    int f=statesFleets(i);
    if(conf.corFlag(f)>9.5){
      if(f!=lastSuchFleet){
        lastSuchFleet=f;
	whichPar++; 
      }
      //std::cout<<statesFleets<<"\n";
      //std::cout<<i<<" "<<f<<" "<<lastSuchFleet<<" "<<whichPar<<"\n";
      phi(i)=toInterval((Type)par.itrans_rhoTime(whichPar),Type(-1.0),Type(1.0),Type(2.0));
    }
  }
  return phi;
  })

SAM_SPECIALIZATION(vector<double> get_fphi(dataSet<double>&, confSet&, paraSet<double>&, array<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> get_fphi(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&));
  
template <class Type>
Type nllF(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, forecastSet<Type>& forecast, array<Type> &logF, data_indicator<vector<Type>,Type> &keep, objective_function<Type> *of)SOURCE({
  Type nll=0; 
  int stateDimF=logF.dim[0];
  int timeSteps=logF.dim[1];
  // int stateDimN=conf.keyLogFsta.dim[1];
  vector<Type> sdLogFsta=exp(par.logSdLogFsta);
  array<Type> resF(stateDimF,timeSteps-1);

  
  if(conf.corFlag(0)==3){
    SAM_ASSERT(getCatchFleets(dat.fleetTypes).size() == 1,"separable F correlation structure is only implemented for a single catch fleet.");
    // Only works for one catch fleet!
    return(f_fun::nllFseparable(dat, conf, par, forecast, logF, keep ,of));
  }
 
  //density::MVNORM_t<Type> neg_log_densityF(fvar);
  matrix<Type> fvar = get_fvar(dat, conf, par, logF);
  vector<Type> phi = get_fphi(dat, conf, par, logF);
  MVMIX_t<Type> neg_log_densityF(fvar,Type(conf.fracMixF));
  Eigen::LLT< Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> > lltCovF(fvar);
  matrix<Type> LF = lltCovF.matrixL();
  matrix<Type> LinvF = LF.inverse();

  if(conf.initState){
    resF.col(0) = LinvF*(vector<Type>(logF.col(0)-par.initF));
    nll+=neg_log_densityF(logF.col(0)-par.initF);//density::MVNORM(diagonalMatrix(Type(0.1),stateDimF))((vector<Type>)(logF.col(0)-par.initF));
    SIMULATE_F(of){
	if(conf.simFlag(0)==0){
	  logF.col(0)=par.initF+neg_log_densityF.simulate();
	}
      }
  }
  
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
      nll+=neg_log_densityF(logF.col(i)-phi*logF.col(i-1)); // F-Process likelihood
      SIMULATE_F(of){
	if(conf.simFlag(0)==0){
	  logF.col(i)=logF.col(i-1)+neg_log_densityF.simulate();
	}
      }
    }
  }

  if(CppAD::Variable(keep.sum()) && conf.initState == 0){ // add wide prior for first state, but _only_ when computing ooa residuals
    Type huge = 10;
    for (int i = 0; i < stateDimF; i++) nll -= dnorm(logF(i, 0), Type(0), huge, true);  
  } 

  if(conf.resFlag==1){
    ADREPORT_F(resF,of);
  }
  REPORT_F(fvar,of);
  REPORT_F(phi,of);
  return nll;
  }
  )



SAM_SPECIALIZATION(double nllF(dataSet<double>&, confSet&, paraSet<double>&, forecastSet<double>&, array<double>&, data_indicator<vector<double>,double>&, objective_function<double>*));
SAM_SPECIALIZATION(TMBad::ad_aug nllF(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, forecastSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&, data_indicator<vector<TMBad::ad_aug>,TMBad::ad_aug>&, objective_function<TMBad::ad_aug>*));
