template <class Type>
Type trans(Type x){
  return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);
}

template <class Type>
Type nllF(array<Type> logF, 
          int timeSteps,
          int corFlag,
          int simFlag,
          int resFlag,
          int stateDimN,
          array<int> keyLogFsta,
          array<int> keyVarF,
          vector<Type> itrans_rho,
          vector<Type> sdLogFsta,
          data_indicator<vector<Type>,Type> keep, 
          objective_function<Type> *of
	  ){
  Type nll=0; 
  int stateDimF=logF.dim[0];
  array<Type> resF(logF.dim[0],logF.dim[1]-1);
  matrix<Type> fvar(stateDimF,stateDimF);
  matrix<Type> fcor(stateDimF,stateDimF);
  vector<Type> fsd(stateDimF);  

  if(corFlag==0){
    fcor.setZero();
  }

  for(int i=0; i<stateDimF; ++i){
    fcor(i,i)=1.0;
  }

  if(corFlag==1){
    for(int i=0; i<stateDimF; ++i){
      for(int j=0; j<i; ++j){
        fcor(i,j)=trans(itrans_rho(0));
        fcor(j,i)=fcor(i,j);
      }
    } 
  }

  if(corFlag==2){
    for(int i=0; i<stateDimF; ++i){
      for(int j=0; j<i; ++j){
        fcor(i,j)=pow(trans(itrans_rho(0)),abs(Type(i-j)));
        fcor(j,i)=fcor(i,j);
      }
    } 
  }

  int i,j;
  for(i=0; i<stateDimF; ++i){
    for(j=0; j<stateDimN; ++j){
      if(keyLogFsta(0,j)==i)break;
    }
    fsd(i)=sdLogFsta(keyVarF(0,j));
  }
 
  for(i=0; i<stateDimF; ++i){
    for(j=0; j<stateDimF; ++j){
      fvar(i,j)=fsd(i)*fsd(j)*fcor(i,j);
    }
  }
  density::MVNORM_t<Type> neg_log_densityF(fvar);
  Eigen::LLT< Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> > lltCovF(fvar);
  matrix<Type> LF = lltCovF.matrixL();
  matrix<Type> LinvF = LF.inverse();

  for(int i=1;i<timeSteps;i++){
    resF.col(i-1) = LinvF*(vector<Type>(logF.col(i)-logF.col(i-1)));    
    nll+=neg_log_densityF(logF.col(i)-logF.col(i-1)); // F-Process likelihood
    SIMULATE_F(of){
      if(simFlag==0){
        logF.col(i)=logF.col(i-1)+neg_log_densityF.simulate();
      }
    }
  }

  if(CppAD::Variable(keep.sum())){ // add wide prior for first state, but _only_ when computing ooa residuals
    Type huge = 10;
    for (int i = 0; i < stateDimF; i++) nll -= dnorm(logF(i, 0), Type(0), huge, true);  
  } 

  if(resFlag==1){
    ADREPORT_F(resF,of);
  }
  return nll;
}
