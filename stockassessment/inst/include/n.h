template <class Type>
Type nllN(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, array<Type> &logF, data_indicator<vector<Type>,Type> &keep, objective_function<Type> *of){ 
  Type nll=0;
  int stateDimN=logN.dim[0];
  int timeSteps=logN.dim[1];
  array<Type> resN(stateDimN,timeSteps-1);
  matrix<Type> nvar(stateDimN,stateDimN);
  vector<Type> varLogN=exp(par.logSdLogN*Type(2.0));
  for(int i=0; i<stateDimN; ++i){
    for(int j=0; j<stateDimN; ++j){
      if(i!=j){nvar(i,j)=0.0;}else{nvar(i,j)=varLogN(conf.keyVarLogN(i));}
    }
  }
  density::MVNORM_t<Type> neg_log_densityN(nvar);
  Eigen::LLT< Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> > lltCovN(nvar);
  matrix<Type> LN = lltCovN.matrixL();
  matrix<Type> LinvN = LN.inverse();

  for(int i = 1; i < timeSteps; ++i){ 
    vector<Type> predN = predNFun(dat,conf,par,logN,logF,i); 

    resN.col(i-1) = LinvN*(vector<Type>(logN.col(i)-predN));    
    nll+=neg_log_densityN(logN.col(i)-predN); // N-Process likelihood 
    SIMULATE_F(of){
      if(conf.simFlag==0){
        logN.col(i) = predN + neg_log_densityN.simulate();
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
