template <class Type>
Type nllP(confSet &conf, paraSet<Type> &par, array<Type> &logP, data_indicator<vector<Type>,Type> &keep, objective_function<Type> *of){

 Type nll=0; 
  int stateDimP=logP.dim[0];
  vector<Type> varLogP=exp(par.logSdLogP*Type(2.0));
  int timeSteps=logP.dim[1];
  array<Type> resP(stateDimP,timeSteps-1);
  //Now take care of P

  matrix<Type> pvar(stateDimP,stateDimP);
  matrix<Type> pcor(stateDimP,stateDimP);
  vector<Type> psd(stateDimP);
  
  //pvar similar to F setup
  for(int i=0; i<stateDimP; ++i){
    for(int j=0; j<stateDimP; ++j){
      if(i!=j){pcor(i,j)=par.rhop(0);}else{pcor(i,j)=1.0;}
    }
    psd(i)=varLogP(conf.keyVarLogP(i));
  }
  for(int i=0; i<stateDimP; ++i){
    for(int j=0; j<stateDimP; ++j){
      pvar(i,j)=psd(i)*psd(j)*pcor(i,j);
    }
  }
  
//  MVNORM_t<Type> neg_log_densityP(pvar);
//  matrix<Type> logpredSCB(stateDimP,noYearsLAI);
//  for(int i=0;i<noYearsLAI;i++){
//    for(int j=0;j<stateDimP;j++){
//      logpredSCB(j,i) = logP(j,i) + logssb(i+noYears-noYearsLAI);
//    }
//  }
//  for(int i=1;i<noYearsLAI;i++){ //Start from year 1972 onwards
//    ans+=neg_log_densityP(logpredSCB.col(i)-logpredSCB.col(i-1)); //P-Process likelihood
//  }



  density::MVNORM_t<Type> neg_log_densityP(pvar);
  Eigen::LLT< Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> > lltCovP(pvar);
  matrix<Type> LP = lltCovP.matrixL();
  matrix<Type> LinvP = LP.inverse();

  for(int i=1;i<timeSteps;i++){
    resP.col(i-1) = LinvP*(vector<Type>(logP.col(i)-logP.col(i-1)));    
    nll+=neg_log_densityP(logP.col(i)-logP.col(i-1)); // P-Process likelihood
    SIMULATE_F(of){
      if(conf.simFlag==0){
        logP.col(i)=logP.col(i-1)+neg_log_densityP.simulate();
      }
    }
  }
  
  if(CppAD::Variable(keep.sum())){ // add wide prior for first state, but _only_ when computing ooa residuals
    Type huge = 10;
    for (int i = 0; i < stateDimP; i++) nll -= dnorm(logP(i, 0), Type(0), huge, true);  
  } 

  if(conf.resFlag==1){
    ADREPORT_F(resP,of);
  }
  REPORT_F(pvar,of);
  
  return nll;
}




  
