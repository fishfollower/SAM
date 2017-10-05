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
 
  vector<Type> predN(stateDimN); 
  Type thisSSB=Type(0); 
  for(int i=1;i<timeSteps;i++){ 
    if(conf.stockRecruitmentModelCode==0){ // straight RW 
      predN(0)=logN(0,i-1);
    }else{
      if((i-conf.minAge)>=0){
        thisSSB=ssbi(dat,conf,logN,logF,i-conf.minAge);
      }else{
        thisSSB=ssbi(dat,conf,logN,logF,0); // use first in beginning       
      } 
      if(conf.stockRecruitmentModelCode==1){//ricker
        predN(0)=par.rec_loga(0)+log(thisSSB)-exp(par.rec_logb(0))*thisSSB;
      }else{
        if(conf.stockRecruitmentModelCode==2){//BH
          predN(0)=par.rec_loga(0)+log(thisSSB)-log(1.0+exp(par.rec_logb(0))*thisSSB); 
        }else{
          error("SR model code not recognized");
        }
      }
    }
  
    for(int j=1; j<stateDimN; ++j){
      if(conf.keyLogFsta(0,j-1)>(-1)){
        predN(j)=logN(j-1,i-1)-exp(logF(conf.keyLogFsta(0,j-1),i-1))-dat.natMor(i-1,j-1); 
      }else{
        predN(j)=logN(j-1,i-1)-dat.natMor(i-1,j-1); 
      }
    }  
    if(conf.maxAgePlusGroup==1){
      predN(stateDimN-1)=log(exp(logN(stateDimN-2,i-1)-exp(logF(conf.keyLogFsta(0,stateDimN-2),i-1))-dat.natMor(i-1,stateDimN-2))+
                             exp(logN(stateDimN-1,i-1)-exp(logF(conf.keyLogFsta(0,stateDimN-1),i-1))-dat.natMor(i-1,stateDimN-1))); 
    }
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
