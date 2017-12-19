template <class Type>
int yearsPFun(confSet &conf, dataSet<Type> &dat){
  
  int noFleets=conf.keyLogFsta.dim[0];
  Type minYear = dat.years(dat.noYears-1);
  Type maxYear = 0;
  int noYearsLAI;
  for(int f=0;f<noFleets;f++){
  	if(dat.fleetTypes(f)==6){
      for(int y=0;y<dat.noYears;y++){
        if(!isNAINT(dat.idx1(f,y))){
          if(dat.years(y)<minYear){
	  	    minYear = dat.years(y);
		    break;	
		  }
	    }
	  }
	  for(int y=0;y<dat.noYears;y++){
        if(!isNAINT(dat.idx1(f,y))){
          if(dat.years(y)>maxYear){
	  	    maxYear = dat.years(y);
	  	  }
        }
	  }	  
    }  	
    noYearsLAI = CppAD::Integer(maxYear - minYear + 1);
  }
  return noYearsLAI;
}


template <class Type>
array<Type> scalePFun(confSet &conf, dataSet<Type> &dat, array<Type> &logP){
  
  int noFleets=conf.keyLogFsta.dim[0];
  int noYearsLAI = yearsPFun(conf,dat);
  int nlogP = logP.dim[0]+1;
  array<Type> logPS(nlogP,logP.dim[1]);
    
  Type totProp;
  for(int j=0;j<noYearsLAI;j++){
    totProp=0;
    for(int i=0;i<(nlogP-1);i++){
      totProp += exp(logP(i,j));
    }
    for(int i=0; i<(nlogP-1);i++){
      logPS(i+1,j) = log(exp(logP(i,j)) / (1+totProp));
    }
    logPS(0,j) = log(1 - totProp / (1+totProp));
  }      
  return logPS;
}

template <class Type>
vector<Type> scaleWeekFun(paraSet<Type> &par, dataSet<Type> &dat, array<Type> &logP){

  int nlogP = logP.dim[0]+1;
  int maxLAIsurv = par.logAlphaSCB.size()+nlogP;
  vector<Type> varAlphaSCB(maxLAIsurv);
  //Take contribution of each survey to component and scale to 1
  int indx; 
  for(int i=0; i<nlogP;i++){
    Type totProp_alpha = 0;
    int idxmin=0; int idxmax=0;
    idxmin = dat.minWeek(i);
    idxmax = dat.maxWeek(i);
    
    for(int j=(idxmin+1);j<=idxmax;j++){
      // Substract i because aSCB is only 7 long but I'm estimating 11 valus
      indx = j - 1 - i;
      totProp_alpha += exp(par.logAlphaSCB(indx));
    }
    for(int j=(idxmin+1); j<=idxmax; ++j){
      indx = j -1 - i;
      varAlphaSCB(j) = log(exp(par.logAlphaSCB(indx)) / (1+totProp_alpha));
    }
    varAlphaSCB(idxmin) = log(1 - totProp_alpha / (1+totProp_alpha));
  }
  return varAlphaSCB;
}


template <class Type>
array<Type> totFFun(confSet &conf, array<Type> &logF){
  int noFleets=conf.keyLogFsta.dim[0];
  int stateDimN=conf.keyLogFsta.dim[1];
  int timeSteps=logF.dim[1];
  array<Type> totF(stateDimN,timeSteps);
  totF.setZero();  
  for(int i=0;i<timeSteps;i++){ 
    for(int j=0; j<stateDimN; ++j){
      for(int f=0; f<noFleets; ++f){
        if(conf.keyLogFsta(f,j)>(-1)){
          totF(j,i) += exp(logF(conf.keyLogFsta(f,j),i));
        }
      }
    }
  }
  return totF;
}

template <class Type>
Type ssbi(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &totF, int i, array<Type> &logF){
  int stateDimN=logN.dim[0];
  int noFleets=conf.keyLogFsta.dim[0];
  Type ssb=0;
  vector<Type> Z(stateDimN);
  Z.setZero();
  for(int j=0; j<stateDimN; ++j){
    Z(j) += dat.natMor(i,j)*dat.propM(i,j);
    for(int f=0; f<noFleets;f++){
      if(conf.keyLogFsta(f,j)>(-1)){
        Z(j) += exp(logF(conf.keyLogFsta(f,j),i))*dat.propF(i,j,f);
      }
    }
    ssb+=exp(logN(j,i))*exp(-Z(j))*dat.propMat(i,j)*dat.stockMeanWeight(i,j);
  }
  return ssb;
}

template <class Type>
vector<Type> ssbFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF){
  int timeSteps=logF.dim[1];
  array<Type> totF=totFFun(conf, logF);
  vector<Type> ssb(timeSteps);
  ssb.setZero();
  for(int i=0;i<timeSteps;i++){
    ssb(i)=ssbi(dat,conf,logN,totF,i,logF);
  }
  return ssb;
}

template <class Type>
vector<Type> catchFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF){
  int len=dat.catchMeanWeight.dim(0);
  int noFleets=conf.keyLogFsta.dim[0];
  array<Type> totF=totFFun(conf, logF);
  vector<Type> cat(len);
  cat.setZero();
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      Type z=dat.natMor(y,a-conf.minAge);
      z+=totF(a-conf.minAge,y);
      for(int f=0; f<noFleets;f++){
        if(conf.keyLogFsta(f,a-conf.minAge)>(-1)){
          cat(y)+=exp(logF(conf.keyLogFsta(f,a-conf.minAge),y))/z*exp(logN(a-conf.minAge,y))*(Type(1.0)-exp(-z))*dat.catchMeanWeight(y,a-conf.minAge,f);
        }
      }
    }
  }
  return cat;
}

template <class Type>
array<Type> catchByFleetFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF){
  int len=dat.catchMeanWeight.dim(0);
  int noFleets=conf.keyLogFsta.dim[0];
  array<Type> totF=totFFun(conf, logF);
  array<Type> cat(len,noFleets);
  cat.setZero();
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      Type z=dat.natMor(y,a-conf.minAge);
      z+=totF(a-conf.minAge,y);
      for(int f=0; f<noFleets;f++){
        if(conf.keyLogFsta(f,a-conf.minAge)>(-1)){
          cat(y,f)+=exp(logF(conf.keyLogFsta(f,a-conf.minAge),y))/z*exp(logN(a-conf.minAge,y))*(Type(1.0)-exp(-z))*dat.catchMeanWeight(y,a-conf.minAge,f);
        }
      }
    }
  }
  return cat;
}

template <class Type>
vector<Type> varLogCatchFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, paraSet<Type> par){
  int len=dat.catchMeanWeight.dim(0);
  int noFleets=conf.keyLogFsta.dim[0];
  array<Type> totF=totFFun(conf, logF);
  vector<Type> cat(len);
  cat.setZero();
  vector<Type> varLogCat(len);
  varLogCat.setZero();
  Type CW=0;
  Type Ca=0;
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      Type z=dat.natMor(y,a-conf.minAge);
      z+=totF(a-conf.minAge,y);
      for(int f=0; f<noFleets;f++){
        if(conf.keyLogFsta(f,a-conf.minAge)>(-1)){
          CW=dat.catchMeanWeight(y,a-conf.minAge,f);
          Ca=exp(logF(conf.keyLogFsta(f,a-conf.minAge),y))/z*exp(logN(a-conf.minAge,y))*(Type(1.0)-exp(-z));
          cat(y)+=Ca*CW;
          varLogCat(y)+=exp(2.0*par.logSdLogObs(conf.keyVarObs(f,a-conf.minAge)))*CW*CW*Ca*Ca; 
    	}
      }
    }
    varLogCat(y)/=cat(y)*cat(y);
  }
  return varLogCat;
}

template <class Type>
vector<Type> fsbFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF){
  int len=dat.catchMeanWeight.dim(0);
  int noFleets=conf.keyLogFsta.dim[0];
  array<Type> totF=totFFun(conf, logF);
  vector<Type> fsb(len);
  fsb.setZero();
  Type sumF;
  for(int y=0;y<len;y++){  // calc logfsb
    sumF=Type(0);
    for(int a=conf.minAge;a<=conf.maxAge;a++){
	  for(int f=0; f<noFleets;f++){  
        if(conf.keyLogFsta(f,a-conf.minAge)>(-1)){
          sumF+=exp(logF(conf.keyLogFsta(f,a-conf.minAge),y));
        }
      }
    }
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      Type z=dat.natMor(y,a-conf.minAge);
      z+=totF(a-conf.minAge,y);
      for(int f=0; f<noFleets;f++){  
        if(conf.keyLogFsta(f,a-conf.minAge)>(-1)){
          fsb(y)+=(exp(logF(conf.keyLogFsta(f,a-conf.minAge),y))/sumF)*exp(logN(a-conf.minAge,y))*exp(-Type(0.5)*z)*dat.catchMeanWeight(y,a-conf.minAge,f);
        }
      }
    }
  }
  return fsb;
}

template <class Type>
vector<Type> tsbFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN){
  int timeSteps=logN.dim[1];
  vector<Type> tsb(timeSteps);
  tsb.setZero();  
  for(int y=0;y<timeSteps;y++){  
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      tsb(y)+=exp(logN(a-conf.minAge,y))*dat.stockMeanWeight(y,a-conf.minAge);
    }
  }
  return tsb;
}

template <class Type>
vector<Type> rFun(array<Type> &logN){
  int timeSteps=logN.dim[1];
  vector<Type> R(timeSteps);
  R.setZero();
  for(int y=0;y<timeSteps;y++){  
    R(y)=exp(logN(0,y));
  }
  return R;
}

template <class Type>
vector<Type> fbarFun(confSet &conf, array<Type> &logF){
  int timeSteps=logF.dim[1];
  array<Type> totF=totFFun(conf, logF);
  vector<Type> fbar(timeSteps);
  fbar.setZero();
  for(int y=0;y<timeSteps;y++){  
    for(int a=conf.fbarRange(0);a<=conf.fbarRange(1);a++){  
      fbar(y)+=totF(a-conf.minAge,y);
    }
    fbar(y)/=Type(conf.fbarRange(1)-conf.fbarRange(0)+1);
  }
  return fbar;
}
