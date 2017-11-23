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
  Type ssb=0;
  //for(int j=0; j<stateDimN; ++j){
  //  ssb+=exp(logN(j,i))*exp(-totF(j,i)*dat.propF(i,j)-dat.natMor(i,j)*dat.propM(i,j))*dat.propMat(i,j)*dat.stockMeanWeight(i,j);
  //}
  int noFleets=conf.keyLogFsta.dim[0];
  for(int f=0; f<noFleets;f++){
	for(int j=0; j<stateDimN; ++j){	
	  if(conf.keyLogFsta(f,j)>(-1)){
	  	ssb+=exp(logN(j,i))*exp(-exp(logF(conf.keyLogFsta(f,j),i))*dat.propF(i,j,f)-dat.natMor(i,j)*dat.propM(i,j))*dat.propMat(i,j)*dat.stockMeanWeight(i,j);
	  }
	}
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
