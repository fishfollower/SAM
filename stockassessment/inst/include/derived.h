template <class Type>
vector<Type> ssbFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF){
  int timeSteps=logF.dim[1];
  int stateDimN=logN.dim[0];
  vector<Type> ssb(timeSteps);
  ssb.setZero();
  for(int i=0;i<timeSteps;i++){
    for(int j=0; j<stateDimN; ++j){
      if(conf.keyLogFsta(0,j)>(-1)){
        ssb(i)+=exp(logN(j,i))*exp(-exp(logF(conf.keyLogFsta(0,j),i))*dat.propF(i,j)-dat.natMor(i,j)*dat.propM(i,j))*dat.propMat(i,j)*dat.stockMeanWeight(i,j);
      }else{
        ssb(i)+=exp(logN(j,i))*exp(-dat.natMor(i,j)*dat.propM(i,j))*dat.propMat(i,j)*dat.stockMeanWeight(i,j);
      }
    }
  }
  return ssb;
}

template <class Type>
vector<Type> catchFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF){
  int len=dat.catchMeanWeight.dim(0);
  vector<Type> cat(len);
  cat.setZero();
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      Type z=dat.natMor(y,a-conf.minAge);
      if(conf.keyLogFsta(0,a-conf.minAge)>(-1)){
        z+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y));
        cat(y)+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y))/z*exp(logN(a-conf.minAge,y))*(Type(1.0)-exp(-z))*dat.catchMeanWeight(y,a-conf.minAge);
      }
    }
  }
  return cat;
}

template <class Type>
vector<Type> fsbFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF){
  int len=dat.catchMeanWeight.dim(0);
  vector<Type> fsb(len);
  fsb.setZero();
  Type sumF;
  for(int y=0;y<len;y++){  // calc logfsb
    sumF=Type(0);
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      if(conf.keyLogFsta(0,a-conf.minAge)>(-1)){
        sumF+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y));
      }
    }
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      Type z=dat.natMor(y,a-conf.minAge);
      if(conf.keyLogFsta(0,a-conf.minAge)>(-1)){
        z+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y));
        fsb(y)+=(exp(logF(conf.keyLogFsta(0,a-conf.minAge),y))/sumF)*exp(logN(a-conf.minAge,y))*exp(-Type(0.5)*z)*dat.catchMeanWeight(y,a-conf.minAge);
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
  vector<Type> fbar(timeSteps);
  fbar.setZero();
  for(int y=0;y<timeSteps;y++){  
    for(int a=conf.fbarRange(0);a<=conf.fbarRange(1);a++){  
      fbar(y)+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y));
    }
    fbar(y)/=Type(conf.fbarRange(1)-conf.fbarRange(0)+1);
  }
  return fbar;
}
