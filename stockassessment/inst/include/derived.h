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
vector<Type> catchFun(array<Type> &logF,
                      array<Type> &logN,
                      int minAge,
                      int maxAge,
                      array<int> &keyLogFsta,
                      array<Type> &catchMeanWeight,
                      array<Type> &natMor
		      ){
  int len=catchMeanWeight.dim(0);
  vector<Type> cat(len);
  cat.setZero();
  for(int y=0;y<len;y++){
    for(int a=minAge;a<=maxAge;a++){  
      Type z=natMor(y,a-minAge);
      if(keyLogFsta(0,a-minAge)>(-1)){
        z+=exp(logF(keyLogFsta(0,a-minAge),y));
        cat(y)+=exp(logF(keyLogFsta(0,a-minAge),y))/z*exp(logN(a-minAge,y))*(Type(1.0)-exp(-z))*catchMeanWeight(y,a-minAge);
      }
    }
  }
  return cat;
}

template <class Type>
vector<Type> fsbFun(array<Type> &logF,
                    array<Type> &logN,
                    int minAge,
                    int maxAge,
                    array<int> &keyLogFsta,
                    array<Type> &catchMeanWeight,
                    array<Type> &natMor
		    ){
  int len=catchMeanWeight.dim(0);
  vector<Type> fsb(len);
  fsb.setZero();
  Type sumF;
  for(int y=0;y<len;y++){  // calc logfsb
    sumF=Type(0);
    for(int a=minAge;a<=maxAge;a++){  
      if(keyLogFsta(0,a-minAge)>(-1)){
        sumF+=exp(logF(keyLogFsta(0,a-minAge),y));
      }
    }
    for(int a=minAge;a<=maxAge;a++){  
      Type z=natMor(y,a-minAge);
      if(keyLogFsta(0,a-minAge)>(-1)){
        z+=exp(logF(keyLogFsta(0,a-minAge),y));
        fsb(y)+=(exp(logF(keyLogFsta(0,a-minAge),y))/sumF)*exp(logN(a-minAge,y))*exp(-Type(0.5)*z)*catchMeanWeight(y,a-minAge);
      }
    }
  }
  return fsb;
}

template <class Type>
vector<Type> tsbFun(array<Type> &logN,
                    int minAge,
                    int maxAge,
                    int timeSteps,
                    array<Type> stockMeanWeight){
  vector<Type> tsb(timeSteps);
  tsb.setZero();  
  for(int y=0;y<timeSteps;y++){  
    for(int a=minAge;a<=maxAge;a++){  
      tsb(y)+=exp(logN(a-minAge,y))*stockMeanWeight(y,a-minAge);
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
vector<Type> fbarFun(array<Type> &logF,
                     int minAge, 
                     int timeSteps, 
                     vector<int> &fbarRange,
                     array<int> &keyLogFsta){
  vector<Type> fbar(timeSteps);
  fbar.setZero();
  for(int y=0;y<timeSteps;y++){  
    for(int a=fbarRange(0);a<=fbarRange(1);a++){  
      fbar(y)+=exp(logF(keyLogFsta(0,a-minAge),y));
    }
    fbar(y)/=Type(fbarRange(1)-fbarRange(0)+1);
  }
  return fbar;
}
