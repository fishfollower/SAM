#ifndef SAM_DERIVED_HPP
#define SAM_DERIVED_HPP

// template <class Type>
// Type ssbi(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, int i){
//   int stateDimN=logN.dim[0];
//   Type ssb=0;
//   for(int j=0; j<stateDimN; ++j){
//     if(conf.keyLogFsta(0,j)>(-1)){
//       ssb += exp(logN(j,i))*exp(-exp(logF(conf.keyLogFsta(0,j),i))*dat.propF(i,j)-dat.natMor(i,j)*dat.propM(i,j))*dat.propMat(i,j)*dat.stockMeanWeight(i,j);
//     }else{
//       ssb += exp(logN(j,i))*exp(-dat.natMor(i,j)*dat.propM(i,j))*dat.propMat(i,j)*dat.stockMeanWeight(i,j);
//     }
//   }
//   return ssb;
// }
template <class Type>
Type ssbi(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, int i, bool give_log = false){
  int stateDimN=logN.dim[0];
  Type logssb = R_NegInf;
  for(int j=0; j<stateDimN; ++j){
    if(dat.propMat(i,j) > 0){
      Type lssbNew = R_NegInf;
      if(conf.keyLogFsta(0,j)>(-1)){
        lssbNew = logN(j,i) - exp(logF(conf.keyLogFsta(0,j),i)) * (dat.propF(i,j)) - dat.natMor(i,j)*dat.propM(i,j) + log(dat.propMat(i,j)) + log(dat.stockMeanWeight(i,j));	
      }else{
        lssbNew = logN(j,i) - dat.natMor(i,j)*dat.propM(i,j) + log(dat.propMat(i,j)) + log(dat.stockMeanWeight(i,j));
      }
      logssb = logspace_add2(logssb, lssbNew);
    }
  }
  if(give_log)
    return logssb;
  return exp(logssb);
}


template <class Type>
vector<Type> ssbFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, bool give_log = false){
  int timeSteps=logF.dim[1];
  vector<Type> ssb(timeSteps);
  ssb.setConstant(R_NegInf);
  for(int i=0;i<timeSteps;i++){
    ssb(i)=ssbi(dat,conf,logN,logF,i, give_log);
  }
  return ssb;
}

// template <class Type>
// matrix<Type> catchFunAge(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF){
//   int len=dat.catchMeanWeight.dim(0);
//   matrix<Type> cat(conf.maxAge - conf.minAge + 1, len);
//   cat.setZero();
//   for(int y=0;y<len;y++){
//     for(int a=conf.minAge;a<=conf.maxAge;a++){  
//       Type z=dat.natMor(y,a-conf.minAge);
//       if(conf.keyLogFsta(0,a-conf.minAge)>(-1)){
//         z+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y));
//         cat(a-conf.minAge, y)+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y))/z*exp(logN(a-conf.minAge,y))*(Type(1.0)-exp(-z))*dat.catchMeanWeight(y,a-conf.minAge);
//       }
//     }
//   }
//   return cat;
// }

template <class Type>
matrix<Type> catchFunAge(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, bool give_log = false){
  int len=dat.catchMeanWeight.dim(0);
  matrix<Type> cat(conf.maxAge - conf.minAge + 1, len);
  cat.setConstant(R_NegInf);
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      Type logz=log(dat.natMor(y,a-conf.minAge) + 1e-12);
      if(conf.keyLogFsta(0,a-conf.minAge)>(-1)){
        logz = logspace_add2(logz, logF(conf.keyLogFsta(0,a-conf.minAge),y));
	Type tmp = logF(conf.keyLogFsta(0,a-conf.minAge),y) - logz + logN(a-conf.minAge,y) + logspace_sub2(Type(0.0),-exp(logz)) + log(dat.catchMeanWeight(y,a-conf.minAge));
	cat(a-conf.minAge, y) = logspace_add2(cat(a-conf.minAge, y), tmp);
      }
    }
  }
  if(give_log)
    return cat;
  return cat.array().exp().matrix();
}


// template <class Type>
// vector<Type> catchFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF){
//   matrix<Type> cat = catchFunAge(dat,conf,logN,logF);  
//   return (vector<Type>)cat.colwise().sum();
// }

template <class Type>
vector<Type> catchFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, bool give_log = false){
  matrix<Type> cat = catchFunAge(dat,conf,logN,logF, true);
  vector<Type> catY(cat.cols());
  catY.setConstant(R_NegInf);
  for(int i = 0; i < cat.cols(); ++i)
    for(int j = 0; j < cat.rows(); ++j)
      catY(i) = logspace_add2(catY(i), cat(j,i));
  if(give_log)
    return catY;
  return exp(catY);
}


template <class Type>
vector<Type> varLogCatchFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, paraSet<Type> par){
  int len=dat.landFrac.dim(0);
  vector<Type> cat(len);
  cat.setZero();
  vector<Type> varLogCat(len);
  varLogCat.setZero();
  Type CW=0;
  Type Ca=0;
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      Type z=dat.natMor(y,a-conf.minAge);
      if(conf.keyLogFsta(0,a-conf.minAge)>(-1)){
        z+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y));
        CW=dat.catchMeanWeight(y,a-conf.minAge);
        Ca=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y))/z*exp(logN(a-conf.minAge,y))*(Type(1.0)-exp(-z));
        cat(y)+=Ca*CW;
        varLogCat(y)+=exp(2.0*par.logSdLogObs(conf.keyVarObs(0,a-conf.minAge)))*CW*CW*Ca*Ca; 
      }
    }
    varLogCat(y)/=cat(y)*cat(y);
  }
  return varLogCat;
}

template <class Type>
vector<Type> landFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF){
  int len=dat.landMeanWeight.dim(0);
  vector<Type> land(len);
  land.setZero();
  Type LW=0;
  Type LF=0;
  Type LWLF=0;
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){
      Type z=dat.natMor(y,a-conf.minAge);
      if(conf.keyLogFsta(0,a-conf.minAge)>(-1)){
        z+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y));
        LW=dat.landMeanWeight(y,a-conf.minAge);
        if(LW<0){LW=0;}
	LF=dat.landFrac(y,a-conf.minAge);
        if(LF<0){LF=0;}
        LWLF=LW*LF;
        land(y)+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y))/z*exp(logN(a-conf.minAge,y))*(Type(1.0)-exp(-z))*LWLF;
      }
    }
  }
  return land;
}

template <class Type>
vector<Type> varLogLandFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, paraSet<Type> par){
  int len=dat.landMeanWeight.dim(0);
  vector<Type> land(len);
  land.setZero();
  vector<Type> varLogLand(len);
  varLogLand.setZero();
  Type LW=0;
  Type LF=0;
  Type LWLF=0;
  Type Ca=0;
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){
      Type z=dat.natMor(y,a-conf.minAge);
      if(conf.keyLogFsta(0,a-conf.minAge)>(-1)){
        z+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y));
        LW=dat.landMeanWeight(y,a-conf.minAge);
        if(LW<0){LW=0;}
        LF=dat.landFrac(y,a-conf.minAge);
        if(LF<0){LF=0;}
        LWLF=LW*LF;
        Ca=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y))/z*exp(logN(a-conf.minAge,y))*(Type(1.0)-exp(-z));
        land(y)+=Ca*LWLF;
        varLogLand(y)+=exp(2.0*par.logSdLogObs(conf.keyVarObs(0,a-conf.minAge)))*LWLF*LWLF*Ca*Ca;
      }
    }
    varLogLand(y)/=land(y)*land(y);
  }
  return varLogLand;
}


template <class Type>
vector<Type> disFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF){
  int len=dat.disMeanWeight.dim(0);
  vector<Type> dis(len);
  dis.setZero();
  Type DW=0;
  Type DF=0;
  Type DWDF=0;
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){
      Type z=dat.natMor(y,a-conf.minAge);
      if(conf.keyLogFsta(0,a-conf.minAge)>(-1)){
        z+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y));
        DW=dat.disMeanWeight(y,a-conf.minAge);
        if(DW<0){DW=0;}
	DF=1.0 - (dat.landFrac(y,a-conf.minAge) - 1e-8);
        if(DF<0){DF=0;}
        DWDF=DW*DF;
        dis(y)+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y))/z*exp(logN(a-conf.minAge,y))*(Type(1.0)-exp(-z))*DWDF;
      }
    }    
  }
  return dis;
}


template <class Type>
vector<Type> fsbFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF){
  int len=dat.landFrac.dim(0);
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

// template <class Type>
// vector<Type> fbarFun(confSet &conf, array<Type> &logF){
//   int timeSteps=logF.dim[1];
//   vector<Type> fbar(timeSteps);
//   fbar.setZero();
//   for(int y=0;y<timeSteps;y++){  
//     for(int a=conf.fbarRange(0);a<=conf.fbarRange(1);a++){  
//       fbar(y)+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y));
//     }
//     fbar(y)/=Type(conf.fbarRange(1)-conf.fbarRange(0)+1);
//   }
//   return fbar;
// }


template <class Type>
vector<Type> fbarFun(confSet &conf, array<Type> &logF, bool give_log = false){
  int timeSteps=logF.dim[1];
  vector<Type> fbar(timeSteps);
  fbar.setConstant(R_NegInf);
  for(int y=0;y<timeSteps;y++){  
    for(int a=conf.fbarRange(0);a<=conf.fbarRange(1);a++){
      fbar(y) = logspace_add2(fbar(y), logF(conf.keyLogFsta(0,a-conf.minAge),y));
    }
    fbar(y) -= log(Type(conf.fbarRange(1)-conf.fbarRange(0)+1));
  }
  if(give_log)
    return fbar;
  return exp(fbar);
}


template <class Type>
vector<Type> landFbarFun(dataSet<Type> &dat, confSet &conf, array<Type> &logF){
  int timeSteps=dat.landFrac.dim(0);
  vector<Type> fbar(timeSteps);
  fbar.setZero();
  for(int y=0;y<timeSteps;y++){  
    for(int a=conf.fbarRange(0);a<=conf.fbarRange(1);a++){  
      fbar(y)+= dat.landFrac(y,a-conf.minAge) * exp(logF(conf.keyLogFsta(0,a-conf.minAge),y));
    }
    fbar(y)/=Type(conf.fbarRange(1)-conf.fbarRange(0)+1);
  }
  return fbar;
}

template <class Type>
vector<Type> disFbarFun(dataSet<Type> &dat, confSet &conf, array<Type> &logF){
  int timeSteps=dat.landFrac.dim(0);
  vector<Type> fbar(timeSteps);
  fbar.setZero();
  for(int y=0;y<timeSteps;y++){  
    for(int a=conf.fbarRange(0);a<=conf.fbarRange(1);a++){  
      fbar(y)+= (1.0 - (dat.landFrac(y,a-conf.minAge) - 1e-8)) * exp(logF(conf.keyLogFsta(0,a-conf.minAge),y));
    }
    fbar(y)/=Type(conf.fbarRange(1)-conf.fbarRange(0)+1);
  }
  return fbar;
}

#endif
