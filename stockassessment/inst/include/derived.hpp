#pragma once
#ifndef SAM_DERIVED_HPP
#define SAM_DERIVED_HPP

// #ifdef TMBAD_FRAMEWORK

template <class Type>
array<Type> totFFun(dataSet<Type>& dat, confSet &conf, array<Type> &logF, int Ftype = 0){
  int noFleets=conf.keyLogFsta.dim[0];
  int stateDimN=conf.keyLogFsta.dim[1];
  int timeSteps=logF.dim[1];
  if(Ftype > 0)
    timeSteps = dat.landFrac.dim[0];
  array<Type> totF(stateDimN,timeSteps);
  totF.setZero();  
  for(int i=0;i<timeSteps;i++){ 
    for(int j=0; j<stateDimN; ++j){
      for(int f=0; f<noFleets; ++f){
        if(conf.keyLogFsta(f,j)>(-1)){
	  Type Fval = exp(logF(conf.keyLogFsta(f,j),i));
	  if(Ftype == 1){
	    Fval *= dat.landFrac(i,j,f);
	  }else if(Ftype == 2){
	    Fval *= (1.0 - dat.landFrac(i,j,f));
	  }
          totF(j,i) += Fval;
        }
      }
    }
  }
  return totF;
}

template <class Type>
Type ssbi(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, MortalitySet<Type>& mort, int i, bool give_log = false){
  int stateDimN=logN.dim[0];
  // int noFleets=conf.keyLogFsta.dim[0];
  Type logssb = R_NegInf;
  for(int j=0; j<stateDimN; ++j){
    if(dat.propMat(i,j) > 0){
      //Type lssbNew = logN(j,i) - dat.natMor(i,j)*dat.propM(i,j) + log(dat.propMat(i,j)) + log(dat.stockMeanWeight(i,j));
      // for(int f=0; f<noFleets;f++){
      // 	if(conf.keyLogFsta(f,j)>(-1)){
      // 	  lssbNew -= exp(logF(conf.keyLogFsta(f,j),i)) * (dat.propF(i,j,f));
      // 	}
      // }
      Type lssbNew = logN(j,i) + log(mort.ssbSurvival_before(j,i)) + log(dat.propMat(i,j)) + log(dat.stockMeanWeight(i,j));
      logssb = logspace_add_SAM(logssb, lssbNew);
    }
  }
  if(give_log)
    return logssb;
  return exp(logssb);
}

template <class Type>
vector<Type> ssbFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF,MortalitySet<Type>& mort, bool give_log = false){
  int timeSteps=logF.dim[1];
  array<Type> totF=totFFun(dat, conf, logF);
  vector<Type> ssb(timeSteps);
  ssb.setConstant(R_NegInf);
  for(int i=0;i<timeSteps;i++){
    ssb(i)=ssbi(dat,conf,logN,logF,mort,i, give_log);
  }
  return ssb;
}

template <class Type>
matrix<Type> catchFunAge(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, MortalitySet<Type>& mort, bool give_log = false){
  int len=dat.landFrac.dim(0);
  int noFleets=conf.keyLogFsta.dim(0);
  array<Type> totF=totFFun(dat, conf, logF);
  matrix<Type> cat(conf.maxAge - conf.minAge + 1, len);
  cat.setZero();
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      // Type z=dat.natMor(y,a-conf.minAge);
      // z+=totF(a-conf.minAge,y);
      for(int f=0; f<noFleets;f++){
	if(dat.fleetTypes(f) == 0){ // Only catch fleets
	// 	if(conf.keyLogFsta(f,a-conf.minAge)>(-1)){
	// 	  cat(a-conf.minAge, y)+=exp(logF(conf.keyLogFsta(f,a-conf.minAge),y))/z*exp(logN(a-conf.minAge,y))*(Type(1.0)-exp(-z))*dat.catchMeanWeight(y,a-conf.minAge);
	// 	}
	// }
	  cat(a-conf.minAge, y) += exp(logN(a-conf.minAge, y)) * mort.fleetCumulativeIncidence(a-conf.minAge,y, f) * dat.catchMeanWeight(y, a-conf.minAge, f);
	}
      }
    }
  }
  if(give_log)
    return cat.array().log().matrix();  
  return cat;
}

template <class Type>
array<Type> catchByFleetFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, MortalitySet<Type>& mort){
  int len=dat.catchMeanWeight.dim(0);
  int noFleets=conf.keyLogFsta.dim[0];
  // array<Type> totF=totFFun(dat, conf, logF);
  array<Type> cat(len,noFleets);
  cat.setZero();
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      // Type z=dat.natMor(y,a-conf.minAge);
      // z+=totF(a-conf.minAge,y);
      for(int f=0; f<noFleets;f++){
	if(dat.fleetTypes(f) == 0){ // Only catch fleets
      //   if(conf.keyLogFsta(f,a-conf.minAge)>(-1)){
      //     cat(y,f)+=exp(logF(conf.keyLogFsta(f,a-conf.minAge),y))/z*exp(logN(a-conf.minAge,y))*(Type(1.0)-exp(-z))*dat.catchMeanWeight(y,a-conf.minAge,f);
      //   }
	  cat(y,f) += exp(logN(a-conf.minAge,y)) * mort.fleetCumulativeIncidence(a-conf.minAge,y,f) * dat.catchMeanWeight(y,a-conf.minAge,f);
	}
      }
    }
  }
  return cat;
}


template <class Type>
vector<Type> catchFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF,MortalitySet<Type>& mort, bool give_log = false){
  matrix<Type> cat = catchFunAge(dat,conf,logN,logF,mort, give_log);
  if(give_log){
    vector<Type> catY(cat.cols());
    catY.setConstant(R_NegInf);
    for(int i = 0; i < cat.cols(); ++i)
      for(int j = 0; j < cat.rows(); ++j)
	catY(i) = logspace_add_SAM(catY(i), cat(j,i));
    return catY;
  }
  vector<Type> catY(cat.cols());
  catY.setZero();
  for(int i = 0; i < cat.cols(); ++i)
    for(int j = 0; j < cat.rows(); ++j)
      catY(i) += cat(j,i);
  return catY;
}
// template <class Type>
// vector<Type> catchFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF){
//   matrix<Type> cat = catchFunAge(dat,conf,logN,logF);  
//   return (vector<Type>)cat.colwise().sum();
// }
// #endif


template <class Type>
  vector<Type> varLogCatchFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, paraSet<Type> par, MortalitySet<Type>& mort){
  int len=dat.catchMeanWeight.dim(0);
  int noFleets=conf.keyLogFsta.dim(0);
  // array<Type> totF=totFFun(dat, conf, logF);
  vector<Type> cat(len);
  cat.setZero();
  vector<Type> varLogCat(len);
  varLogCat.setZero();
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      // Type z=dat.natMor(y,a-conf.minAge);
      // z+=totF(a-conf.minAge,y);
      for(int f=0; f<noFleets;f++){
	if(dat.fleetTypes(f) == 0 && conf.keyVarObs(f,a-conf.minAge) > (-1)){ // Only catch fleets
        // if(conf.keyLogFsta(f,a-conf.minAge)>(-1)){
          Type CW=dat.catchMeanWeight(y,a-conf.minAge,f);
          // Ca=exp(logF(conf.keyLogFsta(f,a-conf.minAge),y))/z*exp(logN(a-conf.minAge,y))*(Type(1.0)-exp(-z));
	  Type Ca = exp(logN(a-conf.minAge,y)) * mort.fleetCumulativeIncidence(a-conf.minAge,y,f);
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
vector<Type> landFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, MortalitySet<Type>& mort){
  int len=dat.landMeanWeight.dim(0);
  int noFleets=conf.keyLogFsta.dim(0);
  // array<Type> totF=totFFun(dat, conf, logF);
  vector<Type> land(len);
  land.setZero();
  Type LW=0;
  Type LF=0;
  Type LWLF=0;
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){
      // Type z=dat.natMor(y,a-conf.minAge);
      // z+=totF(a-conf.minAge,y);
      for(int f=0; f<noFleets;f++){
	if(dat.fleetTypes(f) == 0){ // Only catch fleets
        // if(conf.keyLogFsta(f,a-conf.minAge)>(-1)){
          LW=dat.landMeanWeight(y,a-conf.minAge,f);
          if(LW<0){LW=0;}
          LF=dat.landFrac(y,a-conf.minAge,f);
          if(LF<0){LF=0;}
          LWLF=LW*LF;
          //land(y)+=exp(logF(conf.keyLogFsta(f,a-conf.minAge),y))/z*exp(logN(a-conf.minAge,y))*(Type(1.0)-exp(-z))*LWLF;
	  land(y) += exp(logN(a-conf.minAge,y)) * mort.fleetCumulativeIncidence(a-conf.minAge,y,f) * LWLF;
	}
      }
    }
  }
  return land;
}

template <class Type>
vector<Type> varLogLandFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, paraSet<Type> par, MortalitySet<Type>& mort){
  int len=dat.landMeanWeight.dim(0);
  int noFleets=conf.keyLogFsta.dim[0];
  // array<Type> totF=totFFun(conf, logF);
  vector<Type> land(len);
  land.setZero();
  vector<Type> varLogLand(len);
  varLogLand.setZero();
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){
      // Type z=dat.natMor(y,a-conf.minAge);
      // z+=totF(a-conf.minAge,y);
      for(int f=0; f<noFleets;f++){
	if(dat.fleetTypes(f) == 0 && conf.keyVarObs(f,a-conf.minAge) > (-1)){ // Only catch fleets
         // if(conf.keyLogFsta(f,a-conf.minAge)>(-1)){
          Type LW=dat.landMeanWeight(y,a-conf.minAge,f);
          if(LW<0){LW=0;}
          Type LF=dat.landFrac(y,a-conf.minAge,f);
          if(LF<0){LF=0;}
          Type LWLF=LW*LF;
          //Ca=exp(logF(conf.keyLogFsta(f,a-conf.minAge),y))/z*exp(logN(a-conf.minAge,y))*(Type(1.0)-exp(-z));
	  Type Ca = exp(logN(a-conf.minAge,y)) * mort.fleetCumulativeIncidence(a-conf.minAge,y,f);
          land(y)+=Ca*LWLF;
          varLogLand(y)+=exp(2.0*par.logSdLogObs(conf.keyVarObs(f,a-conf.minAge)))*LWLF*LWLF*Ca*Ca;
	}
      }
    }
    varLogLand(y)/=land(y)*land(y);
  }
  return varLogLand;
}


template <class Type>
vector<Type> disFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, MortalitySet<Type>& mort){
  int len=dat.disMeanWeight.dim(0);
  int noFleets=conf.keyLogFsta.dim[0];
  // array<Type> totF=totFFun(conf, logF);
  vector<Type> dis(len);
  dis.setZero();
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){
      // Type z=dat.natMor(y,a-conf.minAge);
      // z+=totF(a-conf.minAge,y);
      for(int f=0; f<noFleets;f++){
	if(dat.fleetTypes(f) == 0){ // Only catch fleets
         // if(conf.keyLogFsta(f,a-conf.minAge)>(-1)){
	  Type DW=dat.disMeanWeight(y,a-conf.minAge, f);
	  if(DW<0){DW=0;}
	  Type DF=1.0 - (dat.landFrac(y,a-conf.minAge, f) - 1e-8);
	  if(DF<0){DF=0;}
	  Type DWDF=DW*DF;
	  // dis(y)+=exp(logF(conf.keyLogFsta(f,a-conf.minAge),y))/z*exp(logN(a-conf.minAge,y))*(Type(1.0)-exp(-z))*DWDF;
	  dis(y) += exp(logN(a-conf.minAge,y)) * mort.fleetCumulativeIncidence(a-conf.minAge,y,f) * DWDF;
	}
      }    
    }
  }  
  return dis;
}


template <class Type>
vector<Type> fsbFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, MortalitySet<Type>& mort){
  // TODO: Update to use mort
  int len=dat.catchMeanWeight.dim(0);
  int noFleets=conf.keyLogFsta.dim[0];
  array<Type> totF=totFFun(dat, conf, logF);
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

template<class Type>
Type fbari(dataSet<Type>& dat, confSet &conf, array<Type> &logF, int i, bool give_log = false){
  Type fbar = 0.0;
  array<Type> totF=totFFun(dat, conf, logF);
  for(int a=conf.fbarRange(0);a<=conf.fbarRange(1);a++){  
    fbar+=totF(a-conf.minAge,i);
  }
  fbar/=Type(conf.fbarRange(1)-conf.fbarRange(0)+1.0);
  if(give_log)
    return log(fbar);
  return fbar;
}
template <class Type>
vector<Type> fbarFun(dataSet<Type>& dat, confSet &conf, array<Type> &logF, bool give_log = false){
  int timeSteps=logF.dim[1];
  vector<Type> res(timeSteps);
  res.setZero();
  for(int y=0;y<timeSteps;y++){  
    res(y) = fbari(dat, conf, logF, y, give_log);
  }
  return res;
}
// #endif

template<class Type>
vector<Type> fbarByFleeti(confSet &conf, array<Type> &logF, int i, bool give_log = false){
  vector<Type> fbar(conf.keyLogFsta.dim[0]);
  fbar.setZero();
  for(int f = 0; f < fbar.size(); ++f){
    for(int a=conf.fbarRange(0);a<=conf.fbarRange(1);a++){  
      if(conf.keyLogFsta(f,a-conf.minAge)>(-1)){
	fbar(f) += exp(logF(conf.keyLogFsta(f,a-conf.minAge),i));
      }
    }
    fbar(f)/=Type(conf.fbarRange(1)-conf.fbarRange(0)+1.0);
  }
  if(give_log)
    return (vector<Type>)fbar.log();
  return fbar;
}
template <class Type>
matrix<Type> fbarByFleet(confSet &conf, array<Type> &logF, bool give_log = false){
  int timeSteps=logF.dim[1];
  matrix<Type> res(conf.keyLogFsta.dim[0], timeSteps);
  res.setZero();
  for(int y=0;y<timeSteps;y++){  
    res.col(y) = fbarByFleeti(conf, logF, y, give_log);
  }
  return res;
}
// #endif


template <class Type>
vector<Type> landFbarFun(dataSet<Type> &dat, confSet &conf, array<Type> &logF){
  int timeSteps=dat.landFrac.dim(0);
  array<Type> totF=totFFun(dat, conf, logF, 1);
  vector<Type> fbar(timeSteps);
  fbar.setZero();
  for(int y=0;y<timeSteps;y++){  
    for(int a=conf.fbarRange(0);a<=conf.fbarRange(1);a++){  
      fbar(y)+= totF(a-conf.minAge,y);
    }
    fbar(y)/=Type(conf.fbarRange(1)-conf.fbarRange(0)+1);
  }
  return fbar;
}

template <class Type>
vector<Type> disFbarFun(dataSet<Type> &dat, confSet &conf, array<Type> &logF){
  int timeSteps=dat.landFrac.dim(0);
  array<Type> totF=totFFun(dat, conf, logF, 2);
  vector<Type> fbar(timeSteps);
  fbar.setZero();
  for(int y=0;y<timeSteps;y++){  
    for(int a=conf.fbarRange(0);a<=conf.fbarRange(1);a++){  
      fbar(y)+= totF(a-conf.minAge, y);
    }
    fbar(y)/=Type(conf.fbarRange(1)-conf.fbarRange(0)+1);
  }
  return fbar;
}



#endif
