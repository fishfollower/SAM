SAM_DEPENDS(convenience)
SAM_DEPENDS(logspace)
SAM_DEPENDS(define)
SAM_DEPENDS(incidence)


template <class Type>
int yearsPFun(confSet &conf, dataSet<Type> &dat) SOURCE({
  
  int noFleets=conf.keyLogFsta.dim[0];
  Type minYear = dat.years(dat.noYears-1);
  Type maxYear = 0;
  int noYearsLAI = -1;
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
  })


SAM_SPECIALIZATION(int yearsPFun<double>(confSet&, dataSet<double>&));
SAM_SPECIALIZATION(int yearsPFun<TMBad::ad_aug>(confSet&, dataSet<TMBad::ad_aug>&));

template <class Type>
array<Type> scalePFun(confSet &conf, dataSet<Type> &dat, array<Type> &logP)SOURCE({
    // int noYearsLAI = yearsPFun(conf,dat);
    int nlogP = logP.dim[0]+1;
    array<Type> logPS(nlogP,logP.dim[1]);
    logPS.setZero();
    if(logPS.dim[0] == 1 || logPS.dim[1] == 0)
      return logPS;
    
    Type totProp;
    for(int j=0;j<logP.dim[1];j++){
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
  })

SAM_SPECIALIZATION(array<double> scalePFun(confSet&, dataSet<double>&, array<double>&));
SAM_SPECIALIZATION(array<TMBad::ad_aug> scalePFun(confSet&, dataSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&));


template <class Type>
vector<Type> scaleWeekFun(paraSet<Type> &par, dataSet<Type> &dat, array<Type> &logP)SOURCE({

  int nlogP = logP.dim[0]+1;
  int maxLAIsurv = par.logAlphaSCB.size()+nlogP;
  vector<Type> varAlphaSCB(maxLAIsurv);
  varAlphaSCB.setZero();
  if(logP.dim[0] == 0 || logP.dim[1] == 0)
      return varAlphaSCB;
  
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
  })

SAM_SPECIALIZATION(vector<double> scaleWeekFun(paraSet<double>&, dataSet<double>&, array<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> scaleWeekFun(paraSet<TMBad::ad_aug>&, dataSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&));



template <class Type>
Type ssbi(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, MortalitySet<Type>& mort, int i, bool give_log DEFARG(= false))SOURCE({
  int stateDimN=logN.dim[0];
  Type logssb = R_NegInf;
  for(int j=0; j<stateDimN; ++j){
    if(dat.propMat(i,j) > 0){
      Type lssbNew = logN(j,i) + log(mort.ssbSurvival_before(j,i)) + log(dat.propMat(i,j)) + log(dat.stockMeanWeight(i,j));
      logssb = logspace_add_SAM(logssb, lssbNew);
    }
  }
  if(give_log)
    return logssb;
  return exp(logssb);
  })

SAM_SPECIALIZATION(double ssbi(dataSet<double>&, confSet&, array<double>&, array<double>&, MortalitySet<double>&, int, bool));
SAM_SPECIALIZATION(TMBad::ad_aug ssbi(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>&, int, bool));

template <class Type>
vector<Type> ssbFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF,MortalitySet<Type>& mort, bool give_log DEFARG(= false))SOURCE({
  int timeSteps=logF.dim[1];
  array<Type> totF=totFFun(dat, conf, logF);
  vector<Type> ssb(timeSteps);
  ssb.setConstant(R_NegInf);
  for(int i=0;i<timeSteps;i++){
    ssb(i)=ssbi(dat,conf,logN,logF,mort,i, give_log);
  }
  return ssb;
  });

SAM_SPECIALIZATION(vector<double> ssbFun(dataSet<double>&, confSet&, array<double>&, array<double>&, MortalitySet<double>&, bool));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> ssbFun(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>&, bool));


template <class Type>
matrix<Type> catchFunAge(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, MortalitySet<Type>& mort, bool give_log DEFARG(= false))SOURCE({
  int len=dat.landFrac.dim(0);
  int noFleets=conf.keyLogFsta.dim(0);
  array<Type> totF=totFFun(dat, conf, logF);
  matrix<Type> cat(conf.maxAge - conf.minAge + 1, len);
  cat.setZero();
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      for(int f=0; f<noFleets;f++){
	if(dat.fleetTypes(f) == 0){ // Only catch fleets
	  // cat(a-conf.minAge, y) += exp(logN(a-conf.minAge, y)) * mort.fleetCumulativeIncidence(a-conf.minAge,y, f) * dat.catchMeanWeight(y, a-conf.minAge, f);
	  cat(a-conf.minAge, y) += exp(logN(a-conf.minAge,y) + mort.logFleetSurvival_before(a-conf.minAge,y,f) + log(mort.fleetCumulativeIncidence(a,y,f))) * dat.catchMeanWeight(y, a-conf.minAge, f);
	}
      }
    }
  }
  if(give_log)
    return cat.array().log().matrix();  
  return cat;
  });

SAM_SPECIALIZATION(matrix<double> catchFunAge(dataSet<double>&, confSet&, array<double>&, array<double>&, MortalitySet<double>&, bool));
SAM_SPECIALIZATION(matrix<TMBad::ad_aug> catchFunAge(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>&, bool));


template <class Type>
array<Type> catchByFleetFunAge(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, MortalitySet<Type>& mort)SOURCE({
  int len=dat.landFrac.dim(0);
  int noFleets=conf.keyLogFsta.dim(0);
  array<Type> totF=totFFun(dat, conf, logF);
  array<Type> cat(conf.maxAge - conf.minAge + 1, len, noFleets);
  cat.setZero();
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      for(int f=0; f<noFleets;f++){
	if(dat.fleetTypes(f) == 0){ // Only catch fleets
	  // cat(a-conf.minAge, y) += exp(logN(a-conf.minAge, y)) * mort.fleetCumulativeIncidence(a-conf.minAge,y, f) * dat.catchMeanWeight(y, a-conf.minAge, f);
	  cat(a-conf.minAge, y, f) += exp(logN(a-conf.minAge,y) + mort.logFleetSurvival_before(a-conf.minAge,y,f) + log(mort.fleetCumulativeIncidence(a,y,f))) * dat.catchMeanWeight(y, a-conf.minAge, f);
	}
      }
    }
  }
  return cat;
  });

SAM_SPECIALIZATION(array<double> catchByFleetFunAge(dataSet<double>&, confSet&, array<double>&, array<double>&, MortalitySet<double>&));
SAM_SPECIALIZATION(array<TMBad::ad_aug> catchByFleetFunAge(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>&));



template <class Type>
array<Type> catchByFleetFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, MortalitySet<Type>& mort)SOURCE({
  int len=dat.catchMeanWeight.dim(0);
  int noFleets=conf.keyLogFsta.dim[0];
  array<Type> cat(len,noFleets);
  cat.setZero();
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      for(int f=0; f<noFleets;f++){
	if(dat.fleetTypes(f) == 0){ // Only catch fleets
  	  cat(y,f) += exp(logN(a-conf.minAge,y)) * mort.fleetCumulativeIncidence(a-conf.minAge,y,f) * dat.catchMeanWeight(y,a-conf.minAge,f);
	}
      }
    }
  }
  return cat;
  })

SAM_SPECIALIZATION(array<double> catchByFleetFun(dataSet<double>&, confSet&, array<double>&, array<double>&, MortalitySet<double>&));
SAM_SPECIALIZATION(array<TMBad::ad_aug> catchByFleetFun(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>&));


template <class Type>
vector<Type> catchFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF,MortalitySet<Type>& mort, bool give_log DEFARG(= false))SOURCE({
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
  });

SAM_SPECIALIZATION(vector<double> catchFun(dataSet<double>&, confSet&, array<double>&, array<double>&, MortalitySet<double>&, bool));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> catchFun(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>&, bool));


template <class Type>
vector<Type> varLogCatchFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, paraSet<Type>& par, MortalitySet<Type>& mort)SOURCE({
  int len=dat.catchMeanWeight.dim(0);
  int noFleets=conf.keyLogFsta.dim(0);
  vector<Type> cat(len);
  cat.setZero();
  vector<Type> varLogCat(len);
  varLogCat.setZero();
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      for(int f=0; f<noFleets;f++){
	if(dat.fleetTypes(f) == 0 && conf.keyVarObs(f,a-conf.minAge) > (-1)){ // Only catch fleets
          Type CW=dat.catchMeanWeight(y,a-conf.minAge,f);
  	  Type Ca = exp(logN(a-conf.minAge,y)) * mort.fleetCumulativeIncidence(a-conf.minAge,y,f);
          cat(y)+=Ca*CW;
          varLogCat(y)+=exp(2.0*par.logSdLogObs(conf.keyVarObs(f,a-conf.minAge)))*CW*CW*Ca*Ca;
	}
      }
    }
    varLogCat(y)/=cat(y)*cat(y);
  }
  return varLogCat;
  })

SAM_SPECIALIZATION(vector<double> varLogCatchFun(dataSet<double>&, confSet&, array<double>&, array<double>&, paraSet<double>&, MortalitySet<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> varLogCatchFun(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&,paraSet<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>&));


template <class Type>
vector<Type> landFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, MortalitySet<Type>& mort)SOURCE({
  int len=dat.landMeanWeight.dim(0);
  int noFleets=conf.keyLogFsta.dim(0);
  vector<Type> land(len);
  land.setZero();
  Type LW=0;
  Type LF=0;
  Type LWLF=0;
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){
      for(int f=0; f<noFleets;f++){
	if(dat.fleetTypes(f) == 0){ // Only catch fleets
          LW=dat.landMeanWeight(y,a-conf.minAge,f);
          if(LW<0){LW=0;}
          LF=dat.landFrac(y,a-conf.minAge,f);
          if(LF<0){LF=0;}
          LWLF=LW*LF;
	  land(y) += exp(logN(a-conf.minAge,y)) * mort.fleetCumulativeIncidence(a-conf.minAge,y,f) * LWLF;
	}
      }
    }
  }
  return land;
  });


SAM_SPECIALIZATION(vector<double> landFun(dataSet<double>&, confSet&, array<double>&, array<double>&, MortalitySet<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> landFun(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>&));


template <class Type>
vector<Type> varLogLandFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, paraSet<Type>& par, MortalitySet<Type>& mort)SOURCE({
  int len=dat.landMeanWeight.dim(0);
  int noFleets=conf.keyLogFsta.dim[0];
  vector<Type> land(len);
  land.setZero();
  vector<Type> varLogLand(len);
  varLogLand.setZero();
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){
      for(int f=0; f<noFleets;f++){
	if(dat.fleetTypes(f) == 0 && conf.keyVarObs(f,a-conf.minAge) > (-1)){ // Only catch fleets
          Type LW=dat.landMeanWeight(y,a-conf.minAge,f);
          if(LW<0){LW=0;}
          Type LF=dat.landFrac(y,a-conf.minAge,f);
          if(LF<0){LF=0;}
          Type LWLF=LW*LF;
	  Type Ca = exp(logN(a-conf.minAge,y)) * mort.fleetCumulativeIncidence(a-conf.minAge,y,f);
          land(y)+=Ca*LWLF;
          varLogLand(y)+=exp(2.0*par.logSdLogObs(conf.keyVarObs(f,a-conf.minAge)))*LWLF*LWLF*Ca*Ca;
	}
      }
    }
    varLogLand(y)/=land(y)*land(y);
  }
  return varLogLand;
  });


SAM_SPECIALIZATION(vector<double> varLogLandFun(dataSet<double>&, confSet&, array<double>&, array<double>&, paraSet<double>&, MortalitySet<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> varLogLandFun(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, paraSet<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>&));



template <class Type>
vector<Type> disFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, MortalitySet<Type>& mort) SOURCE({
  int len=dat.disMeanWeight.dim(0);
  int noFleets=conf.keyLogFsta.dim[0];
  vector<Type> dis(len);
  dis.setZero();
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){
      for(int f=0; f<noFleets;f++){
	if(dat.fleetTypes(f) == 0){ // Only catch fleets
	  Type DW=dat.disMeanWeight(y,a-conf.minAge, f);
	  if(DW<0){DW=0;}
	  Type DF=1.0 - (dat.landFrac(y,a-conf.minAge, f) - 1e-8);
	  if(DF<0){DF=0;}
	  Type DWDF=DW*DF;
	  dis(y) += exp(logN(a-conf.minAge,y)) * mort.fleetCumulativeIncidence(a-conf.minAge,y,f) * DWDF;
	}
      }    
    }
  }  
  return dis;
  });


SAM_SPECIALIZATION(vector<double> disFun(dataSet<double>&, confSet&, array<double>&, array<double>&, MortalitySet<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> disFun(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>&));



template <class Type>
vector<Type> fsbFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, MortalitySet<Type>& mort)SOURCE({
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
  })


SAM_SPECIALIZATION(vector<double> fsbFun(dataSet<double>&, confSet&, array<double>&, array<double>&, MortalitySet<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> fsbFun(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>&));


template <class Type>
vector<Type> tsbFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN)SOURCE({
  int timeSteps=logN.dim[1];
  vector<Type> tsb(timeSteps);
  tsb.setZero();  
  for(int y=0;y<timeSteps;y++){  
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      tsb(y)+=exp(logN(a-conf.minAge,y))*dat.stockMeanWeight(y,a-conf.minAge);
    }
  }
  return tsb;
  })

  
SAM_SPECIALIZATION(vector<double> tsbFun(dataSet<double>&, confSet&, array<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> tsbFun(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&));


template <class Type>
vector<Type> rFun(array<Type> &logN)SOURCE({
  int timeSteps=logN.dim[1];
  vector<Type> R(timeSteps);
  R.setZero();
  for(int y=0;y<timeSteps;y++){  
    R(y)=exp(logN(0,y));
  }
  return R;
  });

SAM_SPECIALIZATION(vector<double> rFun(array<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> rFun(array<TMBad::ad_aug>&));

template<class Type>
Type fbari(dataSet<Type>& dat, confSet &conf, array<Type> &logF, int i, bool give_log DEFARG(= false))SOURCE({
  Type fbar = 0.0;
  array<Type> totF=totFFun(dat, conf, logF);
  for(int a=conf.fbarRange(0);a<=conf.fbarRange(1);a++){  
    fbar+=totF(a-conf.minAge,i);
  }
  fbar/=Type(conf.fbarRange(1)-conf.fbarRange(0)+1.0);
  if(give_log)
    return log(fbar);
  return fbar;
  });

SAM_SPECIALIZATION(double fbari(dataSet<double>&, confSet&, array<double>&, int, bool));
SAM_SPECIALIZATION(TMBad::ad_aug fbari(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, int, bool));


template <class Type>
vector<Type> fbarFun(dataSet<Type>& dat, confSet &conf, array<Type> &logF, bool give_log DEFARG(= false))SOURCE({
  int timeSteps=logF.dim[1];
  vector<Type> res(timeSteps);
  res.setZero();
  for(int y=0;y<timeSteps;y++){  
    res(y) = fbari(dat, conf, logF, y, give_log);
  }
  return res;
    });

SAM_SPECIALIZATION(vector<double> fbarFun(dataSet<double>&, confSet&, array<double>&, bool));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> fbarFun(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, bool));


template<class Type>
vector<Type> fbarByFleeti(confSet &conf, array<Type> &logF, int i, bool give_log DEFARG(= false))SOURCE({
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
    })

SAM_SPECIALIZATION(vector<double> fbarByFleeti(confSet&, array<double>&, int, bool));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> fbarByFleeti(confSet&, array<TMBad::ad_aug>&, int, bool));



template <class Type>
matrix<Type> fbarByFleet(confSet &conf, array<Type> &logF, bool give_log DEFARG(= false))SOURCE({
  int timeSteps=logF.dim[1];
  matrix<Type> res(conf.keyLogFsta.dim[0], timeSteps);
  res.setZero();
  for(int y=0;y<timeSteps;y++){  
    res.col(y) = fbarByFleeti(conf, logF, y, give_log);
  }
  return res;
  });

  
SAM_SPECIALIZATION(matrix<double> fbarByFleet(confSet&, array<double>&, bool));
SAM_SPECIALIZATION(matrix<TMBad::ad_aug> fbarByFleet(confSet&, array<TMBad::ad_aug>&, bool));



template <class Type>
  vector<Type> landFbarFun(dataSet<Type> &dat, confSet &conf, array<Type> &logF)SOURCE({
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
    })

  
SAM_SPECIALIZATION(vector<double> landFbarFun(dataSet<double>&, confSet&, array<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> landFbarFun(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&));


template <class Type>
  vector<Type> disFbarFun(dataSet<Type> &dat, confSet &conf, array<Type> &logF)SOURCE({
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
    })

  
SAM_SPECIALIZATION(vector<double> disFbarFun(dataSet<double>&, confSet&, array<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> disFbarFun(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&));


