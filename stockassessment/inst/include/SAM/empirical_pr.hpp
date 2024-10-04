SAM_DEPENDS(convenience)
SAM_DEPENDS(logspace)
SAM_DEPENDS(define)
SAM_DEPENDS(incidence)


template<class Type>
Type empiricalSPR_i(dataSet<Type> &dat, confSet &conf, array<Type>& logN, MortalitySet<Type>& mort, int i, bool give_log DEFARG(= false))SOURCE({
    int stateDimN = logN.dim(0);
    int timeSteps = logN.dim(1);
    if(i+stateDimN > timeSteps)
      return R_NaReal;
    Type logR = logN(0,i);
    //Type logssb = R_NegInf;
    Type v = 0.0;
    for(int j = 0; j < stateDimN - 1; ++j){
      if(dat.propMat(i,j) > 0){
	Type lssbNew = logN(j,i+j) + log(mort.ssbSurvival_before(j,i+j)) + log(dat.propMat(i+j,j)) + log(dat.stockMeanWeight(i+j,j));
	v += exp(lssbNew); //logssb = loglogspace_add_SAM(logssb, lssbNew);
      }
    }
    Type lNPlusLast = R_NegInf;
    if(conf.maxAgePlusGroup(0)==1 && stateDimN > 1){
      int j = stateDimN - 1;
      Type v1 = logN(j-1,i+j-1) - mort.cumulativeHazard(j-1,i+j-1);
      Type v2 = logN(j,i+j-1) - mort.cumulativeHazard(j,i+j-1);
      // Proportion of logN plus group that "came from" cohort:
      Type lN = logN(j,i+j) - (v1 - logspace_add(v1,v2));
      Type lssbNew = lN + log(mort.ssbSurvival_before(j,i+j)) + log(dat.propMat(i+j,j)) + log(dat.stockMeanWeight(i+j,j));
      //logssb = logspace_add_SAM(logssb, lssbNew);
      v += exp(lssbNew);
      lNPlusLast = lN;
    }else{
      int j = stateDimN - 1;
      Type lssbNew = logN(j,i+j) + log(mort.ssbSurvival_before(j,i+j)) + log(dat.propMat(i+j,j)) + log(dat.stockMeanWeight(i+j,j));
      //logssb = logspace_add_SAM(logssb, lssbNew);
      v += exp(lssbNew);
      lNPlusLast = logN(j,i+j-1);
    }
    // Loop ahead in time
    for(int q = 1; q < 30; ++q){
      int j = stateDimN - 1;
      int indx = std::min(i+j-1 + q, (int)mort.cumulativeHazard.cols()-1);
      lNPlusLast -= mort.cumulativeHazard(j-1,indx-1);
      Type lssbNew = lNPlusLast + log(mort.ssbSurvival_before(j,indx)) + log(dat.propMat(indx,j)) + log(dat.stockMeanWeight(indx,j));
      //logssb = logspace_add_SAM(logssb, lssbNew);
      v += exp(lssbNew);
    }
    //Type logSPR = logssb - logR;
    Type logSPR = log(v) - logR;
    if(give_log)
      return logSPR;
    return exp(logSPR);  
  })

SAM_SPECIALIZATION(double empiricalSPR_i(dataSet<double>&, confSet&, array<double>&, MortalitySet<double>&, int, bool));
SAM_SPECIALIZATION(TMBad::ad_aug empiricalSPR_i(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>&, int, bool));

template<class Type>
vector<Type> empiricalSPR(dataSet<Type> &dat, confSet &conf, array<Type>& logN, MortalitySet<Type>& mort, bool give_log)SOURCE({
    int timeSteps = logN.dim(1);
    vector<Type> r(timeSteps);
    r.setConstant(R_NaReal);
    for(int i = 0; i < timeSteps; ++i)
      r(i) = empiricalSPR_i(dat, conf, logN, mort, i, give_log);
    return r;
  })

  SAM_SPECIALIZATION(vector<double> empiricalSPR(dataSet<double>&, confSet&, array<double>&, MortalitySet<double>&, bool));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> empiricalSPR(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>&, bool));


template<class Type>
Type empiricalYPR_i(dataSet<Type> &dat, confSet &conf, array<Type>& logN, MortalitySet<Type>& mort, int CT, int i, bool give_log DEFARG(= false))SOURCE({
  int stateDimN = logN.dim(0);
  int timeSteps = dat.landFrac.dim(0); //logN.dim(1);
  int noFleets=conf.keyLogFsta.dim(0);
  if(i+stateDimN > timeSteps)
    return R_NaReal;
  Type logR = logN(0,i);
  //Type logCatch = R_NegInf;
  Type Catch = 0.0;
  for(int j = 0; j < stateDimN - 1; ++j){
    for(int f = 0; f < noFleets; ++f){
      if(dat.fleetTypes(f) == 0){ // Only catch fleets
	Type logC = logN(j,i+j) + mort.logFleetSurvival_before(j,i+j,f) + log(mort.fleetCumulativeIncidence(j,i+j,f));
	if(CT == 0){		// Catch
	  logC += log(dat.catchMeanWeight(i+j, j, f));
	}else if(CT == 1){		// Landing
	  Type LW=dat.landMeanWeight(i+j,j,f);
	  Type LF=dat.landFrac(i+j,j,f);
	  if(LW > 0 && LF > 0){
	    logC += log(LW) + log(LF);
	  }else{
	    logC = R_NegInf;
	  }
	}else if(CT == 2){		// Discard
	  Type DW=dat.disMeanWeight(i+j,j, f);
	  Type DF=1.0 - (dat.landFrac(i+j,j, f) - 1e-8);
	  if(DW > 0 && DF > 0){
	    logC += log(DW) + log(DF);
	  }else{
	    logC = R_NegInf;
	  }
	}
	//logCatch = logspace_add_SAM(logCatch,logC);
	Catch += exp(logC);
      }
    }
  }
  int j = stateDimN - 1;
  Type lN = logN(j,i+j);
  if(conf.maxAgePlusGroup(0)==1 && stateDimN > 1){
    Type v1 = logN(j-1,i+j-1) - mort.cumulativeHazard(j-1,i+j-1);
    Type v2 = logN(j,i+j-1) - mort.cumulativeHazard(j,i+j-1);
    // Proportion of logN plus group that "came from" cohort:
    lN -= v1 - logspace_add(v1,v2);
  }
  for(int f = 0; f < noFleets; ++f){
    if(dat.fleetTypes(f) == 0){ // Only catch fleets
      Type logC = lN + mort.logFleetSurvival_before(j,i+j,f) + log(mort.fleetCumulativeIncidence(j,i+j,f));
      if(CT == 0){		// Catch
	logC += log(dat.catchMeanWeight(i+j, j, f));
      }else if(CT == 1){		// Landing
	Type LW=dat.landMeanWeight(i+j,j,f);
	Type LF=dat.landFrac(i+j,j,f);
	if(LW > 0 && LF > 0){
	  logC += log(LW) + log(LF);
	}else{
	  logC = R_NegInf;
	}
      }else if(CT == 2){		// Discard
	Type DW=dat.disMeanWeight(i+j,j, f);
	Type DF=1.0 - (dat.landFrac(i+j,j, f) - 1e-8);
	if(DW > 0 && DF > 0){
	  logC += log(DW) + log(DF);
	}else{
	  logC = R_NegInf;
	}
      }
      //logCatch = logspace_add_SAM(logCatch,logC);
      Catch += exp(logC);
    }
  }
  // Loop ahead in time
  Type lNPlusLast = lN;
  for(int q = 1; q < 30; ++q){
    int j = stateDimN - 1;
    int indx = std::min(i+j-1 + q, timeSteps-1);
    lNPlusLast -= mort.cumulativeHazard(j-1,indx-1);
    for(int f = 0; f < noFleets; ++f){
      if(dat.fleetTypes(f) == 0){ // Only catch fleets
	Type logC = lNPlusLast + mort.logFleetSurvival_before(j,indx,f) + log(mort.fleetCumulativeIncidence(j,indx,f));
	if(CT == 0){		// Catch
	  logC += log(dat.catchMeanWeight(indx, j, f));
	}else if(CT == 1){		// Landing
	  Type LW=dat.landMeanWeight(indx,j,f);
	  Type LF=dat.landFrac(indx,j,f);
	  if(LW > 0 && LF > 0){
	    logC += log(LW) + log(LF);
	  }else{
	    logC = R_NegInf;
	  }
	}else if(CT == 2){		// Discard
	  Type DW=dat.disMeanWeight(indx,j, f);
	  Type DF=1.0 - (dat.landFrac(indx,j, f) - 1e-8);
	  if(DW > 0 && DF > 0){
	    logC += log(DW) + log(DF);
	  }else{
	    logC = R_NegInf;
	  }
	}
	//logCatch = logspace_add_SAM(logCatch,logC);
	Catch += exp(logC);
      }
    }
  }
  //Type logYPR = logCatch - logR;
  Type logYPR = log(Catch) - logR;
  if(give_log)
    return logYPR;
  return exp(logYPR);  
  })


SAM_SPECIALIZATION(double empiricalYPR_i(dataSet<double>&, confSet&, array<double>&, MortalitySet<double>&, int, int, bool));
SAM_SPECIALIZATION(TMBad::ad_aug empiricalYPR_i(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>&, int, int, bool));



template<class Type>
vector<Type> empiricalYPR(dataSet<Type> &dat, confSet &conf, array<Type>& logN, MortalitySet<Type>& mort, int CT, bool give_log)SOURCE({
    int timeSteps = dat.landFrac.dim(0); //logN.dim(1);
    vector<Type> r(timeSteps);
    r.setConstant(R_NaReal);
    for(int i = 0; i < timeSteps; ++i)
      r(i) = empiricalYPR_i(dat, conf, logN, mort, CT, i, give_log);
    return r;
  })

  SAM_SPECIALIZATION(vector<double> empiricalYPR(dataSet<double>&, confSet&, array<double>&, MortalitySet<double>&, int, bool));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> empiricalYPR(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>&, int, bool));


