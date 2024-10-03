SAM_DEPENDS(define)
SAM_DEPENDS(incidence)
SAM_DEPENDS(recruitment)
SAM_DEPENDS(derived)
SAM_DEPENDS(convenience)

template <class Type>
vector<Type> predNFun(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, array<Type> &logF, Recruitment<Type> &recruit, MortalitySet<Type>& mort, int i)SOURCE({
  int stateDimN=logN.dim[0];
  //array<Type> totF=totFFun(dat,conf, logF);
  vector<Type> predN(stateDimN);
  //predN.setZero();
  predN.setConstant(R_NegInf);
  Type logThisSSB=Type(R_NegInf);

  if((i-conf.minAge)>=0){
    //logThisSSB=ssbi(dat,conf,logN,logF,mort,i-conf.minAge, true);
    logThisSSB=spawningPotentiali(dat,conf,par,logN,logF,mort,i-conf.minAge, true);
  }else{
    //logThisSSB=ssbi(dat,conf,logN,logF,mort,0, true); // use first in beginning
    logThisSSB=spawningPotentiali(dat,conf,par,logN,logF,mort,0, true); // use first in beginning       
  }

  Type lastLogR = R_NaReal;
  if(i > 0)
    lastLogR = logN(0,i-1);    
  predN(0) = recruit(logThisSSB, lastLogR, dat.years(0) + i); // dat.years(0) + i is needed for forecast

  if(par.rec_transphi.size() > 0){
    vector<Type> phi = logitroots2ARpar(par.rec_transphi);
    for(int j = 1; j <= par.rec_transphi.size(); ++j){
      //Type logThatSSB = ssbi(dat,conf,logN,logF,mort,std::max(i-conf.minAge-j,0), true);
      Type logThatSSB = spawningPotentiali(dat,conf,par,logN,logF,mort,std::max(i-conf.minAge-j,0), true);
      Type thatMu = recruit(logThatSSB, logN(0,std::max(i-1-j,0)), dat.years(0) + i - j);
      predN(0) += phi(j-1) * (logN(0,std::max(i-j,0)) - thatMu);
    }
  }
  
  switch(conf.logNMeanAssumption(0)){
  case 0:			// Median on natural scale
    predN(0) += 0.0;
    break;
  case 1:			// Mean on natural scale
    predN(0) -= 0.5 * exp(2.0 * par.logSdLogN(conf.keyVarLogN(0)));
    break;
  case 2:			// Mode on natural scale
    predN(0) += exp(2.0 * par.logSdLogN(conf.keyVarLogN(0)));
    break;
  default:
      Rf_error("logNMeanCorrection not implemented.");
    break;    
  }
		       
  for(int j=1; j<stateDimN; ++j){
    predN(j)=logN(j-1,i-1) - mort.cumulativeHazard(j-1,i-1); //totF(j-1,i-1)-dat.natMor(i-1,j-1); 
  }
  if(conf.maxAgePlusGroup(0)==1){// plusgroup adjustment if catches need them 
    Type v1 = predN(stateDimN-1); // Already updated above
    Type v2 = logN(stateDimN-1,i-1) - mort.cumulativeHazard(stateDimN-1,i-1); //totF(stateDimN-1,i-1) - dat.natMor(i-1,stateDimN-1); // Remaining in plus group from last year
    predN(stateDimN-1) = logspace_add_SAM(v1,v2);
  }

  for(int j=1; j<stateDimN; ++j){
    switch(conf.logNMeanAssumption(1)){
    case 0:			// Median on natural scale
      predN(j) += 0.0;
      break;
    case 1:			// Mean on natural scale
      predN(j) -= 0.5 * exp(2.0 * par.logSdLogN(conf.keyVarLogN(j)));
      break;
    case 2:			// Mode on natural scale
      predN(j) += exp(2.0 * par.logSdLogN(conf.keyVarLogN(j)));
      break;
    default:
      Rf_error("logNMeanCorrection not implemented.");
      break;    
    }
  }
  
  return predN;  
  });

SAM_SPECIALIZATION(vector<double> predNFun(dataSet<double>&, confSet&, paraSet<double>&, array<double>&, array<double>&, Recruitment<double>&, MortalitySet<double>&, int));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> predNFun(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, Recruitment<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>&, int));
