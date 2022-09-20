#pragma once
#ifndef SAM_PREDN_HPP
#define SAM_PREDN_HPP

template <class Type>
vector<Type> predNFun(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, array<Type> &logF, Recruitment<Type> &recruit, MortalitySet<Type>& mort, int i){
  int stateDimN=logN.dim[0];
  //array<Type> totF=totFFun(dat,conf, logF);
  vector<Type> predN(stateDimN);
  //predN.setZero();
  predN.setConstant(R_NegInf);
  Type logThisSSB=Type(R_NegInf);

  if((i-conf.minAge)>=0){
    logThisSSB=ssbi(dat,conf,logN,logF,mort,i-conf.minAge, true);
  }else{
    logThisSSB=ssbi(dat,conf,logN,logF,mort,0, true); // use first in beginning       
  }

  Type lastLogR = R_NaReal;
  if(i > 0)
    lastLogR = logN(0,i-1);    
  predN(0) = recruit(logThisSSB, lastLogR, dat.years(0) + i); // dat.years(0) + i is needed for forecast

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
    predN(j)=logN(j-1,i-1) - mort.totalZ(j-1,i-1); //totF(j-1,i-1)-dat.natMor(i-1,j-1); 
  }
  if(conf.maxAgePlusGroup(0)==1){// plusgroup adjustment if catches need them 
    Type v1 = predN(stateDimN-1); // Already updated above
    Type v2 = logN(stateDimN-1,i-1) - mort.totalZ(stateDimN-1,i-1); //totF(stateDimN-1,i-1) - dat.natMor(i-1,stateDimN-1); // Remaining in plus group from last year
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
}

#endif

