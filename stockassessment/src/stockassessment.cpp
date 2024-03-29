//  --------------------------------------------------------------------------
// Copyright (c) 2014, Anders Nielsen <an@aqua.dtu.dk>,    
// Casper Berg <cbe@aqua.dtu.dk>, Kasper Kristensen <kkr@aqua.dtu.dk>,
// Mollie Brooks <molbr@aqua.dtu.dk>,
// and Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//   * Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//   * Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in the
//     documentation and/or other materials provided with the distribution.
//   * Neither the name of the assessment tool SAM nor the
//     names of its contributors may be used to endorse or promote products
//     derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
// ARE DISCLAIMED. IN NO EVENT SHALL ANDERS NIELSEN, CASPER BERG OR KASPER 
// KRISTENSEN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//  --------------------------------------------------------------------------

// R_init_stockassessment is now defined in main.cpp
//#define TMB_LIB_INIT R_init_stockassessment
//#define TMB_SAFEBOUNDS
// #define TMB_MAX_ORDER 4
// #include "TMB.h"

#define WITH_SAM_LIB
#include "SAM.h"


template<class Type>
Type objective_function<Type>::operator() ()
{

  dataSet<Type> dataset;
  DATA_INTEGER(noFleets); dataset.noFleets=noFleets; 
  DATA_IVECTOR(fleetTypes); dataset.fleetTypes=fleetTypes;    
  DATA_VECTOR(sampleTimes); dataset.sampleTimes=sampleTimes; 
  DATA_INTEGER(noYears); dataset.noYears=noYears; 
  DATA_VECTOR(years); dataset.years=years; 
  DATA_IVECTOR(minAgePerFleet); dataset.minAgePerFleet=minAgePerFleet; 
  DATA_IVECTOR(maxAgePerFleet); dataset.maxAgePerFleet=maxAgePerFleet; 
  DATA_INTEGER(nobs); dataset.nobs=nobs; 
  DATA_IARRAY(idx1); dataset.idx1=idx1;     // minimum index of obs by fleet x year
  DATA_IARRAY(idx2); dataset.idx2=idx2;     // maximum index of obs by fleet x year
  DATA_IARRAY(idxCor); dataset.idxCor=idxCor;    
  DATA_IVECTOR(minWeek); dataset.minWeek=minWeek;
  DATA_IVECTOR(maxWeek); dataset.maxWeek=maxWeek;
  DATA_IARRAY(aux); dataset.aux=aux; 
  DATA_VECTOR(logobs); dataset.logobs=logobs; 
  DATA_VECTOR(weight); dataset.weight=weight; 
  DATA_VECTOR_INDICATOR(keep, logobs); //dataset.keep=keep; 
  DATA_ARRAY(propMat); dataset.propMat=propMat; 
  DATA_ARRAY(stockMeanWeight); dataset.stockMeanWeight=stockMeanWeight;  
  DATA_ARRAY(catchMeanWeight); dataset.catchMeanWeight=catchMeanWeight; 
  DATA_ARRAY(natMor); dataset.natMor=natMor; 
  DATA_ARRAY(landFrac); dataset.landFrac=landFrac; 
  DATA_ARRAY(disMeanWeight); dataset.disMeanWeight=disMeanWeight; 
  DATA_ARRAY(landMeanWeight); dataset.landMeanWeight=landMeanWeight; 
  DATA_ARRAY(propF); dataset.propF=propF; 
  DATA_ARRAY(propM); dataset.propM=propM; 
  DATA_STRUCT(corList,listMatrixFromR); dataset.corList=corList; //Include correlation structures
  DATA_IARRAY(sumKey); dataset.sumKey=sumKey; 

  DATA_STRUCT(forecast, forecastSet);
  DATA_STRUCT(referencepoints, referencepointList);

  confSet confset;
  DATA_INTEGER(minAge); confset.minAge=minAge; 
  DATA_INTEGER(maxAge); confset.maxAge=maxAge; 
  DATA_IVECTOR(maxAgePlusGroup); confset.maxAgePlusGroup=maxAgePlusGroup; 
  DATA_IARRAY(keyLogFsta); confset.keyLogFsta=keyLogFsta; 
  DATA_IVECTOR(corFlag); confset.corFlag=corFlag; 
  DATA_IARRAY(keyLogFpar); confset.keyLogFpar=keyLogFpar; 
  DATA_IARRAY(keyQpow); confset.keyQpow=keyQpow; 
  DATA_IARRAY(keyVarF); confset.keyVarF=keyVarF; 
  DATA_IVECTOR(keyVarLogN); confset.keyVarLogN=keyVarLogN;  
  DATA_IVECTOR(keyVarLogP); confset.keyVarLogP=keyVarLogP; //coupling of SD of RWs of logP
  DATA_IARRAY(keyVarObs); confset.keyVarObs=keyVarObs; 
  DATA_FACTOR(obsCorStruct); confset.obsCorStruct=obsCorStruct;  
  DATA_IARRAY(keyCorObs); confset.keyCorObs=keyCorObs; 
  DATA_INTEGER(stockRecruitmentModelCode); confset.stockRecruitmentModelCode=stockRecruitmentModelCode; 
  DATA_VECTOR(constRecBreaks); vector<double> constRecBreaksDouble(constRecBreaks.size()); for(int i=0; i<constRecBreaks.size(); ++i){constRecBreaksDouble(i)=asDouble(constRecBreaks(i));} confset.constRecBreaks=constRecBreaksDouble; 
  DATA_INTEGER(noScaledYears); confset.noScaledYears=noScaledYears; 
  DATA_IVECTOR(keyScaledYears); confset.keyScaledYears=keyScaledYears; 
  DATA_IMATRIX(keyParScaledYA); confset.keyParScaledYA=keyParScaledYA; 
  DATA_IVECTOR(fbarRange); confset.fbarRange=fbarRange; 
  DATA_IVECTOR(keyBiomassTreat); confset.keyBiomassTreat=keyBiomassTreat;   
  DATA_IVECTOR(simFlag); confset.simFlag=simFlag;  //1 means simulations should not redo F and N
  DATA_INTEGER(resFlag); confset.resFlag=resFlag;  
  DATA_FACTOR(obsLikelihoodFlag); confset.obsLikelihoodFlag=obsLikelihoodFlag; 
  DATA_IVECTOR(fixVarToWeight);
  vector<int> fixVarToWeightInt(noFleets);
  if(fixVarToWeight.size()==1){
    for(int i=0; i<fixVarToWeightInt.size(); ++i){
      fixVarToWeightInt(i)=fixVarToWeight(0);  
    }
  }else{
    fixVarToWeightInt=fixVarToWeight;  
  }
  confset.fixVarToWeight=fixVarToWeightInt;
  
  DATA_SCALAR(fracMixF); confset.fracMixF=asDouble(fracMixF); 
  DATA_VECTOR(fracMixN);
  vector<double> fracMixNDouble(maxAge-minAge+1);
  for(int i=0; i<fracMixNDouble.size(); ++i){
    if(fracMixN.size()==1){
      fracMixNDouble(i)=asDouble(fracMixN(0));
    }else{
      fracMixNDouble(i)=asDouble(fracMixN(i));
    }
  }
  confset.fracMixN=fracMixNDouble;   
  DATA_VECTOR(fracMixObs); vector<double> fracMixObsDouble(fracMixObs.size()); for(int i=0; i<fracMixObs.size(); ++i){fracMixObsDouble(i)=asDouble(fracMixObs(i));} confset.fracMixObs=fracMixObsDouble; 
  DATA_IARRAY(predVarObsLink); confset.predVarObsLink=predVarObsLink;
  DATA_INTEGER(stockWeightModel); confset.stockWeightModel=stockWeightModel;
  DATA_IVECTOR(keyStockWeightMean); confset.keyStockWeightMean=keyStockWeightMean;
  DATA_IVECTOR(keyStockWeightObsVar); confset.keyStockWeightObsVar=keyStockWeightObsVar; 
  DATA_INTEGER(catchWeightModel); confset.catchWeightModel=catchWeightModel;
  DATA_IMATRIX(keyCatchWeightMean); confset.keyCatchWeightMean=keyCatchWeightMean;
  DATA_IMATRIX(keyCatchWeightObsVar); confset.keyCatchWeightObsVar=keyCatchWeightObsVar; 
  DATA_INTEGER(matureModel); confset.matureModel=matureModel;
  DATA_IVECTOR(keyMatureMean); confset.keyMatureMean=keyMatureMean;
  DATA_INTEGER(mortalityModel); confset.mortalityModel=mortalityModel;
  DATA_IVECTOR(keyMortalityMean); confset.keyMortalityMean=keyMortalityMean;
  DATA_IVECTOR(keyMortalityObsVar); confset.keyMortalityObsVar=keyMortalityObsVar; 
  DATA_IMATRIX(keyXtraSd); confset.keyXtraSd=keyXtraSd; 
  DATA_IVECTOR(logNMeanAssumption); confset.logNMeanAssumption=logNMeanAssumption;
  DATA_INTEGER(initState); confset.initState = initState;

  DATA_INTEGER(reportingLevel);

  paraSet<Type> paraset;
  PARAMETER_VECTOR(logFpar); paraset.logFpar=logFpar;  
  PARAMETER_VECTOR(logQpow); paraset.logQpow=logQpow;  
  PARAMETER_VECTOR(logSdLogFsta); paraset.logSdLogFsta=logSdLogFsta;  
  PARAMETER_VECTOR(logSdLogN); paraset.logSdLogN=logSdLogN;  
  PARAMETER_VECTOR(logSdLogP); paraset.logSdLogP=logSdLogP;       //Beta random walk components var
  PARAMETER_VECTOR(logSdLogObs); paraset.logSdLogObs=logSdLogObs; 
  PARAMETER_VECTOR(logSdLogTotalObs); paraset.logSdLogTotalObs=logSdLogTotalObs; 
  PARAMETER_VECTOR(transfIRARdist); paraset.transfIRARdist=transfIRARdist; //transformed distances for IRAR cor obs structure
  PARAMETER_VECTOR(sigmaObsParUS); paraset.sigmaObsParUS=sigmaObsParUS; //choleski elements for unstructured cor obs structure
  PARAMETER_VECTOR(rec_pars); paraset.rec_pars=rec_pars;  
  PARAMETER_VECTOR(itrans_rho); paraset.itrans_rho=itrans_rho;  
  PARAMETER_VECTOR(rhop); paraset.rhop=rhop;                       //Correlation of beta RW components
  PARAMETER_VECTOR(logScale); paraset.logScale=logScale; 
  PARAMETER_VECTOR(logitReleaseSurvival); paraset.logitReleaseSurvival=logitReleaseSurvival;    
  PARAMETER_VECTOR(logitRecapturePhi); paraset.logitRecapturePhi=logitRecapturePhi;    
  PARAMETER_VECTOR(logAlphaSCB); paraset.logAlphaSCB=logAlphaSCB;

  PARAMETER_VECTOR(sepFalpha); paraset.sepFalpha=sepFalpha;    
  PARAMETER_VECTOR(sepFlogitRho); paraset.sepFlogitRho=sepFlogitRho;    
  PARAMETER_VECTOR(sepFlogSd); paraset.sepFlogSd=sepFlogSd;    
  PARAMETER_VECTOR(predVarObs); paraset.predVarObs=predVarObs;
  PARAMETER_VECTOR(logPhiSW); paraset.logPhiSW=logPhiSW;
  PARAMETER_VECTOR(logSdProcLogSW); paraset.logSdProcLogSW=logSdProcLogSW;
  PARAMETER_VECTOR(meanLogSW); paraset.meanLogSW=meanLogSW;
  PARAMETER_VECTOR(logSdLogSW); paraset.logSdLogSW=logSdLogSW;
  PARAMETER_VECTOR(logPhiCW); paraset.logPhiCW=logPhiCW;
  PARAMETER_VECTOR(logSdProcLogCW); paraset.logSdProcLogCW=logSdProcLogCW;
  PARAMETER_VECTOR(meanLogCW); paraset.meanLogCW=meanLogCW;
  PARAMETER_VECTOR(logSdLogCW); paraset.logSdLogCW=logSdLogCW;
  PARAMETER_VECTOR(logPhiMO); paraset.logPhiMO=logPhiMO;
  PARAMETER_VECTOR(logSdProcLogitMO); paraset.logSdProcLogitMO=logSdProcLogitMO;
  PARAMETER_VECTOR(meanLogitMO); paraset.meanLogitMO=meanLogitMO;
  PARAMETER_VECTOR(logSdMO); paraset.logSdMO=logSdMO;
  PARAMETER_VECTOR(logPhiNM); paraset.logPhiNM=logPhiNM;
  PARAMETER_VECTOR(logSdProcLogNM); paraset.logSdProcLogNM=logSdProcLogNM;
  PARAMETER_VECTOR(meanLogNM); paraset.meanLogNM=meanLogNM;
  PARAMETER_VECTOR(logSdLogNM); paraset.logSdLogNM=logSdLogNM;
  PARAMETER_VECTOR(logXtraSd); paraset.logXtraSd=logXtraSd;
  PARAMETER_VECTOR(initF); paraset.initF = initF;
  PARAMETER_VECTOR(initN); paraset.initN = initN;
  
  
  // Forecast FMSY
  PARAMETER(logFScaleMSY); paraset.logFScaleMSY = logFScaleMSY;
  PARAMETER(implicitFunctionDelta); paraset.implicitFunctionDelta = implicitFunctionDelta;

  // YPR reference points
  // PARAMETER(logScaleFmsy); paraset.logScaleFmsy = logScaleFmsy;
  // PARAMETER(logScaleFmypyl); paraset.logScaleFmypyl = logScaleFmypyl;
  // PARAMETER(logScaleFmdy); paraset.logScaleFmdy = logScaleFmdy;
  // PARAMETER(logScaleFmax); paraset.logScaleFmax = logScaleFmax;
  // PARAMETER_VECTOR(logScaleFxdYPR); paraset.logScaleFxdYPR = logScaleFxdYPR;
  // PARAMETER_VECTOR(logScaleFxB0); paraset.logScaleFxB0 = logScaleFxB0;
  // PARAMETER(logScaleFcrash); paraset.logScaleFcrash = logScaleFcrash;
  // PARAMETER(logScaleFext); paraset.logScaleFext = logScaleFext;
  // PARAMETER_VECTOR(logScaleFxPercent); paraset.logScaleFxPercent = logScaleFxPercent;
  // PARAMETER(logScaleFlim); paraset.logScaleFlim = logScaleFlim;
  // PARAMETER_MATRIX(logScaleFmsyRange); paraset.logScaleFmsyRange = logScaleFmsyRange;

  PARAMETER(splinePenalty); paraset.splinePenalty = splinePenalty;
  
  PARAMETER_ARRAY(logF); 
  PARAMETER_ARRAY(logN);

  PARAMETER_ARRAY(logSW);
  PARAMETER_ARRAY(logCW);  
  PARAMETER_ARRAY(logitMO);
  PARAMETER_ARRAY(logNM);    

  PARAMETER_ARRAY(logP);
  
  
  PARAMETER_VECTOR(missing);

    // patch missing 
  int idxmis=0; 
  for(int i=0;i<nobs;i++){
    if(isNA(dataset.logobs(i))){
      dataset.logobs(i)=missing(idxmis++);
    }    
  }

  
  


  Recruitment<Type> recruit = makeRecruitmentFunction(confset, paraset);

  prepareForForecast(forecast, dataset, confset, paraset, logF, logN, recruit);

  Type ans=0.0; //negative log-likelihood

  if(CppAD::Variable(keep.sum())){ // add wide prior for first state, but _only_ when computing ooa residuals
    Type huge = 10.0;
    for (int i = 0; i < missing.size(); i++) ans -= dnorm(missing(i), Type(0.0), huge, true);  
  }
  ans += nllSplinePenalty(dataset, confset, paraset, this);

  
  ans += nllSW(logSW, dataset, confset, paraset, this);
  ans += nllCW(logCW, dataset, confset, paraset, this);
  ans += nllMO(logitMO, dataset, confset, paraset, this);
  ans += nllNM(logNM, dataset, confset, paraset, this);
  
  MortalitySet<Type> mort(dataset, confset, paraset, logF);
  forecast.calculateForecast(logF,logN, dataset, confset, paraset, recruit, mort);    

  ans += nllP(confset, paraset, logP, keep, this);
  
  ans += nllF(dataset, confset, paraset, forecast, logF, keep, this);
  ans += nllN(dataset, confset, paraset, forecast, logN, logF, recruit, mort, keep, this);
  forecastSimulation(dataset, confset, paraset, forecast, logN, logF, recruit,mort, this);

  ans += nllObs(dataset, confset, paraset, forecast, logN, logF, logP, recruit, mort, keep,reportingLevel, this);

  reportDeterministicReferencePoints(dataset, confset, paraset, logN, logF, recruit, referencepoints, this);
  

  return ans;
}
