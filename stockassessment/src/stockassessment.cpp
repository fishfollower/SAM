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
#include <TMB.hpp>
#include "../inst/include/SAM.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  using CppAD::abs;
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
  DATA_STRUCT(forecast, forecastSet); dataset.forecast = forecast;
  DATA_STRUCT(referencepoint, referencepointSet); dataset.referencepoint = referencepoint;

  prepareForForecast(dataset);
    
  confSet confset;
  DATA_INTEGER(minAge); confset.minAge=minAge; 
  DATA_INTEGER(maxAge); confset.maxAge=maxAge; 
  DATA_IVECTOR(maxAgePlusGroup); confset.maxAgePlusGroup=maxAgePlusGroup; 
  DATA_IARRAY(keyLogFsta); confset.keyLogFsta=keyLogFsta; 
  DATA_INTEGER(corFlag); confset.corFlag=corFlag; 
  DATA_IARRAY(keyLogFpar); confset.keyLogFpar=keyLogFpar; 
  DATA_IARRAY(keyQpow); confset.keyQpow=keyQpow; 
  DATA_IARRAY(keyVarF); confset.keyVarF=keyVarF; 
  DATA_IVECTOR(keyVarLogN); confset.keyVarLogN=keyVarLogN;  
  DATA_IARRAY(keyVarObs); confset.keyVarObs=keyVarObs; 
  DATA_FACTOR(obsCorStruct); confset.obsCorStruct=obsCorStruct;  
  DATA_IARRAY(keyCorObs); confset.keyCorObs=keyCorObs; 
  DATA_INTEGER(stockRecruitmentModelCode); confset.stockRecruitmentModelCode=stockRecruitmentModelCode; 
  DATA_VECTOR(constRecBreaks); vector<double> constRecBreaksDouble(constRecBreaks.size()); for(int i=0; i<constRecBreaks.size(); ++i){constRecBreaksDouble(i)=asDouble(constRecBreaks(i));} confset.constRecBreaks=constRecBreaksDouble; 
  DATA_INTEGER(noScaledYears); confset.noScaledYears=noScaledYears; 
  DATA_IVECTOR(keyScaledYears); confset.keyScaledYears=keyScaledYears; 
  DATA_IMATRIX(keyParScaledYA); 
  confset.keyParScaledYA=keyParScaledYA; 
  DATA_IVECTOR(fbarRange); confset.fbarRange=fbarRange; 
  DATA_IVECTOR(keyBiomassTreat); confset.keyBiomassTreat=keyBiomassTreat;   
  DATA_IVECTOR(simFlag); confset.simFlag=simFlag;  //1 means simulations should not redo F and N
  DATA_INTEGER(resFlag); confset.resFlag=resFlag;  
  DATA_FACTOR(obsLikelihoodFlag); confset.obsLikelihoodFlag=obsLikelihoodFlag; 
  DATA_INTEGER(fixVarToWeight); confset.fixVarToWeight=fixVarToWeight; 
  DATA_SCALAR(fracMixF); confset.fracMixF=asDouble(fracMixF); 
  DATA_SCALAR(fracMixN); confset.fracMixN=asDouble(fracMixN); 
  DATA_VECTOR(fracMixObs); vector<double> fracMixObsDouble(fracMixObs.size()); for(int i=0; i<fracMixObs.size(); ++i){fracMixObsDouble(i)=asDouble(fracMixObs(i));} confset.fracMixObs=fracMixObsDouble; 
  DATA_IARRAY(predVarObsLink); confset.predVarObsLink=predVarObsLink;

  paraSet<Type> paraset;
  PARAMETER_VECTOR(logFpar); paraset.logFpar=logFpar;  
  PARAMETER_VECTOR(logQpow); paraset.logQpow=logQpow;  
  PARAMETER_VECTOR(logSdLogFsta); paraset.logSdLogFsta=logSdLogFsta;  
  PARAMETER_VECTOR(logSdLogN); paraset.logSdLogN=logSdLogN;  
  PARAMETER_VECTOR(logSdLogObs); paraset.logSdLogObs=logSdLogObs; 
  PARAMETER_VECTOR(logSdLogTotalObs); paraset.logSdLogTotalObs=logSdLogTotalObs; 
  PARAMETER_VECTOR(transfIRARdist); paraset.transfIRARdist=transfIRARdist; //transformed distances for IRAR cor obs structure
  PARAMETER_VECTOR(sigmaObsParUS); paraset.sigmaObsParUS=sigmaObsParUS; //choleski elements for unstructured cor obs structure
  PARAMETER_VECTOR(rec_pars); paraset.rec_pars=rec_pars;  
  PARAMETER_VECTOR(itrans_rho); paraset.itrans_rho=itrans_rho;  
  PARAMETER_VECTOR(logScale); paraset.logScale=logScale; 
  PARAMETER_VECTOR(logitReleaseSurvival); paraset.logitReleaseSurvival=logitReleaseSurvival;    
  PARAMETER_VECTOR(logitRecapturePhi); paraset.logitRecapturePhi=logitRecapturePhi;

  PARAMETER_VECTOR(sepFalpha); paraset.sepFalpha=sepFalpha;    
  PARAMETER_VECTOR(sepFlogitRho); paraset.sepFlogitRho=sepFlogitRho;    
  PARAMETER_VECTOR(sepFlogSd); paraset.sepFlogSd=sepFlogSd;    
  PARAMETER_VECTOR(predVarObs); paraset.predVarObs=predVarObs;

  // Forecast FMSY
  PARAMETER(logFScaleMSY); paraset.logFScaleMSY = logFScaleMSY;
  PARAMETER(implicitFunctionDelta); paraset.implicitFunctionDelta = implicitFunctionDelta;

  // YPR reference points
  PARAMETER(logScaleFmsy); paraset.logScaleFmsy = logScaleFmsy;
  PARAMETER(logScaleFmax); paraset.logScaleFmax = logScaleFmax;
  PARAMETER(logScaleF01); paraset.logScaleF01 = logScaleF01;
  PARAMETER(logScaleFcrash); paraset.logScaleFcrash = logScaleFcrash;
  PARAMETER_VECTOR(logScaleFxPercent); paraset.logScaleFxPercent = logScaleFxPercent;
  PARAMETER(logScaleFlim); paraset.logScaleFlim = logScaleFlim;
  
  PARAMETER_ARRAY(logF); 
  PARAMETER_ARRAY(logN);
  PARAMETER_VECTOR(missing);

  // patch missing 
  int idxmis=0; 
  for(int i=0;i<nobs;i++){
    if(isNA(dataset.logobs(i))){
      dataset.logobs(i)=missing(idxmis++);
    }    
  }
  
  Type ans=0; //negative log-likelihood

  if(CppAD::Variable(keep.sum())){ // add wide prior for first state, but _only_ when computing ooa residuals
    Type huge = 10;
    for (int i = 0; i < missing.size(); i++) ans -= dnorm(missing(i), Type(0), huge, true);  
  } 

  dataset.forecast.calculateForecast(logF,logN, dataset, confset, paraset);

  ans += nllF(dataset, confset, paraset, logF, keep, this);
  ans += nllN(dataset, confset, paraset, logN, logF, keep, this);
  ans += nllObs(dataset, confset, paraset, logN, logF, keep,  this);

  ans += nllReferencepoints(dataset, confset, paraset, logN, logF, this);
  
  return ans;
}
