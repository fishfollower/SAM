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

#define TMB_LIB_INIT R_init_stockassessment
#include <TMB.hpp>
#include "../inst/include/macros.h"
#include "../inst/include/f.h"
#include "../inst/include/n.h"
#include "../inst/include/predobs.h"
#include "../inst/include/obs.h"
#include "../inst/include/derived.h"

template<class Type>
Type objective_function<Type>::operator() ()
{
  using CppAD::abs;
  DATA_INTEGER(noFleets);
  DATA_IVECTOR(fleetTypes); 
  DATA_VECTOR(sampleTimes);
  DATA_INTEGER(noYears);
  DATA_VECTOR(years);
  DATA_IVECTOR(minAgePerFleet);
  DATA_IVECTOR(maxAgePerFleet);
  DATA_INTEGER(nobs);
  DATA_IARRAY(idx1);    // minimum index of obs by fleet x year
  DATA_IARRAY(idx2);    // maximum index of obs by fleet x year
  DATA_IARRAY(aux);
  DATA_VECTOR(logobs);
  DATA_VECTOR(weight);
  DATA_VECTOR_INDICATOR(keep, logobs);
  DATA_ARRAY(propMat);
  DATA_ARRAY(stockMeanWeight); 
  DATA_ARRAY(catchMeanWeight);
  DATA_ARRAY(natMor);
  DATA_ARRAY(landFrac);
  DATA_ARRAY(disMeanWeight);
  DATA_ARRAY(landMeanWeight);
  DATA_ARRAY(propF);
  DATA_ARRAY(propM);
  DATA_INTEGER(minAge);
  DATA_INTEGER(maxAge);
  DATA_INTEGER(maxAgePlusGroup);
  DATA_IARRAY(keyLogFsta);
  DATA_INTEGER(corFlag);
  DATA_IARRAY(keyLogFpar);
  DATA_IARRAY(keyQpow);
  DATA_IARRAY(keyVarF);
  DATA_IVECTOR(keyVarLogN); 
  DATA_IARRAY(keyVarObs);
  DATA_FACTOR(obsCorStruct); 
  DATA_IARRAY(keyCorObs);
  DATA_INTEGER(stockRecruitmentModelCode);
  DATA_INTEGER(noScaledYears);
  DATA_IVECTOR(keyScaledYears);
  DATA_IARRAY(keyParScaledYA);
  DATA_IVECTOR(fbarRange);
  DATA_IVECTOR(keyBiomassTreat)
  DATA_INTEGER(simFlag); //1 means simulations should not redo F and N
  DATA_INTEGER(resFlag); 
  DATA_FACTOR(obsLikelihoodFlag);
  DATA_INTEGER(fixVarToWeight);

  PARAMETER_VECTOR(logFpar); 
  PARAMETER_VECTOR(logQpow); 
  PARAMETER_VECTOR(logSdLogFsta); 
  PARAMETER_VECTOR(logSdLogN); 
  PARAMETER_VECTOR(logSdLogObs);
  PARAMETER_VECTOR(logSdLogTotalObs);
  PARAMETER_VECTOR(transfIRARdist);//transformed distances for IRAR cor obs structure
  PARAMETER_VECTOR(sigmaObsParUS);//choleski elements for unstructured cor obs structure
  PARAMETER_VECTOR(rec_loga); 
  PARAMETER_VECTOR(rec_logb); 
  PARAMETER_VECTOR(itrans_rho); 
  PARAMETER_VECTOR(logScale);
  PARAMETER_VECTOR(logitReleaseSurvival);   
  PARAMETER_VECTOR(logitRecapturePhi);   
  PARAMETER_ARRAY(logF); 
  PARAMETER_ARRAY(logN);
  PARAMETER_VECTOR(missing);
  int timeSteps=logF.dim[1];
  int stateDimN=logN.dim[0];

  // patch missing 
  int idxmis=0; 
  for(int i=0;i<nobs;i++){
    if(isNA(logobs(i))){
      logobs(i)=missing(idxmis++);
    }    
  }

  vector<Type> R = rFun(logN, timeSteps);
  vector<Type> logR = log(R);  

  testdat<Type> dat;
  dat.logF=logF;
  dat.logN=logN;

  //vector<Type> ssb = ssbFun(logF, logN, timeSteps, stateDimN, keyLogFsta, natMor, propM, propF, propMat, stockMeanWeight);
  vector<Type> ssb = ssbFun(dat, timeSteps, stateDimN, keyLogFsta, natMor, propM, propF, propMat, stockMeanWeight);
  vector<Type> logssb = log(ssb);

  vector<Type> fbar = fbarFun(logF, minAge, timeSteps, fbarRange, keyLogFsta);
  vector<Type> logfbar = log(fbar);

  vector<Type> cat = catchFun(logF, logN, minAge, maxAge, keyLogFsta, catchMeanWeight, natMor);
  vector<Type> logCatch = log(cat);

  vector<Type> fsb = fsbFun(logF, logN, minAge, maxAge, keyLogFsta, catchMeanWeight, natMor);
  vector<Type> logfsb = log(fsb);

  vector<Type> tsb = tsbFun(logN, minAge, maxAge, timeSteps, stockMeanWeight);
  vector<Type> logtsb = log(tsb);

  vector<Type> predObs=predObsFun(logF, logN, logFpar, logScale,
                                  logQpow, nobs, minAge, maxAge,
                                  noScaledYears, fleetTypes,
                                  keyScaledYears, keyQpow,
                                  keyBiomassTreat, aux, keyLogFsta,
                                  keyLogFpar, keyParScaledYA, natMor,
                                  sampleTimes, logssb, logfsb,
                                  logCatch, logitReleaseSurvival);

  Type ans=0; //negative log-likelihood

  ans+=nllF(logF, timeSteps, corFlag, simFlag, resFlag, stateDimN,
            keyLogFsta, keyVarF, itrans_rho, logSdLogFsta, keep, this );

  ans += nllN(logN, logF, timeSteps, stateDimN, minAge,
              maxAgePlusGroup, simFlag, resFlag, keyVarLogN,
              keyLogFsta, stockRecruitmentModelCode, ssb, natMor,
              logSdLogN, rec_loga, rec_logb, keep, this);

  ans += nllObs(nobs, noFleets, noYears, fleetTypes, minAgePerFleet,
                maxAgePerFleet, minAge, maxAge, obsCorStruct,
                logSdLogObs, logSdLogTotalObs, transfIRARdist,
                sigmaObsParUS, logitRecapturePhi, keyVarObs, keyCorObs, 
                aux, obsLikelihoodFlag, idx1, idx2, weight, fixVarToWeight,
                logobs, predObs, keep, this);
  

  if(CppAD::Variable(keep.sum())){ // add wide prior for first state, but _only_ when computing ooa residuals
    Type huge = 10;
    for (int i = 0; i < missing.size(); i++) ans -= dnorm(missing(i), Type(0), huge, true);  
  } 

  SIMULATE {
    REPORT(logF);
    REPORT(logN);
    REPORT(logobs);
  }
  REPORT(predObs);
  ADREPORT(logssb);
  ADREPORT(logfbar);
  ADREPORT(logCatch);
  ADREPORT(logtsb);
  ADREPORT(logR);

  vector<Type> lastLogN = logN.col(timeSteps-1);
  ADREPORT(lastLogN);
  vector<Type> lastLogF = logF.col(timeSteps-1);
  ADREPORT(lastLogF);  

  vector<Type> beforeLastLogN = logN.col(timeSteps-2);
  ADREPORT(beforeLastLogN);
  vector<Type> beforeLastLogF = logF.col(timeSteps-2);
  ADREPORT(beforeLastLogF);  
  return ans;
}
