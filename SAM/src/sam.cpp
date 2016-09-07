//  --------------------------------------------------------------------------
// Copyright (c) 2014, Anders Nielsen <an@aqua.dtu.dk>,    
// Casper Berg <cbe@aqua.dtu.dk>, and Kasper Kristensen <kkr@aqua.dtu.dk>.
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
 
#include <TMB.hpp>
#include <iostream>

/* Parameter transform */
template <class Type>
Type f(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

template <class Type> 
Type square(Type x){return x*x;}

 template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(noFleets);
  DATA_VECTOR(fleetTypes); 
  DATA_VECTOR(sampleTimes);
  DATA_INTEGER(noYears);
  DATA_VECTOR(years);
  DATA_INTEGER(nobs);
  DATA_VECTOR(idx1);
  DATA_VECTOR(idx2);
  DATA_ARRAY(obs);
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
  DATA_IARRAY(keyVarLogN); 
  DATA_IARRAY(keyVarObs); 
  DATA_INTEGER(stockRecruitmentModelCode);
  DATA_INTEGER(noScaledYears);
  DATA_IVECTOR(keyScaledYears);
  DATA_IMATRIX(keyParScaledYA);
  DATA_IVECTOR(fbarRange);

  PARAMETER_VECTOR(logFpar); 
  PARAMETER_VECTOR(logQpow); 
  PARAMETER_VECTOR(logSdLogFsta); 
  PARAMETER_VECTOR(logSdLogN); 
  PARAMETER_VECTOR(logSdLogObs); 
  PARAMETER_VECTOR(rec_loga); 
  PARAMETER_VECTOR(rec_logb); 
  PARAMETER(logit_rho); 
  PARAMETER_VECTOR(logScale); 
  PARAMETER_VECTOR(logScaleSSB); 
  PARAMETER_VECTOR(logPowSSB); 
  PARAMETER_VECTOR(logSdSSB); 
  PARAMETER_ARRAY(logF); 
  PARAMETER_ARRAY(logN);
  int timeSteps=logF.dim[1];
  int stateDimF=logF.dim[0];
  int stateDimN=logN.dim[0];
  Type rho=f(logit_rho);
  vector<Type> sdLogFsta=exp(logSdLogFsta);
  vector<Type> varLogN=exp(logSdLogN*Type(2.0));
  vector<Type> varLogObs=exp(logSdLogObs*Type(2.0));
  vector<Type> ssb(timeSteps);
  vector<Type> logssb(timeSteps);
  vector<Type> fbar(timeSteps);
  vector<Type> logfbar(timeSteps);
  vector<Type> cat(catchMeanWeight.dim(0));
  vector<Type> logCatch(catchMeanWeight.dim(0));
  vector<Type> tsb(timeSteps);
  vector<Type> logtsb(timeSteps);
  vector<Type> logR(timeSteps);
  vector<Type> R(timeSteps);

  Type ans=0; //negative log-likelihood
  
  //First take care of F
  matrix<Type> fvar(stateDimF,stateDimF);
  matrix<Type> fcor(stateDimF,stateDimF);
  vector<Type> fsd(stateDimF);
  
  for(int i=0; i<stateDimF; ++i){
    fcor(i,i)=1.0;
  }

  if(corFlag==1){
    for(int i=0; i<stateDimF; ++i){
      for(int j=0; j<i; ++j){
        fcor(i,j)=rho;
        fcor(j,i)=fcor(i,j);
      }
    } 
  }

  if(corFlag==2){
    for(int i=0; i<stateDimF; ++i){
      for(int j=0; j<i; ++j){
        fcor(i,j)=pow(rho,abs(Type(i-j)));
        fcor(j,i)=fcor(i,j);
      }
    } 
  }

  int i,j;
  for(i=0; i<stateDimF; ++i){
    for(j=0; j<stateDimN; ++j){
      if(keyLogFsta(0,j)==i)break;
    }
    fsd(i)=sdLogFsta(keyVarF(0,j));
  }
 
  for(i=0; i<stateDimF; ++i){
    for(j=0; j<stateDimF; ++j){
      fvar(i,j)=fsd(i)*fsd(j)*fcor(i,j);
    }
  }

  using namespace density;
  MVNORM_t<Type> neg_log_densityF(fvar);
  for(int i=1;i<timeSteps;i++){    
     ans+=neg_log_densityF(logF.col(i)-logF.col(i-1)); // F-Process likelihood 
  }
  
  for(int i=0;i<timeSteps;i++){ // calc ssb
    ssb(i)=0.0;    
    for(int j=0; j<stateDimN; ++j){
      if(keyLogFsta(0,j)>(-1)){
        ssb(i)+=exp(logN(j,i))*exp(-exp(logF(keyLogFsta(0,j),i))*propF(i,j)-natMor(i,j)*propM(i,j))*propMat(i,j)*stockMeanWeight(i,j);
      }else{
        ssb(i)+=exp(logN(j,i))*exp(-natMor(i,j)*propM(i,j))*propMat(i,j)*stockMeanWeight(i,j);
      }
    }
    logssb(i)=log(ssb(i));
  }
  
  //Now take care of N
  matrix<Type> nvar(stateDimN,stateDimN);
  for(int i=0; i<stateDimN; ++i){
    for(int j=0; j<stateDimN; ++j){
      if(i!=j){nvar(i,j)=0.0;}else{nvar(i,j)=varLogN(keyVarLogN(0,i));}
    }
  }
  MVNORM_t<Type> neg_log_densityN(nvar);
  vector<Type> predN(stateDimN); 
  for(int i=1;i<timeSteps;i++){ 
    if(stockRecruitmentModelCode==0){ // straight RW 
      predN(0)=logN(0,i-1);
    }else{
      if(stockRecruitmentModelCode==1){//ricker
        predN(0)=rec_loga(0)+log(ssb(i-1))-exp(rec_logb(0))*ssb(i-1); 
      }else{
        if(stockRecruitmentModelCode==2){//BH
          predN(0)=rec_loga(0)+log(ssb(i-1))-log(1.0+exp(rec_logb(0))*ssb(i-1)); 
        }else{
          error("SR model code not recognized");
        }
      }
    }
  
    for(int j=1; j<stateDimN; ++j){
      if(keyLogFsta(0,j-1)>(-1)){
        predN(j)=logN(j-1,i-1)-exp(logF(keyLogFsta(0,j-1),i-1))-natMor(i-1,j-1); 
      }else{
        predN(j)=logN(j-1,i-1)-natMor(i-1,j-1); 
      }
    }  
    if(maxAgePlusGroup==1){
      predN(stateDimN-1)=log(exp(logN(stateDimN-2,i-1)-exp(logF(keyLogFsta(0,stateDimN-2),i-1))-natMor(i-1,stateDimN-2))+
                             exp(logN(stateDimN-1,i-1)-exp(logF(keyLogFsta(0,stateDimN-1),i-1))-natMor(i-1,stateDimN-1))); 
    }
    ans+=neg_log_densityN(logN.col(i)-predN); // N-Process likelihood 
  }

  // Now finally match to observations
  int f, ft, a, y,yy, scaleIdx;  // a is no longer just ages, but an attribute (e.g. age or length) 
  int minYear=CppAD::Integer((obs(0,0)));
  Type zz;
  vector<Type> predObs(nobs);
  vector<Type> predSd(nobs);
  for(int i=0;i<nobs;i++){
    y=CppAD::Integer(obs(i,0))-minYear;
    f=CppAD::Integer(obs(i,1));
    ft=CppAD::Integer(fleetTypes(f-1));
    a=CppAD::Integer(obs(i,2))-minAge;
    zz=exp(logF(keyLogFsta(0,a),y))+natMor(y,a);
    
    switch(ft){
      case 0:
        predObs(i)=logN(a,y)-log(zz)+log(1-exp(-zz));
        if(keyLogFsta(f-1,a)>(-1)){
          predObs(i)+=logF(keyLogFsta(0,a),y);
        }
        scaleIdx=-1;
        yy=CppAD::Integer(obs(i,0));
        for(int j=0; j<noScaledYears; ++j){
          if(yy==keyScaledYears(j)){
            scaleIdx=keyParScaledYA(j,a);
            if(scaleIdx>=0){
              predObs(i)-=logScale(scaleIdx);
            }
            break;
          }
        }
      break;
  
      case 1:
  	std::cerr<<"Unknown fleet code: "<<ft<<std::endl;
        return(0);
      break;
  
      case 2:
        predObs(i)=logN(a,y)-zz*sampleTimes(f-1);
        if(keyQpow(f-1,a)>(-1)){
          predObs(i)*=exp(logQpow(keyQpow(f-1,a))); 
        }
        if(keyLogFpar(f-1,a)>(-1)){
          predObs(i)+=logFpar(keyLogFpar(f-1,a));
        }
        
      break;
  
      case 3:
  	std::cerr<<"Unknown fleet code: "<<ft<<std::endl;
        return 0;
      break;
  
      case 4:
  	std::cerr<<"Unknown fleet code: "<<ft<<std::endl;
        return 0;
      break;
  
      case 5:
  	std::cerr<<"Unknown fleet code: "<<ft<<std::endl;
        return 0;
      break;
  
      case 6:
  	std::cerr<<"Unknown fleet code: "<<ft<<std::endl;
        return 0;
      break;
  
      case 7:
  	std::cerr<<"Unknown fleet code: "<<ft<<std::endl;
        return 0;
      break;
  
      default:
  	std::cerr<<"Unknown fleet code: "<<ft<<std::endl;
        return 0 ;
      break;
    }    
    predSd(i)=sqrt(varLogObs(keyVarObs(f-1,a)));
    ans+=-dnorm(log(obs(i,3)),predObs(i),predSd(i),true);
  }

  for(int y=0;y<timeSteps;y++){  
    fbar(y)=Type(0);
    for(int a=fbarRange(0);a<=fbarRange(1);a++){  
      fbar(y)+=exp(logF(keyLogFsta(0,a-minAge),y));
    }
    fbar(y)/=Type(fbarRange(1)-fbarRange(0)+1);
    logfbar(y)=log(fbar(y));
  }

  for(int y=0;y<catchMeanWeight.dim(0);y++){  
    cat(y)=Type(0);
    for(int a=minAge;a<=maxAge;a++){  
      Type z=exp(logF(keyLogFsta(0,a-minAge),y))+natMor(y,a-minAge);
      cat(y)+=exp(logF(keyLogFsta(0,a-minAge),y))/z*exp(logN(a-minAge,y))*(Type(1.0)-exp(-z))*catchMeanWeight(y,a-minAge);
    }
    logCatch(y)=log(cat(y));
  }

  for(int y=0;y<timeSteps;y++){  
    tsb(y)=Type(0);
    for(int a=minAge;a<=maxAge;a++){  
      tsb(y)+=exp(logN(a-minAge,y))*stockMeanWeight(y,a-minAge);
    }
    logtsb(y)=log(tsb(y));
  }

  for(int y=0;y<timeSteps;y++){  
    logR(y)=logN(0,y);
    R(y)=exp(logR(y));
  }
  
  REPORT(predObs);
  REPORT(predSd);
  ADREPORT(ssb);
  ADREPORT(logssb);
  ADREPORT(fbar);
  ADREPORT(logfbar);
  ADREPORT(cat);
  ADREPORT(logCatch);
  ADREPORT(tsb);
  ADREPORT(logtsb);
  ADREPORT(R);
  ADREPORT(logR);

  return ans;
}
