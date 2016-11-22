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

/* Parameter transform */
template <class Type>
Type trans(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

bool isNAINT(int x){
  return NA_INTEGER==x;
}

template <class Type>
matrix<Type> setupVarCovMatrix(int minAge, int maxAge, int minAgeFleet, int maxAgeFleet, vector<int> rhoMap, vector<Type> rhoVec, vector<int> sdMap, vector<Type> sdVec){

  int dim = maxAgeFleet-minAgeFleet+1;
  int offset = minAgeFleet-minAge;
  matrix<Type> ret(dim,dim);
  ret.setZero();

  Type rho0 = Type(0.5);
  vector<Type> xvec(dim);
  xvec(0)=Type(0);
  int maxrm=-1;
  if(rhoVec.size()>0){
    for(int i=1; i<xvec.size(); i++) { 
      if(rhoMap(i-1+offset)>=0)
	xvec(i) = xvec(i-1)+rhoVec(rhoMap(i-1+offset)); 
      if(rhoMap(i-1)>maxrm) maxrm=rhoMap(i-1);
    } 
  }
   
  for(int i=0; i<dim; i++)
    for(int j=0; j<dim; j++){
      if(i!=j && maxrm>=0){	
	Type dist = abs(xvec(i)-xvec(j));
     	ret(i,j)=pow( rho0,dist)*sdVec( sdMap(i+offset) )*sdVec( sdMap(j+offset));
      } else if(i==j) ret(i,j) = sdVec( sdMap(i+offset) )*sdVec( sdMap(j+offset));
    }
  
  return ret;
}


template <class Type>
vector<Type> addLogratio(vector<Type> logx){
  int n = logx.size();
  vector<Type> res(n-1);
  for(int i = 0; i < res.size(); ++i)
    res(i) = logx(i) - logx(n-1);
  return res;//log(x.head(x.size()-1)/x.tail(1));
};

template<class Type>
vector<Type> multLogratio(vector<Type> logx){
  vector<Type> res(logx.size()-1);
  for(int i = 0; i < res.size(); ++i)
    res(i) = logx(i)-log(Type(1.0)-exp(logExpSum(logx.head(i+1))));
  return res;
}


template <class Type>
Type log2expsum(vector<Type> x){
  return exp(x).sum();
}

template<class Type>
Type logExpSum(vector<Type> x){
  Type m = max(x);
  return m + log(exp(x-m).sum());
}

template <class Type>
vector<Type> log2proportion(vector<Type> x){
  return exp(x) / log2expsum(x);
}


template<class Type>
matrix<Type> buildJac(vector<Type> x, vector<Type> w){
  matrix<Type> res(x.size(),x.size()); 
  Type xs = x.sum();
  Type xsp = pow(xs,Type(2));
  for(int i = 0; i < res.rows(); ++i){
    for(int j = 0; j < res.cols(); ++j){
      if(i == j){
	res(i,j) = Type(1.0)/xs-x(i)/xsp;
      }else{
	res(i,j) = -x(i)/xsp;
      }
    }
  }
  for(int j = 0; j < res.cols(); ++j){
    res(res.rows()-1,j) = w(j);
  }
  return res;
}


template <class Type>
Type jacobianDet(vector<Type> x){
  vector<Type> w(x.size());
  w.fill(Type(1.0));
  return buildJac(x,w).determinant();
}
template <class Type>
Type jacobianDet(vector<Type> x,vector<Type> w){
  return buildJac(x,w).determinant();
}




template <class Type> 
density::UNSTRUCTURED_CORR_t<Type> getCorrObj(vector<Type> params){
  density::UNSTRUCTURED_CORR_t<Type> ret(params);
  return ret;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
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
  DATA_INTEGER(simFlag); //1 means simulations should not redo F and N
  DATA_FACTOR(obsLikelihoodFlag);

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
  PARAMETER_VECTOR(logScaleSSB); 
  PARAMETER_VECTOR(logPowSSB); 
  PARAMETER_VECTOR(logSdSSB);
  PARAMETER_ARRAY(logF); 
  PARAMETER_ARRAY(logN);
  PARAMETER_VECTOR(missing);
  int timeSteps=logF.dim[1];
  int stateDimF=logF.dim[0];
  int stateDimN=logN.dim[0];
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
  
  vector<Type> IRARdist(transfIRARdist.size()); //[ d_1, d_2, ...,d_N-1 ]
  if(transfIRARdist.size()>0) IRARdist=exp(transfIRARdist);

  vector< vector<Type> > sigmaObsParVec(noFleets);
  int nfleet = maxAgePerFleet(0)-minAgePerFleet(0)+1;
  int dn=nfleet*(nfleet-1)/2;
  int from=-dn, to=-1; 
  for(int f=0; f<noFleets; f++){
    if(obsCorStruct(f)!=2) continue;
    nfleet = maxAgePerFleet(f)-minAgePerFleet(f)+1;
    dn = nfleet*(nfleet-1)/2;
    from=to+1;
    to=to+dn;
    vector<Type> tmp(dn);
    for(int i=from; i<=to; i++) tmp(i-from) = sigmaObsParUS(i);
    sigmaObsParVec(f) = tmp; 
  }

  Type ans=0; //negative log-likelihood


  // patch missing 
  int idxmis=0; 
  for(int i=0;i<nobs;i++){
    if(isNA(logobs(i))){
      logobs(i)=missing(idxmis++);
    }    
  }
  
  //First take care of F
  matrix<Type> fvar(stateDimF,stateDimF);
  matrix<Type> fcor(stateDimF,stateDimF);
  vector<Type> fsd(stateDimF);  

  if(corFlag==0){
    fcor.setZero();
  }

  for(int i=0; i<stateDimF; ++i){
    fcor(i,i)=1.0;
  }

  if(corFlag==1){
    for(int i=0; i<stateDimF; ++i){
      for(int j=0; j<i; ++j){
        fcor(i,j)=trans(itrans_rho(0));
        fcor(j,i)=fcor(i,j);
      }
    } 
  }

  if(corFlag==2){
    for(int i=0; i<stateDimF; ++i){
      for(int j=0; j<i; ++j){
        fcor(i,j)=pow(trans(itrans_rho(0)),abs(Type(i-j)));
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
    SIMULATE{
      if(simFlag==0){
        logF.col(i)=logF.col(i-1)+neg_log_densityF.simulate();
      }
    }
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
      if(i!=j){nvar(i,j)=0.0;}else{nvar(i,j)=varLogN(keyVarLogN(i));}
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
    SIMULATE{
      if(simFlag==0){
        logN.col(i) = predN + neg_log_densityN.simulate();
      }
    }
  }
  

  // Calculate predicted observations
  int f, ft, a, y,yy, scaleIdx;  // a is no longer just ages, but an attribute (e.g. age or length) 
  int minYear=aux(0,0);
  Type zz;
  vector<Type> predObs(nobs);
  vector<Type> predSd(nobs);
  for(int i=0;i<nobs;i++){
    y=aux(i,0)-minYear;
    f=aux(i,1);
    ft=fleetTypes(f-1);
    a=aux(i,2)-minAge;
    zz=exp(logF(keyLogFsta(0,a),y))+natMor(y,a);
    
    switch(ft){
      case 0:
        predObs(i)=logN(a,y)-log(zz)+log(1-exp(-zz));
        if(keyLogFsta(f-1,a)>(-1)){
          predObs(i)+=logF(keyLogFsta(0,a),y);
        }
        scaleIdx=-1;
        yy=aux(i,0);
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
  	error("Unknown fleet code");
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
  	error("Unknown fleet code");
        return 0;
      break;
  
      case 4:
  	error("Unknown fleet code");
        return 0;
      break;
  
      case 5:
  	error("Unknown fleet code");
        return 0;
      break;
  
      case 6:
  	error("Unknown fleet code");
        return 0;
      break;
  
      case 7:
  	error("Unknown fleet code");
        return 0;
      break;
  
      default:
  	error("Unknown fleet code");
        return 0 ;
      break;
    }    
    //predSd(i)=sqrt(varLogObs(keyVarObs(f-1,a)));
    //ans+=-keep(i)*dnorm(logobs(i),predObs(i),predSd(i),true);
    //SIMULATE{
    //  logobs(i)=rnorm(predObs(i), predSd(i));
    //}
  }

  // setup obs likelihoods
  vector< density::MVNORM_t<Type> >  nllVec(noFleets);
  vector< density::UNSTRUCTURED_CORR_t<Type> > neg_log_densityObsUnstruc(noFleets);
  vector< vector<Type> > obsCovScaleVec(noFleets);
  for(int f=0; f<noFleets; ++f){
    int thisdim=maxAgePerFleet(f)-minAgePerFleet(f)+1;
    matrix<Type> cov(thisdim,thisdim);
    cov.setZero();
    
    if(obsCorStruct(f)==0){//ID (independent)  
      for(int i=0; i<thisdim; ++i){
	cov(i,i)=varLogObs(keyVarObs(f,i+minAgePerFleet(f)-minAge));
      }
    } else if(obsCorStruct(f)==1){//(AR) irregular lattice AR
      cov = setupVarCovMatrix(minAge, maxAge, minAgePerFleet(f), maxAgePerFleet(f), keyCorObs.transpose().col(f), IRARdist, keyVarObs.transpose().col(f) , exp(logSdLogObs) );
    } else if(obsCorStruct(f)==2){//(US) unstructured
      neg_log_densityObsUnstruc(f) = getCorrObj(sigmaObsParVec(f));  
      matrix<Type> tmp = neg_log_densityObsUnstruc(f).cov();
      
      tmp.setZero();
      int offset = minAgePerFleet(f)-minAge;
      obsCovScaleVec(f).resize(tmp.rows());
      for(int i=0; i<tmp.rows(); i++) {
	tmp(i,i) = sqrt(varLogObs(keyVarObs(f,i+offset)));
	obsCovScaleVec(f)(i) = tmp(i,i);
      }
      cov  = tmp*matrix<Type>(neg_log_densityObsUnstruc(f).cov()*tmp);

    } else { error("Unkown obsCorStruct code"); }
    if(obsLikelihoodFlag(f) == 1){ // Additive logistic normal needs smaller covariance matrix
      nllVec(f).setSigma(cov.block(0,0,thisdim-1,thisdim-1));
    }else{
      nllVec(f).setSigma(cov);
    }
  }
  
  //eval likelihood 
  for(int y=0;y<noYears;y++){
    int totalParKey = 0;
    for(int f=0;f<noFleets;f++){
      if(!isNAINT(idx1(f,y))){
        int idxfrom=idx1(f,y);
        int idxlength=idx2(f,y)-idx1(f,y)+1;
	switch(obsLikelihoodFlag(f)){
	case 0: // (LN) log-Normal distribution
	  ans += nllVec(f)(logobs.segment(idxfrom,idxlength)-predObs.segment(idxfrom,idxlength));
	  // += ((vector<Type>)logobs.segment(idxfrom,idxlength)).sum();
	  SIMULATE{
	    logobs.segment(idxfrom,idxlength) = predObs.segment(idxfrom,idxlength) + nllVec(f).simulate();
	  }
	  break;
	case 1: // (ALN) Additive logistic-normal proportions + log-normal total numbers
	  ans +=  nllVec(f)(addLogratio((vector<Type>)logobs.segment(idxfrom,idxlength))-addLogratio((vector<Type>)predObs.segment(idxfrom,idxlength)));
	   ans += log(log2proportion((vector<Type>)logobs.segment(idxfrom,idxlength))).sum();
	   ans -= dnorm(log(log2expsum((vector<Type>)logobs.segment(idxfrom,idxlength))),
	   	       log(log2expsum((vector<Type>)predObs.segment(idxfrom,idxlength))),
	   	       exp(logSdLogTotalObs(totalParKey++)),true);
	   ans += log(log2expsum((vector<Type>)logobs.segment(idxfrom,idxlength)));
	   ans -= log(abs(jacobianDet((vector<Type>)logobs.segment(idxfrom,idxlength).exp())));
	  SIMULATE{
	    vector<Type> logProb(idxlength);
	    logProb.setZero();
	    logProb.segment(0,idxlength-1) = addLogratio(((vector<Type>)predObs.segment(idxfrom,idxlength))) + nllVec(f).simulate();
	    Type logDenom = logExpSum(logProb);
	    logProb -= logDenom;
	    Type logTotal = rnorm(log(log2expsum((vector<Type>)predObs.segment(idxfrom,idxlength))),
				  exp(logSdLogTotalObs(totalParKey++)));
	    logobs.segment(idxfrom,idxlength) = logProb + logTotal; 
	  }
	  break;
	default:
	  error("Unknown obsLikelihoodFlag");
	}
      }  
    }  
  }
  
  // derived quantities for ADreport
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

  if(CppAD::Variable(keep.sum())){ // add wide prior for first state, but _only_ when computing ooa residuals
    Type huge = 10;
    for (int i = 0; i < stateDimN; i++) ans -= dnorm(logN(i, 0), Type(0), huge, true);  
    for (int i = 0; i < stateDimF; i++) ans -= dnorm(logF(i, 0), Type(0), huge, true);  
  } 
  
  SIMULATE {
    REPORT(logF);
    REPORT(logN);
    REPORT(logobs);
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
