#define REPORT_F(name,F)					\
if(isDouble<Type>::value && F->current_parallel_region<0) {     \
  defineVar(install(#name),                                     \
            asSEXP_protect(name),F->report);                    \
  UNPROTECT(1);                                                 \
}

#define SIMULATE_F(F)						\
if(isDouble<Type>::value && F->do_simulate)

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

bool isNAINT(int x){
  return NA_INTEGER==x;
}

template <class Type>
matrix<Type> setupVarCovMatrix(int minAge, int maxAge, int minAgeFleet, int maxAgeFleet, vector<int> rhoMap, vector<Type> rhoVec, vector<int> sdMap, vector<Type> sdVec){

  using CppAD::abs;
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
density::UNSTRUCTURED_CORR_t<Type> getCorrObj(vector<Type> params){
  density::UNSTRUCTURED_CORR_t<Type> ret(params);
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
Type nllObs(int noFleets, 
            int noYears,
            vector<int> fleetTypes, 
            vector<int> minAgePerFleet, 
            vector<int> maxAgePerFleet, 
            int minAge,
            int maxAge,
            vector<int> obsCorStruct,
            vector<Type> logSdLogObs,
            vector<Type> logSdLogTotalObs,
            vector<Type> transfIRARdist,
            vector<Type> sigmaObsParUS,
            vector<Type> recapturePhiVec,
            array<int> keyVarObs,
            array<int> keyCorObs,
            vector<int> obsLikelihoodFlag,
            array<int> idx1, 
            array<int> idx2,
            vector<Type> weight, 
            int fixVarToWeight,
            vector<Type> &logobs,  
            vector<Type> predObs,
            data_indicator<vector<Type>,Type> keep, 
            objective_function<Type> *of
           ){
  Type nll=0; 
  // setup obs likelihoods
  vector< density::MVNORM_t<Type> >  nllVec(noFleets);
  vector< density::UNSTRUCTURED_CORR_t<Type> > neg_log_densityObsUnstruc(noFleets);
  vector< vector<Type> > obsCovScaleVec(noFleets);
  vector<Type> varLogObs=exp(logSdLogObs*Type(2.0));
  vector<Type> IRARdist(transfIRARdist.size()); //[ d_1, d_2, ...,d_N-1 ]
  if(transfIRARdist.size()>0) IRARdist=exp(transfIRARdist);
  vector< vector<Type> > sigmaObsParVec(noFleets);
  int aidx;
  vector< matrix<Type> > obsCov(noFleets); // for reporting

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

  for(int f=0; f<noFleets; ++f){
    if(fleetTypes(f)!=5){ 
      int thisdim=maxAgePerFleet(f)-minAgePerFleet(f)+1;
      matrix<Type> cov(thisdim,thisdim);
      cov.setZero();
      if(obsCorStruct(f)==0){//ID (independent)  
        for(int i=0; i<thisdim; ++i){
          if(fleetTypes(f)!=3){
            aidx = i+minAgePerFleet(f)-minAge;
          }else{
            aidx = 0;
          }
  	  cov(i,i)=varLogObs(keyVarObs(f,aidx));
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
        obsCov(f) = cov.block(0,0,thisdim-1,thisdim-1);
      }else{
        nllVec(f).setSigma(cov);
        obsCov(f) = cov;
      }
    }
  }
  //eval likelihood 
  for(int y=0;y<noYears;y++){
    int totalParKey = 0;
    for(int f=0;f<noFleets;f++){
      if(fleetTypes(f)!=5){
        if(!isNAINT(idx1(f,y))){
          int idxfrom=idx1(f,y);
          int idxlength=idx2(f,y)-idx1(f,y)+1;

          vector<Type> currentVar=nllVec(f).cov().diagonal();
          vector<Type> sqrtW(currentVar.size());

  	  switch(obsLikelihoodFlag(f)){
	  case 0: // (LN) log-Normal distribution
            
            for(int idxV=0; idxV<currentVar.size(); ++idxV){
              if(isNA(weight(idxfrom+idxV))){
                sqrtW(idxV)=Type(1.0);
              }else{
                if(fixVarToWeight==1){ 
                  sqrtW(idxV)=sqrt(weight(idxfrom+idxV)/currentVar(idxV));
                }else{
                  sqrtW(idxV)=sqrt(Type(1)/weight(idxfrom+idxV));
                }
              }
            }

	    nll += nllVec(f)((logobs.segment(idxfrom,idxlength)-predObs.segment(idxfrom,idxlength))/sqrtW,keep.segment(idxfrom,idxlength));
            nll += (log(sqrtW)*keep.segment(idxfrom,idxlength)).sum();
	    SIMULATE_F(of){
	      logobs.segment(idxfrom,idxlength) = predObs.segment(idxfrom,idxlength) + (nllVec(f).simulate()*sqrtW);
	    }
	    break;
	  case 1: // (ALN) Additive logistic-normal proportions + log-normal total numbers
	    nll +=  nllVec(f)(addLogratio((vector<Type>)logobs.segment(idxfrom,idxlength))-addLogratio((vector<Type>)predObs.segment(idxfrom,idxlength)));
	    nll += log(log2proportion((vector<Type>)logobs.segment(idxfrom,idxlength))).sum();
	    nll -= dnorm(log(log2expsum((vector<Type>)logobs.segment(idxfrom,idxlength))),
	     	         log(log2expsum((vector<Type>)predObs.segment(idxfrom,idxlength))),
	   	         exp(logSdLogTotalObs(totalParKey++)),true);
	    nll += log(log2expsum((vector<Type>)logobs.segment(idxfrom,idxlength)));
	    nll -= log(abs(jacobianDet((vector<Type>)logobs.segment(idxfrom,idxlength).exp())));
            nll -= logobs.segment(idxfrom,idxlength).sum();
	    SIMULATE_F(of){
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
      }else{ //fleetTypes(f)==5     
        if(!isNAINT(idx1(f,y))){    
          for(int i=idx1(f,y); i<=idx2(f,y); ++i){
            nll += -keep(i)*dnbinom(logobs(i),predObs(i)*recapturePhiVec(i)/(Type(1.0)-recapturePhiVec(i)),recapturePhiVec(i),true);
            SIMULATE_F(of){
	      logobs(i) = rnbinom(predObs(i)*recapturePhiVec(i)/(Type(1.0)-recapturePhiVec(i)),recapturePhiVec(i));
            }
          }
        }   
      }   
    }  
  }
  //if(isDouble<Type>::value && of->current_parallel_region<0) {  
  //  defineVar(install("obsCov"),asSEXP_protect(obsCov),of->report);
  //  UNPROTECT(1);
  //}
  REPORT_F(obsCov,of)
  return nll;
}
