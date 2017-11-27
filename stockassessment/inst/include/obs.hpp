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
      if(rhoMap(i-1+offset)>maxrm) maxrm=rhoMap(i-1+offset);
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
}

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
Type nllObs(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, array<Type> &logF, vector<Type> &predObs, vector<Type> &varLogCatch, data_indicator<vector<Type>,Type> &keep, objective_function<Type> *of){
  using CppAD::abs;
  Type nll=0; 
  // setup obs likelihoods
  vector< density::MVNORM_t<Type> >  nllVec(dat.noFleets);
  vector< density::UNSTRUCTURED_CORR_t<Type> > neg_log_densityObsUnstruc(dat.noFleets);
  vector< vector<Type> > obsCovScaleVec(dat.noFleets);
  vector<Type> varLogObs=exp(par.logSdLogObs*Type(2.0));
  vector<Type> IRARdist(par.transfIRARdist.size()); //[ d_1, d_2, ...,d_N-1 ]
  if(par.transfIRARdist.size()>0) IRARdist=exp(par.transfIRARdist);
  vector< vector<Type> > sigmaObsParVec(dat.noFleets);
  int aidx;
  vector< matrix<Type> > obsCov(dat.noFleets); // for reporting
 
  vector<Type> recapturePhi(par.logitRecapturePhi.size());
  vector<Type> recapturePhiVec(dat.nobs);
  if(par.logitRecapturePhi.size()>0){
    recapturePhi=invlogit(par.logitRecapturePhi);
    for(int j=0; j<dat.nobs; ++j){
      if(!isNAINT(dat.aux(j,7))){
        recapturePhiVec(j)=recapturePhi(dat.aux(j,7)-1);
      }
    }
  }

  int nfleet = dat.maxAgePerFleet(0)-dat.minAgePerFleet(0)+1;
  int dn=nfleet*(nfleet-1)/2;
  int from=-dn, to=-1; 
  for(int f=0; f<dat.noFleets; f++){
    if(conf.obsCorStruct(f)!=2) continue; // skip if not US 
    nfleet = dat.maxAgePerFleet(f)-dat.minAgePerFleet(f)+1;
    if(conf.obsLikelihoodFlag(f) == 1) nfleet-=1; // ALN has dim-1
    dn = nfleet*(nfleet-1)/2;
    from=to+1;
    to=to+dn;
    vector<Type> tmp(dn);
    for(int i=from; i<=to; i++) tmp(i-from) = par.sigmaObsParUS(i);
    sigmaObsParVec(f) = tmp; 
  }

  for(int f=0; f<dat.noFleets; ++f){
    if(!((dat.fleetTypes(f)==5)||(dat.fleetTypes(f)==3)||(dat.fleetTypes(f)==7))){ 
      int thisdim=dat.maxAgePerFleet(f)-dat.minAgePerFleet(f)+1;
      if(conf.obsLikelihoodFlag(f) == 1) thisdim-=1; // ALN has dim-1
      matrix<Type> cov(thisdim,thisdim);
      cov.setZero();
      if(conf.obsCorStruct(f)==0){//ID (independent)  
        for(int i=0; i<thisdim; ++i){
          aidx = i+dat.minAgePerFleet(f)-conf.minAge;
  	  cov(i,i)=varLogObs(conf.keyVarObs(f,aidx));
        }
      } else if(conf.obsCorStruct(f)==1){//(AR) irregular lattice AR
        cov = setupVarCovMatrix(conf.minAge, conf.maxAge, dat.minAgePerFleet(f), dat.maxAgePerFleet(f), conf.keyCorObs.transpose().col(f), IRARdist, conf.keyVarObs.transpose().col(f) , exp(par.logSdLogObs) );
	if(conf.obsLikelihoodFlag(f) == 1){ // ALN has dim-1
	  cov.conservativeResize(thisdim,thisdim); // resize, keep contents but drop last row/col
	}

      } else if(conf.obsCorStruct(f)==2){//(US) unstructured
        neg_log_densityObsUnstruc(f) = getCorrObj(sigmaObsParVec(f));  
        matrix<Type> tmp = neg_log_densityObsUnstruc(f).cov();
  
        tmp.setZero();
        int offset = dat.minAgePerFleet(f)-conf.minAge;
        obsCovScaleVec(f).resize(tmp.rows());
        for(int i=0; i<tmp.rows(); i++) {
	  tmp(i,i) = sqrt(varLogObs(conf.keyVarObs(f,i+offset)));
	  obsCovScaleVec(f)(i) = tmp(i,i);
        }
        cov  = tmp*matrix<Type>(neg_log_densityObsUnstruc(f).cov()*tmp);
      } else { error("Unkown obsCorStruct code"); }
        nllVec(f).setSigma(cov);
        obsCov(f) = cov;
    }else{
      matrix<Type> dummy(1,1);
      dummy(0,0) = R_NaReal;
      obsCov(f) = dummy;
    }
  }
  for(int f=0; f<dat.noFleets; ++f){ // separate loop to insure other fleet's cov has been setup 
    if(dat.fleetTypes(f)==7){ 
      int thisdim=dat.maxAgePerFleet(f)-dat.minAgePerFleet(f)+1;
      matrix<Type> cov(thisdim,thisdim);
      cov.setIdentity(); // place holder
      nllVec(f).setSigma(cov);
      obsCov(f) = cov;
    }
  }
  //eval likelihood 
  for(int y=0;y<dat.noYears;y++){
    int totalParKey = 0;
    for(int f=0;f<dat.noFleets;f++){
      if(!((dat.fleetTypes(f)==5)||(dat.fleetTypes(f)==3))){ 
        if(!isNAINT(dat.idx1(f,y))){
          int idxfrom=dat.idx1(f,y);
          int idxlength=dat.idx2(f,y)-dat.idx1(f,y)+1;

          // ----------------if sum fleet need to update covariance matrix
          if(dat.fleetTypes(f)==7){
            array<Type> totF=totFFun(conf, logF);             
            Type zz=0;
            int thisDim=dat.maxAgePerFleet(f)-dat.minAgePerFleet(f)+1;
            int Nparts=0;
            int offset=0;
            for(int ff=0; ff<dat.noFleets; ++ff){
              if(dat.sumKey(f,ff)==1){++Nparts;}
            }
            matrix<Type> muMat(thisDim,Nparts);
            vector<Type> muSum(thisDim);
            muMat.setZero();
            muSum.setZero();
            matrix<Type> VV(thisDim*Nparts,thisDim*Nparts);
            VV.setZero();
            matrix<Type> G(thisDim*Nparts,thisDim); // Gradient 
            G.setZero();
            matrix<Type> combiCov(thisDim,thisDim); 
            combiCov.setZero();
            int element=-1;
            for(int ff=0; ff<dat.noFleets; ++ff){
              if(dat.sumKey(f,ff)==1){
                ++element; 
                for(int aa=dat.minAgePerFleet(ff); aa<=dat.maxAgePerFleet(ff); ++aa){
                  zz = dat.natMor(y,aa-conf.minAge)+totF(aa-conf.minAge,y);
                  muMat(aa-dat.minAgePerFleet(f),element)=exp(logN(aa-conf.minAge,y)-log(zz)+log(1-exp(-zz))+logF(conf.keyLogFsta(ff,aa-conf.minAge),y));
                  muSum(aa-dat.minAgePerFleet(f)) += muMat(aa-dat.minAgePerFleet(f),element);
                }
                offset=dat.minAgePerFleet(ff)-dat.minAgePerFleet(f);
                VV.block(thisDim*element+offset,thisDim*element+offset,thisDim-offset,thisDim-offset)=nllVec(ff).cov(); // possibly too simple 
              }
            }
            element=-1;
            for(int ff=0; ff<dat.noFleets; ++ff){
              if(dat.sumKey(f,ff)==1){
                ++element;
                for(int aa=0; aa<thisDim; ++aa){
                  G(aa+element*thisDim,aa) = muMat(aa,element)/muSum(aa);
                } 
              }
            }
            combiCov = G.transpose()*VV*G;
            nllVec(f).setSigma(combiCov);
          }
          // ----------------updating of covariance matrix done

          vector<Type> currentVar=nllVec(f).cov().diagonal();
          vector<Type> sqrtW(currentVar.size());

  	  switch(conf.obsLikelihoodFlag(f)){
	  case 0: // (LN) log-Normal distribution
            
            for(int idxV=0; idxV<currentVar.size(); ++idxV){
              if(isNA(dat.weight(idxfrom+idxV))){
                sqrtW(idxV)=Type(1.0);
              }else{
                if(conf.fixVarToWeight==1){ 
                  sqrtW(idxV)=sqrt(dat.weight(idxfrom+idxV)/currentVar(idxV));
                }else{
                  sqrtW(idxV)=sqrt(Type(1)/dat.weight(idxfrom+idxV));
                }
              }
            }

	    nll += nllVec(f)((dat.logobs.segment(idxfrom,idxlength)-predObs.segment(idxfrom,idxlength))/sqrtW,keep.segment(idxfrom,idxlength));
            nll += (log(sqrtW)*keep.segment(idxfrom,idxlength)).sum();
	    SIMULATE_F(of){
	      dat.logobs.segment(idxfrom,idxlength) = predObs.segment(idxfrom,idxlength) + (nllVec(f).simulate()*sqrtW);
	    }
	    break;
	  case 1: // (ALN) Additive logistic-normal proportions + log-normal total numbers
	    nll +=  nllVec(f)(addLogratio((vector<Type>)dat.logobs.segment(idxfrom,idxlength))-addLogratio((vector<Type>)predObs.segment(idxfrom,idxlength)));
	    nll += log(log2proportion((vector<Type>)dat.logobs.segment(idxfrom,idxlength))).sum();
	    nll -= dnorm(log(log2expsum((vector<Type>)dat.logobs.segment(idxfrom,idxlength))),
	     	         log(log2expsum((vector<Type>)predObs.segment(idxfrom,idxlength))),
	   	         exp(par.logSdLogTotalObs(totalParKey++)),true);
	    nll += log(log2expsum((vector<Type>)dat.logobs.segment(idxfrom,idxlength)));
	    nll -= log(abs(jacobianDet((vector<Type>)dat.logobs.segment(idxfrom,idxlength).exp())));
            nll -= dat.logobs.segment(idxfrom,idxlength).sum();
	    SIMULATE_F(of){
	      vector<Type> logProb(idxlength);
	      logProb.setZero();
	      logProb.segment(0,idxlength-1) = addLogratio(((vector<Type>)predObs.segment(idxfrom,idxlength))) + nllVec(f).simulate();
	      Type logDenom = logExpSum(logProb);
	      logProb -= logDenom;
	      Type logTotal = rnorm(log(log2expsum((vector<Type>)predObs.segment(idxfrom,idxlength))),
				    exp(par.logSdLogTotalObs(totalParKey++)));
	      dat.logobs.segment(idxfrom,idxlength) = logProb + logTotal; 
	    }
	    break;
  	    default:
	      error("Unknown obsLikelihoodFlag");
	  }
        }
      }else{ //dat.fleetTypes(f)==5
        if(dat.fleetTypes(f)==5){
          if(!isNAINT(dat.idx1(f,y))){    
            for(int i=dat.idx1(f,y); i<=dat.idx2(f,y); ++i){
              nll += -keep(i)*dnbinom(dat.logobs(i),predObs(i)*recapturePhiVec(i)/(Type(1.0)-recapturePhiVec(i)),recapturePhiVec(i),true);
              SIMULATE_F(of){
	        dat.logobs(i) = rnbinom(predObs(i)*recapturePhiVec(i)/(Type(1.0)-recapturePhiVec(i)),recapturePhiVec(i));
              }
            }
          }
        }else{
          if(dat.fleetTypes(f)==3){
            Type sd=0;
            if(!isNAINT(dat.idx1(f,y))){
              for(int i=dat.idx1(f,y); i<=dat.idx2(f,y); ++i){
                if(conf.keyBiomassTreat(f)==3){
                  sd = sqrt(varLogCatch(y));
                }else{
                  sd = exp(par.logSdLogObs(conf.keyVarObs(f,0)));
                }  
                nll += -keep(i)*dnorm(dat.logobs(i),predObs(i),sd,true);
                SIMULATE_F(of){
  	          dat.logobs(i) = rnorm(predObs(i),sd);
                }
              }
            }    
          }
        }   
      }   
    }  
  }
  REPORT_F(obsCov,of)
  return nll;
}
