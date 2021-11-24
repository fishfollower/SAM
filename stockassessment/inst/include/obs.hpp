#include <algorithm>

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
Type nllObs(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, array<Type> &logF,
	    //vector<Type> &predObs, vector<Type> &varLogCatch,
	    data_indicator<vector<Type>,Type> &keep, objective_function<Type> *of){
  using CppAD::abs;
  Type nll=0;

  // Calculate values to report
  vector<Type> ssb = ssbFun(dat, conf, logN, logF);
  vector<Type> logssb = log(ssb);

  vector<Type> fsb = fsbFun(dat, conf, logN, logF);
  vector<Type> logfsb = log(fsb);

  vector<Type> cat = catchFun(dat, conf, logN, logF);
  vector<Type> logCatch = log(cat);

  matrix<Type> catAge = catchFunAge(dat, conf, logN, logF);
  matrix<Type> logCatchAge = catAge.array().log().matrix();

  vector<Type> land = landFun(dat, conf, logN, logF);
  vector<Type> logLand = log(land);

  vector<Type> varLogCatch = varLogCatchFun(dat, conf, logN, logF, par);

  vector<Type> varLogLand = varLogLandFun(dat, conf, logN, logF, par);

  vector<Type> tsb = tsbFun(dat, conf, logN);
  vector<Type> logtsb = log(tsb);

  vector<Type> R = rFun(logN);
  vector<Type> logR = log(R);  

  vector<Type> fbar = fbarFun(conf, logF);
  vector<Type> logfbar = log(fbar);

  vector<Type> fbarL = landFbarFun(dat, conf, logF);
  vector<Type> logfbarL = log(fbarL);

  vector<Type> predObs = predObsFun(dat, conf, par, logN, logF, logssb, logtsb, logfsb, logCatch, logLand);
  
  // setup obs likelihoods
  vector< MVMIX_t<Type> >  nllVec(dat.noFleets);
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
  vector<Type> logitRecapturePhiVec(dat.nobs);
  if(par.logitRecapturePhi.size()>0){
    recapturePhi=invlogit(par.logitRecapturePhi);
    for(int j=0; j<dat.nobs; ++j){
      if(!isNAINT(dat.aux(j,7))){
        recapturePhiVec(j)=recapturePhi(dat.aux(j,7)-1);
	logitRecapturePhiVec(j) = par.logitRecapturePhi(dat.aux(j,7)-1);
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
    if(!((dat.fleetTypes(f)==5)||(dat.fleetTypes(f)==3))){ 
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
      } else { 
        Rf_error("Unkown obsCorStruct code"); 
      }
      nllVec(f).setSigma(cov, Type(conf.fracMixObs(f)));
      obsCov(f) = cov;
    }else{
      matrix<Type> dummy(1,1);
      dummy(0,0) = R_NaReal;
      obsCov(f) = dummy;
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

          vector<Type> currentVar=nllVec(f).cov().diagonal();
          vector<Type> sqrtW(currentVar.size());

  	      switch(conf.obsLikelihoodFlag(f)){
    	      case 0: // (LN) log-Normal distribution
              
              for(int idxV=0; idxV<currentVar.size(); ++idxV){
                if(isNA(dat.weight(idxfrom+idxV))){
                  sqrtW(idxV)=Type(1.0);
                  int a = dat.aux(idxfrom+idxV,2)-conf.minAge;
                  if(conf.predVarObsLink(f,a)>(-1)){
                    sqrtW(idxV) = sqrt(findLinkV(par.logSdLogObs(conf.keyVarObs(f,a))+(exp(par.predVarObs(conf.predVarObsLink(f,a))) -Type(1))*predObs(idxfrom+idxV))/currentVar(idxV));
                  }
		  for(int idxXtraSd=0; idxXtraSd<(conf.keyXtraSd).rows(); ++idxXtraSd){
		    int realfleet=f+1;
		    int realyear=y+CppAD::Integer(min(dat.years));
		    int realage=dat.aux(idxfrom+idxV,2);		    
		    vector<int> fyac=conf.keyXtraSd.row(idxXtraSd);
		    if((realfleet==fyac(0))&&(realyear==fyac(1))&&(realage==fyac(2))){
                      sqrtW(idxV)=exp(par.logXtraSd(fyac(3)));
		      break;
		    }
		  }
                }else{
                  if(conf.fixVarToWeight==1){
                    sqrtW(idxV)=sqrt(dat.weight(idxfrom+idxV)/currentVar(idxV));
                  }else{
                    sqrtW(idxV)=sqrt(Type(1)/dat.weight(idxfrom+idxV));
                  }
                }
              }
              if(isNAINT(dat.idxCor(f,y))){
      	        nll += nllVec(f)((dat.logobs.segment(idxfrom,idxlength)-predObs.segment(idxfrom,idxlength))/sqrtW,keep.segment(idxfrom,idxlength));
                nll += (log(sqrtW)*keep.segment(idxfrom,idxlength)).sum();
      	        SIMULATE_F(of){
    	            dat.logobs.segment(idxfrom,idxlength) = predObs.segment(idxfrom,idxlength) + (nllVec(f).simulate()*sqrtW);
  	            }
  	          }else{
  	            int thisdim=currentVar.size();
  	            matrix<Type> thiscor=dat.corList(dat.idxCor(f,y));
  	            matrix<Type> thiscov(thisdim,thisdim);
  	            for(int r=0;r<thisdim;++r){
  	              for(int c=0;c<thisdim;++c){
  	                thiscov(r,c)=thiscor(r,c)*sqrt(currentVar(r)*currentVar(c));
  	              }
	      } 
	      MVMIX_t<Type> thisnll(thiscov,Type(conf.fracMixObs(f)));
	      nll+= thisnll((dat.logobs.segment(idxfrom,idxlength)-predObs.segment(idxfrom,idxlength))/sqrtW, keep.segment(idxfrom,idxlength));              
	      nll+= (log(sqrtW)*keep.segment(idxfrom,idxlength)).sum();
	      SIMULATE_F(of){
	        dat.logobs.segment(idxfrom,idxlength) = predObs.segment(idxfrom,idxlength) + thisnll.simulate()*sqrtW;
	      }
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
	    Rf_error("Unknown obsLikelihoodFlag");
	  }
        }
      }else{ //dat.fleetTypes(f)==5
        if(dat.fleetTypes(f)==5){
          if(!isNAINT(dat.idx1(f,y))){    
            for(int i=dat.idx1(f,y); i<=dat.idx2(f,y); ++i){
              //nll += -keep(i)*dnbinom(dat.logobs(i),predObs(i)*recapturePhiVec(i)/(Type(1.0)-recapturePhiVec(i)),recapturePhiVec(i),true);
	      Type log_mu = log(predObs(i));
	      Type log_var_minus_mu = log_mu - logitRecapturePhiVec(i);	      
	      nll += -keep(i)*dnbinom_robust(dat.logobs(i),log_mu,log_var_minus_mu,true);
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
                  if(conf.keyBiomassTreat(f)==4){
                    sd = sqrt(varLogLand(y));
                  }else{
                    sd = exp(par.logSdLogObs(conf.keyVarObs(f,0)));
                  }
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

  SIMULATE_F(of) {
    REPORT_F(logF,of);
    REPORT_F(logN,of);
    vector<Type> logobs=dat.logobs; 
    REPORT_F(logobs,of);
  }

  REPORT_F(obsCov,of);
  REPORT_F(predObs,of);
  ADREPORT_F(logssb,of);
  ADREPORT_F(logfbar,of);
  ADREPORT_F(logCatch,of);
  ADREPORT_F(logLand,of);
  ADREPORT_F(logtsb,of);

  // Additional forecast quantities
  if(dat.forecast.nYears > 0){
    vector<Type> dis = disFun(dat, conf, logN, logF);
    vector<Type> logDis = log(dis);
    
    vector<Type> landFbar = landFbarFun(dat, conf, logF);
    vector<Type> loglandfbar = log(landFbar);
    
    vector<Type> disFbar = disFbarFun(dat, conf, logF);
    vector<Type> logdisfbar = log(disFbar);

    ADREPORT_F(logDis,of);
    ADREPORT_F(loglandfbar,of);
    ADREPORT_F(logdisfbar,of);      
    SIMULATE_F(of) {
      //if(dat.forecast.simFlag[0] == 0 || dat.forecast.simFlag[1] == 0){
	REPORT_F(logssb,of);
	REPORT_F(logfbar,of);
	REPORT_F(logfbarL,of);
	REPORT_F(logCatch,of);
	REPORT_F(logCatchAge,of);
	REPORT_F(logLand,of);
	REPORT_F(logtsb,of);      
	//}
    }
  }

  vector<Type> logLagR(logR.size());
  for(int i=0; i<logR.size(); ++i){
    logLagR(i) = logN(1,i);
  }
  ADREPORT_F(logR,of);
  ADREPORT_F(logLagR,of);

  int timeSteps=logF.dim[1];
  
  vector<Type> lastLogN = logN.col(timeSteps-1);
  ADREPORT_F(lastLogN,of);
  vector<Type> lastLogF = logF.col(timeSteps-1);
  ADREPORT_F(lastLogF,of);  

  vector<Type> beforeLastLogN = logN.col(timeSteps-2);
  ADREPORT_F(beforeLastLogN,of);
  vector<Type> beforeLastLogF = logF.col(timeSteps-2);
  ADREPORT_F(beforeLastLogF,of);  

  if(dat.forecast.nYears > 0 && dat.forecast.FModel(dat.forecast.FModel.size()-1) == dat.forecast.findMSY){
    
    int catchYears = std::min((int)asDouble(dat.forecast.nYears),dat.forecast.nCatchAverageYears);
    Type catchSum = sum((vector<Type>)cat.tail(catchYears)) / (Type)catchYears;
    nll -= par.implicitFunctionDelta * log(catchSum);    
    // Calculate Fmsy
    Type logFMSY = par.logFScaleMSY + logfbar(timeSteps - dat.forecast.nYears - 1);
    ADREPORT_F(logFMSY, of);

    // Output stock status - Positive is good for the stock
    Type logFstatus = logFMSY - logfbar(timeSteps - dat.forecast.nYears - 1);
    Type logSSBstatus = logssb(timeSteps - dat.forecast.nYears - 1) - logssb(timeSteps - 1);
    ADREPORT_F(logFstatus, of);
    ADREPORT_F(logSSBstatus, of);    
  }

  
  return nll;
}
