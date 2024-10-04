SAM_DEPENDS(convenience)
SAM_DEPENDS(define)
SAM_DEPENDS(recruitment)
SAM_DEPENDS(incidence)
SAM_DEPENDS(forecast)
SAM_DEPENDS(mvmix)
SAM_DEPENDS(derived)
SAM_DEPENDS(survival)
SAM_DEPENDS(predobs)
SAM_DEPENDS(reproductive)
SAM_DEPENDS(equilibrium)
SAM_DEPENDS(empirical_pr)
SAM_DEPENDS(dirichlet)

template <class Type>
matrix<Type> setupVarCovMatrix(int dim, int offset, vector<int> rhoMap, vector<Type> rhoVec, vector<int> sdMap, vector<Type> sdVec)SOURCE({

    // int dim = maxAgeFleet-minAgeFleet+1;
    // int offset = minAgeFleet-minAge;
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
	  Type dist = fabs(xvec(i)-xvec(j));
	  ret(i,j)=pow( rho0,dist)*sdVec( sdMap(i+offset) )*sdVec( sdMap(j+offset));
	} else if(i==j) ret(i,j) = sdVec( sdMap(i+offset) )*sdVec( sdMap(j+offset));
      }
    return ret;
  });

SAM_SPECIALIZATION(matrix<double> setupVarCovMatrix(int, int, vector<int>, vector<double>, vector<int>, vector<double>));
SAM_SPECIALIZATION(matrix<TMBad::ad_aug> setupVarCovMatrix(int, int,vector<int>, vector<TMBad::ad_aug>, vector<int>, vector<TMBad::ad_aug>));



template <class Type> 
density::UNSTRUCTURED_CORR_t<Type> getCorrObj(vector<Type> params)SOURCE({
    density::UNSTRUCTURED_CORR_t<Type> ret(params);
    return ret;
  });

SAM_SPECIALIZATION(density::UNSTRUCTURED_CORR_t<double> getCorrObj(vector<double>));
SAM_SPECIALIZATION(density::UNSTRUCTURED_CORR_t<TMBad::ad_aug> getCorrObj(vector<TMBad::ad_aug>));

  


template <class Type>
vector< MVMIX_t<Type> > getnllVec(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par,
				  objective_function<Type> *of)SOURCE({
				      // setup obs likelihoods
				      vector< MVMIX_t<Type> >  nllVec(dat.noFleets);
				      vector< density::UNSTRUCTURED_CORR_t<Type> > neg_log_densityObsUnstruc(dat.noFleets);
				      vector< vector<Type> > obsCovScaleVec(dat.noFleets);
				      vector<Type> varLogObs=exp(par.logSdLogObs*Type(2.0));
				      vector<Type> IRARdist(par.transfIRARdist.size()); //[ d_1, d_2, ...,d_N-1 ]
				      if(par.transfIRARdist.size()>0) IRARdist=exp(par.transfIRARdist);
				      vector< vector<Type> > sigmaObsParVec(dat.noFleets);
				      vector< matrix<Type> > obsCov(dat.noFleets); // for reporting

				      int nfleet = dat.fleetCovarianceSize(0); //dat.maxAgePerFleet(0)-dat.minAgePerFleet(0)+1;
				      int dn=nfleet*(nfleet-1)/2;
				      int from=-dn, to=-1; 
				      for(int f=0; f<dat.noFleets; f++){
					if(conf.obsCorStruct(f)!=2) continue; // skip if not US 
					nfleet = dat.fleetCovarianceSize(f); //dat.maxAgePerFleet(f)-dat.minAgePerFleet(f)+1;
					if(conf.obsLikelihoodFlag(f) == 1 && dat.fleetTypes(f) < 80) nfleet-=1; // ALN has dim-1 (but not for proportion fleets)
					dn = nfleet*(nfleet-1)/2;
					from=to+1;
					to=to+dn;
					vector<Type> tmp(dn);
					for(int i=from; i<=to; i++) tmp(i-from) = par.sigmaObsParUS(i);
					sigmaObsParVec(f) = tmp; 
				      }

				      for(int f=0; f<dat.noFleets; ++f){
					if(!((dat.fleetTypes(f)==5)||(dat.fleetTypes(f)==3)||(dat.fleetTypes(f)==7)||(dat.fleetTypes(f)==6))){
					  //||(dat.fleetTypes(f)==80)||(dat.fleetTypes(f)==90)||(dat.fleetTypes(f)==92))){ // Maybe easier to switch to inclusive instead of exclusive test 
					  int thisdim= dat.fleetCovarianceSize(f); //dat.maxAgePerFleet(f)-dat.minAgePerFleet(f)+1;
					  int offset = 0;
					  if(dat.minAgePerFleet(f) > conf.minAge)
					    offset = dat.minAgePerFleet(f) - conf.minAge;
					  if(conf.obsLikelihoodFlag(f) == 1 && dat.fleetTypes(f) < 80) thisdim-=1; // ALN has dim-1 (but not for proportion fleets)
					  if(conf.obsLikelihoodFlag(f) == 2){
					    matrix<Type> dummy(1,1);
					    dummy(0,0) = R_NaReal;
					    obsCov(f) = dummy;
					    continue;
					  }
					  matrix<Type> cov(thisdim,thisdim);
					  cov.setZero();
					  if(conf.obsCorStruct(f)==0){//ID (independent)  
					    for(int i=0; i<thisdim; ++i){
					      int aidx = i+offset;
					      cov(i,i)=varLogObs(conf.keyVarObs(f,aidx));
					    }
					  }else if(conf.obsCorStruct(f)==1){//(AR) irregular lattice AR
					     cov = setupVarCovMatrix(dat.fleetCovarianceSize(f), offset, conf.keyCorObs.transpose().col(f), IRARdist, conf.keyVarObs.transpose().col(f) , exp(par.logSdLogObs) );
					    if(conf.obsLikelihoodFlag(f) == 1){ // ALN has dim-1
					      cov.conservativeResize(thisdim,thisdim); // resize, keep contents but drop last row/col
					    }

					  } else if(conf.obsCorStruct(f)==2){//(US) unstructured
					    neg_log_densityObsUnstruc(f) = getCorrObj(sigmaObsParVec(f));  
					    matrix<Type> tmp = neg_log_densityObsUnstruc(f).cov();
  
					    tmp.setZero();					   
					    obsCovScaleVec(f).resize(tmp.rows());
					    for(int i=0; i<tmp.rows(); i++) {
					      tmp(i,i) = sqrt(varLogObs(conf.keyVarObs(f,i+offset)));
					      obsCovScaleVec(f)(i) = tmp(i,i);
					    }
					    cov  = tmp*matrix<Type>(neg_log_densityObsUnstruc(f).cov()*tmp);
					  } else { 
					    Rf_error("Unknown obsCorStruct code"); 
					  }
					  nllVec(f).setSigma(cov, Type(conf.fracMixObs(f)));
					  obsCov(f) = cov;
					}else{
					  matrix<Type> dummy(1,1);
					  dummy(0,0) = R_NaReal;
					  obsCov(f) = dummy;
					}
					
				      }
				      REPORT_F(obsCov,of);
				      return nllVec;
				    }
				    );

SAM_SPECIALIZATION(vector< MVMIX_t<double> > getnllVec(dataSet<double>&, confSet&, paraSet<double> &,objective_function<double>*));
SAM_SPECIALIZATION(vector< MVMIX_t<TMBad::ad_aug> > getnllVec(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug> &,objective_function<TMBad::ad_aug>*));

#ifndef WITH_LIB_SAM
namespace obs_fun {

  template <class Type>
  Type log2expsum(vector<Type> x)SOURCE({
    return exp(x).sum();
    })

  SAM_SPECIALIZATION(double log2expsum(vector<double>));
  SAM_SPECIALIZATION(TMBad::ad_aug log2expsum(vector<TMBad::ad_aug>));

  template<class Type>
  Type logExpSum(vector<Type> x)SOURCE({
    Type m = max(x);
    return m + log(exp(x-m).sum());
    })

  SAM_SPECIALIZATION(double logExpSum(vector<double>));
  SAM_SPECIALIZATION(TMBad::ad_aug logExpSum(vector<TMBad::ad_aug>));

  
  template <class Type>
  vector<Type> log2proportion(vector<Type> x)SOURCE({
    return exp(x) / log2expsum(x);
    })

  SAM_SPECIALIZATION(vector<double> log2proportion(vector<double>));
  SAM_SPECIALIZATION(vector<TMBad::ad_aug> log2proportion(vector<TMBad::ad_aug>));


  template <class Type>
  vector<Type> addLogratio(vector<Type> logx)SOURCE({
    int n = logx.size();
    vector<Type> res(n-1);
    for(int i = 0; i < res.size(); ++i)
      res(i) = logx(i) - logx(n-1);
    return res;//log(x.head(x.size()-1)/x.tail(1));
    })

  SAM_SPECIALIZATION(vector<double> addLogratio(vector<double>));
  SAM_SPECIALIZATION(vector<TMBad::ad_aug> addLogratio(vector<TMBad::ad_aug>));

  template<class Type>
  vector<Type> multLogratio(vector<Type> logx)SOURCE({
    vector<Type> res(logx.size()-1);
    for(int i = 0; i < res.size(); ++i)
      res(i) = logx(i)-log(Type(1.0)-exp(logExpSum((vector<Type>)logx.head(i+1))));
    return res;
    })
  
  SAM_SPECIALIZATION(vector<double> multLogratio(vector<double>));
  SAM_SPECIALIZATION(vector<TMBad::ad_aug> multLogratio(vector<TMBad::ad_aug>));

 
  template<class Type>
  matrix<Type> buildJac(vector<Type> x, vector<Type> w)SOURCE({
    matrix<Type> res(x.size(),x.size()); 
    Type xs = x.sum();
    Type xsp = pow(xs,2);
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
    })


  SAM_SPECIALIZATION(matrix<double> buildJac(vector<double>,vector<double>));
  SAM_SPECIALIZATION(matrix<TMBad::ad_aug> buildJac(vector<TMBad::ad_aug>,vector<TMBad::ad_aug>));

  
  template <class Type>
  Type jacobianDet(vector<Type> x)SOURCE({
    vector<Type> w(x.size());
    w.fill(Type(1.0));
    return buildJac(x,w).determinant();
    })

  SAM_SPECIALIZATION(double jacobianDet(vector<double>));
  SAM_SPECIALIZATION(TMBad::ad_aug jacobianDet(vector<TMBad::ad_aug>));

  
  template <class Type>
  Type jacobianDet(vector<Type> x,vector<Type> w)SOURCE({
    return buildJac(x,w).determinant();
    })

  SAM_SPECIALIZATION(double jacobianDet(vector<double>,vector<double>));
  SAM_SPECIALIZATION(TMBad::ad_aug jacobianDet(vector<TMBad::ad_aug>,vector<TMBad::ad_aug>));


 // template<class Type>
 //  matrix<Type> buildJacProportions(vector<Type> x)SOURCE({
 //    matrix<Type> res(x.size(),x.size()); 
 //    for(int i = 0; i < res.cols(); ++i)
 //      res(i,i) = Type(1.0)/x(i);
 //    for(int j = 0; j < res.cols(); ++j){
 //      res(res.rows()-1,j) = Type(-1.0)/x(x.size()-1);
 //    }
 //    return res;
 //    })


 //  SAM_SPECIALIZATION(matrix<double> buildJacProportions(vector<double>));
 //  SAM_SPECIALIZATION(matrix<TMBad::ad_aug> buildJacProportions(vector<TMBad::ad_aug>));

  
 //  template <class Type>
 //  Type jacobianDetProportions(vector<Type> x)SOURCE({
 //    return buildJacProportions(x).determinant();
 //    })

 //  SAM_SPECIALIZATION(double jacobianDetProportions(vector<double>));
 //  SAM_SPECIALIZATION(TMBad::ad_aug jacobianDetProportions(vector<TMBad::ad_aug>));

  

  template <class Type>
  Type findLinkV(Type k, int n DEFARG(=0))SOURCE({
    // small helper function to solve exp(v)-exp(k-v/2)-1=0 for v
    Type v = log(exp(k)+Type(1));
    for(int i=0; i<n; ++i){
      v -= (exp(v)-exp(k-0.5*v)-1.0)/(exp(v)+.5*exp(k-0.5*v));
    }
    return v;
    })

  SAM_SPECIALIZATION(double findLinkV(double,int));
  SAM_SPECIALIZATION(TMBad::ad_aug findLinkV(TMBad::ad_aug,int));

}
#endif

template <class Type>
Type nllObs(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, forecastSet<Type>& forecast, array<Type> &logN, array<Type> &logF, array<Type>& logP, array<Type>& logitFseason,
	    Recruitment<Type> &recruit,
	    MortalitySet<Type>& mort,
	    data_indicator<vector<Type>,Type> &keep,
	    int reportingLevel,
	    objective_function<Type> *of)
  SOURCE({
      Rcout << "From nllObs 1: " << dat.idxCor(0,0) << "\n";
      Rcout << "A\n";
      Type nll=0.0;
      // Calculate values to report
      vector<Type> logssb = ssbFun(dat, conf, logN, logF,mort, true);
      vector<Type> ssb = exp(logssb);
      Rcout << "From nllObs 2: " << dat.idxCor(0,0) << "\n";

      vector<Type> logerb = erbFun(dat, conf, par, logN, logF,mort, true);
      vector<Type> erb = exp(logerb);
      Rcout << "From nllObs 3: " << dat.idxCor(0,0) << "\n";
      vector<Type> logrelativeerb = logerb - logssb;
      vector<Type> relativeerb = exp(logrelativeerb);
      

      
      vector<Type> fsb = fsbFun(dat, conf, logN, logF,mort);
      vector<Type> logfsb = log(fsb);
      Rcout << "From nllObs 4: " << dat.idxCor(0,0) << "\n";

      vector<Type> logCatch = catchFun(dat, conf, logN, logF,mort, true);
      vector<Type> cat = exp(logCatch);
      Rcout << "From nllObs 5: " << dat.idxCor(0,0) << "\n";

      matrix<Type> logCatchAge = catchFunAge(dat, conf, logN, logF,mort, true);
      matrix<Type> catAge = logCatchAge.array().exp().matrix();
      Rcout << "From nllObs 6: " << dat.idxCor(0,0) << "\n";

      array<Type> catchByFleet = catchByFleetFun(dat, conf, logN, logF, mort);
      array<Type> logCatchByFleet = catchByFleet;
      for(int i=0; i<logCatchByFleet.dim(0); ++i){
	for(int j=0; j<logCatchByFleet.dim(1); ++j){
	  logCatchByFleet(i,j)=log(catchByFleet(i,j));
	}
      }
      Rcout << "From nllObs 7: " << dat.idxCor(0,0) << "\n";

      // For caytable
      array<Type> catchByFleetAge = catchByFleetFunAgeNum(dat, conf, logN, logF, mort);
      array<Type> logCatchByFleetAge = catchByFleetAge;
      for(int i=0; i<logCatchByFleetAge.dim(0); ++i){
	for(int j=0; j<logCatchByFleetAge.dim(1); ++j){
	  for(int k=0; k<logCatchByFleetAge.dim(2); ++k){
	    logCatchByFleetAge(i,j,k)=log(logCatchByFleetAge(i,j,k));
	  }
	}
      }
      Rcout << "From nllObs 8: " << dat.idxCor(0,0) << "\n";


      vector<Type> land = landFun(dat, conf, logN, logF, mort);
      vector<Type> logLand = log(land);
      Rcout << "From nllObs 9: " << dat.idxCor(0,0) << "\n";

      vector<Type> varLogCatch = varLogCatchFun(dat, conf, logN, logF, par, mort);
      Rcout << "From nllObs 10: " << dat.idxCor(0,0) << "\n";

      vector<Type> varLogLand = varLogLandFun(dat, conf, logN, logF, par, mort);
      Rcout << "From nllObs 11: " << dat.idxCor(0,0) << "\n";

      vector<Type> tsb = tsbFun(dat, conf, logN);
      vector<Type> logtsb = log(tsb);
      Rcout << "From nllObs 12: " << dat.idxCor(0,0) << "\n";

      vector<Type> R = rFun(logN);
      vector<Type> logR = log(R);  
      Rcout << "From nllObs 12: " << dat.idxCor(0,0) << "\n";

      vector<Type> fbar = fbarFun(dat,conf, logF);
      vector<Type> logfbar = log(fbar);
      vector<Type> logfbar_Effective = Effective_fbar(dat,conf,mort,true);
      Rcout << "From nllObs 13: " << dat.idxCor(0,0) << "\n";

      vector<Type> fbarL = landFbarFun(dat, conf, logF);
      vector<Type> logfbarL = log(fbarL);
      Rcout << "From nllObs 14: " << dat.idxCor(0,0) << "\n";

      array<Type> comps = scalePFun(conf, dat, logP);
      vector<Type> weekContrib = scaleWeekFun(par, dat, logP);
      int noYearsLAI = yearsPFun(conf,dat);
      Rcout << "B\n";
      Rcout << "From nllObs 15: " << dat.idxCor(0,0) << "\n";

      if(reportingLevel > 0){
	NOT_SIMULATE_F(of){  
	  vector<Type> logLifeExpectancy = log(lifeexpectancy(dat, conf, logF));
	  matrix<Type> logLifeExpectancyAge = lifeexpectancyAge(dat, conf, logF).array().log().matrix();
	  vector<Type> logLifeExpectancyRec = log(lifeexpectancyRec(dat, conf, logF));
	  ADREPORT_F(logLifeExpectancy,of);
	  ADREPORT_F(logLifeExpectancyRec,of);
	  ADREPORT_F(logLifeExpectancyAge,of);
	  REPORT_F(logLifeExpectancy,of);
	  REPORT_F(logLifeExpectancyRec,of);
	  REPORT_F(logLifeExpectancyAge,of);

      
	  vector<Type> logYLTF = log(yearsLostFishing(dat, conf, logF));
	  matrix<Type> logYLTFf = yearsLostFishingFleet(dat, conf, logF).array().log().matrix();
	  vector<Type> logYLTM = log(yearsLostOther(dat, conf, logF));
	  vector<Type> logYNL = log(temporaryLifeExpectancy(dat, conf, logF));
	  ADREPORT_F(logYLTF, of);
	  ADREPORT_F(logYLTFf, of);
	  ADREPORT_F(logYLTM, of);
	  ADREPORT_F(logYNL, of);	
	  REPORT_F(logYLTF, of);
	  REPORT_F(logYLTFf, of);
	  REPORT_F(logYLTM, of);
	  REPORT_F(logYNL, of);	
	  vector<Type> logrmax = log(rmax(dat,conf,par,recruit));
	  vector<Type> logGenerationLength = log(generationLength(dat,conf,par));
	  ADREPORT_F(logrmax, of);
	  ADREPORT_F(logGenerationLength, of);
 
	  vector<Type> logYPR = yieldPerRecruit(dat,conf,par,logF, true);
	  vector<Type> logSPR = spawnersPerRecruit(dat,conf,par,logF, true);
	  vector<Type> logSe = equilibriumBiomass(dat,conf,par,logF, true);
	  vector<Type> logB0 = B0(dat,conf,par,logF, true);
	  ADREPORT_F(logYPR, of);
	  ADREPORT_F(logSPR, of);
	  ADREPORT_F(logSe, of);
	  ADREPORT_F(logB0, of);

	  vector<Type> logEmpiricalSPR = empiricalSPR(dat, conf, logN, mort, true);
	  ADREPORT_F(logEmpiricalSPR,of);
	  vector<Type> logEmpiricalYPR = empiricalYPR(dat, conf, logN, mort, 0, true);
	  ADREPORT_F(logEmpiricalYPR,of);

	  
	}
      }
      Rcout << "C\n";

      vector<Type> predObs = predObsFun(dat, conf, par, logN, logF, comps, logitFseason, weekContrib, mort, logssb, logtsb, logfsb, logCatch, logLand, logfbar, noYearsLAI);
      vector< MVMIX_t<Type> > nllVec = getnllVec(dat, conf, par, of);
      Rcout << "D\n";
      Rcout << "From nllObs 16: " << dat.idxCor(0,0) << "\n";

      vector<Type> recapturePhi(par.logitRecapturePhi.size());
      vector<Type> recapturePhiVec(dat.nobs);
      vector<Type> logitRecapturePhiVec(dat.nobs);
      if(par.logitRecapturePhi.size()>0){
	recapturePhi=invlogit(par.logitRecapturePhi);
	for(int j=0; j<dat.nobs; ++j){
	  int indx = CppAD::Integer(dat.auxData(j,4));
	  if(!isNAINT(indx) && dat.fleetTypes(dat.aux(j,1)-1) == 5){
	    recapturePhiVec(j)=recapturePhi(indx-1);
	    logitRecapturePhiVec(j) = par.logitRecapturePhi(indx-1);
	  }
	}
      }
      Rcout << "E\n";
      Rcout << "From nllObs 17: " << dat.idxCor(0,0) << "\n";

      //eval likelihood
      int noYears = dat.idx1.dim(1); //dat.noYears; Also works when forecast has new data
      for(int y=0;y<noYears;y++){
	Rcout << "DATA year" << y << "\n";
	int totalParKey = 0;
	for(int f=0;f<dat.noFleets;f++){
	Rcout << "\tDATA fleet" << f << "\n";
	Rcout << "\tFrom nllObs 18: " << dat.idxCor(f,y) << "\n";
	  if(!((dat.fleetTypes(f)==5)||(dat.fleetTypes(f)==3)||(dat.fleetTypes(f)==6)||(dat.fleetTypes(f)==80)||(dat.fleetTypes(f)==90)||(dat.fleetTypes(f)==92))){
	    Rcout << "\t\t1\n";
	    if(!isNAINT(dat.idx1(f,y))){
	      Rcout << "\t\t2\n";
	      int idxfrom=dat.idx1(f,y);
	      int idxlength=dat.idx2(f,y)-dat.idx1(f,y)+1;
	      Rcout << "\t\t3\n";
	    
	      // ----------------if sum fleet need to update covariance matrix
	      if(dat.fleetTypes(f)==7){
		//array<Type> totF=totFFun(dat,conf, logF);             
		//Type zz=0;
		int thisDim=dat.fleetCovarianceSize(f); //dat.maxAgePerFleet(f)-dat.minAgePerFleet(f)+1;
		int Nparts=0;
		int offset=0;
		int setoff=0;
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
		      //zz = dat.natMor(y,aa-conf.minAge)+totF(aa-conf.minAge,y);
		      Type ci = mort.fleetCumulativeIncidence(aa-conf.minAge,y,ff);
		      muMat(aa-dat.minAgePerFleet(f),element)= exp(logN(aa-conf.minAge,y)) * ci; //exp(logN(aa-conf.minAge,y)-log(zz)+log(1-exp(-zz))+logF(conf.keyLogFsta(ff,aa-conf.minAge),y));
		      muSum(aa-dat.minAgePerFleet(f)) += muMat(aa-dat.minAgePerFleet(f),element);
		    }
		    offset=dat.minAgePerFleet(ff)-dat.minAgePerFleet(f);
		    setoff=dat.maxAgePerFleet(f)-dat.maxAgePerFleet(ff);
		    VV.block(thisDim*element+offset,thisDim*element+offset,thisDim-offset-setoff,thisDim-offset-setoff)=nllVec(ff).cov(); // possibly too simple 
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
		REPORT_F(VV,of);
		combiCov = G.transpose()*VV*G;
		REPORT_F(combiCov,of);
		nllVec(f).setSigma(combiCov);
	      }
	      Rcout << "\t\t4\n";
	      // ----------------updating of covariance matrix done
	      vector<Type> currentVar=nllVec(f).cov().diagonal();
	      vector<Type> sqrtW(currentVar.size());
	      Rcout << "\t\t5\n";
	    
	      switch(conf.obsLikelihoodFlag(f)){
	      case 0: // (LN) log-Normal distribution
		Rcout << "\t\t6\n";
		for(int idxV=0; idxV<currentVar.size(); ++idxV){
		  if(isNA(dat.weight(idxfrom+idxV))){
		    Rcout << "\t\t7-A\n";
		    sqrtW(idxV)=Type(1.0);
		    int a = dat.aux(idxfrom+idxV,2)-conf.minAge;
		    if(conf.predVarObsLink(f,a)>(-1)){
		      sqrtW(idxV) = sqrt(obs_fun::findLinkV(par.logSdLogObs(conf.keyVarObs(f,a))+(exp(par.predVarObs(conf.predVarObsLink(f,a))) -Type(1))*predObs(idxfrom+idxV),0)/currentVar(idxV));
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
		    Rcout << "\t\t7-B\n";
		    if(conf.fixVarToWeight(f)==1){
		      sqrtW(idxV)=sqrt(dat.weight(idxfrom+idxV)/currentVar(idxV));
		    }else{
		      sqrtW(idxV)=sqrt(Type(1)/dat.weight(idxfrom+idxV));
		    }
		  }
		}
		Rcout << "\t\t8\n";
		Rcout << dat.idxCor(f,y) << ", " << isNAINT(dat.idxCor(f,y)) << ", " << isNA((double)dat.idxCor(f,y)) << "\n";
		Rcout << R_NaInt << ", " << NA_INTEGER << ", " << INT_MIN << "\n";
		Rcout << (dat.idxCor(f,y) == R_NaInt) << ", " << (dat.idxCor(f,y) == NA_INTEGER) << ", " << (dat.idxCor(f,y) == INT_MIN) << "\n";
		Rcout << R_NaReal << ", " << (int)R_NaReal << "\n";
		if(isNAINT(dat.idxCor(f,y))){
		  Rcout << "\t\t9-A\n";
		  nll += nllVec(f)((dat.logobs.segment(idxfrom,idxlength)-predObs.segment(idxfrom,idxlength))/sqrtW,keep.segment(idxfrom,idxlength));
		  nll += (log(sqrtW)*keep.segment(idxfrom,idxlength)).sum();
		  SIMULATE_F(of){
		    if((conf.simFlag(2)==0 && y < dat.noYears) ||
		       (forecast.nYears > 0 &&
			forecast.forecastYear(y) > 0 &&
			forecast.simFlag(2)==0 &&
			y >= dat.noYears)){
		      dat.logobs.segment(idxfrom,idxlength) = predObs.segment(idxfrom,idxlength) + (nllVec(f).simulate()*sqrtW);
		    }
		  }
		}else{
		  Rcout << "\t\t9-B\n";
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
		    if((conf.simFlag(2)==0 && y < dat.noYears) ||
		       (forecast.nYears > 0 &&
			forecast.forecastYear(y) > 0 &&
			forecast.simFlag(2)==0 &&
			y >= dat.noYears)){
		      dat.logobs.segment(idxfrom,idxlength) = predObs.segment(idxfrom,idxlength) + thisnll.simulate()*sqrtW;
		    }
		  }
		}
		Rcout << "\t\t10\n";
		break;
	      case 1: // (ALN) Additive logistic-normal proportions + log-normal total numbers
		nll +=  nllVec(f)(obs_fun::addLogratio((vector<Type>)dat.logobs.segment(idxfrom,idxlength))-obs_fun::addLogratio((vector<Type>)predObs.segment(idxfrom,idxlength)));
		nll += log(obs_fun::log2proportion((vector<Type>)dat.logobs.segment(idxfrom,idxlength))).sum();
		nll -= dnorm(log(obs_fun::log2expsum((vector<Type>)dat.logobs.segment(idxfrom,idxlength))),
			     log(obs_fun::log2expsum((vector<Type>)predObs.segment(idxfrom,idxlength))),
			     exp(par.logSdLogTotalObs(totalParKey)),true);
		nll += log(obs_fun::log2expsum((vector<Type>)dat.logobs.segment(idxfrom,idxlength)));
		nll -= log(fabs(obs_fun::jacobianDet((vector<Type>)dat.logobs.segment(idxfrom,idxlength).exp())));
		nll -= dat.logobs.segment(idxfrom,idxlength).sum();
		SIMULATE_F(of){
		  if((conf.simFlag(2)==0 && y < dat.noYears) ||
		     (forecast.nYears > 0 &&
		      forecast.forecastYear(y) > 0 &&
		      forecast.simFlag(2)==0 &&
		      y >= dat.noYears)){
		    vector<Type> logProb(idxlength);
		    logProb.setZero();
		    logProb.segment(0,idxlength-1) = obs_fun::addLogratio(((vector<Type>)predObs.segment(idxfrom,idxlength))) + nllVec(f).simulate();
		    Type logDenom = obs_fun::logExpSum(logProb);
		    logProb -= logDenom;
		    Type logTotal = rnorm(log(obs_fun::log2expsum((vector<Type>)predObs.segment(idxfrom,idxlength))),
					  exp(par.logSdLogTotalObs(totalParKey)));
		    dat.logobs.segment(idxfrom,idxlength) = logProb + logTotal;
		  }
		}
		totalParKey++;
		break;
	      case 2: //Dirichlet
		Rf_error("Dirichlet distribution can only be used for fleet types 80, 90, 92");
	      default:
		Rf_error("Unknown obsLikelihoodFlag");
	      }
	    }
	  }else if(dat.fleetTypes(f)==5){
	    if(!isNAINT(dat.idx1(f,y))){    
	      for(int i=dat.idx1(f,y); i<=dat.idx2(f,y); ++i){
		//nll += -keep(i)*dnbinom(dat.logobs(i),predObs(i)*recapturePhiVec(i)/(Type(1.0)-recapturePhiVec(i)),recapturePhiVec(i),true);
		Type log_mu = log(predObs(i));
		Type log_var_minus_mu = log_mu - logitRecapturePhiVec(i);	      
		nll += -keep(i)*dnbinom_robust(dat.logobs(i),log_mu,log_var_minus_mu,true);
		SIMULATE_F(of){
		  if((conf.simFlag(2)==0 && y < dat.noYears) ||
		     (forecast.nYears > 0 &&
		      forecast.forecastYear(y) > 0 &&
		      forecast.simFlag(2)==0 &&
		      y >= dat.noYears)){
		    dat.logobs(i) = rnbinom(predObs(i)*recapturePhiVec(i)/(Type(1.0)-recapturePhiVec(i)),recapturePhiVec(i));
		  }
		}
	      }
	    }
	  }else if(dat.fleetTypes(f)==3||(dat.fleetTypes(f)==6)){
	    Type sd=0;
	    if(!isNAINT(dat.idx1(f,y))){
	      for(int i=dat.idx1(f,y); i<=dat.idx2(f,y); ++i){
		if(conf.keyBiomassTreat(f)==3){
		  sd = sqrt(varLogCatch(y));
		}else{
		  if(conf.keyBiomassTreat(f)==4){
		    sd = sqrt(varLogLand(y));
		  }else{
		    if(isNA(dat.weight(i))){
		      sd = exp(par.logSdLogObs(conf.keyVarObs(f,0)));
		    }else{
		      if(conf.fixVarToWeight(f)==1){
                        sd = sqrt(dat.weight(i));
		      }else{
                        sd = exp(par.logSdLogObs(conf.keyVarObs(f,0)))/sqrt(dat.weight(i));
		      }
		    }
		  }
		}
		nll += -keep(i)*dnorm(dat.logobs(i),predObs(i),sd,true);
		SIMULATE_F(of){
		  if((conf.simFlag(2)==0 && y < dat.noYears) ||
		     (forecast.nYears > 0 &&
		      forecast.forecastYear(y) > 0 &&
		      forecast.simFlag(2)==0 &&
		      y >= dat.noYears)){
		    dat.logobs(i) = rnorm(predObs(i),sd);
		  }
		}
	      }
	    }
	  }else if(dat.fleetTypes(f)==80){
	    if(!isNAINT(dat.idx1(f,y))){	    
	      int iMin = dat.idx1(f,y);
	      int nSeasons = CppAD::Integer(dat.auxData(iMin,5));
	      // Observation on additive logistic scale
	      vector<Type> log_X(nSeasons);
	      log_X.setConstant(R_NegInf);
	      // Prediction on log-proportion scale
	      vector<Type> log_P(nSeasons);
	      log_P.setConstant(R_NegInf);
	      data_indicator<vector<Type>, Type> K2 = keep.segment(dat.idx1(f,y),dat.idx2(f,y)-dat.idx1(f,y)+1);
	      Type xs = 1.0;
	      Type ps = 0.0;
	      for(int i=dat.idx1(f,y); i<=dat.idx2(f,y); ++i){
		log_X(CppAD::Integer(dat.auxData(i,4))-1) = dat.logobs(i);
		xs += exp(dat.logobs(i));
		log_P(CppAD::Integer(dat.auxData(i,4))-1) = predObs(i);
		ps += exp(predObs(i));
	      }
	      log_X(nSeasons-1) = 0.0;
	      log_P(nSeasons-1) = log(1.0 - squeeze(ps));
	      if(conf.obsLikelihoodFlag(f) == 1){ // Additive logistic normal
		// Keep free observations, already on correct scale
		vector<Type> logXuse = log_X.segment(0,nSeasons-1);
		// Additive logistic transformation of predicted proportions
		vector<Type> logPuse = log_P.segment(0,nSeasons-1) - log_P(nSeasons-1);
		// Quick fix for now
		for(int i = 0; i < logXuse.size(); ++i)
		  nll += nllVec(f)((vector<Type>)(logXuse - logPuse), K2);	       
	      }else if(conf.obsLikelihoodFlag(f) == 2){ // Dirichlet
		Type log_alpha = par.logSdLogObs(conf.keyVarObs(f,0));
		// Transform log_X to log proportions
		for(int i = 0; i < log_X.size(); ++i)
		  log_X(i) -= log(xs);
		// Type d = obs_fun::jacobianDetProportions((vector<Type>)log_X.exp());
		// nll -= ddirichlet(log_X,log_P,-log_alpha,K2,true) + log(fabs(d));
		nll -= ddirichlet(log_X,log_P,-log_alpha,K2,true);
	      }else{
		Rf_error("Fleet type 80 must use obsLikelihoodFlag ALN or Dirichlet");
	      }
	    }
	  }else if(dat.fleetTypes(f) == 90){
	    // Do nothing
	  }else if(dat.fleetTypes(f) == 92){
	    // Do nothing
	  }else{
	    Rf_error("Fleet type not implemented");
	  }   
	}   
      }
      Rcout << "F\n";
      SIMULATE_F(of) {
	REPORT_F(logF,of);
	REPORT_F(logN,of);
	vector<Type> logobs=dat.logobs; 
	REPORT_F(logobs,of);

	vector<Type> logEmpiricalSPR = empiricalSPR(dat, conf, logN, mort, true);
	REPORT_F(logEmpiricalSPR,of);
	vector<Type> logEmpiricalYPR = empiricalYPR(dat, conf, logN, mort, 0, true);
	REPORT_F(logEmpiricalYPR,of);
	vector<Type> logEmpiricalYPR_L = empiricalYPR(dat, conf, logN, mort, 1, true);
	REPORT_F(logEmpiricalYPR_L,of);
	vector<Type> logEmpiricalYPR_D = empiricalYPR(dat, conf, logN, mort, 2, true);
	REPORT_F(logEmpiricalYPR_D,of);
      }
      Rcout << "G\n";
      // REPORT_F(obsCov,of);
      REPORT_F(predObs,of);
      if(reportingLevel >= 0){
	ADREPORT_F(logssb,of);
	ADREPORT_F(logerb,of);
	ADREPORT_F(logrelativeerb,of);
	ADREPORT_F(logfbar,of);
	ADREPORT_F(logfbar_Effective,of);
	ADREPORT_F(logCatch,of);
	ADREPORT_F(logCatchByFleet,of);
	ADREPORT_F(logLand,of);
	ADREPORT_F(logtsb,of);
      
	REPORT_F(logCatchByFleetAge,of);

	REPORT_F(comps, of);
	ADREPORT_F(comps, of);
	REPORT_F(weekContrib, of);
      }
      Rcout << "H\n";

      // Additional forecast quantities
      if(forecast.nYears > 0){
	Rcout << "I\n";
	vector<Type> dis = disFun(dat, conf, logN, logF, mort);
	vector<Type> logDis = log(dis);
    
	vector<Type> landFbar = landFbarFun(dat, conf, logF);
	vector<Type> loglandfbar = log(landFbar);
    
	vector<Type> disFbar = disFbarFun(dat, conf, logF);
	vector<Type> logdisfbar = log(disFbar);

	matrix<Type> logFbarByFleet = fbarByFleet(conf, logF, true);
    
	ADREPORT_F(logDis,of);
	ADREPORT_F(loglandfbar,of);
	ADREPORT_F(logdisfbar,of);
	ADREPORT_F(logCatchAge, of);
	ADREPORT_F(logFbarByFleet, of);
    
	// SIMULATE_F(of) {
	//if(dat.forecast.simFlag[0] == 0 || dat.forecast.simFlag[1] == 0){
	REPORT_F(logssb,of);
	REPORT_F(logerb,of);
	REPORT_F(logrelativeerb,of);
	REPORT_F(logfbar,of);
	REPORT_F(logfbarL,of);
	REPORT_F(logCatch,of);
	REPORT_F(logCatchAge,of);
	REPORT_F(logCatchByFleet,of);
	REPORT_F(logCatchByFleetAge,of);
	REPORT_F(logFbarByFleet, of);
	REPORT_F(logLand,of);
	REPORT_F(logtsb,of);      
	//}
	// }
      }
      Rcout << "J\n";

      int timeSteps=logF.dim[1];
      if(reportingLevel >= 0){
	vector<Type> logLagR(logR.size());
	for(int i=0; i<logR.size(); ++i){
	  logLagR(i) = logN(1,i);
	}
	ADREPORT_F(logR,of);
	ADREPORT_F(logLagR,of);
  
	vector<Type> lastLogN = logN.col(timeSteps-1);
	ADREPORT_F(lastLogN,of);
	vector<Type> lastLogF = logF.col(timeSteps-1);
	ADREPORT_F(lastLogF,of);  

	vector<Type> beforeLastLogN = logN.col(timeSteps-2);
	ADREPORT_F(beforeLastLogN,of);
	vector<Type> beforeLastLogF = logF.col(timeSteps-2);
	ADREPORT_F(beforeLastLogF,of);
      }
      Rcout << "K\n";
      if(forecast.nYears > 0 && forecast.FModel(forecast.FModel.size()-1) == forecast.findMSY){
	Rcout << "L\n";

	int catchYears = std::min((int)asDouble(forecast.nYears),forecast.nCatchAverageYears);
	Type catchSum = sum((vector<Type>)cat.tail(catchYears)) / (Type)catchYears;
	nll -= par.implicitFunctionDelta * log(catchSum);    
	// Calculate Fmsy
	Type logFMSY = par.logFScaleMSY + logfbar(timeSteps - forecast.nYears - 1);
	ADREPORT_F(logFMSY, of);

	// Output stock status - Positive is good for the stock
	Type logFstatus = logFMSY - logfbar(timeSteps - forecast.nYears - 1);
	Type logSSBstatus = logssb(timeSteps - forecast.nYears - 1) - logssb(timeSteps - 1);
	ADREPORT_F(logFstatus, of);
	ADREPORT_F(logSSBstatus, of);    
      }
      Rcout << "Done with nllObs\n";
      return nll;
    }
    )

SAM_SPECIALIZATION(double nllObs(dataSet<double>&, confSet&, paraSet<double>&, forecastSet<double>&, array<double>&, array<double>&, array<double>&,  array<double>&, Recruitment<double>&, MortalitySet<double>&, data_indicator<vector<double>,double>&, int, objective_function<double>*));
SAM_SPECIALIZATION(TMBad::ad_aug nllObs(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, forecastSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, Recruitment<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>&, data_indicator<vector<TMBad::ad_aug>,TMBad::ad_aug>&, int, objective_function<TMBad::ad_aug>*));
