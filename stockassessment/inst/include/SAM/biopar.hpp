SAM_DEPENDS(convenience)
SAM_DEPENDS(define)
SAM_DEPENDS(forecast)

HEADER(
template <class Type>
struct bioResult {
  Type nll;

  matrix<Type> Sigma_GMRF;
  matrix<Type> mu;
  int nYears;
  int nAges;
 
  bioResult(Type nll_, matrix<Type> Sigma_GMRF_, matrix<Type> mu_, int nYears_, int nAges_) :
    nll(nll_), Sigma_GMRF(Sigma_GMRF_), mu(mu_), nYears(nYears_), nAges(nAges_) {};

  bioResult() :
    nll(), Sigma_GMRF(), mu(), nYears(), nAges() {};
  

 matrix<Type> getSigmaMarginal(int y){
    if(y < 0 || y >= mu.rows())
      Rf_error("Year is out of bounds in biopar getSigmaSubsets");
    matrix<Type> Sigma_11(nAges,nAges); // (Ages year y) x (Ages year y)
    Sigma_11.setZero();
    
    for(int j1 = 0; j1 < nAges; ++j1){
      for(int j2 = 0; j2 < nAges; ++j2){
	Sigma_11(j1,j2) = Sigma_GMRF(y + j1 * nYears,y + j2 * nYears);
      }
    }   
    return Sigma_11;
  }
  
  vector<matrix<Type> > getSigmaSubsets(int y){
    if(y <= 0 || y >= mu.rows())
      Rf_error("Year is out of bounds in biopar getSigmaSubsets");
    matrix<Type> Sigma_11(nAges,nAges); // (Ages year y) x (Ages year y)
    Sigma_11.setZero();
    matrix<Type> Sigma_12(nAges,nAges); // (Ages year y) x (Ages year y-1) 
    Sigma_12.setZero();
    matrix<Type> Sigma_21(nAges,nAges); // (Ages year y-1) x (Ages year y) 
    Sigma_21.setZero();
    matrix<Type> Sigma_22(nAges,nAges); // (Ages year y-1) x (Ages year y-1) 
    Sigma_22.setZero();
    
    for(int j1 = 0; j1 < nAges; ++j1){
      for(int j2 = 0; j2 < nAges; ++j2){
	Sigma_11(j1,j2) = Sigma_GMRF(y + j1 * nYears,y + j2 * nYears);
	Sigma_12(j1,j2) = Sigma_GMRF(y + j1 * nYears,(y-1) + j2 * nYears);
	Sigma_21(j1,j2) = Sigma_GMRF((y-1) + j1 * nYears,y + j2 * nYears);
	Sigma_22(j1,j2) = Sigma_GMRF((y-1) + j1 * nYears,(y-1) + j2 * nYears);
      }
    }
    vector<matrix<Type> > res(4);
    res(0) = Sigma_11;
    res(1) = Sigma_12;
    res(2) = Sigma_21;
    res(3) = Sigma_22;
    return res;
  }

  vector<Type> simulate(vector<Type>& x, int y){
    if(y < 0 || y > mu.rows())
      Rf_error("Year is out of bounds in biopar simulation");
    if(y == 0){
      matrix<Type> Sigma = getSigmaMarginal(y);
      density::MVNORM_t<Type> dens(Sigma);
      vector<Type> mu1 = mu.row(y);
      return mu1 + dens.simulate();
    }
    vector<matrix<Type> > Sigmas = getSigmaSubsets(y);
    matrix<Type> Prec_22 = atomic::matinv(Sigmas(3));
    matrix<Type> MeanAdj = Sigmas(1) * Prec_22;
    matrix<Type> Cov_1given2 = Sigmas(0) - Sigmas(1) * Prec_22 * Sigmas(2);
    density::MVNORM_t<Type> dens(Cov_1given2);
    vector<Type> mu1 = mu.row(y);
    vector<Type> mu2 = mu.row(y-1);
    vector<Type> tmp = x - mu2;
    vector<Type> pred = mu1 + MeanAdj * tmp;
    return pred + dens.simulate();
  }
  

}
       )

template <class Type>
bioResult<Type> nllBioProcess(array<Type> P, vector<Type> meanVec, vector<int> keyMeanVec, vector<Type> logPhi, Type logSdP)SOURCE({
    int nrow=P.dim[0];		// Years
    int ncol=P.dim[1];		// Ages
    int n=nrow*ncol;
    vector<int> r(n); r.setZero();   
    vector<int> c(n); c.setZero();
    int idx=0;
    for(int j=0; j<ncol; ++j){
      for(int i=0; i<nrow; ++i){
        r(idx)=i;
	c(idx)=j;
	++idx;
      }
    }
    matrix<Type> Wc(n,n); Wc.setZero();
    matrix<Type> Wd(n,n); Wd.setZero();
    matrix<Type> Wp(n,n); Wp.setZero();    
    for(int i=0; i<n; ++i){
      for(int j=0; j<n; ++j){
	if((c(i)==c(j))&&(abs(r(i)-r(j))==1)){
	  Wc(i,j)=1;
	  Wc(i,i)-=1;
        }   
 	if( (((r(i)-r(j))==1)&&((c(i)-c(j))==1))||(((r(i)-r(j)==(-1)))&&((c(i)-c(j))==(-1))) ){
       	  Wd(i,j)=1;
	  Wd(i,i)-=1;
	}
	if(logPhi.size()==3){
          if ((c(i)==(ncol-1)) && (c(j)==(ncol-1)) && (abs(r(i)-r(j))==1) ){
            Wp(i,j)=1;
	    Wp(i,i)-=1;
	  }
	}

      }
    }

    array<Type> mP(nrow,ncol);
    for(int i=0; i<nrow; ++i){
      for(int j=0; j<ncol; ++j){
	mP(i,j)=meanVec(keyMeanVec(j));
      }
    }

    vector<Type> phi=exp(logPhi);
    matrix<Type> I(Wc.rows(),Wc.cols());
    I.setIdentity();
    matrix<Type> Q=I-phi(0)*Wc-phi(1)*Wd;
    if(logPhi.size()==3){
      Q-=phi(2)*Wp;
    }

    matrix<Type> Sigma_GMRF = atomic::matinv(Q) * exp(2.0*logSdP);

    Type nll = density::SCALE(density::GMRF(asSparseMatrix(Q)),exp(logSdP))((P-mP).vec());

    bioResult<Type> res(nll, Sigma_GMRF, mP.matrix(), P.dim[0], P.dim[1]);

    return res;
  });

SAM_SPECIALIZATION(bioResult<double> nllBioProcess(array<double>, vector<double>, vector<int>, vector<double>, double));
SAM_SPECIALIZATION(bioResult<TMBad::ad_aug> nllBioProcess(array<TMBad::ad_aug>, vector<TMBad::ad_aug>, vector<int>, vector<TMBad::ad_aug>, TMBad::ad_aug));


template <class Type>
Type nllSW(array<Type> &logSW, dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, forecastSet<Type>& forecast, objective_function<Type> *of) SOURCE({
    if(conf.stockWeightModel>=1){
      Type nll=0;
      array<Type> sw=dat.stockMeanWeight;
      bioResult<Type> br = nllBioProcess(logSW, par.meanLogSW, conf.keyStockWeightMean, par.logPhiSW, par.logSdProcLogSW(0));
      nll += br.nll;
      for(int i=0; i<sw.dim[0]; ++i){
	for(int j=0; j<sw.dim[1]; ++j){
	  if(!isNA(sw(i,j))){
	    nll += -dnorm(log(sw(i,j)),logSW(i,j),exp(par.logSdLogSW(conf.keyStockWeightObsVar(j))),true);	   
	  }	  
	  dat.stockMeanWeight(i,j)=exp(logSW(i,j));
	}
	SIMULATE_F(of){
	  if((forecast.nYears > 0 && forecast.simFlag(3) == 0 && ((!forecast.useModelLastN && forecast.forecastYear(i) >= 1) || forecast.forecastYear(i) > 1)) || (conf.simFlag(3)==0 && i > 0)){
	    vector<Type> v = logSW.matrix().row(i-1);
	    vector<Type> p = br.simulate(v,i);
	    for(int j=0; j<sw.dim[1]; ++j){
	      logSW(i,j) = p(j);
	      dat.stockMeanWeight(i,j)=exp(logSW(i,j));
	    }
	  }else if((forecast.nYears > 0 && forecast.simFlag(3) == 1 && ((!forecast.useModelLastN && forecast.forecastYear(i) >= 1) || forecast.forecastYear(i) > 1))){
	    for(int j=0; j<sw.dim[1]; ++j){
	      logSW(i,j) = log(sw(i,j));
	      dat.stockMeanWeight(i,j) = sw(i,j);
	    }
	  }
	}
      }
      array<Type> stockMeanWeight = dat.stockMeanWeight;
      REPORT_F(stockMeanWeight,of);
      ADREPORT_F(logSW,of);	// Needed for R based forecast
      int timeSteps = dat.years.size();
      vector<Type> lastLogSW = logSW.matrix().row(timeSteps-1);
      ADREPORT_F(lastLogSW,of);
      vector<Type> beforeLastLogSW = logSW.matrix().row(timeSteps-2);
      ADREPORT_F(beforeLastLogSW,of);
      return nll;
    }
    return Type(0);
  } );

SAM_SPECIALIZATION(double nllSW(array<double>&, dataSet<double>&, confSet&, paraSet<double>&, forecastSet<double>&, objective_function<double>*));
SAM_SPECIALIZATION(TMBad::ad_aug nllSW(array<TMBad::ad_aug>&, dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, forecastSet<TMBad::ad_aug>&, objective_function<TMBad::ad_aug>*));

template <class Type>
Type nllCW(array<Type> &logCW, dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, forecastSet<Type>& forecast, objective_function<Type> *of) SOURCE( {
  if(conf.catchWeightModel>=1){

    Type nll=0;
    array<Type> cw=dat.catchMeanWeight;
    for(int k = 0; k < logCW.dim[2]; ++k){
      array<Type> lcw = logCW.col(k);
      //nll += nllBioProcess(lcw, par.meanLogCW, (vector<int>)conf.keyCatchWeightMean.row(k), (vector<Type>)par.logPhiCW.col(k), (Type)par.logSdProcLogCW(k));
      bioResult<Type> br = nllBioProcess(lcw, par.meanLogCW, (vector<int>)conf.keyCatchWeightMean.row(k), (vector<Type>)par.logPhiCW.col(k), (Type)par.logSdProcLogCW(k));
      nll += br.nll;
    // }
    // for(int k=0; k<cw.dim[2]; ++k){
      for(int i=0; i<cw.dim[0]; ++i){
	for(int j=0; j<cw.dim[1]; ++j){
	  if(!isNA(cw(i,j,k))){
	    nll += -dnorm(log(cw(i,j,k)),logCW(i,j,k),exp(par.logSdLogCW(conf.keyCatchWeightObsVar(k,j))),true);
	  }
	  dat.catchMeanWeight(i,j,k)=exp(logCW(i,j,k));
	}
	SIMULATE_F(of){
	  if((forecast.nYears > 0 && forecast.simFlag(3) == 0 && ((!forecast.useModelLastN && forecast.forecastYear(i) >= 1) || forecast.forecastYear(i) > 1)) || (conf.simFlag(3)==0 && i > 0)){
	    vector<Type> v = logCW.col(k).matrix().row(i-1);
	    vector<Type> p = br.simulate(v,i);
	    for(int j=0; j<cw.dim[1]; ++j){
	      logCW(i,j,k) = p(j);
	      dat.catchMeanWeight(i,j,k)=exp(logCW(i,j,k));
	    }
	  }else if((forecast.nYears > 0 && forecast.simFlag(3) == 1 && ((!forecast.useModelLastN && forecast.forecastYear(i) >= 1) || forecast.forecastYear(i) > 1))){
	    for(int j=0; j<cw.dim[1]; ++j){
	      logCW(i,j,k) = log(cw(i,j,k));
	      dat.catchMeanWeight(i,j,k) = cw(i,j,k);
	    }
	  }
	}
      }
    }
    array<Type> catchMeanWeight = dat.catchMeanWeight;
    REPORT_F(catchMeanWeight,of);
    ADREPORT_F(logCW,of);	// Needed for R based forecast
    int timeSteps = dat.years.size();
    vector<Type> lastLogCW = logCW.matrix().row(timeSteps-1);
    ADREPORT_F(lastLogCW,of);
    vector<Type> beforeLastLogCW = logCW.matrix().row(timeSteps-2);
    ADREPORT_F(beforeLastLogCW,of);
    return nll;
  }
  return Type(0);
  } );


SAM_SPECIALIZATION(double nllCW(array<double>&, dataSet<double>&, confSet&, paraSet<double>&, forecastSet<double>&, objective_function<double>*));
SAM_SPECIALIZATION(TMBad::ad_aug nllCW(array<TMBad::ad_aug>&, dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, forecastSet<TMBad::ad_aug>&, objective_function<TMBad::ad_aug>*));


template <class Type>
  Type nllMO(array<Type> &logitMO, dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, forecastSet<Type>& forecast, objective_function<Type> *of) SOURCE({
  if(conf.matureModel>=1){
    Type nll=0;
    array<Type> mo=dat.propMat;
    //nll += nllBioProcess(logitMO, par.meanLogitMO, conf.keyMatureMean, par.logPhiMO, par.logSdProcLogitMO(0));
    bioResult<Type> br = nllBioProcess(logitMO, par.meanLogitMO, conf.keyMatureMean, par.logPhiMO, par.logSdProcLogitMO(0));
    nll += br.nll;
    Type m,a,b, prec;
    for(int i=0; i<mo.dim[0]; ++i){
      for(int j=0; j<mo.dim[1]; ++j){
	m = invlogit(logitMO(i,j));
        if(!isNA(mo(i,j))){

 	  //Francisco Cribari-Neto, Achim Zeileis (2010). "Beta Regression in R."
	  prec=invlogit(par.logSdMO(conf.keyMatureObsVar(j)))*1.0e3;
	  a = m*prec;                 //
          b = (Type(1)-m)*prec;       //v=mu*(1-mu)/(1+precision)
	  //v=m*(1-m)/(1+exp(logPrecision));
	
	  if(mo(i,j)>0.5){
            nll += -dbeta(Type(squash(1-mo(i,j))),b,a,true);
	  }else{
	    nll += -dbeta(Type(squash(mo(i,j))),a,b,true);
	  }
        }
	dat.propMat(i,j)=m;
      }
      SIMULATE_F(of){
	if((forecast.nYears > 0 && forecast.simFlag(3) == 0 && ((!forecast.useModelLastN && forecast.forecastYear(i) >= 1) || forecast.forecastYear(i) > 1)) || (conf.simFlag(3)==0 && i > 0)){
	  vector<Type> v = logitMO.matrix().row(i-1);
	  vector<Type> p = br.simulate(v,i);
	  for(int j=0; j<mo.dim[1]; ++j){
	    logitMO(i,j) = p(j);
	    dat.propMat(i,j)=invlogit(logitMO(i,j));
	  }
	}else if((forecast.nYears > 0 && forecast.simFlag(3) == 1 && ((!forecast.useModelLastN && forecast.forecastYear(i) >= 1) || forecast.forecastYear(i) > 1))){
	    for(int j=0; j<mo.dim[1]; ++j){
	      logitMO(i,j) = logit(mo(i,j));
		dat.propMat(i,j) = invlogit(logitMO(i,j));
	    }
	}
      }
    }
    array<Type> propMat = dat.propMat;
    REPORT_F(propMat, of);
    ADREPORT_F(logitMO,of);	// Needed for R based forecast
    int timeSteps = dat.years.size();
    vector<Type> lastLogitMO = logitMO.matrix().row(timeSteps-1);
    ADREPORT_F(lastLogitMO,of);
    vector<Type> beforeLastLogitMO = logitMO.matrix().row(timeSteps-2);
    ADREPORT_F(beforeLastLogitMO,of);
    return nll;
  }
  return Type(0);
    } )


SAM_SPECIALIZATION(double nllMO(array<double>&, dataSet<double>&, confSet&, paraSet<double>&, forecastSet<double>&, objective_function<double>*));
SAM_SPECIALIZATION(TMBad::ad_aug nllMO(array<TMBad::ad_aug>&, dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, forecastSet<TMBad::ad_aug>&, objective_function<TMBad::ad_aug>*));


template <class Type>
Type nllNM(array<Type> &logNM, dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, forecastSet<Type>& forecast, objective_function<Type> *of) SOURCE({
    Type nll=0;
    array<Type> nm=dat.natMor;
    if(conf.mortalityModel>=1 || conf.mortalityModelMeanStructure >= 0){
      //array<Type> nm=dat.natMor;
      //nll += nllBioProcess(logNM, par.meanLogNM, conf.keyMortalityMean, par.logPhiNM, par.logSdProcLogNM(0));
      vector<int> fakeKeyMM(conf.keyMortalityMean.size());
      for(int i = 0; i < fakeKeyMM.size(); ++i)
	fakeKeyMM(i) = i;
      vector<Type> mlnm(conf.keyMortalityMean.size());
      mlnm.setZero();
      /*
	conf.mortalityModelMeanStructure
	0: Free parameters by age, combined with keyMortalityMean
	1: Decreasing parameters by age, combined with keyMortalityMean
	2: Weibull
	3: Lorenzen allometric and length-inverse
       */
      if(conf.mortalityModelMeanStructure == 0){
	for(int i = 0; i < fakeKeyMM.size(); ++i)
	  mlnm(i) = par.meanLogNM(conf.keyMortalityMean(i));
      }else if(conf.mortalityModelMeanStructure == 1){
	vector<Type> mmX(par.meanLogNM.size());
	mmX(par.meanLogNM.size()-1) = par.meanLogNM(par.meanLogNM.size()-1);
	for(int qq = par.meanLogNM.size()-2; qq >= 0; qq--)
	  mmX(qq) = mmX(qq+1) + exp(par.meanLogNM(qq));
	for(int i = 0; i < fakeKeyMM.size(); ++i)
	  mlnm(i) = mmX(conf.keyMortalityMean(i));
      }else if(conf.mortalityModelMeanStructure == 2){
	// Weibull (2 par)
	// Hazard: p * lambda ^ p * a ^(p-1)
	// Cumulative hazard: (lambda * a) ^ p
	// Cumulative / average hazard between age a and a+1: ((a+1)^p-a^p) * lambda^p
	Type p = 1.0 / (1.0 + exp(par.meanLogNM(0)));
	Type logLambda = par.meanLogNM(1);
	for(int i = 0; i < mlnm.size(); ++i){
	  Type a = conf.minAge + i;
	  mlnm(i) = logspace_sub_SAM(p*log(a+1),p*log(a+Type(1e-8))) + p * logLambda;
	}
      }else if(conf.mortalityModelMeanStructure == 3){
	// Lorenzen allometric and length-inverse (4 par)
	Type logMLinf = par.meanLogNM(0);
	Type K = exp(par.meanLogNM(1));
	Type a0 = -exp(par.meanLogNM(2));
	Type c = -exp(par.meanLogNM(3));
	for(int i = 0; i < mlnm.size(); ++i){
	  Type a = conf.minAge + i + 0.5;
	  mlnm(i) = logMLinf + c * logspace_sub_SAM(Type(0.0), -K * (a-a0));
	}
      }else if(conf.mortalityModelMeanStructure == 4){
	// Charnov (2 par)
	Type logK = par.meanLogNM(0);
	Type K = exp(par.meanLogNM(0));
	Type a0 = -exp(par.meanLogNM(1));
	for(int i = 0; i < mlnm.size(); ++i){
	  Type a = conf.minAge + i + 0.5;
	  mlnm(i) = logK + (-1.5) * logspace_sub_SAM(Type(0.0), -K * (a-a0));
	}
      }else if(conf.mortalityModelMeanStructure == 5){
	// Caddy (2 par)
	Type logMtau = par.meanLogNM(0);
	Type logB = par.meanLogNM(1);
	for(int i = 0; i < mlnm.size(); ++i){
	  Type a = conf.minAge + i + 0.5;
	  mlnm(i) = logspace_add_SAM(logMtau, logB - log(a+Type(1e-8)));
	}
      }else if(conf.mortalityModelMeanStructure == 6){
	// Siler (5 par)
	Type logMc = par.meanLogNM(0);
	Type logMd = par.meanLogNM(1);
	Type logMs = par.meanLogNM(2);
	Type bd = exp(par.meanLogNM(3));
	Type bs = exp(par.meanLogNM(4));
	for(int i = 0; i < mlnm.size(); ++i){
	  Type a = conf.minAge + i + 0.5;
	  mlnm(i) = logspace_add_SAM(logMc, logspace_add_SAM(logMd - bd * a, logMs - bs * a));
	}
      }else{
	Rf_error("Wrong mortalityModelMeanStructure");
      }
      REPORT_F(mlnm,of);
      bioResult<Type> br;
      if(conf.mortalityModel >= 1){
        br = nllBioProcess(logNM, mlnm, fakeKeyMM, par.logPhiNM, par.logSdProcLogNM(0));
	// bioResult<Type> br = nllBioProcess(logNM, par.meanLogNM, conf.keyMortalityMean, par.logPhiNM, par.logSdProcLogNM(0));
	nll += br.nll;
      }
      array<Type> lX(dat.natMor.dim[0],dat.natMor.dim[1]);
      lX.setZero();
      if(dat.Mcovariate.dim[1] > 0){
	// Proportional hazard covariates
	Type dt = dat.Mcovariate.dim[1];
	for(int i = 0; i < lX.dim[0]; ++i){ // Year
	  for(int j = 0; j < lX.dim[1]; ++j){ // Age
	    for(int k = 0; k < dat.Mcovariate.dim[2]; ++k){ // Covariate
	      for(int t = 0; t < dat.Mcovariate.dim[1]; ++t){ // Covariate time point
		if(!isNAINT(conf.keyMortalityCovariate(j)))
		  lX(i,j) += par.Mbeta(conf.keyMortalityCovariate(j),k) * dat.Mcovariate(i,t,k) / dt;	
	      }
	    }
	  }
	}
      }
      for(int i=0; i<nm.dim[0]; ++i){
	for(int j=0; j<nm.dim[1]; ++j){
	  Type eLNM = 0.0;
	   if(conf.mortalityModel >= 1){
	     eLNM = logNM(i,j);
	   }else{
	     eLNM = mlnm(j);
	   }
	  if(!isNA(nm(i,j))){
	    nll += -dnorm(log(nm(i,j)),eLNM + lX(i,j),exp(par.logSdLogNM(conf.keyMortalityObsVar(j))),true);
	  }
	  dat.natMor(i,j)=exp(eLNM + lX(i,j));
	}
	SIMULATE_F(of){
	  if(conf.mortalityModel >= 1){
	    if((forecast.nYears > 0 && forecast.simFlag(3) == 0 && ((!forecast.useModelLastN && forecast.forecastYear(i) >= 1) || forecast.forecastYear(i) > 1)) || (conf.simFlag(3)==0 && i > 0)){
	      vector<Type> v = logNM.matrix().row(i-1);
	      vector<Type> p = br.simulate(v,i);
	
	      for(int j=0; j<nm.dim[1]; ++j){
		logNM(i,j) = p(j);
		dat.natMor(i,j)=exp(logNM(i,j));
	      }
	    }else if((forecast.nYears > 0 && forecast.simFlag(3) == 1 && ((!forecast.useModelLastN && forecast.forecastYear(i) >= 1) || forecast.forecastYear(i) > 1))){
	      for(int j=0; j<nm.dim[1]; ++j){
		logNM(i,j) = log(nm(i,j));
		dat.natMor(i,j) = nm(i,j);
	      }
	    }
	  }
	}
      }
      array<Type> natMor = dat.natMor;
      REPORT_F(natMor, of);
      if(conf.mortalityModel >= 1){
	ADREPORT_F(logNM,of);	// Needed for R based forecast
      }
      int timeSteps = dat.years.size();
      if(conf.mortalityModel >= 1){
	vector<Type> lastLogNM = logNM.matrix().row(timeSteps-1);
	ADREPORT_F(lastLogNM,of);
	vector<Type> beforeLastLogNM = logNM.matrix().row(timeSteps-2);
	ADREPORT_F(beforeLastLogNM,of);
      }
      //return nll;
    }
    if(conf.keyScaleMModel > 0){
      vector<Type> m = dat.natMor.matrix().colwise().mean();
      for(int i=0; i<dat.natMor.dim[0]; ++i){
	for(int j=0; j<dat.natMor.dim[1]; ++j){
	  if(conf.keyScaleMModel == 1){ // Scale all M with one number
	    dat.natMor(i,j) = exp(log(dat.natMor(i,j)) + par.scaleMpars(0));
	  }else if(conf.keyScaleMModel == 2){ // Scale each age
	    dat.natMor(i,j) = exp(log(dat.natMor(i,j)) + par.scaleMpars(j));
	  }else if(conf.keyScaleMModel == 3){ // Double Gompertz-Makeham	    
	    Type lambda = exp(par.scaleMpars(0));
	    Type a1 = exp(par.scaleMpars(1));
	    Type a2 = exp(par.scaleMpars(2));
	    Type b1 = exp(par.scaleMpars(3));
	    Type b2 = exp(par.scaleMpars(4));
	    Type ag = (Type)j / (Type)dat.natMor.dim[1];
	    if(conf.mortalityModel>=1){
	      dat.natMor(i,j) = dat.natMor(i,j)/nm(i,j) * (a1*exp(-b1*ag) + a2*exp(b2*ag) + lambda);
	    }else{
	      dat.natMor(i,j) = dat.natMor(i,j)/m(j) * (a1*exp(-b1*ag) + a2*exp(b2*ag) + lambda);
	    }
	  // }else if(conf.keyScaleMModel == 4){ // Double Gompertz-Makeham + covariates
	  //   Type lambda = exp(par.scaleMpars(0));
	  //   Type a1 = exp(par.scaleMpars(1));
	  //   Type a2 = exp(par.scaleMpars(2));
	  //   Type b1 = exp(par.scaleMpars(3));
	  //   Type b2 = exp(par.scaleMpars(4));
	  //   Type ag = (Type)j / (Type)dat.natMor.dim[1];
	  //   Type hm = m(j);
	  //   if(conf.mortalityModel>=1){
	  //     hm = nm(i,j);
	  //   }
	  //   Type h0 = (a1*exp(-b1*ag) + a2*exp(b2*ag) + lambda) * dat.natMor(i,j)/hm;
	  //   Type logPH = 0.0;
	  //   for(int k = 0; k < dat.Mcovariate.dim[1]; ++k)
	  //     logPH += par.scaleMpars(5+k) * dat.Mcovariate(i,j);
	  //   dat.natMor(i,j) =  h0 * exp(logPH);
	  }
	}
      }
      array<Type> natMor = dat.natMor;
      REPORT_F(natMor, of);
    }

    //return Type(0);
    return nll;
    })


SAM_SPECIALIZATION(double nllNM(array<double>&, dataSet<double>&, confSet&, paraSet<double>&, forecastSet<double>&, objective_function<double>*));
SAM_SPECIALIZATION(TMBad::ad_aug nllNM(array<TMBad::ad_aug>&, dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, forecastSet<TMBad::ad_aug>&, objective_function<TMBad::ad_aug>*));
