SAM_DEPENDS(convenience)
SAM_DEPENDS(define)
SAM_DEPENDS(forecast)

HEADER(
template <class Type>
struct bioResult {
  Type nll;
  matrix<Type> Sigma_RR;
  matrix<Type> MeanAdjScale;
  vector<Type> mu;

  density::MVNORM_t<Type> dens;

  bioResult(Type nll_, matrix<Type> Sigma_RR_, matrix<Type> MeanAdjScale_, vector<Type> mu_) :
    nll(nll_), Sigma_RR(Sigma_RR_), MeanAdjScale(MeanAdjScale_), mu(mu_), dens(Sigma_RR_) {};
  
  vector<Type> pred(vector<Type>& x){
    return mu + MeanAdjScale * (vector<Type>)(x - mu);
  }
  vector<Type> simulate(vector<Type>& x){
    return pred(x) + dens.simulate();
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
    
    // Conditional distribution of last time given one before
    matrix<Type> Qlast(2*ncol,2*ncol); Qlast.setZero();
    idx = 0;
    for(int i = 0; i < n; ++i)
      if(r(i) > 0 && r(i) < 3)
	for(int j = 0; j < n; ++j) // Find ages
	  if(r(j) > 0 && r(j) < 3){	// last two years
	    int ki = c(i) + ncol * (2 - r(i));
	    int kj = c(j) + ncol * (2 - r(j));
	    Qlast(ki,kj) = Q(i,j);
	  }
    Qlast *= exp(-2.0 * logSdP);

    // matrix<Type> SigmaLast = atomic::matinv(Qlast);
    // // Conditioning  on first column
    // matrix<Type> SigmaLast_Ran_Cond = SigmaLast.block(ncol,0,ncol,ncol);
    // matrix<Type> SigmaLast_Ran_Ran = SigmaLast.block(ncol,ncol,ncol,ncol);
    // matrix<Type> Q_Cond_Cond = Qlast.block(0,0,ncol,ncol);
    // matrix<Type> Sigma_Ran = SigmaLast_Ran_Ran - SigmaLast_Ran_Cond * Q_Cond_Cond * (matrix<Type>)SigmaLast_Ran_Cond.transpose();
    matrix<Type> Q_RR = Qlast.block(ncol,ncol,ncol,ncol);
    matrix<Type> Q_CR = Qlast.block(ncol,0,ncol,ncol);
    matrix<Type> Sigma_RR = atomic::matinv(Q_RR);
    matrix<Type> MeanAdjScale = Sigma_RR * Q_CR;
    vector<Type> mu(mP.dim(1));
    for(int j=0; j<ncol; ++j){
	mu(j)=meanVec(keyMeanVec(j));
    }
    
    Type nll = density::SCALE(density::GMRF(asSparseMatrix(Q)),exp(logSdP))((P-mP).vec());

    bioResult<Type> res(nll, Sigma_RR, MeanAdjScale, mu);

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
	  if((forecast.nYears > 0 && forecast.forecastYear(i) > 0) || (conf.simFlag(1)==0 && i > 0)){
	    vector<Type> v = logSW.matrix().row(i-1);
	    vector<Type> p = br.simulate(v);
	    for(int j=0; j<sw.dim[1]; ++j){
	      logSW(i,j) = p(j);
	      dat.stockMeanWeight(i,j)=exp(logSW(i,j));
	    }
	  }
	}
      }
      array<Type> stockMeanWeight = dat.stockMeanWeight;
      REPORT_F(stockMeanWeight,of);
      ADREPORT_F(logSW,of);	// Needed for R based forecast
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
	  if((forecast.nYears > 0 && forecast.forecastYear(i) > 0) || (conf.simFlag(1)==0 && i > 0)){
	    vector<Type> v = logCW.col(k).matrix().row(i-1);
	    vector<Type> p = br.simulate(v);
	    for(int j=0; j<cw.dim[1]; ++j){
	      logCW(i,j,k) = p(j);
	      dat.catchMeanWeight(i,j,k)=exp(logCW(i,j,k));
	    }
	  }
	}
      }
    }
    array<Type> catchMeanWeight = dat.catchMeanWeight;
    REPORT_F(catchMeanWeight,of);
    ADREPORT_F(logCW,of);	// Needed for R based forecast
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
	  prec=invlogit(par.logSdMO(0))*1.0e3;
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
	if((forecast.nYears > 0 && forecast.forecastYear(i) > 0) || (conf.simFlag(1)==0 && i > 0)){
	  vector<Type> v = logitMO.matrix().row(i-1);
	  vector<Type> p = br.simulate(v);
	  for(int j=0; j<mo.dim[1]; ++j){
	    logitMO(i,j) = p(j);
	    dat.propMat(i,j)=invlogit(logitMO(i,j));
	  }
	}
      }
    }
    array<Type> propMat = dat.propMat;
    REPORT_F(propMat, of);
    ADREPORT_F(logitMO,of);	// Needed for R based forecast
    return nll;
  }
  return Type(0);
    } )


SAM_SPECIALIZATION(double nllMO(array<double>&, dataSet<double>&, confSet&, paraSet<double>&, forecastSet<double>&, objective_function<double>*));
SAM_SPECIALIZATION(TMBad::ad_aug nllMO(array<TMBad::ad_aug>&, dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, forecastSet<TMBad::ad_aug>&, objective_function<TMBad::ad_aug>*));


template <class Type>
Type nllNM(array<Type> &logNM, dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, forecastSet<Type>& forecast, objective_function<Type> *of) SOURCE({
  if(conf.mortalityModel>=1){
    Type nll=0;
    array<Type> nm=dat.natMor;
    //nll += nllBioProcess(logNM, par.meanLogNM, conf.keyMortalityMean, par.logPhiNM, par.logSdProcLogNM(0));
    bioResult<Type> br = nllBioProcess(logNM, par.meanLogNM, conf.keyMortalityMean, par.logPhiNM, par.logSdProcLogNM(0));
    nll += br.nll;
    for(int i=0; i<nm.dim[0]; ++i){
      for(int j=0; j<nm.dim[1]; ++j){
        if(!isNA(nm(i,j))){
          nll += -dnorm(log(nm(i,j)),logNM(i,j),exp(par.logSdLogNM(conf.keyMortalityObsVar(j))),true);
        }
	dat.natMor(i,j)=exp(logNM(i,j));
      }
      SIMULATE_F(of){
	if((forecast.nYears > 0 && forecast.forecastYear(i) > 0) || (conf.simFlag(1)==0 && i > 0)){
	  vector<Type> v = logNM.matrix().row(i-1);
	  vector<Type> p = br.simulate(v);
	  for(int j=0; j<nm.dim[1]; ++j){
	    logNM(i,j) = p(j);
	    dat.natMor(i,j)=exp(logNM(i,j));
	  }
	}
      }
    }
    array<Type> natMor = dat.natMor;
    REPORT_F(natMor, of);
    ADREPORT_F(logNM,of);	// Needed for R based forecast
    return nll;
  }
  return Type(0);
  })


SAM_SPECIALIZATION(double nllNM(array<double>&, dataSet<double>&, confSet&, paraSet<double>&, forecastSet<double>&, objective_function<double>*));
SAM_SPECIALIZATION(TMBad::ad_aug nllNM(array<TMBad::ad_aug>&, dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, forecastSet<TMBad::ad_aug>&, objective_function<TMBad::ad_aug>*));
