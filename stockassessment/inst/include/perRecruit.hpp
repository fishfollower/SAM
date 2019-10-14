
#define SAM_NegInf -100.0
#define SAM_Zero exp(-20.0)
// exp(SAM_NegInf)

template<class Type>
Type dFunctionalSR(Type ssb, vector<Type> rp, int srmc){
  vector<AD<Type> > rp2(rp.size());
  rp2 = rp.template cast<AD<Type> >();
  CppAD::vector<AD<Type> > x( 1 );
  x[0] = ssb;
  CppAD::Independent(x);
  CppAD::vector<AD<Type> > y( 1 );
  y[0] = exp(functionalStockRecruitment(x[0], rp2, srmc));
  CppAD::ADFun<Type> F(x, y);
  CppAD::vector<Type> x_eval( 1 );
  x_eval[0] = ssb;
  vector<Type> r = F.Jacobian(x_eval);
  return r[0];
}


template<class Type>
struct PERREC_t {
  Type logFbar;
  Type logYPR;
  Type logSPR;
  Type logSe;
  Type logRe;
  Type logYe;  
};


template<class Type, class T>
PERREC_t<T> perRecruit(T Fbar, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, vector<Type>& sel, vector<int> aveYears, int nYears = 300){
  
  // Prepare data
  dataSet<T> newDat = dat.template cast<T>();
  int nMYears = dat.noYears;
  // propMat
  extendArray(newDat.propMat, nMYears, nYears, aveYears, false);
   // stockMeanWeight
  extendArray(newDat.stockMeanWeight, nMYears, nYears, aveYears, false);
  // catchMeanWeight
  extendArray(newDat.catchMeanWeight, nMYears, nYears, aveYears, false);
  // natMor
  extendArray(newDat.natMor, nMYears, nYears, aveYears, false);
  // landFrac
  extendArray(newDat.landFrac, nMYears, nYears, aveYears, false);
  // disMeanWeight
  extendArray(newDat.disMeanWeight, nMYears, nYears, aveYears, false);
  // landMeanWeight
  extendArray(newDat.landMeanWeight, nMYears, nYears, aveYears, false);
  // propF
  extendArray(newDat.propF, nMYears, nYears, aveYears, false);
  // propM
  extendArray(newDat.propM, nMYears, nYears, aveYears, false);
  newDat.noYears = nYears;

  // Prepare conf
  confSet newConf = conf;
  //// Set random walk recruitment
  newConf.stockRecruitmentModelCode = 0;

  // Prepare parameters
  paraSet<T> newPar = par.template cast<T>();

  vector<T> selT = sel.template cast<T>();
  // Make logF array
  array<T> logF(selT.size(), nYears);
  logF.setZero();
  for(int i = 0; i < nYears; ++i)
    logF.col(i) += log(selT) + log(Fbar);

  // Make logN array - start with one recruit
  int nAge = conf.maxAge - conf.minAge + 1;
  array<T> logN(nAge, nYears);
  logN.setZero();
  logN += SAM_NegInf;
  logN(0,0) = 0.0;

  // Run loop over years
  for(int i = 1; i < nYears; ++i){
    // predN
    logN.col(i) = predNFun(newDat, newConf, newPar, logN, logF, i);
    // remove recruitment
    logN(0,i) = SAM_NegInf;
  }

  // Calculate yield
  vector<T> cat(nYears);
  switch(newDat.referencepoint.catchType){
  case referencepointSet<T>::CatchType::totalCatch:
    cat = catchFun(newDat, newConf, logN, logF);
    break;
  case referencepointSet<T>::CatchType::landings:
    cat = landFun(newDat, newConf, logN, logF);
    break;
  case referencepointSet<T>::CatchType::discard:
    cat = disFun(newDat, newConf, logN, logF);
    break;
  default:
    Rf_error("Unknown reference point catch type.");
      break;
  }
  T logYPR = log(sum(cat));
  // Calculate spawners
  vector<T> ssb = ssbFun(newDat, newConf, logN, logF);
  T logSPR = log(sum(ssb));
  T lambda = exp(logSPR);


  if(conf.stockRecruitmentModelCode == 0 ||
     conf.stockRecruitmentModelCode == 3){
    PERREC_t<T> res = {log(Fbar), // logFbar
		       logYPR,	   // logYPR
		       logSPR,	   // logSPR
		       R_NaReal,	   // logSe
		       R_NaReal,		 // logRe
		       R_NaReal}; // logYe
    return res;
  }
  



  // Calculate Se
  T Se = 0.0; //R_NegInf;
  
  T dsr0 = 10000.0;
  if(conf.stockRecruitmentModelCode != 62)
    dsr0 = dFunctionalSR(T(SAM_Zero), newPar.rec_pars, conf.stockRecruitmentModelCode);
  
  switch(conf.stockRecruitmentModelCode){
  case 0: // straight RW 
    Rf_error("Equilibrium SSB not implemented");
    break;
  case 1: //ricker
    Se = exp(-newPar.rec_pars(1)) * log((exp(newPar.rec_pars(0)) * lambda));
    break;
  case 2:  //BH
    Se = (exp(newPar.rec_pars(0)) * lambda - 1.0) * exp(-newPar.rec_pars(1));
    break;
  case 3: //Constant mean
    Rf_error("Equilibrium SSB not implemented");
    break;
  case 61: // Hockey stick
    Se = lambda * exp(newPar.rec_pars(0));
    break;
  case 62: // AR1 (on log-scale)
    Se = lambda * exp(newPar.rec_pars(0));
    break;
  case 63: //Bent hyperbola / Hockey-stick-like
    Se = (2.0 * sqrt(exp(2.0 * newPar.rec_pars(0)) + exp(2.0 * newPar.rec_pars(2)) / 4.0) / (lambda * exp(newPar.rec_pars(1))) - 2.0 * exp(newPar.rec_pars(0)) - 2.0 * sqrt(exp(2.0 * newPar.rec_pars(0)) + exp(2.0 * newPar.rec_pars(2)) / 4.0)) / ( 1.0 / ((lambda * lambda * exp(2.0 * newPar.rec_pars(1)))) - 2.0 / (lambda * exp(newPar.rec_pars(1)))  );  
    break;
  case 64: // Cushing
    Se = exp(1.0 / (1 - exp(newPar.rec_pars(1))) * (newPar.rec_pars(0) + log(lambda)));
    break;
  case 65: // Shepherd
    Se = exp( newPar.rec_pars(1) + 1.0 / exp(newPar.rec_pars(2)) * log(exp(newPar.rec_pars(0)) * lambda - 1.0));
    break;
    default:
      error("SR model code not recognized");
    break;   
  }

  // Type logSe = log(Se);
  // Type logSe = CppAD::CondExpGt(exp(-logSPR), dsr0, Type(SAM_NegInf), log(Se));
  T logSe = CppAD::CondExpGt(exp(-logSPR), dsr0, log(Se) - 10.0 * (exp(-logSPR) - dsr0), log(Se));

  // Return
  PERREC_t<T> res = {log(Fbar), // logFbar
			logYPR,	   // logYPR
			logSPR,	   // logSPR
			logSe,	   // logSe
			logSe - logSPR,		 // logRe
			logSe - logSPR + logYPR}; // logYe

  return res;
}

  

template<class Type>
struct REFERENCE_POINTS {


  

  // Input data
  dataSet<Type>& dat;
  confSet& conf;
  paraSet<Type>& par;

  int nYears;
  vector<int> aveYears;
  vector<int> selYears;

  array<Type>& logN;
  array<Type>& logF;
  
  // Input F

  Type logFsq;			// Status quo
  Type logF0;			// "No" fishing
  Type logFmsy; 		// Maximizes yield
  Type logFmax;			// Maximizes yield per recruit
  Type logF01;			// F such that YPR'(0) = 0.1 * YPR'(F)
  Type logFcrash;		// F such that 1/SPR(f) = SR'(0) (i.e. stock crashes if slope of spawner-per-recruit in origin is less than slope of stock-recruitment model in origin)
  Type logF35;			// F such that SSB is reduced to 35% of unfished stock (Se(F) = 0.35 * Se(0) )
  //Type logFmed;			// Fishing  mortality  rate  F  corresponding  to  a  SSB/R  equal  to  the  inverse  of  the  50th  percentile of the observed R/SSB
  Type logFlim;			// F such that Se(F) = Blim (for hockey-stick-like stock recruitment only)
  //Type logFpa;			// F that corresponds to logBpa
   

  // Corresponding SSB??
  Type logBsq;
  Type logB0;
  Type logBmsy;
  Type logBmax;
  Type logB01;
  Type logBcrash;
  Type logB35;
  Type logBlim;			// Known from model parameters (for hockey-stick-like stock recruitment only)
  //Type logBpa;			// Ba = Blim * exp(1.645 * sigma) where sigma is the standard deviation of log(SSB) at the start of the year following the terminal year of the assessment if sigma is unknown, 0.2 can be used as default.

  // Corresponding Yield??
  Type logYsq;
  Type logY0;
  Type logYmsy;
  Type logYmax;
  Type logY01;
  Type logYcrash;
  Type logY35;
  Type logYlim;			// Known from model parameters (for hockey-stick-like stock recruit

  Type logYPRsq;
  Type logYPR0;
  Type logYPRmsy;
  Type logYPRmax;
  Type logYPR01;
  Type logYPRcrash;
  Type logYPR35;
  Type logYPRlim;

  Type logSPRsq;
  Type logSPR0;
  Type logSPRmsy;
  Type logSPRmax;
  Type logSPR01;
  Type logSPRcrash;
  Type logSPR35;
  Type logSPRlim;

  
  // Derived values
  vector<Type> sel;

  CppAD::ADFun<Type> FSR;
  CppAD::ADFun<Type> FYPR;
  CppAD::ADFun<Type> FSPR;
  
  REFERENCE_POINTS(){}
  REFERENCE_POINTS(dataSet<Type>& dat_,
		   confSet& conf_,
		   paraSet<Type>& par_,
		   array<Type>& logN_,
		   array<Type>& logF_	        
		   ): dat(dat_), conf(conf_), par(par_),
		      logN(logN_),
		      logF(logF_)
  {
    nYears = dat.referencepoint.nYears;
    aveYears = dat.referencepoint.aveYears;
    selYears = dat.referencepoint.selYears;

    // Calculate current selectivity and status quo
    vector<Type> fbartmp(selYears.size());
    fbartmp.setZero();
    sel = vector<Type>(logF.rows());
    sel.setZero();

    for(int y = 0; y < selYears.size(); ++y){
      sel += exp(logF.col(selYears(y)));
      for(int a = conf.fbarRange(0); a <= conf.fbarRange(1); a++){  
	fbartmp(y) += exp(logF(conf.keyLogFsta(0,a-conf.minAge),selYears(y)));
      }
      fbartmp(y) /= Type(conf.fbarRange(1)-conf.fbarRange(0)+1);
    }

    logFsq = log(sum(fbartmp)) - log(fbartmp.size());
    sel *= 1.0 / sum(fbartmp);
    
    logBsq = log(Se(exp(logFsq)));
    logYsq = log(yield(exp(logFsq)));
    logYPRsq = log(YPR(exp(logFsq)));
    logSPRsq = log(SPR(exp(logFsq)));

    logF0 = SAM_NegInf;
    logB0 = log(Se(exp(logF0)));
    logY0 = log(yield(exp(logF0)));
    logYPR0 = log(YPR(exp(logF0)));
    logSPR0 = log(SPR(exp(logF0)));

    // Calculate actual F values
    if(CppAD::Variable(par.logScaleFmsy)){
      logFmsy = logFsq + par.logScaleFmsy;
      logBmsy = log(Se(exp(logFmsy)));
      logYmsy = log(yield(exp(logFmsy)));
      logYPRmsy = log(YPR(exp(logFmsy)));
      logSPRmsy = log(SPR(exp(logFmsy)));
    }else{
      logFmsy = R_NaReal;// R_NaReal;
      logBmsy = R_NaReal;
      logYmsy = R_NaReal;
      logYPRmsy = R_NaReal;
      logSPRmsy = R_NaReal;
    }

    
    if(CppAD::Variable(par.logScaleFmax)){
      logFmax = logFsq + par.logScaleFmax;
      logBmax = log(Se(exp(logFmax)));
      logYmax = log(yield(exp(logFmax)));
      logYPRmax = log(YPR(exp(logFmax)));
      logSPRmax = log(SPR(exp(logFmax)));
    }else{
      logFmax = R_NaReal;//R_NaReal;
      logBmax = R_NaReal;
      logYmax = R_NaReal;
      logYPRmax = R_NaReal;
      logSPRmax = R_NaReal;
    }


    if(CppAD::Variable(par.logScaleF01)){
      logF01 = logFsq + par.logScaleF01;
      logB01 = log(Se(exp(logF01)));
      logY01 = log(yield(exp(logF01)));
      logYPR01 = log(YPR(exp(logF01)));
      logSPR01 = log(SPR(exp(logF01)));
    }else{
      logF01 = R_NaReal;//R_NaReal;
      logB01 = R_NaReal;
      logY01 = R_NaReal;
      logYPR01 = R_NaReal;
      logSPR01 = R_NaReal;
    }


    if(CppAD::Variable(par.logScaleFcrash)){
      logFcrash = logFsq + par.logScaleFcrash;
      logBcrash = log(Se(exp(logFcrash)));
      logYcrash = log(yield(exp(logFcrash)));
      logYPRcrash = log(YPR(exp(logFcrash)));
      logSPRcrash = log(SPR(exp(logFcrash)));
    }else{
      logFcrash = R_NaReal;//R_NaReal;
      logBcrash = R_NaReal;
      logYcrash = R_NaReal;
      logYPRcrash = R_NaReal;
      logSPRcrash = R_NaReal;
    }


    if(CppAD::Variable(par.logScaleF35)){
      logF35 = logFsq + par.logScaleF35;
      logB35 = log(Se(exp(logF35)));
      logY35 = log(yield(exp(logF35)));
      logYPR35 = log(YPR(exp(logF35)));
      logSPR35 = log(SPR(exp(logF35)));
    }else{
      logF35 = R_NaReal;//R_NaReal;
      logB35 = R_NaReal;
      logY35 = R_NaReal;
      logYPR35 = R_NaReal;
      logSPR35 = R_NaReal;
    }
    if(CppAD::Variable(par.logScaleFlim) &&
       (conf.stockRecruitmentModelCode == 61 ||
	conf.stockRecruitmentModelCode == 63)){    
      logFlim = logFsq + par.logScaleFlim;
      logYlim = log(yield(exp(logFlim)));
      if(conf.stockRecruitmentModelCode == 61){
	logBlim = par.rec_pars(1);
      }else if(conf.stockRecruitmentModelCode == 63){
	logBlim = par.rec_pars(0);      
      }
      logYPRlim = log(YPR(exp(logFlim)));
      logSPRlim = log(SPR(exp(logFlim)));
    }else{
      logFlim = R_NaReal;//R_NaReal;
      if(conf.stockRecruitmentModelCode == 61){
	logBlim = par.rec_pars(1);
      }else if(conf.stockRecruitmentModelCode == 63){
	logBlim = par.rec_pars(0);      
      }else{
	logBlim = R_NaReal;
      }
      logYlim = R_NaReal;
      logYPRlim = R_NaReal;
      logSPRlim = R_NaReal;
    }
 
 
    
    // Prepare AD
    vector<Type> Fsqvec(1);
    Fsqvec(0) = exp(logFsq);
    CppAD::vector<AD<Type> > x1( Fsqvec );
    CppAD::vector<AD<Type> > y1( 1 );
    CppAD::Independent(x1);
    y1[0] = YPR(x1);
    FYPR = CppAD::ADFun<Type>(x1, y1);
    
    CppAD::vector<AD<Type> > x2( Fsqvec );
    CppAD::vector<AD<Type> > y2( 1 );
    CppAD::Independent(x2);
    y2[0] = SPR(x2);
    FSPR = CppAD::ADFun<Type>(x2, y2);

    vector<Type> Bsqvec(1);
    Bsqvec(0) = exp(logBsq);
    CppAD::vector<AD<Type> > x3( Bsqvec );
    CppAD::vector<AD<Type> > y3( 1 );
    CppAD::Independent(x3);
    y3[0] = SR(x3);
    FSR = CppAD::ADFun<Type>(x3, y3);

  }

  Type YPR(Type Fbar){
    PERREC_t<Type> r = perRecruit<Type, Type>(Fbar, dat, conf, par, sel, aveYears, nYears);
    return exp(r.logYPR);
  }

  AD<Type> YPR(CppAD::vector<AD<Type> > Fbar){
    PERREC_t<AD<Type> > r = perRecruit<Type, AD<Type> >(Fbar[0], dat, conf, par, sel, aveYears, nYears);
    return exp(r.logYPR);
  }
  
  Type dYPR(Type Fbar){
      vector<Type> Fv(1);
      Fv(0) = Fbar;
      CppAD::vector<Type> x_eval( Fv );
      return FYPR.Jacobian(x_eval)[0];
  }

  Type SPR(Type Fbar){
    PERREC_t<Type> r = perRecruit<Type, Type>(Fbar, dat, conf, par, sel, aveYears, nYears);
    return exp(r.logSPR);
  }

  AD<Type> SPR(CppAD::vector<AD<Type> > Fbar){
    PERREC_t<AD<Type> > r = perRecruit<Type, AD<Type> >(Fbar[0], dat, conf, par, sel, aveYears, nYears);
    return exp(r.logSPR);
  }

  Type dSPR(Type Fbar){
    vector<Type> Fv(1);
    Fv(0) = Fbar;
    CppAD::vector<Type> x_eval( Fv );
    return FSPR.Jacobian(x_eval)[0];
  }

  Type Se(Type Fbar){
    PERREC_t<Type> r = perRecruit<Type, Type>(Fbar, dat, conf, par, sel, aveYears, nYears);
    return exp(r.logSe);
  }

  Type yield(Type Fbar){
    PERREC_t<Type> r = perRecruit<Type, Type>(Fbar, dat, conf, par, sel, aveYears, nYears);
    return exp(r.logYe);
  }

  Type SR(Type ssb){
    return exp(functionalStockRecruitment(ssb, par.rec_pars, conf.stockRecruitmentModelCode));
  }

  AD<Type> SR(CppAD::vector<AD<Type> > ssb){
    AD<Type> s0 = ssb[0];
    vector<AD<Type> > rp(par.rec_pars.size());
    rp = par.rec_pars.template cast<AD<Type> >();
    if(conf.stockRecruitmentModelCode == 0)
      return 0.0;
    return exp(functionalStockRecruitment(s0, rp, conf.stockRecruitmentModelCode));
  }


  Type dSR(Type ssb){
    vector<Type> Fv(1);
    Fv(0) = ssb;
    CppAD::vector<Type> x_eval( Fv );
    return FSR.Jacobian(x_eval)[0];
  }

  // vector<Type> dSR(Type ssb, vector<Type> rp){
  // CppAD::vector<AD<Type> > x( rp );
  //     CppAD::vector<AD<Type> > y( 1 );
  //     CppAD::Independent(x);
  //     y[0] = functionalStockRecruitment(ssb, rp, conf.stockRecruitmentModelCode);
  //     CppAD::ADFun<Type> F(x, y);
  //     CppAD::vector<Type> x_eval( rp );
  //     return F.Jacobian(x_eval);
  // }


  
  // Calculate "likelihood" contribution to estimate reference points
  Type operator()() {

    Type nll = 0.0;

    if(CppAD::Variable(par.logScaleFmsy)){
      nll -= log(yield(exp(logFmsy)));
    }

    if(CppAD::Variable(par.logScaleFmax)){
      nll -= log(YPR(exp(logFmax)));
    }

    if(CppAD::Variable(par.logScaleF01)){
     Type tmp = 0.1 * dYPR(Type(SAM_Zero)) - dYPR(exp(logF01));
      nll += tmp * tmp;
    }

    if(CppAD::Variable(par.logScaleFcrash)){
      Type tmp = dSR(Type(SAM_Zero)) - (1.0 / SPR(exp(logFcrash)));
      nll += tmp * tmp;
    }
    
    if(CppAD::Variable(par.logScaleF35)){
      Type tmp = 0.35 * SPR(Type(SAM_Zero)) - SPR(exp(logF35));
      nll += tmp * tmp;
    }

    if(!isNA(logBlim) && CppAD::Variable(par.logScaleFlim)){
      Type tmp = logBlim - log(Se(exp(logFlim)));
      nll += tmp * tmp;
    }
    
    return nll;
  }

};

template<class Type>
Type nllReferencepoints(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, array<Type> &logF, objective_function<Type> *of){

  if(dat.referencepoint.nYears == 0)
    return 0.0;

  REFERENCE_POINTS<Type> referencepoint(dat, conf, par, logN, logF);


  if(dat.referencepoint.Fsequence.size() > 0){
    vector<Type> Fseq = dat.referencepoint.Fsequence;
    vector<Type> logYPR(Fseq.size());
    logYPR.setZero();
    vector<Type> logSPR(Fseq.size());
    logSPR.setZero();
    vector<Type> logSe(Fseq.size());
    logSe.setZero();
    vector<Type> logYe(Fseq.size());
    logYe.setZero();
    vector<Type> logRe(Fseq.size());
    logRe.setZero();

    for(int i = 0; i < Fseq.size(); ++i){
      PERREC_t<Type> v = perRecruit<Type, Type>(Fseq(i),
						referencepoint.dat,
						referencepoint.conf,
						referencepoint.par,
						referencepoint.sel,
						referencepoint.aveYears,
						referencepoint.nYears);
      logYPR(i) = v.logYPR;
      logSPR(i) = v.logSPR;
      logSe(i) = v.logSe;
      logYe(i) = v.logYe;
      logRe(i) = v.logRe;
	
    }
    ADREPORT_F(logYPR, of);
    ADREPORT_F(logSPR, of);
    ADREPORT_F(logSe, of);
    ADREPORT_F(logYe, of);
    ADREPORT_F(logRe, of);
  }
  
  ADREPORT_F(referencepoint.logFsq,of);
  ADREPORT_F(referencepoint.logBsq,of);
  ADREPORT_F(referencepoint.logYsq,of);
  ADREPORT_F(referencepoint.logYPRsq,of);
  ADREPORT_F(referencepoint.logSPRsq,of);

  ADREPORT_F(referencepoint.logF0,of);
  ADREPORT_F(referencepoint.logB0,of);
  ADREPORT_F(referencepoint.logY0,of);
  ADREPORT_F(referencepoint.logYPR0,of);
  ADREPORT_F(referencepoint.logSPR0,of);
  
  ADREPORT_F(referencepoint.logFmsy,of);
  ADREPORT_F(referencepoint.logBmsy,of);
  ADREPORT_F(referencepoint.logYmsy,of);
  ADREPORT_F(referencepoint.logYPRmsy,of);
  ADREPORT_F(referencepoint.logSPRmsy,of);

  ADREPORT_F(referencepoint.logFmax,of);
  ADREPORT_F(referencepoint.logBmax,of);
  ADREPORT_F(referencepoint.logYmax,of);
  ADREPORT_F(referencepoint.logYPRmax,of);
  ADREPORT_F(referencepoint.logSPRmax,of);
 
  ADREPORT_F(referencepoint.logF01,of);
  ADREPORT_F(referencepoint.logB01,of);
  ADREPORT_F(referencepoint.logY01,of);
  ADREPORT_F(referencepoint.logYPR01,of);
  ADREPORT_F(referencepoint.logSPR01,of);

  ADREPORT_F(referencepoint.logFcrash,of);
  ADREPORT_F(referencepoint.logBcrash,of);
  ADREPORT_F(referencepoint.logYcrash,of);
  ADREPORT_F(referencepoint.logYPRcrash,of);
  ADREPORT_F(referencepoint.logSPRcrash,of);

  ADREPORT_F(referencepoint.logF35,of);
  ADREPORT_F(referencepoint.logB35,of);
  ADREPORT_F(referencepoint.logY35,of);
  ADREPORT_F(referencepoint.logYPR35,of);
  ADREPORT_F(referencepoint.logSPR35,of);

  ADREPORT_F(referencepoint.logFlim,of);
  ADREPORT_F(referencepoint.logBlim,of);
  ADREPORT_F(referencepoint.logYlim,of);
  ADREPORT_F(referencepoint.logYPRlim,of);
  ADREPORT_F(referencepoint.logSPRlim,of);

  Type ans = referencepoint();
  return ans;
}




// R functions

extern "C" {

  SEXP perRecruitR(SEXP Fbar, SEXP dat, SEXP conf, SEXP pl, SEXP sel, SEXP aveYears, SEXP nYears){
    dataSet<double> d0(dat);
    confSet c0(conf);
    paraSet<double> p0(pl);
    vector<double> s0 = asVector<double>(sel);
    vector<int> a0 = asVector<int>(aveYears);
    double Fbar0 = Rf_asReal(Fbar);
    int nY0 = Rf_asInteger(nYears);
 
    PERREC_t<double> y = perRecruit<double, double>(Fbar0, d0, c0, p0, s0, a0, nY0);
    const char *resNms[] = {"logF", "logYPR", "logSPR", "logSe", "logRe", "logYe", ""}; // Must end with ""
    SEXP res;
    PROTECT(res = Rf_mkNamed(VECSXP, resNms));
    SET_VECTOR_ELT(res, 0, asSEXP(y.logFbar));
    SET_VECTOR_ELT(res, 1, asSEXP(y.logYPR));
    SET_VECTOR_ELT(res, 2, asSEXP(y.logSPR));
    SET_VECTOR_ELT(res, 3, asSEXP(y.logSe));
    SET_VECTOR_ELT(res, 4, asSEXP(y.logRe));
    SET_VECTOR_ELT(res, 5, asSEXP(y.logYe));

    UNPROTECT(1);    
    return res;

  }


  SEXP stockRecruitmentModelR(SEXP ssb, SEXP rec_pars, SEXP code){
    double b = Rf_asReal(ssb);
    vector<double> rp = asVector<double>(rec_pars);
    int srmc = Rf_asInteger(code);
	
    double v = exp(functionalStockRecruitment(b, rp, srmc));

    vector<AD<double> > rp2( rp.size() );
    for(int i = 0; i < rp2.size(); ++i)
      rp2(i) = AD<double>(rp(i));
    CppAD::Independent(rp2);

    vector<AD<double> > x( 1 );
    x[0] = b;

    vector<AD<double> > y( 1 );
    y[0] = exp(functionalStockRecruitment(x[0], rp2, srmc));
    CppAD::ADFun<double> F(x, y);
    vector<double> x_eval( rp );
    vector<double> r = F.Jacobian(x_eval);

    const char *resNms[] = {"Recruits", "Gradient", ""}; // Must end with ""
    SEXP res;
    PROTECT(res = Rf_mkNamed(VECSXP, resNms));
    SET_VECTOR_ELT(res, 0, asSEXP(v));
    SET_VECTOR_ELT(res, 1, asSEXP(r));
 
    UNPROTECT(1);    
    return res;
      
  }
  
}
