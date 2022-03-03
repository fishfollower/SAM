#ifndef SAM_PERRECRUIT_STOCHASTIC_HPP
#define SAM_PERRECRUIT_STOCHASTIC_HPP

#define SAM_NegInf -20.0
#define SAM_NIZero -10.0

#define SAM_Zero exp(SAM_NIZero)

// Functions to simulate population dynamics and reference points

// 1) Spawners per recruit
// 2) Yield per recruit
// 3) Equilibrium Biomass
// 4) Equilibrium yield
// (life years lost is not random for fixed F and M)
// 5) calculate reference points for a given simulation

// Not to be used with CppAD
template<class Type>
PERREC_t<Type> stochPerRecruit(Type logFbar, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, vector<Type>& logSel, vector<int> aveYears, int nYears, int CT, vector<Type> logNinit){

  // Extend arrays
  dataSet<Type> newDat = dat;
  int nMYears = dat.noYears;
  // propMat
  bool guessNY = nYears < 1;
  if(guessNY){		// Calculate nYears from M
    array<Type> pm = newDat.propMat;
    extendArray(pm, nMYears, 1, aveYears, false);
    Type Mlt = pm(pm.dim[0]-1,pm.dim[1]-1); // Last year last age
    Type Flt = exp(logFbar + logSel(logSel.size()-1)); // Last age
    Type Zlt = Mlt + Flt;
    nYears = std::max(pm.dim[1] - (-60) / Zlt, (Type)pm.dim[1] * 20.0);       
  }
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


  // Make logF array
  array<Type> logF(logSel.size(), nYears);
  logF.setZero();
  for(int i = 0; i < nYears; ++i)
    logF.col(i) = logSel + logFbar;

  // Make logN array - start with one recruit
  int nAge = conf.maxAge - conf.minAge + 1;
  array<Type> logN(nAge, nYears);
  logN.setConstant(R_NegInf);
  logN(0,0) = 0.0;
  

  matrix<Type> nvar = get_nvar(newDat, conf, par, logN, logF);
  MVMIX_t<Type> neg_log_densityN(nvar,Type(conf.fracMixN),false);
  for(int i = 1; i < nYears; ++i){
    vector<Type> predN = predNFun(newDat, conf, par, logN, logF, i);
    logN.col(i) = predN + neg_log_densityN.simulate();
    // Remove recruitment
    logN(0,i) = R_NegInf;    
  }
  
  vector<Type> cat(nYears);
  cat.setZero();
  typename referencepointSet<Type>::CatchType catchType = static_cast<typename referencepointSet<Type>::CatchType>(CT);
  switch(catchType){
  case referencepointSet<Type>::totalCatch:
    cat = catchFun(newDat, conf, logN, logF);
    break;
  case referencepointSet<Type>::landings:
    cat = landFun(newDat, conf, logN, logF);
    break;
  case referencepointSet<Type>::discard:
    cat = disFun(newDat, conf, logN, logF);
    break;
  default:
    Rf_error("Unknown reference point catch type.");
    break;
  }
  Type logYPR = log(sum(cat));//

  Type logYLTF = log(yearsLostFishing_i(newDat, conf, logF, newDat.natMor.dim(0)-1, conf.minAge, conf.maxAge));
  Type logLifeExpectancy = log(temporaryLifeExpectancy_i(newDat, conf, logF, newDat.natMor.dim(0)-1, conf.minAge, 10 * conf.maxAge) + (Type)conf.minAge);
  vector<Type> ssb = ssbFun(newDat, conf, logN, logF);
  Type logSPR = log(sum(ssb)); //log(sum(ssb)); log(sum(ssb) + (T)exp(-12.0));


  ////////////////////////////////////////////////////////////////////////////////
  // Survival calculations                                                      //
  ////////////////////////////////////////////////////////////////////////////////
  Type discYPR = 0.0;
  // Custom recursive calculation inspired by survival.hpp
  Type yl_logp = 0.0;
  Type yl_q = 0.0;
  Type yl = 0.0;
 
  for(int aa = 0; aa < cat.size(); ++aa){
    int j = std::min(std::max(aa-conf.minAge,0), newDat.natMor.cols()-1);		// Cohort age index
    int i = std::min(aa, newDat.natMor.rows()-1);
    Type M = newDat.natMor(i,j);
    Type F = 0.0;    
    if(conf.keyLogFsta(0,j)>(-1) && aa >= conf.minAge)
      F += exp(logF(conf.keyLogFsta(0,j),i));
    Type Z = M + F;
    yl += yl_q + exp(yl_logp) * F / Z * (1.0 - 1.0 / Z * (1.0 - exp(-Z)));    
    discYPR += cat(aa) * exp(-yl);
    // discYe += catYe(aa) * exp(-yl);
    yl_q += exp(yl_logp) * F / Z * (1.0 - exp(-Z));
    yl_logp += -Z;
  }
  Type logDiscYPR = log(discYPR);
  ////////////////////////////////////////////////////////////////////////////////
  // End survival calculations                                                  //
  ////////////////////////////////////////////////////////////////////////////////


  
  if(conf.stockRecruitmentModelCode == 0){//  ||
    // conf.stockRecruitmentModelCode == 3){
    PERREC_t<Type> res = {logFbar, // logFbar
			  logYPR,	   // logYPR
			  logSPR,	   // logSPR
			  R_NaReal,	   // logSe
			  R_NaReal,		 // logRe
			  R_NaReal,// logYe
			  R_NaReal,// dSR0
			  logLifeExpectancy, // logLifeExpectancy
			  logYLTF, // logYearsLost
			  logDiscYPR, // logDiscYPR
			  R_NaReal	 // logDiscYe
    }; 
    return res;
  }

 

  // Simulate biomass
  //  int nYearsEqMax = (guessNY) ? 1000 : nYears;
  array<Type> logNeq(nAge, nYears);
  logNeq.setConstant(R_NegInf);
  // logNeq.col(0) = logNinit;
 
  // array<Type> logFeq0(logSel.size(), nYearsEqMax);
  // logFeq0.setZero();
  // for(int i = 0; i < nYearsEqMax; ++i)
  //   logFeq0.col(i) = logSel + logFbar;

  
  if(guessNY){
    // Improve Initial value by starting stochastic simulation at deterministic equilibrium, but avoid a crashed equilibrium.
    vector<Type> logNinit2 = logNinit;
    array<Type> logFInit(logF.rows(), 2);
    for(int i = 0; i < logFInit.cols(); ++i)
      logF.col(i) = logSel + logFbar;
    // array<Type> noF(logF.rows(),1);
    // noF.setConstant(R_NegInf);
    Type imp = 100.0;
    while(imp > 1e-6){
      array<Type> tmp(nAge,2);
      tmp.setZero();
      tmp.col(0) = logNinit2;
      // Try unfished equilibrium? Or adaptively reduce F if crashing?
      vector<Type> lnitmp = predNFun(newDat, conf, par, tmp, logFInit, 1);
      if(conf.minAge == 0){
	tmp.col(1) = lnitmp;
	lnitmp(0) = predNFun(newDat, conf, par, tmp, logFInit, 1)(0);
      }
      if(lnitmp(0) < logNinit(0) + log(0.01)){ // Probably heading for crash
	imp = 0.0;
      }else{
	imp = (logNinit2 - lnitmp).abs().maxCoeff();
      }
      logNinit2 = lnitmp;    
    }
    logNeq.col(0) = logNinit2;
    // for(int i = 0; i < logNinit2.size(); ++i) Rprintf("%d: %.03f\n",i,logNinit2(i));
    // Rprintf("\n");
  }
   
  matrix<Type> nvareq = get_nvar(newDat, conf, par, logNeq, logF);
  MVMIX_t<Type> neg_log_densityNeq(nvareq,Type(conf.fracMixN));
  
  // vector<Type> runMean = exp(logNeq0.col(0));
  // runMean.setZero();
  // int nYearsEq = 1;
  for(int i = 1; i < nYears; ++i){
    //Harvest control rule?
    vector<Type> predN = predNFun(newDat, conf, par, logNeq, logF, i);
    vector<Type> noiseN = neg_log_densityNeq.simulate();
    logNeq.col(i) = predN + noiseN;
    if(conf.minAge == 0){
      // In this case, predicted recruitment was wrong since it depends on SSB the same year
      // Predict recruitment again with updated ssb (assuming maturity at age 0 is 0):
      Type predRec = predNFun(newDat, conf, par, logNeq, logF, i)(0);
      // Overwrite recruitment, but keep simulated noise to retain correlation:
      logNeq(0,i) = predRec + noiseN(0);      
    }
  }

  vector<Type> catYe(nYears);
  catYe.setZero();
  switch(catchType){
  case referencepointSet<Type>::totalCatch:
    catYe = catchFun(newDat, conf, logNeq, logF);
    break;
  case referencepointSet<Type>::landings:
    catYe = landFun(newDat, conf, logNeq, logF);
    break;
  case referencepointSet<Type>::discard:
    catYe = disFun(newDat, conf, logNeq, logF);
    break;
  default:
    Rf_error("Unknown reference point catch type.");
    break;
  }

  vector<Type> ssbeq = ssbFun(newDat, conf, logNeq, logF);
  Type logSe = log(ssbeq(nYears-1));
  Type logYe = log(catYe(nYears-1));
  Type logRe = logNeq(0,nYears-1);
  Type logDiscYe = logDiscYPR + logRe;


  Type dsr0 = 10000.0;
  if(conf.stockRecruitmentModelCode != 0 &&
     conf.stockRecruitmentModelCode != 3 &&
     conf.stockRecruitmentModelCode != 62 &&
     conf.stockRecruitmentModelCode != 65 &&
     (conf.stockRecruitmentModelCode != 68) && // || newPar.rec_pars[2] < 0) &&
     (conf.stockRecruitmentModelCode != 69) &&
     conf.stockRecruitmentModelCode != 90 &&
     conf.stockRecruitmentModelCode != 91 &&
     conf.stockRecruitmentModelCode != 92){ // || newPar.rec_pars[2] < 0 ))
    dsr0 = dFunctionalSR(Type(SAM_Zero), par.rec_pars, conf.stockRecruitmentModelCode);
  }
  if(conf.stockRecruitmentModelCode == 68 || conf.stockRecruitmentModelCode == 69){
    dsr0 = CppAD::CondExpLt(par.rec_pars[2],
 			    (Type)0.0,
 			    dFunctionalSR(Type(SAM_Zero),
 					  par.rec_pars,
 					  conf.stockRecruitmentModelCode),
 			    dsr0);
  }

  if(conf.stockRecruitmentModelCode == 90 ||
     conf.stockRecruitmentModelCode == 91 ||
     conf.stockRecruitmentModelCode == 92){
    // dSplineSR uses logssb
    dsr0 = dSplineSR((Type)SAM_Zero,
 		     (vector<Type>)conf.constRecBreaks.template cast<Type>(),
 		     par.rec_pars,
 		     conf.stockRecruitmentModelCode);
  }

  
  // Return
  PERREC_t<Type> res = {logFbar, // logFbar
			logYPR,	// logYPR
			logSPR,	// logSPR
			logSe,	// logSe
			logRe,	// logRe
			logYe,	// logYe
			dsr0,
			logLifeExpectancy, // logLifeExpectancy
			logYLTF, // logYearsLost
			logDiscYPR,
			logDiscYe,
  };	// DSR0


  // PERREC_t<Type> res = {R_NaReal, // logFbar
  // 		       R_NaReal,	   // logYPR
  // 		       R_NaReal,	   // logSPR
  // 		       R_NaReal,	   // logSe
  // 		       R_NaReal,		 // logRe
  // 		       R_NaReal,// logYe
  // 		       R_NaReal,// dSR0
  // 		       R_NaReal, // logLifeExpectancy
  // 		       R_NaReal, // logYearsLost
  // 		       R_NaReal, // logDiscYPR
  // 		       R_NaReal	 // logDiscYe
  // }; 

  return res;
}







#endif
