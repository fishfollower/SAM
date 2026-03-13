SAM_DEPENDS(define)
SAM_DEPENDS(convenience)

#ifdef SAM_NegInf
#undef SAM_NegInf
#define SAM_NegInf R_NegInf
#endif

template<class Type>
vector<Type> toLogSeasonEffect(vector<Type> x, vector<int> isFishingSeason, vector<double> seasonTimes)SOURCE({
  vector<Type> r(x.size()+1);
  Type rs = 0.0;
  r.setZero();
  for(int i = 0; i < x.size(); ++i){
    Type v = (x(i));
    r(i) = v;
    rs = logspace_add(rs, v);
  }
  // r(r.size()-1) = 0.0;
  r -= rs;
  vector<Type> res(isFishingSeason.size());
  res.setConstant(SAM_NegInf);
  int indx = 0;
  for(int i = 0; i < res.size(); ++i)
    if(isFishingSeason(i)){
      res(i) = r(indx++) - log(Type(seasonTimes(i+1) - seasonTimes(i)));
    }
  return res;
  })


template <class Type>
vector<Type> totFFun(dataSet<Type>& dat, confSet &conf, array<Type> &logF, int y, int Ftype) SOURCE({
  int noFleets=conf.keyLogFsta.dim[0];
  int stateDimN=conf.keyLogFsta.dim[1];
  vector<Type> totF(stateDimN);
  totF.setZero();  
  for(int j=0; j<stateDimN; ++j){
    for(int f=0; f<noFleets; ++f){
      if(conf.keyLogFsta(f,j)>(-1)){
	Type Fval = exp(logF(conf.keyLogFsta(f,j),y)) * (dat.sampleTimesEnd(f) - dat.sampleTimesStart(f));
	if(Ftype == 1){
	  Fval *= dat.landFrac(y,j,f);
	}else if(Ftype == 2){
	  Fval *= (1.0 - dat.landFrac(y,j,f));
	}
	totF(j) += Fval;
      }
    }
  }
  return totF;
  })

SAM_SPECIALIZATION(vector<double> totFFun(dataSet<double>&, confSet&, array<double>&, int, int));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> totFFun(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, int, int));

template <class Type>
array<Type> totFFun(dataSet<Type>& dat, confSet &conf, array<Type> &logF, int Ftype DEFARG(= 0)) SOURCE({
    int stateDimN=conf.keyLogFsta.dim[1];
    int timeSteps=logF.dim[1];
    if(Ftype > 0)
      timeSteps = dat.landFrac.dim[0];
    array<Type> totF(stateDimN,timeSteps);
    totF.setZero();  
    for(int i=0;i<timeSteps;i++){ 
      totF.col(i) = totFFun(dat,conf,logF,i,Ftype);
    }
    return totF;
  })


SAM_SPECIALIZATION(array<double> totFFun(dataSet<double>&, confSet&, array<double>&, int));
SAM_SPECIALIZATION(array<TMBad::ad_aug> totFFun(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, int));



HEADER(
template<class Type>
struct MortalitySet {
  
  // TODO: Switch to log-scale calculations?
  matrix<Type> logCumulativeHazard;		// age x year
  array<Type> logCumulativeHazard_F;		// age x year x fleet
  array<Type> logCumulativeHazard_M;		// age x year x cause (M1 + risks)
  array<Type> FullYear_logCumulativeIncidence_Fishing;		// age x year x fleets
  array<Type> FullYear_logCumulativeIncidence_Other;		// age x year x cause (M1 + risks)
  matrix<Type> FullYear_logSurvival;		// age x year
  matrix<Type> Effective_logF;
  matrix<Type> Effective_logM;
  // matrix<Type> totalF;
  // array<Type> totalZseason;	// season x age x year
  // array<Type> totalFseason;
  // matrix<Type> totalFCI;
  array<Type> logFleetSurvival_before; // age x year x fleet (including surveys)
  array<Type> fleetLogCumulativeIncidence; // age x year x fleet (including surveys)
  array<Type> otherLogCumulativeIncidence; // age x year x causes (M1 + risks)
  matrix<Type> ssbLogSurvival_before;	     // age x year
  array<Type> Fseason;			     // seasons x year x (processes+1)

  std::vector<Type> activeHazard_breakpoints; // Vector of time points
  vector<int> activeHazard_season; // Key to the fishing season (-1 if none)
  array<int> activeHazard_F;	   // Indicator, 1 if fleet is active (brkpoints x ages x fleets)
  matrix<int> activeHazardMap_risk;	// Indicator, 1 if non-fishing hazard is active. (brkpoints  x causes)

  array<Type> logHazard_breakpoints;	// age x year x activeHazard_breakpoints
  array<Type> logHazard_F_breakpoints;	// age x year x fleet x activeHazard_breakpoints
  array<Type> logHazard_M_breakpoints;	// age x year x causes x activeHazard_breakpoints
  array<Type> logCIF_F_breakpoints;	// age x year x fleet x activeHazard_breakpoints
  array<Type> logCIF_M_breakpoints;	// age x year x causes x activeHazard_breakpoints
  
  inline MortalitySet():
    logCumulativeHazard(),
    logCumulativeHazard_F(),
    logCumulativeHazard_M(),
    FullYear_logCumulativeIncidence_Fishing(),
    FullYear_logCumulativeIncidence_Other(),
    FullYear_logSurvival(),
    Effective_logF(),
    Effective_logM(),
    // totalF(),
    // totalZseason(),
    // totalFseason(),
    logFleetSurvival_before(),
    fleetLogCumulativeIncidence(),
    otherLogCumulativeIncidence(),
    ssbLogSurvival_before(),
    Fseason(),
    activeHazard_breakpoints(),
    activeHazard_season(),
    activeHazard_F(),
    activeHazardMap_risk(),
    logHazard_breakpoints(),
    logHazard_F_breakpoints(),
    logHazard_M_breakpoints(),
    logCIF_F_breakpoints(),
    logCIF_M_breakpoints()
  {}
  
  template<class T>
  inline MortalitySet(const MortalitySet<T> x) : logCumulativeHazard(x.logCumulativeHazard),
						 logCumulativeHazard_F(x.logCumulativeHazard_F,x.logCumulativeHazard_F.dim),
						 logCumulativeHazard_M(x.logCumulativeHazard_M,x.logCumulativeHazard_M.dim),
						 FullYear_logCumulativeIncidence_Fishing(x.FullYear_logCumulativeIncidence_Fishing, x.FullYear_logCumulativeIncidence_Fishing.dim),
						 FullYear_logCumulativeIncidence_Other(x.FullYear_logCumulativeIncidence_Other, x.FullYear_logCumulativeIncidence_Other.dim),
						 FullYear_logSurvival(x.FullYear_logSurvival),
						 Effective_logF(x.Effective_logF),
						 Effective_logM(x.Effective_logM),
						 // totalF(x.totalF),
						 // totalZseason(x.totalZseason),
						 // totalFseason(x.totalFseason),
						 logFleetSurvival_before(x.logFleetSurvival_before,x.logFleetSurvival_before.dim),
						 fleetLogCumulativeIncidence(x.fleetLogCumulativeIncidence,x.fleetLogCumulativeIncidence.dim),
						 otherLogCumulativeIncidence(x.otherLogCumulativeIncidence,x.otherLogCumulativeIncidence.dim),
						 ssbLogSurvival_before(x.ssbLogSurvival_before),
						 Fseason(x.Fseason),
						 activeHazard_breakpoints(x.activeHazard_breakpoints),
						 activeHazard_season(x.activeHazard_season),
						 activeHazard_F(x.activeHazard_F,x.activeHazard_F.dim),
						 activeHazardMap_risk(x.activeHazardMap_risk),
						 logHazard_breakpoints(x.logHazard_breakpoints,x.logHazard_breakpoints.dim),
						 logHazard_F_breakpoints(x.logHazard_F_breakpoints,x.logHazard_F_breakpoints.dim),
						 logHazard_M_breakpoints(x.logHazard_M_breakpoints,x.logHazard_M_breakpoints.dim),
						 logCIF_F_breakpoints(x.logCIF_F_breakpoints,x.logCIF_F_breakpoints.dim),
						 logCIF_M_breakpoints(x.logCIF_M_breakpoints,x.logCIF_M_breakpoints.dim)
						 
  {}
  
  MortalitySet(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, array<Type>& logitFseason);

  void updateSeasons(dataSet<Type>& dat, confSet& conf, array<Type>& logitFseason, int y);
  void updateHazards(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, int y);
  // For simulation based forecast
  void updateYear(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, array<Type>& logitFseason, int y);

  Type logSurvival(int a, int y, Type t0, Type t1);
  Type logCIF(int fleet, int a, int y, Type t0, Type t1);
  Type partialLogCIF(int fleet, int a, int y, Type t0, Type t1);
  Type logCIFr(int risk, int a, int y, Type t0, Type t1);
  Type partialLogCIFr(int risk, int a, int y, Type t0, Type t1);
  
});



template<class Type>
SEXP asSEXP(const MortalitySet<Type> &x) SOURCE({
    const char *resNms[] = {"logCumulativeHazard", "logCumulativeHazard_F", "logCumulativeHazard_M", "FullYear_logCumulativeIncidence_Fishing", "FullYear_logCumulativeIncidence_Other", "FullYear_logSurvival", "Effective_logF", "Effective_logM", "logFleetSurvival_before","fleetLogCumulativeIncidence","otherLogCumulativeIncidence","ssbLogSurvival_before", "Fseason", "activeHazard_breakpoints", "activeHazard_season", "activeHazard_F","activeHazardMap_risk", "logHazard_breakpoints", "logHazard_F_breakpoints", "logHazard_M_breakpoints", "logCIF_F_breakpoints", "logCIF_M_breakpoints", ""}; // Must end with ""
    SEXP res;
    PROTECT(res = Rf_mkNamed(VECSXP, resNms));
    SET_VECTOR_ELT(res, 0, asSEXP(x.logCumulativeHazard));
    SET_VECTOR_ELT(res, 1, asSEXP(x.logCumulativeHazard_F));
    SET_VECTOR_ELT(res, 2, asSEXP(x.logCumulativeHazard_M));
    SET_VECTOR_ELT(res, 3, asSEXP(x.FullYear_logCumulativeIncidence_Fishing));
    SET_VECTOR_ELT(res, 4, asSEXP(x.FullYear_logCumulativeIncidence_Other));
    SET_VECTOR_ELT(res, 5, asSEXP(x.FullYear_logSurvival));
    SET_VECTOR_ELT(res, 6, asSEXP(x.Effective_logF));
    SET_VECTOR_ELT(res, 7, asSEXP(x.Effective_logM));
    SET_VECTOR_ELT(res, 8, asSEXP(x.logFleetSurvival_before));
    SET_VECTOR_ELT(res, 9, asSEXP(x.fleetLogCumulativeIncidence));
    SET_VECTOR_ELT(res, 10, asSEXP(x.otherLogCumulativeIncidence));
    SET_VECTOR_ELT(res, 11, asSEXP(x.ssbLogSurvival_before));
    SET_VECTOR_ELT(res, 12, asSEXP(x.Fseason));
    vector<Type> tmp(x.activeHazard_breakpoints.size());
    for(int i = 0; i < tmp.size(); ++i)
      tmp(i) = x.activeHazard_breakpoints[i];
    SET_VECTOR_ELT(res, 13, asSEXP(tmp));
    SET_VECTOR_ELT(res, 14, asSEXP(x.activeHazard_season));
    SET_VECTOR_ELT(res, 15, asSEXP(x.activeHazard_F));
    SET_VECTOR_ELT(res, 16, asSEXP(x.activeHazardMap_risk));
    SET_VECTOR_ELT(res, 17, asSEXP(x.logHazard_breakpoints));
    SET_VECTOR_ELT(res, 18, asSEXP(x.logHazard_F_breakpoints));
    SET_VECTOR_ELT(res, 19, asSEXP(x.logHazard_M_breakpoints));
    SET_VECTOR_ELT(res, 20, asSEXP(x.logCIF_F_breakpoints));
    SET_VECTOR_ELT(res, 21, asSEXP(x.logCIF_M_breakpoints));

    // Report RiskHazard function over range of input
    
    UNPROTECT(1);    
    return res;
  })


SAM_SPECIALIZATION(SEXP asSEXP(const MortalitySet<double>&));
SAM_SPECIALIZATION(SEXP asSEXP(const MortalitySet<TMBad::ad_aug>&));




SOURCE(
       template<class Type>
       void MortalitySet<Type>::updateSeasons(dataSet<Type>& dat, confSet& conf, array<Type>& logitFseason, int y){
	 //Fseason.setZero();
	 // Type NFseason = conf.isFishingSeason.sum();
	 // If no season info, use constant (i.e. 0)
	 for(int s = 0; s < Fseason.dim(0); ++s){
	   Fseason(s,y,0) = 0.0; //(Type)log(conf.isFishingSeason(s)) - log(NFseason);
	   //Fseason(s,y,0) -= log(conf.seasonTimes(s+1)-conf.seasonTimes(s));    
	 }
	 // Loop over processes
	 for(int p = 0; p < logitFseason.dim(2); ++p){
	   if(y + dat.years(0) <= (Type)conf.seasonFirstYear){
	     for(int s = 0; s < Fseason.dim(0); ++s)
	       Fseason(s,y,p+1) = 0.0;
	   }else{
	     vector<Type> x0(logitFseason.dim(0));
	     for(int qq = 0; qq < x0.size(); ++qq) x0(qq) = logitFseason(qq,y,p);
	     vector<Type> tfs = toLogSeasonEffect(x0, conf.isFishingSeason, conf.seasonTimes);
	     for(int s = 0; s < Fseason.dim(0); ++s)
	       Fseason(s,y,p+1) = tfs(s);
	   }
	 }
	 return;
       }
       )

SOURCE(
template<class Type>
void MortalitySet<Type>::updateHazards(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF,int y){ // Add par and logN
  // Hazard and cumulative incidence per interval
  int nIntervals = logHazard_breakpoints.dim(2);
  int nAges = logHazard_breakpoints.dim(0);
  int nFleet = conf.keyLogFsta.dim(0);
  int nRisk = dat.CompRisk.size() + 1;
  

  vector<Type> tmpCumHaz(nAges);
  tmpCumHaz.setConstant(SAM_NegInf);
  matrix<Type> tmpCumHaz_F(nAges,nFleet);
  tmpCumHaz_F.setConstant(SAM_NegInf);
  matrix<Type> tmpCumHaz_M(nAges,nRisk);
  tmpCumHaz_M.setConstant(SAM_NegInf);

  matrix<Type> tmpFYCIF_F(nAges,nFleet);
  tmpFYCIF_F.setConstant(SAM_NegInf);

  matrix<Type> tmpFYCIF_Other(nAges,nRisk);
  tmpFYCIF_Other.setConstant(SAM_NegInf);

  // logHazard_breakpoints.setConstant(SAM_NegInf);
  // logHazard_M_breakpoints.setConstant(SAM_NegInf);
  // logHazard_F_breakpoints.setConstant(SAM_NegInf);
  
  for(int i = 0; i < nIntervals; ++i){
    int s = activeHazard_season(i);
    Type dt = activeHazard_breakpoints[i+1] - activeHazard_breakpoints[i];
    for(int a = 0; a < nAges; ++a){
      logHazard_breakpoints(a,y,i) = SAM_NegInf;
      // Natural hazards (natMor + risks)
      // M1 (always active
      //if(activeHazard_M(i,a,0)){
      logHazard_breakpoints(a,y,i) = logspace_add_SAM(logHazard_breakpoints(a,y,i),log(dat.natMor(y,a)));
      logHazard_M_breakpoints(a,y,0,i) = log(dat.natMor(y,a)); //logspace_add_SAM(logHazard_M_breakpoints(a,y,0,i),log(dat.natMor(y,a)));
      tmpCumHaz_M(a,0) = logspace_add_SAM(tmpCumHaz_M(a,0),logHazard_M_breakpoints(a,y,0,i) + log(dt));
      // }
      // Competing risks
      if(nRisk > 1){
	for(int r = 1; r < nRisk; ++r)
	  if(activeHazardMap_risk(i,r-1) > -1){
	    int brk = activeHazardMap_risk(i,r-1);
	    int p = conf.keyCompRisk(r-1,a);
	    Type logh = SAM_NegInf;
	    if(p >= 0){
	      int yx = std::min(y,(int)dat.CompRisk(r-1).X.rows()-1);
	      Type mUse = par.cp_m(p);
	      if(!(R_IsNA(dat.CompRisk(r-1).mMin) || R_IsNA(dat.CompRisk(r-1).mMax)))
	        mUse = toInterval((Type)par.cp_m(p), (Type)dat.CompRisk(r-1).mMin, (Type)dat.CompRisk(r-1).mMax, Type(1.0));
	      logh = logRiskHazard(dat.CompRisk(r-1).X(yx,brk),mUse,par.cp_loga(p),par.cp_logb(p), dat.CompRisk(r-1).Model);
	    }
	    logHazard_breakpoints(a,y,i) = logspace_add_SAM(logHazard_breakpoints(a,y,i),logh);
	    logHazard_M_breakpoints(a,y,r,i) = logh; // logspace_add_SAM(logHazard_M_breakpoints(a,y,r,i),logh);
	    tmpCumHaz_M(a,r) = logspace_add_SAM(tmpCumHaz_M(a,r), logHazard_M_breakpoints(a,y,r,i) + log(dt));
	  }
      }
	// Fleets
      for(int f = 0; f < nFleet; ++f){
	logHazard_F_breakpoints(a,y,f,i) = SAM_NegInf;
	int indx = conf.keyLogFsta(f,a);
	if(indx > (-1) && activeHazard_F(i,a,f)){
	  int j = conf.keyLogFseason(f,a);	  
	  Type logFs = logF(indx,y);
	  if(s > (-1))
	    logFs += Fseason(s,y,j+1);
	  logHazard_breakpoints(a,y,i) = logspace_add_SAM(logHazard_breakpoints(a,y,i), logFs);
	  logHazard_F_breakpoints(a,y,f,i) = logFs; //logspace_add_SAM(logHazard_F_breakpoints(a,y,f,i),logFs);
	  tmpCumHaz_F(a,f) = logspace_add_SAM(tmpCumHaz_F(a,f), logHazard_F_breakpoints(a,y,f,i) + log(dt));
	}
      }
      Type Stmp = -exp(tmpCumHaz(a));
      tmpCumHaz(a) = logspace_add_SAM(tmpCumHaz(a),logHazard_breakpoints(a,y,i) + log(dt));
      // Calculate CIF M
      logCIF_M_breakpoints(a,y,0,i) = SAM_NegInf;
      //if(activeHazard_M(i,a,0))
      for(int r = 0; r < nRisk; ++r){
	//CIF_M_breakpoints(a,y,r,i) = Hazard_M_breakpoints(a,y,r,i) / Hazard_breakpoints(a,y,i) * (1.0 - exp(-Hazard_breakpoints(a,y,i) * dt));
	logCIF_M_breakpoints(a,y,r,i) = logHazard_M_breakpoints(a,y,r,i) - logHazard_breakpoints(a,y,i) + logspace_sub_SAM(Type(0.0),-exp(logHazard_breakpoints(a,y,i)) * dt);
	tmpFYCIF_Other(a,r) = logspace_add_SAM(tmpFYCIF_Other(a,r), Stmp + logCIF_M_breakpoints(a,y,r,i));
      }
      // Calculate CIF F
      for(int f = 0; f < nFleet; ++f){
	logCIF_F_breakpoints(a,y,f,i) = SAM_NegInf;
	int indx = conf.keyLogFsta(f,a);
	if(indx > (-1) && activeHazard_F(i,a,f)){
	  int j = conf.keyLogFseason(f,a);
	  // Type logFs = logF(indx,y);
	  // if(s > (-1))
	  //   logFs += Fseason(s,y,j+1);
	  //CIF_F_breakpoints(a,y,f,i) = exp(logFs) / Hazard_breakpoints(a,y,i) * (1.0 - exp(-Hazard_breakpoints(a,y,i) * dt));
	  logCIF_F_breakpoints(a,y,f,i) = logHazard_F_breakpoints(a,y,f,i) - logHazard_breakpoints(a,y,i) + logspace_sub_SAM(Type(0.0),-exp(logHazard_breakpoints(a,y,i)) * dt);
	}
	tmpFYCIF_F(a,f) = logspace_add_SAM(tmpFYCIF_F(a,f), Stmp + logCIF_F_breakpoints(a,y,f,i));
      }
    }
  }
  for(int a = 0; a < nAges; ++a){
    logCumulativeHazard(a,y) = tmpCumHaz(a);
    FullYear_logSurvival(a,y) = -exp(tmpCumHaz(a));
    for(int f = 0; f < nFleet; ++f)
      FullYear_logCumulativeIncidence_Fishing(a,y,f) = tmpFYCIF_F(a,f);
    for(int r = 0; r < nRisk; ++r)
      FullYear_logCumulativeIncidence_Other(a,y,r) = tmpFYCIF_Other(a,r);
    for(int f = 0; f < nFleet; ++f)
      logCumulativeHazard_F(a,y,f) = tmpCumHaz_F(a,f);
    for(int r = 0; r < nRisk; ++r)
      logCumulativeHazard_M(a,y,r) = tmpCumHaz_M(a,r);
    // Update effective F and M
    Type v = exp(FullYear_logSurvival(a,y));
    Type u = exp(FullYear_logCumulativeIncidence_Fishing(a,y,0));
    for(int f = 1; f < nFleet; ++f)
      u += exp(FullYear_logCumulativeIncidence_Fishing(a,y,f));
    Effective_logF(a,y) = log(-log(v)) + log(u) - log(1.0 - v);
    Effective_logM(a,y) = log(-log(v)) + log(1.0 - u - v) - log(1.0 - v);
  }
  return;
}
       )


SOURCE(
       template<class Type>
       Type MortalitySet<Type>::logSurvival(int a, int y, Type t0, Type t1){
	 // Loop over hazard intervals
	 Type logS0 = 0.0;	 
	 for(int i = 0; i < (int)activeHazard_breakpoints.size() - 1; ++i){
	   // Skip interval if it ends before time interval
	   if(activeHazard_breakpoints[i+1] <= t0)
	     continue;
	   // Break loop if interval starts after time interval
	   if(activeHazard_breakpoints[i] >= t1)
	     break;
	   // Otherwise, look at overlap between intervals	   
	   Type Astart = std::max(activeHazard_breakpoints[i],t0);
	   Type Aend =std::min(activeHazard_breakpoints[i+1],t1);
	   // Hazard is constant, so cumulative hazard is hazard times interval length
	   logS0 -= exp(logHazard_breakpoints(a,y,i)) * (Aend - Astart);
	 }
	 return logS0;
       }
       )

SOURCE(
       template<class Type>
       Type MortalitySet<Type>::logCIF(int fleet, int a, int y, Type t0, Type t1){
	 // Cumulative incidence between t0 and t1

	 // Loop over hazard intervals (start from probability one of surviving until t0
	 Type logS0 = 0.0;
	 Type vlogCIF = SAM_NegInf;
	 for(int i = 0; i < (int)activeHazard_breakpoints.size() - 1; ++i){
	   // Skip interval if it ends before time interval
	   if(activeHazard_breakpoints[i+1] <= t0)
	     continue;
	   // Break loop if interval starts after time interval
	   if(activeHazard_breakpoints[i] >= t1)
	     break;
	   
	   // Otherwise, look at overlap between intervals
	   // Full interval
	   Type Astart = std::max(activeHazard_breakpoints[i],t0);
	   Type Aend =std::min(activeHazard_breakpoints[i+1],t1);
	   int f0 = fleet;
	   int f1 = fleet+1;
	   if(fleet < 0){
	     f0 = 0;
	     f1 = logCIF_F_breakpoints.dim(2);
	   }
	   for(int f = f0; f < f1; ++f){	     
	     if(t0 <= activeHazard_breakpoints[i] &&
		t1 >= activeHazard_breakpoints[i+1]){ // Full interval (pre-calculated)	     
	       //vCIF += exp(logS0) * CIF_F_breakpoints(a,y,fleet,i);
	       vlogCIF = logspace_add_SAM(vlogCIF, logS0 + logCIF_F_breakpoints(a,y,fleet,i));
	     }else{
	       //vCIF += exp(logS0) * Hazard_F_breakpoints(a,y,f,i) / Hazard_breakpoints(a,y,i) * (1.0 - exp(-Hazard_breakpoints(a,y,i) * (Aend-Astart)));
	       Type tmp = logHazard_F_breakpoints(a,y,f,i) - logHazard_breakpoints(a,y,i) + logspace_sub_SAM(Type(0.0), -exp(logHazard_breakpoints(a,y,i)) * (Aend-Astart));
	       vlogCIF = logspace_add_SAM(vlogCIF, logS0 + tmp);
	     }
	   }	 
	   // Update survival before next interval: Hazard is constant, so cumulative hazard is hazard times interval length
	   logS0 -= exp(logHazard_breakpoints(a,y,i)) * (Aend - Astart);
	 }
	 return vlogCIF;
       }
       )


SOURCE(
       template<class Type>
       Type MortalitySet<Type>::partialLogCIF(int fleet, int a, int y, Type t0, Type t1){
	 // Probability of survival until t0 and dying from fishing in (t0,t1)
	 Type logS0 = this->logSurvival(a,y,(Type)0.0,t0);
	 Type vlogCIF = this->logCIF(fleet,a,y,t0,t1);
	 return logS0 + vlogCIF;
       }
       )



SOURCE(
       template<class Type>
       Type MortalitySet<Type>::logCIFr(int risk, int a, int y, Type t0, Type t1){
	 // Cumulative incidence between t0 and t1

	 // Loop over hazard intervals (start from probability one of surviving until t0
	 Type logS0 = 0.0;
	 Type vlogCIF = SAM_NegInf;
	 for(int i = 0; i < (int)activeHazard_breakpoints.size() - 1; ++i){
	   // Skip interval if it ends before time interval
	   if(activeHazard_breakpoints[i+1] <= t0)
	     continue;
	   // Break loop if interval starts after time interval
	   if(activeHazard_breakpoints[i] >= t1)
	     break;
	   
	   // Otherwise, look at overlap between intervals
	   // Full interval
	   Type Astart = std::max(activeHazard_breakpoints[i],t0);
	   Type Aend =std::min(activeHazard_breakpoints[i+1],t1);
	   int r0 = risk;
	   int r1 = risk+1;
	   if(risk < 0){
	     r0 = 0;
	     r1 = logCIF_M_breakpoints.dim(2);
	   }
	   for(int r = r0; r < r1; ++r){	     
	     if(t0 <= activeHazard_breakpoints[i] &&
		t1 >= activeHazard_breakpoints[i+1]){ // Full interval (pre-calculated)	     
	       //vCIF += exp(logS0) * CIF_M_breakpoints(a,y,r,i);
	       vlogCIF = logspace_add_SAM(logS0, logCIF_M_breakpoints(a,y,r,i));
	     }else{
	       //vCIF += exp(logS0) * Hazard_M_breakpoints(a,y,r,i) / Hazard_breakpoints(a,y,i) * (1.0 - exp(-Hazard_breakpoints(a,y,i) * (Aend-Astart)));
	       Type tmp = logHazard_M_breakpoints(a,y,r,i) - logHazard_breakpoints(a,y,i) + logspace_sub_SAM(Type(0.0), -exp(logHazard_breakpoints(a,y,i)) * (Aend-Astart));
	       vlogCIF = logspace_add_SAM(vlogCIF, logS0 + tmp);
	     }
	   }	 
	   // Update survival before next interval: Hazard is constant, so cumulative hazard is hazard times interval length
	   logS0 -= exp(logHazard_breakpoints(a,y,i)) * (Aend - Astart);
	 }
	 return vlogCIF;
       }
       )


SOURCE(
       template<class Type>
       Type MortalitySet<Type>::partialLogCIFr(int risk, int a, int y, Type t0, Type t1){
	 // Probability of survival until t0 and dying from fishing in (t0,t1)
	 Type logS0 = this->logSurvival(a,y,(Type)0.0,t0);
	 Type vlogCIF = this->logCIFr(risk,a,y,t0,t1);
	 return logS0 + vlogCIF;
       }
       )





SOURCE(
	 template<class Type>
	 MortalitySet<Type>::MortalitySet(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, array<Type>& logitFseason) :
	 logCumulativeHazard(conf.keyLogFsta.dim(1), dat.natMor.dim(0)),
	 logCumulativeHazard_F(conf.keyLogFsta.dim(1), dat.natMor.dim(0), conf.keyLogFsta.dim(0)),
	 logCumulativeHazard_M(conf.keyLogFsta.dim(1), dat.natMor.dim(0), dat.CompRisk.size()+1),
	 FullYear_logCumulativeIncidence_Fishing(conf.keyLogFsta.dim(1), dat.natMor.dim(0),conf.keyLogFsta.dim(0)),
	 FullYear_logCumulativeIncidence_Other(conf.keyLogFsta.dim(1), dat.natMor.dim(0),dat.CompRisk.size()+1),
	 FullYear_logSurvival(conf.keyLogFsta.dim(1), dat.natMor.dim(0)),
	 Effective_logF(conf.keyLogFsta.dim(1), dat.natMor.dim(0)),
	 Effective_logM(conf.keyLogFsta.dim(1), dat.natMor.dim(0)),

	 // totalF(totFFun(dat,conf,logF).matrix()),
	 // totalZseason(conf.seasonTimes.size()-1, conf.keyLogFsta.dim(1), dat.natMor.dim(0)),
	 // totalFseason(conf.seasonTimes.size()-1, conf.keyLogFsta.dim(1), dat.natMor.dim(0)),
	 logFleetSurvival_before(conf.keyLogFsta.dim(1), dat.natMor.dim(0), conf.keyLogFsta.dim(0)),
	 fleetLogCumulativeIncidence(conf.keyLogFsta.dim(1), dat.natMor.dim(0), conf.keyLogFsta.dim(0)),
	 otherLogCumulativeIncidence(conf.keyLogFsta.dim(1), dat.natMor.dim(0), dat.CompRisk.size()+1),
	 ssbLogSurvival_before(conf.keyLogFsta.dim(1), dat.natMor.dim(0)),
	 Fseason(logitFseason.dim(0)+1,logitFseason.dim(1),logitFseason.dim(2)+1),
	 activeHazard_breakpoints(),	  
	 activeHazard_season(),
	 activeHazard_F(),
	 activeHazardMap_risk(),
	 logHazard_breakpoints(),
	 logHazard_F_breakpoints(),
	 logHazard_M_breakpoints(),
	 logCIF_F_breakpoints(),
	 logCIF_M_breakpoints()

	 {
	   int nYear = dat.natMor.dim(0);
	   logCumulativeHazard.setConstant(SAM_NegInf);
	   logCumulativeHazard_F.setConstant(SAM_NegInf);
	   logCumulativeHazard_M.setConstant(SAM_NegInf);
	   FullYear_logCumulativeIncidence_Fishing.setConstant(SAM_NegInf);
	   FullYear_logCumulativeIncidence_Other.setConstant(SAM_NegInf);
	   FullYear_logSurvival.setConstant(SAM_NegInf);
	   Effective_logF.setConstant(SAM_NegInf);
	   Effective_logM.setConstant(SAM_NegInf);
	   
	   Fseason.setZero();
	   // totalZseason.setZero();
	   // totalFseason.setZero();
	   logFleetSurvival_before.setConstant(SAM_NegInf);
	   fleetLogCumulativeIncidence.setConstant(SAM_NegInf);
	   otherLogCumulativeIncidence.setConstant(SAM_NegInf);
	   ssbLogSurvival_before.setConstant(SAM_NegInf);

	   // activeHazard breakpoints
	   std::vector<Type> ahb_tmp(0);
	   // Natural mortality M1 is constant from start to end
	   ahb_tmp.push_back(0.0);
	   ahb_tmp.push_back(1.0);
	   // Risks
	   for(int r = 0; r < dat.CompRisk.size(); ++r){
	     vector<double> brk = dat.CompRisk(r).breakpoints;
	     for(int i = 0; i < brk.size(); ++i)
	       ahb_tmp.push_back(brk(i));
	   }
	   // Fleets
	   for(int f = 0; f < dat.sampleTimesStart.size(); ++f){
	     // only changes hazards if fleet type is 0
	     if(dat.fleetTypes(f) == 0){
	       ahb_tmp.push_back(dat.sampleTimesStart(f));
	       ahb_tmp.push_back(dat.sampleTimesEnd(f));
	     }
	   }
	   // Seasons
	   for(int s = 0; s < conf.seasonTimes.size(); ++s){
	     ahb_tmp.push_back(conf.seasonTimes(s));
	   }
	   // Sort (no parameters)
	   std::sort(ahb_tmp.begin(), ahb_tmp.end());
	   // Push result
	   activeHazard_breakpoints.push_back(ahb_tmp[0]);
	   for(int i = 1; i < (int)ahb_tmp.size(); ++i)
	     if(ahb_tmp[i-1] != ahb_tmp[i])
	       activeHazard_breakpoints.push_back(ahb_tmp[i]);

	   // Activity indicators
	   int nInterval = activeHazard_breakpoints.size()-1;
	   activeHazard_season = vector<int>(nInterval);
	   activeHazard_season.setConstant(-1);
	   // Loop over seasons
	   for(int s = 0; s < conf.seasonTimes.size()-1; ++s){
	     if(conf.isFishingSeason(s)){
	       // Loop over breakpoints
	       for(int i = 0; i < nInterval; ++i){
		 // If season starts before interval and ends after, it is active
		 if(conf.seasonTimes(s) <= activeHazard_breakpoints[i] &&
		    conf.seasonTimes(s+1) >= activeHazard_breakpoints[i+1]){
		   activeHazard_season(i) = s;
		 }
	       }
	     }
	   }
	   int nRisk = dat.CompRisk.size();
	   activeHazardMap_risk = matrix<int>(nInterval,nRisk);
	   if(nRisk > 0){
	     activeHazardMap_risk.setConstant(-1);
	     // Loop over risks (exclude M1)
	     for(int r = 0; r < nRisk; ++r){
	       int lastBrk = 0;
	       vector<double> brks = dat.CompRisk(r).breakpoints;
	       for(int i = 0; i < nInterval; ++i){
		 for(int q = lastBrk; q < brks.size()-1; ++q){ //Since brks are ordered, the hazard interval cannot be lower than the previous
		   double t0 = brks[q];
		   double t1 = brks[q+1];		   
		   if(t0 <= activeHazard_breakpoints[i] &&
		      t1 >= activeHazard_breakpoints[i+1]){
		     activeHazardMap_risk(i,r) = q;
		     break;
		   }
		 }
	       }
	     }
	   }
	   int nFleets = conf.keyLogFsta.dim(0);
	   int nAges = conf.keyLogFsta.dim(1);
	   
	   activeHazard_F = array<int>(nInterval,nAges,nFleets);
	   activeHazard_F.setZero();
	   // Loop over fleets
	   for(int f = 0; f < nFleets; ++f){
	     // Loop over ages
	     for(int a = 0; a < nAges; ++a){
	       int i = conf.keyLogFsta(f,a);
	       if(i > (-1)){	       
		 // Loop over breakpoints
		 for(int i = 0; i < nInterval; ++i){
		   // If it is an active fishing season, fishing starts before interval and ends after, it is active
		   if(activeHazard_season(i) > (-1) &&
		      dat.sampleTimesStart(f) <= activeHazard_breakpoints[i] &&
		      dat.sampleTimesEnd(f) >= activeHazard_breakpoints[i+1]){
		     activeHazard_F(i,a,f) = 1;
		   }
		 }
	       }
	     }
	   }
	   // Non-fishing hazard. Currently only one that is always active
	   // activeHazard_M = array<int>(nInterval,nAges,1);
	   // activeHazard_M.setConstant(1);

	   
	   logHazard_breakpoints = array<Type>(nAges,nYear,nInterval);
	   logHazard_F_breakpoints = array<Type>(nAges,nYear,nFleets,nInterval);
	   logHazard_M_breakpoints = array<Type>(nAges,nYear,nRisk+1,nInterval);
	   logCIF_F_breakpoints = array<Type>(nAges,nYear,nFleets,nInterval);
	   logCIF_M_breakpoints = array<Type>(nAges,nYear,nRisk+1,nInterval);
	   
	   for(int y = 0; y < nYear; ++y)
	     this->updateYear(dat,conf,par,logF,logitFseason,y);
	   return;
	 }
       )

SOURCE(
       template<class Type>
       void MortalitySet<Type>::updateYear(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, array<Type>& logitFseason, int y){
	 
	 int nFleet = conf.keyLogFsta.dim(0);
	 int nRisk = dat.CompRisk.size() + 1;
	 int nAge = conf.keyLogFsta.dim(1);
	 int nYear = dat.natMor.dim(0);
	 // int nSeason = conf.seasonTimes.size()-1;
	 if(y > nYear || y < 0)
	   Rf_error("MortalitySet.updateYear: Year not in range");

	 // Update Fseason
	 this->updateSeasons(dat, conf,logitFseason,y);
	 // Update hazards
	 this->updateHazards(dat,conf,par, logF,y);

	 // Update everything else
	 for(int a = 0; a < nAge; ++a){
	   Type vssb = SAM_NegInf; //dat.natMor(y,a) * dat.propM(y,a);
	   bool didit = false;
	   for(int r = 0; r < nRisk; ++r){
	     if(dat.propM(y,a) > 0){
	       didit = true;
	       vssb = logspace_add_SAM(vssb, logCumulativeHazard_M(a,y,r) + log(dat.propM(y,a)));
	     }
	     otherLogCumulativeIncidence(a,y,r) = this->logCIFr(r, a, y, (Type)0, (Type)1);
	   }
	   for(int f = 0; f < nFleet; ++f){
	     logFleetSurvival_before(a,y,f) = this->logSurvival(a, y, Type(0.0), (Type)dat.sampleTimesStart(f));
	     fleetLogCumulativeIncidence(a,y,f) = this->logCIF(f, a, y, (Type)0, (Type)1); //(Type)dat.sampleTimesStart(f), (Type)dat.sampleTimesEnd(f));
	     int i = conf.keyLogFsta(f,a);
	     if(i > (-1) && dat.propF(y,a,f) > 0){
	       didit = true;
	     //   vssb += exp(logF(i,y)) * dat.propF(y,a,f);
	       vssb = logspace_add_SAM(vssb, logCumulativeHazard_F(a,y,f) + log(dat.propF(y,a,f)));
	       // vssb += cumulativeHazard_F(a,y,f) * dat.propF(y,a,f);
	     }
	   }
	   int yx = std::min(y, (int)dat.recruitmentTimeOfYear.size()-1);
	   if(dat.recruitmentTimeOfYear(yx) >= 0 && dat.recruitmentTimeOfYear(yx) <= 1){
	     ssbLogSurvival_before(a,y) = this->logSurvival(a, y, Type(0.0), dat.recruitmentTimeOfYear(yx));
	   }else if(!didit){
	     ssbLogSurvival_before(a,y) = 0;
	   }else{
	     ssbLogSurvival_before(a,y) = -exp(vssb);
	   }
	 }
	 return;
       
	 
	 // // Update everything else	 
	 // for(int a = 0; a < nAge; ++a){
	 //   Type newTZ = dat.natMor(y,a);
	 //   Type newTF = 0.0;
	 //   for(int f = 0; f < nFleet; ++f){
	 //     int i = conf.keyLogFsta(f,a);
	 //     if(i > (-1)){
	 //       Type ff = exp(logF(i,y)) * (dat.sampleTimesEnd(f) - dat.sampleTimesStart(f));
	 //       newTZ += ff;
	 //       newTF += ff;
	 //     }
	 //   }
	 //   totalZ(a,y) = newTZ;
	 //   totalF(a,y) = newTF;
	 //   if(dat.natMor(y,a) > 0){
	 //     Type v = (1.0 - exp(-totalZ(a,y))) / totalZ(a,y);
	 //     otherCumulativeIncidence(a,y,0) = dat.natMor(y,a) * v;
	 //     Type vssb = dat.natMor(y,a) * dat.propM(y,a);
	 //     // Calculate totalZ per season
	 //     for(int s = 0; s < nSeason; ++s){
	 //       totalZseason(s,a,y) = dat.natMor(y,a) * (conf.seasonTimes(s+1)-conf.seasonTimes(s));
	 //       totalFseason(s,a,y) = 0.0;
	 //     }
	 //     for(int f = 0; f < nFleet; ++f){
	 //       int i = conf.keyLogFsta(f,a);
	 //       int j = conf.keyLogFseason(f,a);	
	 //       for(int s = 0; s < nSeason; ++s){		     
	 // 	 if(conf.isFishingSeason(s) && i > (-1)){
	 // 	   // If season does not end before fleet starts
	 // 	   // and season does not start before fleet ends
	 // 	   if(!(conf.seasonTimes(s+1) <= dat.sampleTimesStart(f)) &&
	 // 	      !(conf.seasonTimes(s) >= dat.sampleTimesEnd(f))){
	 // 	     // Start time is maximum of season and fleet start (no parameters)
	 // 	     Type As = std::max((Type)conf.seasonTimes(s),(Type)dat.sampleTimesStart(f));
	 // 	     // End time is minimum of season and fleet end (no parameters)
	 // 	     Type Ae = std::min((Type)conf.seasonTimes(s+1),(Type)dat.sampleTimesEnd(f));
	 // 	     Type Fs = exp(logF(i,y) + Fseason(s,y,j+1)) * (Ae - As);
	 // 	     //totalZseason(s,a,y) += exp(logF(i,y)) * Fseason(s) * (conf.seasonTimes(s+1)-conf.seasonTimes(s));
	 // 	     totalFseason(s,a,y) += Fs;
	 // 	     totalZseason(s,a,y) += Fs;
	 // 	   }
	 // 	 }
	 //       }
	 //     }
	 //     // Done updating hazards
	     
	 //     // Calculate survival before fleet starts and CIF
	 //     for(int f = 0; f < nFleet; ++f){
	 //       logFleetSurvival_before(a,y,f) = this->logSurvival(dat, conf, par, logF, a, y, (Type)dat.sampleTimesStart(f));
	 //       fleetCumulativeIncidence(a,y,f) = this->CIF(dat, conf, par, logF, f, a, y, (Type)dat.sampleTimesStart(f), (Type)dat.sampleTimesEnd(f));
	 //       int i = conf.keyLogFsta(f,a);
	 //       if(i > (-1))
	 // 	 vssb += exp(logF(i,y)) * dat.propF(y,a,f);
	 //     }
	 //     ssbSurvival_before(a,y) = exp(-vssb);
	 //   }
	 // }
	 //return;
       }
       );
       
SAM_SPECIALIZATION(struct MortalitySet<double>);
SAM_SPECIALIZATION(struct MortalitySet<TMBad::ad_aug>);
