SAM_DEPENDS(define)


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
  res.setConstant(R_NegInf);
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
  matrix<Type> cumulativeHazard;		// age x year
  // matrix<Type> totalF;
  // array<Type> totalZseason;	// season x age x year
  // array<Type> totalFseason;
  // matrix<Type> totalFCI;
  array<Type> logFleetSurvival_before; // age x year x fleet (including surveys)
  array<Type> fleetCumulativeIncidence; // age x year x fleet (including surveys)
  array<Type> otherCumulativeIncidence; // age x year x causes
  matrix<Type> ssbSurvival_before;	     // age x year
  array<Type> Fseason;			     // seasons x year x (processes+1)

  std::vector<Type> activeHazard_breakpoints; // Vector of time points
  vector<int> activeHazard_season; // Key to the fishing season (-1 if none)
  array<int> activeHazard_F;	   // Indicator, 1 if fleet is active (brkpoints x ages x fleets)
  array<int> activeHazard_M;	// Indicator, 1 if non-fishing hazard is active. (brkpoints x ages x causes)

  array<Type> Hazard_breakpoints;	// age x year x activeHazard_breakpoints
  array<Type> Hazard_F_breakpoints;	// age x year x fleet x activeHazard_breakpoints
  array<Type> CIF_F_breakpoints;	// age x year x fleet x activeHazard_breakpoints
  array<Type> CIF_M_breakpoints;	// age x year x causes x activeHazard_breakpoints
  
  inline MortalitySet():
    cumulativeHazard(),
    // totalF(),
    // totalZseason(),
    // totalFseason(),
    logFleetSurvival_before(),
    fleetCumulativeIncidence(),
    otherCumulativeIncidence(),
    ssbSurvival_before(),
    Fseason(),
    activeHazard_breakpoints(),
    activeHazard_season(),
    activeHazard_F(),
    activeHazard_M(),
    Hazard_breakpoints(),
    Hazard_F_breakpoints(),
    CIF_F_breakpoints(),
    CIF_M_breakpoints()
  {}
  
  template<class T>
  inline MortalitySet(const MortalitySet<T> x) : cumulativeHazard(x.cumulativeHazard),
						 // totalF(x.totalF),
						 // totalZseason(x.totalZseason),
						 // totalFseason(x.totalFseason),
						 logFleetSurvival_before(x.logFleetSurvival_before),
						 fleetCumulativeIncidence(x.fleetCumulativeIncidence),
						 otherCumulativeIncidence(x.otherCumulativeIncidence),
						 ssbSurvival_before(x.ssbSurvival_before),
						 Fseason(x.Fseason),
						 activeHazard_breakpoints(x.activeHazard_breakpoints),
						 activeHazard_season(x.activeHazard_season),
						 activeHazard_F(x.activeHazard_F),
						 activeHazard_M(x.activeHazard_M),
						 Hazard_breakpoints(x.Hazard_breakpoints),
						 Hazard_F_breakpoints(x.Hazard_breakpoints),
						 CIF_F_breakpoints(x.CIF_F_breakpoints),
						 CIF_M_breakpoints(x.CIF_M_breakpoints)
						 
  {}
  
  MortalitySet(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, array<Type>& logitFseason);

  void updateSeasons(confSet& conf, array<Type>& logitFseason, int y);
  void updateHazards(dataSet<Type>& dat, confSet& conf, array<Type>& logF, int y);
  // For simulation based forecast
  void updateYear(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, array<Type>& logitFseason, int y);

  Type logSurvival(int a, int y, Type t0, Type t1);
  Type CIF(int fleet, int a, int y, Type t0, Type t1);
  Type partialCIF(int fleet, int a, int y, Type t0, Type t1);
  
});

SOURCE(
       template<class Type>
       void MortalitySet<Type>::updateSeasons(confSet& conf, array<Type>& logitFseason, int y){
	 //Fseason.setZero();
	 Type NFseason = conf.isFishingSeason.sum();
	 // If no season info, use constant (i.e. 0)
	 for(int s = 0; s < Fseason.dim(0); ++s){
	   Fseason(s,y,0) = 0.0; //(Type)log(conf.isFishingSeason(s)) - log(NFseason);
	   //Fseason(s,y,0) -= log(conf.seasonTimes(s+1)-conf.seasonTimes(s));    
	 }
	 // Loop over processes
	 for(int p = 0; p < logitFseason.dim(2); ++p){
	   vector<Type> x0(logitFseason.dim(0));
	   for(int qq = 0; qq < x0.size(); ++qq) x0(qq) = logitFseason(qq,y,p);
	   vector<Type> tfs = toLogSeasonEffect(x0, conf.isFishingSeason, conf.seasonTimes);
	   for(int s = 0; s < Fseason.dim(0); ++s)
	     Fseason(s,y,p+1) = tfs(s);
	 }
	 return;
       }
       )

SOURCE(
template<class Type>
void MortalitySet<Type>::updateHazards(dataSet<Type>& dat, confSet& conf, array<Type>& logF,int y){
  // Hazard and cumulative incidence per interval
  int nIntervals = Hazard_breakpoints.dim(2);
  int nAges = Hazard_breakpoints.dim(0);
  int nFleet = conf.keyLogFsta.dim(0);	 

  vector<Type> tmpCumHaz(nAges);
  tmpCumHaz.setZero();
  for(int i = 0; i < nIntervals; ++i){
    int s = activeHazard_season(i);
    Type dt = activeHazard_breakpoints[i+1] - activeHazard_breakpoints[i];
    for(int a = 0; a < nAges; ++a){
      Hazard_breakpoints(a,y,i) = 0.0;      
      // Natural hazards (only natMor)
      if(activeHazard_M(i,a,0))
        Hazard_breakpoints(a,y,i) += dat.natMor(y,a);
      // Fleets
      for(int f = 0; f < nFleet; ++f){
	Hazard_F_breakpoints(a,y,f,i) = 0.0;
	int indx = conf.keyLogFsta(f,a);
	if(indx > (-1) && activeHazard_F(i,a,f)){
	  int j = conf.keyLogFseason(f,a);	  
	  Type logFs = logF(indx,y);
	  if(s > (-1))
	    logFs += Fseason(s,y,j+1);
	  Hazard_breakpoints(a,y,i) += exp(logFs);
	  Hazard_F_breakpoints(a,y,f,i) += exp(logFs);
	}
      }
      // Calculate CIF M
      CIF_M_breakpoints(a,y,0,i) = 0.0;
      if(activeHazard_M(i,a,0))
	CIF_M_breakpoints(a,y,0,i) = dat.natMor(y,a) / Hazard_breakpoints(a,y,i) * (1.0 - exp(-Hazard_breakpoints(a,y,i) * dt));
      // Calculate CIF F
      for(int f = 0; f < nFleet; ++f){
	CIF_F_breakpoints(a,y,f,i) = 0.0;
	int indx = conf.keyLogFsta(f,a);
	if(indx > (-1) && activeHazard_F(i,a,f)){
	  int j = conf.keyLogFseason(f,a);
	  Type logFs = logF(indx,y);
	  if(s > (-1))
	    logFs += Fseason(s,y,j+1);
	  CIF_F_breakpoints(a,y,f,i) = exp(logFs) / Hazard_breakpoints(a,y,i) * (1.0 - exp(-Hazard_breakpoints(a,y,i) * dt));
	}
      }
      tmpCumHaz(a) += Hazard_breakpoints(a,y,i) * dt;
    }
  }
  for(int a = 0; a < nAges; ++a){
    cumulativeHazard(a,y) = tmpCumHaz(a);
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
	   logS0 -= Hazard_breakpoints(a,y,i) * (Aend - Astart);
	 }
	 return logS0;
	 // // Survival until t
	 // int nFleet = conf.keyLogFsta.dim(0);
	 // // int nAge = conf.keyLogFsta.dim(1);
	 // // int nYear = dat.natMor.dim(0);
	 // int nSeason = conf.seasonTimes.size()-1;
	 // Type logS0 = 0.0;
	 // for(int s = 0; s < nSeason; ++s){	   
	 //   if(conf.seasonTimes(s+1) < t){
	 //     // If season ends before t0, add full period hazard
	 //     logS0 -= totalZseason(s,a,y);
	 //   }else if(conf.seasonTimes(s) > t){
	 //     // If season starts after t0, do nothing
	 //   }else{
	 //     // If season ends after t0 and starts before t0
	 //     // Natural mortality for period
	 //     //logS0 -= dat.natMor(y,a) * (std::min((Type)t,(Type)conf.seasonTimes(s+1)) - conf.seasonTimes(s));
	 //     logS0 -= dat.natMor(y,a) * (t - conf.seasonTimes(s));
	 //     if(conf.isFishingSeason(s)){
	 //       // Loop over fleets
	 //       for(int f = 0; f < nFleet; ++f){
	 // 	 int i = conf.keyLogFsta(f,a);
	 // 	 if(i > (-1)){
	 // 	   int j = conf.keyLogFseason(f,a);
	 // 	   // If fleet starts before t0
	 // 	   if((Type)dat.sampleTimesStart(f) < t){
	 // 	     // Start time is maximum of season and fleet start (no parameters)
	 // 	     Type As = std::max((Type)conf.seasonTimes(s),(Type)dat.sampleTimesStart(f));
	 // 	     // End time is minimum of fleet end and t0 (no parameters)
	 // 	     Type Ae = std::min((Type)t,(Type)dat.sampleTimesEnd(f));
	 // 	     logS0 -= exp(logF(i,y) + Fseason(s,y,j+1)) * (Ae - As);

	 // 	     /*
	 // 	   // If fleet overlaps season
	 // 	   // If season does not end before fleet starts
	 // 	   // and season does not start before fleet ends
	 // 	   if(!(conf.seasonTimes(s+1) <= (Type)dat.sampleTimesStart(f)) &&
	 // 	      !(conf.seasonTimes(s) >= dat.sampleTimesEnd(f))){
	 // 	     // Fleet/Season Start time is maximum of season and fleet start (no parameters)
	 // 	     Type As = std::max((Type)conf.seasonTimes(s),(Type)dat.sampleTimesStart(f));
	 // 	     // Fleet/Season Start time is maximum of season and fleet start (no parameters)
	 // 	     Type Ae = std::min((Type)conf.seasonTimes(s),(Type)dat.sampleTimesEnd(f));
	 // 	     if(As < t){ // If fleet starts before t
	 // 	       // End time is minimum of fleet/season end and t0 (no parameters)
	 // 	       Ae = std::min((Type)t,(Type)Ae);
	 // 	       logS0 -= exp(logF(i,y) + Fseason(s,y,j+1)) * (Ae - As);
	 // 	     }
	 // 	     */
	 // 	   }
	 // 	 }
	 //       }
	 //     }
	 //   }
	 // }
	 // return logS0;
       }
       )

SOURCE(
       template<class Type>
       Type MortalitySet<Type>::CIF(int fleet, int a, int y, Type t0, Type t1){
	 // Cumulative incidence between t0 and t1

	 // Loop over hazard intervals (start from probability one of surviving until t0
	 Type logS0 = 0.0;
	 Type vCIF = 0.0;
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
	     f1 = CIF_F_breakpoints.dim(2);
	   }
	   for(int f = f0; f < f1; ++f){	     
	     if(t0 <= activeHazard_breakpoints[i] &&
		t1 >= activeHazard_breakpoints[i+1]){ // Full interval (pre-calculated)	     
	       vCIF += exp(logS0) * CIF_F_breakpoints(a,y,fleet,i);
	     }else{
	       vCIF += exp(logS0) * Hazard_F_breakpoints(a,y,f,i) / Hazard_breakpoints(a,y,i) * (1.0 - exp(-Hazard_breakpoints(a,y,i) * (Aend-Astart)));
	     }
	   }	 
	   // Update survival before next interval: Hazard is constant, so cumulative hazard is hazard times interval length
	   logS0 -= Hazard_breakpoints(a,y,i) * (Aend - Astart);
	 }
	 return vCIF;


	 
       // 	 int nFleet = conf.keyLogFsta.dim(0);
       // 	 // int nAge = conf.keyLogFsta.dim(1);
       // 	 // int nYear = dat.natMor.dim(0);
       // 	 int nSeason = conf.seasonTimes.size()-1;
       // 	 Type vCIF = 0.0;
       // 	 Type logPS = 0.0;
       // 	 int i = conf.keyLogFsta(fleet,a);
       // 	 if(i < 0)
       // 	   return 0.0;
       // 	 int j = conf.keyLogFseason(fleet,a);
       // 	 Type STS = std::max(t0,dat.sampleTimesStart(fleet));
       // 	 Type STE = std::min(t1,dat.sampleTimesEnd(fleet));
       // 	 for(int s = 0; s < nSeason; ++s){
       // 	   // Survival * F / Z * (1-exp(-Z*dt))
       // 	   // If season does not end before fleet starts
       // 	   // and season does not start before fleet ends
       // 	   if(!(conf.seasonTimes(s+1) <= STS) &&
       // 	      !(conf.seasonTimes(s) >= STE)){
       // 	     // Start time is maximum of season and fleet start (no parameters)
       // 	     Type As = std::max((Type)conf.seasonTimes(s),(Type)STS);
       // 	     // End time is minimum of seson and fleet end (no parameters)
       // 	     Type Ae = std::min((Type)conf.seasonTimes(s+1),(Type)STE);
       // 	     Type logFs = logF(i,y) + Fseason(s,y,j+1);
       // 	     if((Type)conf.seasonTimes(s) > (Type)STS &&
       // 		(Type)STE > (Type)conf.seasonTimes(s+1)){
       // 	       // Full season (totalZseason is already scaled in time
       // 	       if(conf.isFishingSeason(s) && Ae > As)
       // 		 vCIF += exp(logPS + logFs - log(totalZseason(s,a,y))) * (1.0 - exp(-totalZseason(s,a,y))); // * (Ae-As)));
       // 	       logPS -= totalZseason(s,a,y);// * (Ae-As);
       // 	     }else{
       // 	       // Find Z for this part of season
       // 	       // NOTE: A loop is needed to account for fleets with zero F in part of season!
       // 	       Type dt = Ae - As;
       // 	       Type tmpZ = dat.natMor(y,a) * dt;
       // 	       if(conf.isFishingSeason(s)){
       // 		 for(int ff = 0; ff < nFleet; ++ff){
       // 		   int ii = conf.keyLogFsta(ff,a);
       // 		   Type STSs = std::max(t0,dat.sampleTimesStart(ff));
       // 		   Type STEs = std::min(t1,dat.sampleTimesEnd(ff));		     
       // 		   if(ii > (-1) &&
       // 		      (!(conf.seasonTimes(s+1) <= STS) &&
       // 		       !(conf.seasonTimes(s) >= STE))){
       // 		     int jj = conf.keyLogFseason(ff,a);
       // 		     // Start time is maximum of season and fleet start (no parameters)
       // 		     Type Ass = std::max(As,STSs);
       // 		     // End time is minimum of season and fleet end (no parameters)
       // 		     Type Aee = std::min(Ae,STEs);
       // 		     Type Fs = exp(logF(ii,y) + Fseason(s,y,jj+1)) * (Aee - Ass);
       // 		     if(conf.isFishingSeason(s))
       // 		       tmpZ += Fs;
       // 		   }
       // 		 }
       // 		 //if(dat.natMor(y,a) > 0)
       // 		 vCIF += exp(logPS + logFs - log(tmpZ)) * (1.0 - exp(-tmpZ * (Ae-As)));
       // 	       }
       // 	       logPS -= tmpZ;	       
       // 	     }
       // 	     }
       // 	 }
       // 	 return vCIF;
       }
       )


SOURCE(
       template<class Type>
       Type MortalitySet<Type>::partialCIF(int fleet, int a, int y, Type t0, Type t1){
	 // Probability of survival until t0 and dying from fishing in (t0,t1)
	 // int nFleet = conf.keyLogFsta.dim(0);
	 // int nAge = conf.keyLogFsta.dim(1);
	 // int nYear = dat.natMor.dim(0);
	 // int nSeason = conf.seasonTimes.size()-1;
	 Type logS0 = this->logSurvival(a,y,(Type)0.0,t0);
	 Type vCIF = this->CIF(fleet,a,y,t0,t1);
	 return exp(logS0) * vCIF;
       }
       )


SOURCE(
	 template<class Type>
	 MortalitySet<Type>::MortalitySet(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, array<Type>& logitFseason) :
	 cumulativeHazard(conf.keyLogFsta.dim(1), dat.natMor.dim(0)),
	 // totalF(totFFun(dat,conf,logF).matrix()),
	 // totalZseason(conf.seasonTimes.size()-1, conf.keyLogFsta.dim(1), dat.natMor.dim(0)),
	 // totalFseason(conf.seasonTimes.size()-1, conf.keyLogFsta.dim(1), dat.natMor.dim(0)),
	 logFleetSurvival_before(conf.keyLogFsta.dim(1), dat.natMor.dim(0), conf.keyLogFsta.dim(0)),
	 fleetCumulativeIncidence(conf.keyLogFsta.dim(1), dat.natMor.dim(0), conf.keyLogFsta.dim(0)),
	 otherCumulativeIncidence(conf.keyLogFsta.dim(1), dat.natMor.dim(0), 1),
	 ssbSurvival_before(conf.keyLogFsta.dim(1), dat.natMor.dim(0)),
	 Fseason(logitFseason.dim(0)+1,logitFseason.dim(1),logitFseason.dim(2)+1),
	 activeHazard_breakpoints(),	  
	 activeHazard_season(),
	 activeHazard_F(),
	 activeHazard_M(),
	 Hazard_breakpoints(),
	 Hazard_F_breakpoints(),
	 CIF_F_breakpoints(),
	 CIF_M_breakpoints()

	 {
	   int nYear = dat.natMor.dim(0);
	   cumulativeHazard.setZero();
	   Fseason.setZero();
	   // totalZseason.setZero();
	   // totalFseason.setZero();
	   logFleetSurvival_before.setZero();
	   fleetCumulativeIncidence.setZero();
	   otherCumulativeIncidence.setZero();
	   ssbSurvival_before.setZero();

	   // activeHazard breakpoints
	   std::vector<Type> ahb_tmp(0);
	   // Natural mortality is constant from start to end
	   ahb_tmp.push_back(0.0);
	   ahb_tmp.push_back(1.0);
	   // Fleets
	   for(int f = 0; f < dat.sampleTimesStart.size(); ++f){
	     ahb_tmp.push_back(dat.sampleTimesStart(f));
	     ahb_tmp.push_back(dat.sampleTimesEnd(f));
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
	   activeHazard_M = array<int>(nInterval,nAges,1);
	   activeHazard_M.setConstant(1);

	   
	   Hazard_breakpoints = array<Type>(nAges,nYear,nInterval);
	   Hazard_F_breakpoints = array<Type>(nAges,nYear,nFleets,nInterval);
	   CIF_F_breakpoints = array<Type>(nAges,nYear,nFleets,nInterval);
	   CIF_M_breakpoints = array<Type>(nAges,nYear,1,nInterval);
	   
	   for(int y = 0; y < nYear; ++y)
	     this->updateYear(dat,conf,par,logF,logitFseason,y);
	   return;
	 }
       )

SOURCE(
       template<class Type>
       void MortalitySet<Type>::updateYear(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, array<Type>& logitFseason, int y){
	 
	 int nFleet = conf.keyLogFsta.dim(0);
	 int nAge = conf.keyLogFsta.dim(1);
	 int nYear = dat.natMor.dim(0);
	 int nSeason = conf.seasonTimes.size()-1;
	 if(y > nYear || y < 0)
	   Rf_error("MortalitySet.updateYear: Year not in range");

	 // Update Fseason
	 this->updateSeasons(conf,logitFseason,y);
	 // Update hazards
	 this->updateHazards(dat,conf,logF,y);

	 // Update everything else
	 for(int a = 0; a < nAge; ++a){
	   Type vssb = dat.natMor(y,a) * dat.propM(y,a);
	   for(int f = 0; f < nFleet; ++f){
	     logFleetSurvival_before(a,y,f) = this->logSurvival(a, y, Type(0.0), (Type)dat.sampleTimesStart(f));
	     fleetCumulativeIncidence(a,y,f) = this->CIF(f, a, y, (Type)dat.sampleTimesStart(f), (Type)dat.sampleTimesEnd(f));
	     int i = conf.keyLogFsta(f,a);
	     if(i > (-1))
	       vssb += exp(logF(i,y)) * dat.propF(y,a,f);
	   }
	   ssbSurvival_before(a,y) = exp(-vssb);
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
