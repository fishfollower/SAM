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
      res(i) = r(indx++) / Type(seasonTimes(i+1) - seasonTimes(i));
    }
  return res;
  })


template <class Type>
array<Type> totFFun(dataSet<Type>& dat, confSet &conf, array<Type> &logF, int Ftype DEFARG(= 0)) SOURCE({
  int noFleets=conf.keyLogFsta.dim[0];
  int stateDimN=conf.keyLogFsta.dim[1];
  int timeSteps=logF.dim[1];
  if(Ftype > 0)
    timeSteps = dat.landFrac.dim[0];
  array<Type> totF(stateDimN,timeSteps);
  totF.setZero();  
  for(int i=0;i<timeSteps;i++){ 
    for(int j=0; j<stateDimN; ++j){
      for(int f=0; f<noFleets; ++f){
        if(conf.keyLogFsta(f,j)>(-1)){
	  Type Fval = exp(logF(conf.keyLogFsta(f,j),i)) * (dat.sampleTimesEnd(f) - dat.sampleTimesStart(f));
	  if(Ftype == 1){
	    Fval *= dat.landFrac(i,j,f);
	  }else if(Ftype == 2){
	    Fval *= (1.0 - dat.landFrac(i,j,f));
	  }
          totF(j,i) += Fval;
        }
      }
    }
  }
  return totF;
  })

SAM_SPECIALIZATION(array<double> totFFun(dataSet<double>&, confSet&, array<double>&, int));
SAM_SPECIALIZATION(array<TMBad::ad_aug> totFFun(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, int));



HEADER(
template<class Type>
struct MortalitySet {
  // TODO: Switch to log-scale calculations?
  matrix<Type> totalZ;		// age x year
  array<Type> totalZseason;	// season x age x year
  // matrix<Type> totalFCI;
  array<Type> logFleetSurvival_before; // age x year x fleet (including surveys)
  array<Type> fleetCumulativeIncidence; // age x year x fleet (including surveys)
  array<Type> otherCumulativeIncidence; // age x year x causes
  matrix<Type> ssbSurvival_before;	     // age x year

  inline MortalitySet():
    totalZ(),
    totalZseason(),
    logFleetSurvival_before(),
    fleetCumulativeIncidence(),
    otherCumulativeIncidence(),
    ssbSurvival_before() {}
  // : totalZ(),
  // 		   fleetSurvival_before(),
  // 		   fleetCumulativeIncidence(),
  // 		   otherCumulativeIncidence(),
  // 		   ssbSurvival_before() {}
  
  template<class T>
  inline MortalitySet(const MortalitySet<T> x) : totalZ(x.totalZ),
						 totalZseason(x.totalZseason),
					  logFleetSurvival_before(x.logFleetSurvival_before),
					  fleetCumulativeIncidence(x.fleetCumulativeIncidence),
					  otherCumulativeIncidence(x.otherCumulativeIncidence),
					  ssbSurvival_before(x.ssbSurvival_before) {}
  
  MortalitySet(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, array<Type>& logitFseason);

  // For simulation based forecast
  void updateYear(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, array<Type>& logitFseason, int y);
});

SOURCE(
	 template<class Type>
	 MortalitySet<Type>::MortalitySet(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, array<Type>& logitFseason) :
	 totalZ(totFFun(dat,conf,logF).matrix()),
	 totalZseason(conf.seasonTimes.size()-1, conf.keyLogFsta.dim(1), dat.natMor.dim(0)),
	 logFleetSurvival_before(conf.keyLogFsta.dim(1), dat.natMor.dim(0), conf.keyLogFsta.dim(0)),
	 fleetCumulativeIncidence(conf.keyLogFsta.dim(1), dat.natMor.dim(0), conf.keyLogFsta.dim(0)),
	 otherCumulativeIncidence(conf.keyLogFsta.dim(1), dat.natMor.dim(0), 1),
	 ssbSurvival_before(conf.keyLogFsta.dim(1), dat.natMor.dim(0))
	 {
	   // Only constant mortality over entire year (TODO: implement other options)
	   // int nFleet = conf.keyLogFsta.dim(0);
	   // int nAge = conf.keyLogFsta.dim(1);
	   int nYear = dat.natMor.dim(0);
	   // int nSeason = conf.seasonTimes.size()-1;
	   // totalFCI = matrix(nAge,nYear);
	   // totalFCI.setZero();
	   totalZseason.setZero();
	   logFleetSurvival_before.setZero();
	   fleetCumulativeIncidence.setZero();
	   otherCumulativeIncidence.setZero();
	   ssbSurvival_before.setZero();

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
	 for(int a = 0; a < nAge; ++a){
	   Type newTZ = dat.natMor(y,a);
	   for(int f = 0; f < nFleet; ++f){
	     int i = conf.keyLogFsta(f,a);
	     if(i > (-1))
	       newTZ += exp(logF(i,y)) * (dat.sampleTimesEnd(f) - dat.sampleTimesStart(f));
	   }
	   totalZ(a,y) = newTZ;
	   if(totalZ(a,y) > 0){
	     Type v = (1.0 - exp(-totalZ(a,y))) / totalZ(a,y);
	     otherCumulativeIncidence(a,y,0) = dat.natMor(y,a) * v;
	     Type vssb = dat.natMor(y,a) * dat.propM(y,a);
	     // Calculate totalZ per season
	     for(int s = 0; s < nSeason; ++s){
	       totalZseason(s,a,y) = dat.natMor(y,a) * (conf.seasonTimes(s+1)-conf.seasonTimes(s));
	     }
	     for(int f = 0; f < nFleet; ++f){
	       int i = conf.keyLogFsta(f,a);
	       int j = conf.keyLogFseason(f,a);
	       // TODO: Fseson is probably calculated more times than needed
	       vector<Type> Fseason(totalZseason.dim(0)); Fseason.setZero();
	       if(j > (-1)){
		 vector<Type> x0(logitFseason.dim(0));
		 for(int qq = 0; qq < x0.size(); ++qq) x0(qq) = logitFseason(qq,y,j);
		 Fseason = toLogSeasonEffect(x0, conf.isFishingSeason, conf.seasonTimes);
	       }else{
		 Type NFseason = conf.isFishingSeason.sum();
		 for(int qq = 0; qq < Fseason.size(); ++qq){
		   Fseason(qq) = (Type)log(conf.isFishingSeason(qq)) - log(NFseason);
		   Fseason(qq) -= log(conf.seasonTimes(qq+1)-conf.seasonTimes(qq));
		 }
	       }
	       for(int s = 0; s < nSeason; ++s){		     
		 if(conf.isFishingSeason(s) && i > (-1)){
		   // If season does not end before fleet starts
		   // and season does not start before fleet ends
		   if(!(conf.seasonTimes(s+1) <= dat.sampleTimesStart(f)) &&
		      !(conf.seasonTimes(s) >= dat.sampleTimesEnd(f))){
		     // Start time is maximum of season and fleet start (no parameters)
		     Type As = std::max((Type)conf.seasonTimes(s),(Type)dat.sampleTimesStart(f));
		     // End time is minimum of seson and fleet end (no parameters)
		     Type Ae = std::min((Type)conf.seasonTimes(s+1),(Type)dat.sampleTimesEnd(f));
		     Type Fs = exp(logF(i,y) + Fseason(s));
		     //totalZseason(s,a,y) += exp(logF(i,y)) * Fseason(s) * (conf.seasonTimes(s+1)-conf.seasonTimes(s));
		     totalZseason(s,a,y) += Fs * (Ae - As);
		   }
		 }
	       }
	     }
	     // Calculate survival before fleet starts
	     for(int f = 0; f < nFleet; ++f){
	       logFleetSurvival_before(a,y,f) = 0.0;
	       fleetCumulativeIncidence(a,y,f) = 0.0;
	       //fleetSurvival_before(a,y,f) = exp(-totalZ(a,y) * dat.sampleTimes(f));
	       for(int s = 0; s < nSeason; ++s){
		 if(conf.isFishingSeason(s)){
		   // // If season starts before fleet
		   if(conf.seasonTimes(s) <= dat.sampleTimesStart(f)){
		     // If season ends before fleet starts
		     if(conf.seasonTimes(s+1) <= dat.sampleTimesStart(f)){
		       logFleetSurvival_before(a,y,f) -= totalZseason(s,a,y) * (conf.seasonTimes(s+1)-conf.seasonTimes(s));
		       // Else if season ends after fleet
		     }else{
		       logFleetSurvival_before(a,y,f) -= totalZseason(s,a,y) * (dat.sampleTimesStart(f)-conf.seasonTimes(s));
		     }
		   } // Do nothing if season starts after fleet		    
		 }
	       }
	       // Fleet cumulative incidence during period
	       int i = conf.keyLogFsta(f,a);
	       if(i > (-1)){	// Has catch
		 int j = conf.keyLogFseason(f,a);	
		 vector<Type> Fseason(totalZseason.dim(0)); Fseason.setZero();
		 if(j > (-1)){
		   vector<Type> x0(logitFseason.dim(0));
		   for(int qq = 0; qq < x0.size(); ++qq) x0(qq) = logitFseason(qq,y,j);
		   Fseason = toLogSeasonEffect(x0, conf.isFishingSeason, conf.seasonTimes);
		 }else{
		   Type NFseason = conf.isFishingSeason.sum();
		   for(int qq = 0; qq < Fseason.size(); ++qq){
		     Fseason(qq) = (Type)log(conf.isFishingSeason(qq)) - log(NFseason);
		     Fseason(qq) -= log(conf.seasonTimes(qq+1)-conf.seasonTimes(qq));
		   }
		 }
		 // fleetCumulativeIncidence(a,y,f) = exp(logF(i,y)) * v;
		 Type logPS = 0.0;
		 for(int s = 0; s < nSeason; ++s){
		   // Survival * F / Z * (1-exp(-Z*dt))
		   if(conf.isFishingSeason(s)){
		     // If season does not end before fleet starts
		     // and season does not start before fleet ends
		     if(!(conf.seasonTimes(s+1) <= dat.sampleTimesStart(f)) &&
			!(conf.seasonTimes(s) >= dat.sampleTimesEnd(f))){
		       // Start time is maximum of season and fleet start (no parameters)
		       Type As = std::max((Type)conf.seasonTimes(s),(Type)dat.sampleTimesStart(f));
		       // End time is minimum of seson and fleet end (no parameters)
		       Type Ae = std::min((Type)conf.seasonTimes(s+1),(Type)dat.sampleTimesEnd(f));
		       Type logFs = logF(i,y) + Fseason(s);
		       fleetCumulativeIncidence(a,y,f) += exp(logPS + logFs - log(totalZseason(s,a,y))) * (1.0 - exp(-totalZseason(s,a,y) * (Ae-As)));
		       // Update survival
		       logPS -= totalZseason(s,a,y) * (Ae-As);
		     }
		   }
		     }		    
		     vssb += exp(logF(i,y)) * dat.propF(y,a,f);
		   }// Otherwise, stay 0
		 }
		 ssbSurvival_before(a,y) = exp(-vssb);
	     }
	   }
	   return;
	 }
       );
  
SAM_SPECIALIZATION(struct MortalitySet<double>);
SAM_SPECIALIZATION(struct MortalitySet<TMBad::ad_aug>);
