SAM_DEPENDS(define)

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
	  Type Fval = exp(logF(conf.keyLogFsta(f,j),i));
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
  // matrix<Type> totalFCI;
  array<Type> fleetSurvival_before; // age x year x fleet (including surveys)
  array<Type> fleetCumulativeIncidence; // age x year x fleet (including surveys)
  array<Type> otherCumulativeIncidence; // age x year x causes
  matrix<Type> ssbSurvival_before;	     // age x year

  inline MortalitySet():
    totalZ(),
    fleetSurvival_before(),
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
					  fleetSurvival_before(x.fleetSurvival_before),
					  fleetCumulativeIncidence(x.fleetCumulativeIncidence),
					  otherCumulativeIncidence(x.otherCumulativeIncidence),
					  ssbSurvival_before(x.ssbSurvival_before) {}
  
  MortalitySet(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF);

  // For simulation based forecast
  void updateYear(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, int y);
});

SOURCE(
	 template<class Type>
	 MortalitySet<Type>::MortalitySet(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF) :
	 totalZ(totFFun(dat,conf,logF).matrix()),
	 fleetSurvival_before(conf.keyLogFsta.dim(1), dat.natMor.dim(0), conf.keyLogFsta.dim(0)),
	 fleetCumulativeIncidence(conf.keyLogFsta.dim(1), dat.natMor.dim(0), conf.keyLogFsta.dim(0)),
	 otherCumulativeIncidence(conf.keyLogFsta.dim(1), dat.natMor.dim(0), 1),
	 ssbSurvival_before(conf.keyLogFsta.dim(1), dat.natMor.dim(0))
	 {
	   // Only constant mortality over entire year (TODO: implement other options)
	   int nFleet = conf.keyLogFsta.dim(0);
	   int nAge = conf.keyLogFsta.dim(1);
	   int nYear = dat.natMor.dim(0);
	   // totalFCI = matrix(nAge,nYear);
	   // totalFCI.setZero();
	   fleetSurvival_before.setZero();
	   fleetCumulativeIncidence.setZero();
	   otherCumulativeIncidence.setZero();
	   ssbSurvival_before.setZero();
    
	   for(int a = 0; a < nAge; ++a){
	     for(int y = 0; y < nYear; ++y){
	       totalZ(a,y) += dat.natMor(y,a);
	       if(totalZ(a,y) > 0){
		 Type v = (1.0 - exp(-totalZ(a,y))) / totalZ(a,y);
		 otherCumulativeIncidence(a,y,0) = dat.natMor(y,a) * v;
		 // TODO: Implement option to have time of spawning as fraction of year
		 Type vssb = dat.natMor(y,a) * dat.propM(y,a);
		 for(int f = 0; f < nFleet; ++f){
		   fleetSurvival_before(a,y,f) = exp(-totalZ(a,y) * dat.sampleTimes(f));
		   int i = conf.keyLogFsta(f,a);
		   if(i > (-1)){	// Has catch
		     fleetCumulativeIncidence(a,y,f) = exp(logF(i,y)) * v;
		     // totalFCI(a,y) += fleetCumulativeIncidence(a,y,f);
		     // TODO: Implement option to have time of spawning as fraction of year
		     vssb += exp(logF(i,y)) * dat.propF(y,a,f);
		   }// Otherwise, stay 0
		 }
		 ssbSurvival_before(a,y) = exp(-vssb);
	       }
	     }
	   }
	   return;
	 }
	 )

SOURCE(
	 template<class Type>
	 void MortalitySet<Type>::updateYear(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, int y){
	   int nFleet = conf.keyLogFsta.dim(0);
	   int nAge = conf.keyLogFsta.dim(1);
	   int nYear = dat.natMor.dim(0);
	   if(y > nYear || y < 0)
	     Rf_error("MortalitySet.updateYear: Year not in range");
	   for(int a = 0; a < nAge; ++a){
	     Type newTZ = dat.natMor(y,a);
	     for(int f = 0; f < nFleet; ++f){
	       int i = conf.keyLogFsta(f,a);
	       if(i > (-1))
		 newTZ += exp(logF(i,y));
	     }
	     totalZ(a,y) = newTZ;
	     if(totalZ(a,y) > 0){
	       Type v = (1.0 - exp(-totalZ(a,y))) / totalZ(a,y);	  
	       otherCumulativeIncidence(a,y,0) = dat.natMor(y,a) * v;
	       // TODO: Implement option to have time of spawning as fraction of year
	       Type vssb = dat.natMor(y,a) * dat.propM(y,a);
	       for(int f = 0; f < nFleet; ++f){
		 int i = conf.keyLogFsta(f,a);
		 fleetSurvival_before(a,y,f) = exp(-totalZ(a,y) * dat.sampleTimes(f));
		 if(i > (-1)){	// Has catch
		   fleetCumulativeIncidence(a,y,f) = exp(logF(i,y)) * v;
		   // totalFCI(a,y) += fleetCumulativeIncidence(a,y,f);
		   // TODO: Implement option to have time of spawning as fraction of year
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
