#pragma once
#ifndef SAM_INCIDENCE_HPP
#define SAM_INCIDENCE_HPP


// template<class Type>
// struct Mortality {
//   Type Hazard(int yi, int ai, Type t0, Type t1);
//   Type CumulativeHazard(int yi, int ai, Type t0, Type t1);
// };

// template<class Type>
// struct Constant_F_Mortality : Mortality {
//   array<Type> logF
// }



template<class Type>
struct MortalitySet {
  // TODO: Switch to log-scale calculations?
  matrix<Type> totalZ;		// age x year
  // matrix<Type> totalFCI;
  array<Type> fleetSurvival_before; // age x year x fleet (including surveys)
  array<Type> fleetCumulativeIncidence; // age x year x fleet (including surveys)
  array<Type> otherCumulativeIncidence; // age x year x causes
  matrix<Type> ssbSurvival_before;	     // age x year

  template<class T>
  MortalitySet(const MortalitySet<T> x) : totalZ(x.totalZ),
					  fleetSurvival_before(x.fleetSurvival_before),
					  fleetCumulativeIncidence(x.fleetCumulativeIncidence),
					  otherCumulativeIncidence(x.otherCumulativeIncidence),
					  ssbSurvival_before(x.ssbSurvival_before) {}
  
  MortalitySet(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF) :
    totalZ(totFFun(conf,logF).matrix()),
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

  // For simulation based forecast
  void updateYear(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, int y){
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
  
};

  


#endif
