#pragma once
#ifndef SAM_TOF_HPP
#define SAM_TOF_HPP


enum ConstraintType {
		     Constrain_Fbar,
		     Constrain_Catch,
		     Constrain_SSB,
		     Constrain_TSB,
		     Constrain_Landing,
		     Constrain_KeepRelF,
		     Constrain_NONE = 99
};


template<class Type>
struct FConstraint {
  int Amin;
  int Amax;
  int fleet;
  int relative; 		// -3: Absolute; -2: relative to last year; -1: relative to total; >=0: relative to fleet #
  ConstraintType cstr;
  Type target;

  FConstraint() = default;

  FConstraint(SEXP x){
    Amin = (int)*REAL(getListElement(x,"Amin"));
    Amax = (int)*REAL(getListElement(x,"Amax"));
    fleet = (int)*REAL(getListElement(x,"fleet"));
    relative = (int)*REAL(getListElement(x,"relative"));
    cstr = static_cast<ConstraintType>((int)*REAL(getListElement(x,"cstr")));
    target = (Type)*REAL(getListElement(x,"target"));
  };

  template<class T>
  FConstraint(const FConstraint<T>& x) : Amin(x.Amin),
					 Amax(x.Amax),
					 fleet(x.fleet),
					 relative(x.relative),
					 cstr(x.cstr),
					 target(x.target) {}
};

template<class Type>
struct FConstraintList : vector<FConstraint<Type> > {
  FConstraintList() : vector<FConstraint<Type> >() {};
  FConstraintList(int n) : vector<FConstraint<Type> >(n) {};
  FConstraintList(SEXP x){ 
    (*this).resize(Rf_length(x));
    for(int i=0; i<Rf_length(x); i++){
      (*this)(i) = FConstraint<Type>(VECTOR_ELT(x, i));
    }
  }
  template<class T>
  FConstraintList(const FConstraintList<T>& other) : vector<FConstraint<Type> >(other.size()) {
    for(int i = 0; i < (int)other.size(); ++i)
      (*this)(i) = FConstraint<T>(other(i));
  }
  
};




namespace ConstrainCalculations {



 
  // Newton approach


  vector<int> getCatchFleets(vector<int> fleetTypes){
    std::vector<int> r;
    for(int i = 0; i < (int)fleetTypes.size(); ++i)
      if(fleetTypes(i) == 0)
	r.push_back(i);
    return vector<int>(r);
  };
  
  template<class Type>
  vector<Type> getFleetLogFbar(dataSet<Type>& dat, confSet& conf, vector<int>& cFleets, vector<Type>& logF, int a0, int a1){
    vector<Type> fbar(cFleets.size());
    fbar.setConstant(0.0);
    for(int i = 0; i < cFleets.size(); ++i){
      int f = cFleets(i);
      for(int a=a0; a<=a1; a++){
	if(conf.keyLogFsta(f,a-conf.minAge) > (-1))
	  fbar(i) += exp(logF(conf.keyLogFsta(f,a-conf.minAge)));
      }
      fbar(i) /= Type(a1-a0 + 1.0);
    }
    return log(fbar);
  };

  template<class Type>
  Type getFleetCatch(dataSet<Type>& dat, confSet& conf, vector<int>& cFleets, array<Type>& logN, vector<Type>& logF, int y, int a0, int a1, int fleet){
    Type logCat = R_NegInf;
    // Type cat = 0.0;
    for(int a=a0; a<=a1; a++){
      Type logZa = log(dat.natMor(y, a-conf.minAge));
      // Type Za = dat.natMor(y, a-conf.minAge);
      for(int ii = 0; ii < cFleets.size(); ++ii){
	int f = cFleets(ii);
	if(conf.keyLogFsta(f,a-conf.minAge) > (-1)){
	  // Za += exp(logF(conf.keyLogFsta(f,a-conf.minAge)));
	  logZa = logspace_add_SAM(logZa, logF(conf.keyLogFsta(f,a-conf.minAge)));
	}
      }
      Type logv = logspace_sub_SAM(Type(0.0), -exp(logZa)) - logZa;
      // Type v = (1.0 - exp(-Za)) / Za;
      int f0 = 0;
      int f1 = cFleets.size()-1;
      if(fleet > (-1)){
	f0 = fleet;
	f1 = fleet;
      }
    for(int ii = f0; ii <= f1; ++ii){
      int f = cFleets(ii);
	if(conf.keyLogFsta(f,a-conf.minAge) > (-1)){
	  Type logFI = logv + logF(conf.keyLogFsta(f,a-conf.minAge));
	  // Type FI = v * exp(logF(conf.keyLogFsta(f,a-conf.minAge)));
	  // Type cc = FI * exp(logN(a-conf.minAge,y)) * dat.catchMeanWeight(y, a-conf.minAge, f);
	  Type lc =  logFI + logN(a-conf.minAge,y) + log(dat.catchMeanWeight(y, a-conf.minAge, f));
	  logCat = logspace_add_SAM(logCat, lc);
	  // cat += cc;
	}
      }
    }
    // return log(cat);
    return logCat;
  }
  
  template<class Type>
  Type getFleetLanding(dataSet<Type>& dat, confSet& conf, vector<int>& cFleets, array<Type>& logN, vector<Type>& logF, int y, int a0, int a1, int fleet){
    Type logCat = R_NegInf;
    for(int a=a0; a<=a1; a++){
      Type logZa = log(dat.natMor(y, a-conf.minAge));
      for(int ii = 0; ii < cFleets.size(); ++ii){
	int f = cFleets(ii);
	if(conf.keyLogFsta(f,a-conf.minAge) > (-1))
	  logZa = logspace_add_SAM(logZa, logF(conf.keyLogFsta(f,a-conf.minAge)));
      }
      Type logv = logspace_sub_SAM(Type(0.0), -exp(logZa)) - logZa;
      int f0 = 0;
      int f1 = cFleets.size()-1;
      if(fleet > (-1)){
	f0 = fleet;
	f1 = fleet;
      }
    for(int ii = f0; ii <= f1; ++ii){
      int f = cFleets(ii);	
	if(conf.keyLogFsta(f,a-conf.minAge) > (-1)){
	  Type LW = dat.landMeanWeight(y, a-conf.minAge, f);
	  Type LF = dat.landFrac(y,a-conf.minAge,f);
	  if(LW > 0 && LF > 0){
	    Type logFI = logv + logF(conf.keyLogFsta(f,a-conf.minAge));
	    Type lc =  log(LF) + logFI + logN(a-conf.minAge,y) + log(LW);
	    logCat = logspace_add_SAM(logCat, lc);
	  }
	}
      }
    }
    return logCat;
  }

  // Begining of next year
  template<class Type>
  Type getSSB(dataSet<Type>& dat, confSet& conf, vector<int>& cFleets, Recruitment<Type> &recruit, array<Type>& logN, vector<Type>& logF, int y, int a0, int a1, bool rel = false){
    // Current SSB for recruitment
    int yn = std::min(y+1, dat.propMat.dim[0]-1);
    Type logThisSSB = R_NegInf;
    int ys = std::max(yn-conf.minAge, 0); // Year of birth for predicted logN
    for(int a = 0; a < logN.dim[0]; a++){
      logThisSSB = logspace_add_SAM(logThisSSB, logN(a,ys) + log(dat.propMat(ys,a)) + log(dat.stockMeanWeight(ys,a)));
    }
    // Next N
    vector<Type> logNp(logN.dim[0]);
    logNp.setConstant(R_NegInf);
    //// Recruit
    logNp(0) = recruit(logThisSSB,logN(0,y), dat.years(0)+y+1);
    //// Middle
    for(int a = 1; a < logNp.size(); ++a){
      Type logZa = log(dat.natMor(y, a-1));
      for(int ii = 0; ii < cFleets.size(); ++ii){
	int f = cFleets(ii);
	if(conf.keyLogFsta(f,a-1) > (-1))
	  logZa = logspace_add_SAM(logZa, logF(conf.keyLogFsta(f,a-1)));
      }     
      logNp(a) = logN(a-1,y) - exp(logZa);
    }
    //// Plus group
    if(conf.maxAgePlusGroup(0)==1){
      int a = logNp.size()-1;
      Type v1 = logNp(a);
      Type logZa = log(dat.natMor(y, a));
      for(int ii = 0; ii < cFleets.size(); ++ii){
	int f = cFleets(ii);
	if(conf.keyLogFsta(f,a) > (-1))
	  logZa = logspace_add_SAM(logZa, logF(conf.keyLogFsta(f,a)));
      }
      Type v2 = logN(a, y) - exp(logZa);
      logNp(a) = logspace_add_SAM(v1,v2);
    }    
   // Next SSB
    Type logssb = R_NegInf;
    Type llssb = R_NegInf;
    for(int a = a0; a <= a1; a++){
      if(dat.propMat(yn,a-conf.minAge) > 0)
	logssb = logspace_add_SAM(logssb, logNp(a-conf.minAge) + log(dat.propMat(yn,a-conf.minAge)) + log(dat.stockMeanWeight(yn,a-conf.minAge)));
      if(dat.propMat(yn,a-conf.minAge) > 0)
	llssb = logspace_add_SAM(llssb, logN(a-conf.minAge,y) + log(dat.propMat(y,a-conf.minAge)) + log(dat.stockMeanWeight(y,a-conf.minAge)));
    }
    //
    if(rel)
      return logssb - llssb;
    return logssb;    
  };

  // Begining of next year
  template<class Type>
  Type getTSB(dataSet<Type>& dat, confSet& conf, vector<int>& cFleets, Recruitment<Type> &recruit, array<Type>& logN, vector<Type>& logF, int y, int a0, int a1, bool rel = false){
    // Current SSB for recruitment
    int yn = std::min(y+1, dat.propMat.dim[0]-1);
    Type logThisSSB = R_NegInf;
    int ys = std::max(yn-conf.minAge, 0); // Year of birth for predicted logN
    for(int a = 0; a < logN.dim[0]; a++){
      logThisSSB = logspace_add_SAM(logThisSSB, logN(a,ys) + log(dat.propMat(ys,a)) + log(dat.stockMeanWeight(ys,a)));
    }
    // Next N
    vector<Type> logNp(logN.dim[0]);
    logNp.setConstant(R_NegInf);
    //// Recruit
    logNp(0) = recruit(logThisSSB,logN(0,y), dat.years(0)+y+1);
    //// Middle
    for(int a = 1; a < logNp.size(); ++a){
      Type logZa = log(dat.natMor(y, a-1));
      for(int ii = 0; ii < cFleets.size(); ++ii){
	int f = cFleets(ii);
	if(conf.keyLogFsta(f,a-1) > (-1))
	  logZa = logspace_add_SAM(logZa, logF(conf.keyLogFsta(f,a-1)));
      }     
      logNp(a) = logN(a-1,y) - exp(logZa);
    }
    //// Plus group
    if(conf.maxAgePlusGroup(0)==1){
      int a = logNp.size()-1;
      Type v1 = logNp(a);
      Type logZa = log(dat.natMor(y, a));
      for(int ii = 0; ii < cFleets.size(); ++ii){
	int f = cFleets(ii);
	if(conf.keyLogFsta(f,a) > (-1))
	  logZa = logspace_add_SAM(logZa, logF(conf.keyLogFsta(f,a)));
      }
      Type v2 = logN(a, y) - exp(logZa);
      logNp(a) = logspace_add_SAM(v1,v2);
    }    
    // Next TSB
    Type logtsb = R_NegInf;
    Type lltsb = R_NegInf;
    for(int a = a0; a <= a1; a++){
      logtsb = logspace_add_SAM(logtsb, logNp(a-conf.minAge) + log(dat.stockMeanWeight(yn,a-conf.minAge)));
      lltsb = logspace_add_SAM(lltsb, logN(a-conf.minAge,y) + log(dat.stockMeanWeight(y,a-conf.minAge)));
    }
    //
    if(rel)
      return logtsb - lltsb;
    return logtsb;
  };
  
  typedef TMBad::ad_aug ad;

#define RELATIVE_CONSTRAINTS(LAST,TOTAL,FLEET)				\
  if(cstr.relative == -3){						\
  }else if(cstr.relative == -2){					\
    LAST								\
  }else if(cstr.relative == -1){					\
    TOTAL								\
  }else if(cstr.relative >= 0 && cstr.relative < cFleets.size() && cstr.relative != cstr.fleet){ \
    FLEET								\
      }else{								\
    Rf_error("Wrong relative constraint code");				\
  }							 
  
  struct ForecastF {
    dataSet<ad> dat;
    confSet conf;
    vector<int> cFleets;
    Recruitment<ad> recruit;
    FConstraintList<ad> cstrs;
    vector<ad> lastLogF;
    array<ad> logN;
    int y;

    ad operator()(const vector<ad>& logFs){
      ad kappa = 0.0;
      //
      // matrix<ad> llogF = toFleetMatrix(dat, conf, lastLogF, (vector<ad>)(logFs * ad(0.0)));
      // vector<ad> llfs = (vector<ad>)(lastLogF - lastLFB);
      // matrix<ad> logF = toFleetMatrix(dat, conf, lastLogF, logFs);

      vector<ad> newLogF = lastLogF;
   
      vector<bool> done(newLogF.size());
      done.setConstant(false);
      for(int i = 0; i < cFleets.size(); ++i){
	int f = cFleets(i);
	for(int a = 0; a < conf.keyLogFsta.dim(1); ++a){
	  int indx = conf.keyLogFsta(f,a);
	  if(indx > (-1) && !done(indx)){
	    newLogF(indx) += logFs(i);
	    done(indx) = true;
	  }
	}
      }
      
      for(int i = 0; i < cstrs.size(); ++i){
	FConstraint<ad> cstr = cstrs(i);
	if(cstr.cstr == ConstraintType::Constrain_NONE){
	  continue;
	}
	// Previous F values
	vector<ad> lastFleetLogFbar = getFleetLogFbar(dat, conf, cFleets, lastLogF, cstr.Amin, cstr.Amax);
	ad lastLogFbar = logspace_sum(lastFleetLogFbar);
	// New F values
	vector<ad> fleetLogFbar = getFleetLogFbar(dat, conf, cFleets, newLogF, cstr.Amin, cstr.Amax);
	ad logFbar = logspace_sum(fleetLogFbar);
	// Add constraint

	if(cstr.cstr == ConstraintType::Constrain_Fbar){

	  if(cstr.fleet == (-1)){	// Total F
	    ad trgt = cstr.target;
	    RELATIVE_CONSTRAINTS(
				 // Last year
				 trgt += lastLogFbar;
				 ,
				 // Total
				 Rf_error("A total F constraint can only be relative to last year");
				 ,
				 // Fleet
				 Rf_error("A total F constraint can only be relative to last year");
				 )
	    ad tmp = logFbar - trgt;
	    kappa += tmp * tmp;
	  }else{
	    ad trgt = cstr.target;
	    RELATIVE_CONSTRAINTS(
				 // Last year
				 trgt += lastFleetLogFbar(cstr.fleet);
				 ,
				 // Total
				 trgt += logFbar;
				 ,
				 // Fleet
				 trgt += fleetLogFbar(cstr.relative);
				 )
	      ad tmp = fleetLogFbar(cstr.fleet) - trgt;
	    kappa += tmp * tmp;
	  }

	}else if(cstr.cstr == ConstraintType::Constrain_Catch){

	  ad logC = getFleetCatch(dat, conf, cFleets, logN, newLogF, y, cstr.Amin, cstr.Amax, cstr.fleet);
	  ad trgt = cstr.target;
	  RELATIVE_CONSTRAINTS(
			       // Last year
			       ad logCL = getFleetCatch(dat, conf, cFleets, logN, lastLogF, y-1, cstr.Amin, cstr.Amax, cstr.fleet);
			       trgt += logCL;
			       ,
			       // Total
			       ad logCL = getFleetCatch(dat, conf, cFleets, logN, newLogF, y, cstr.Amin, cstr.Amax, -1);
			       trgt += logCL;
			       ,
			       // Fleet
			       ad logCL = getFleetCatch(dat, conf, cFleets, logN, newLogF, y, cstr.Amin, cstr.Amax, cstr.relative);
			       trgt += logCL;
			       )	
	  ad tmp = logC - trgt;
	  kappa += tmp * tmp;

	}else if(cstr.cstr == ConstraintType::Constrain_SSB){
	  SAM_ASSERT(cstr.relative <= -2, "SSB constraints can only be relative to last year");
	  ad logB = getSSB(dat,conf, cFleets, recruit, logN, newLogF, y, cstr.Amin, cstr.Amax, cstr.relative == (-2));
	  ad tmp = logB - cstr.target;
	  kappa += tmp * tmp;
	  
	}else if(cstr.cstr == ConstraintType::Constrain_TSB){
	  SAM_ASSERT(cstr.relative <= -2, "TSB constraints can only be relative to last year");
	  ad logB = getTSB(dat,conf, cFleets, recruit, logN, newLogF, y, cstr.Amin, cstr.Amax, cstr.relative == (-2));
	  ad tmp = logB - cstr.target;
	  kappa += tmp * tmp;
	  
	}else if(cstr.cstr == ConstraintType::Constrain_Landing){

	  ad logL = getFleetLanding(dat, conf, cFleets, logN, newLogF, y, cstr.Amin, cstr.Amax, cstr.fleet);
	  ad trgt = cstr.target;
	  RELATIVE_CONSTRAINTS(
			       // Last year
			       ad logLL = getFleetLanding(dat, conf, cFleets, logN, lastLogF, y-1, cstr.Amin, cstr.Amax, cstr.fleet);
			       trgt += logLL;
			       ,
			       // Total
			       ad logLL = getFleetLanding(dat, conf, cFleets, logN, newLogF, y, cstr.Amin, cstr.Amax, -1);
			       trgt += logLL;
			       ,
			       // Fleet
			       ad logLL = getFleetLanding(dat, conf, cFleets, logN, newLogF, y, cstr.Amin, cstr.Amax, cstr.relative);
			       trgt += logLL;
			       )	
	  ad tmp = logL - trgt;
	  kappa += tmp * tmp;

	}else if(cstr.cstr == ConstraintType::Constrain_KeepRelF){

	  SAM_ASSERT(cstr.fleet >= 0,"Keep relative Fbar can only be used for fleets, not total");
	  SAM_ASSERT(cstr.relative >= 0 && cstr.fleet != cstr.relative, "Keep relative can only be used relative to other fleets in same year");
	  ad tmp = (lastFleetLogFbar(cstr.fleet) - lastFleetLogFbar(cstr.relative)) - (fleetLogFbar(cstr.fleet) - fleetLogFbar(cstr.relative));
	  kappa += tmp * tmp;	  
	}else{
	  Rf_error("Constraint type not implemented");
	}
      }
      return kappa;
    }
   
  };

};				// End of namespace


// vector<double> tryStart(ConstrainCalculations::ForecastF fc,
// 			vector<double> start,
// 			int n){
//   int i = 0;
//   vector<double> r(start);
//   double v = R_PosInf;
//   while(i++ < n){
//     vector<double> test(start);
//     for(int j = 0; j < test.size(); ++j)
//       test(j) += Rf_runif(-1.0,1.0);
//     vector<TMBad::ad_aug> t2(test);
//     double v2 = asDouble(fc(t2));
//     if(R_finite(v2) && v2 < v){
//       v = v2;
//       r = test;
//     }
//   }
//   return r;
// }

template<class Type>
vector<Type> calculateNewFVec(dataSet<Type>& dat,
			      confSet& conf,
			      paraSet<Type>& par,
			      FConstraintList<Type>& cstrs,
			      vector<Type>& lastLogF,
			      array<Type>& logN,
			      int y,
			      newton::newton_config& cfg){

  paraSet<TMBad::ad_aug> parad(par);
  Recruitment<TMBad::ad_aug> recruit = makeRecruitmentFunction(conf,parad);

  vector<int> cFleets = ConstrainCalculations::getCatchFleets(dat.fleetTypes);

  ConstrainCalculations::ForecastF fc = {dat, conf, cFleets, recruit, cstrs, lastLogF, logN, y};

  
  // vector<double> s0(cFleets.size());
  // s0.setConstant(0);
  // vector<Type> start(tryStart(fc,s0,200));
  vector<Type> start(cFleets.size());
  start.setConstant(0.0);
  
  vector<Type> res = newton::Newton(fc, start, cfg);
  vector<Type> newLogF = lastLogF;// - lastLogFbar;

  vector<bool> done(newLogF.size());
  done.setConstant(false);
  for(int i = 0; i < cFleets.size(); ++i){
    int f = cFleets(i);
    for(int a = 0; a < conf.keyLogFsta.dim(1); ++a){
      int indx = conf.keyLogFsta(f,a);
      if(indx > (-1) && !done(indx)){
	newLogF(indx) += res(i);
	done(indx) = true;
      }
    }
  }
  return newLogF;
}


#endif
