SAM_DEPENDS(define)
SAM_DEPENDS(logspace)
SAM_DEPENDS(recruitment)
SAM_DEPENDS(convenience)
SAM_DEPENDS(newton)
SAM_DEPENDS(hcr)

HEADER(
enum ConstraintType {
		     Constrain_Fbar = 0,
		     Constrain_Catch = 1,
		     Constrain_SSB = 2,
		     Constrain_TSB = 3,
		     Constrain_Landing = 4,
		     Constrain_KeepRelF = 5,
		     Constrain_HCR = 6,
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
  vector<Type> settings;

  inline FConstraint() = default;

  FConstraint(SEXP x);

  template<class T>
  inline FConstraint(const FConstraint<T>& x) : Amin(x.Amin),
						Amax(x.Amax),
						fleet(x.fleet),
						relative(x.relative),
						cstr(x.cstr),
						target(x.target),
						settings(x.settings) {}
};
       )

SOURCE(
	 template<class Type>
	 FConstraint<Type>::FConstraint(SEXP x){
	   Amin = Rf_asInteger(getListElement(x,"Amin",  &isNumericScalar));
	   Amax = Rf_asInteger(getListElement(x,"Amax",  &isNumericScalar));
	   fleet = Rf_asInteger(getListElement(x,"fleet", &isNumericScalar));
	   relative = Rf_asInteger(getListElement(x,"relative",  &isNumericScalar));
	   cstr = static_cast<ConstraintType>(Rf_asInteger(getListElement(x,"cstr",  &isNumericScalar)));
	   target = Rf_asReal(getListElement(x,"target",  &isNumericScalar));
	   settings = asVector<Type>(getListElement(x,"settings", &Rf_isNumeric));
	 }
       );

SAM_SPECIALIZATION(struct FConstraint<double>);
SAM_SPECIALIZATION(struct FConstraint<TMBad::ad_aug>);

SAM_SPECIALIZATION(struct tmbutils::vector<FConstraint<double> >);
SAM_SPECIALIZATION(struct tmbutils::vector<FConstraint<TMBad::ad_aug> >);

HEADER(
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
      (*this)(i) = FConstraint<Type>(other(i));
  }  
};
       )



#ifndef WITH_SAM_LIB
namespace ConstrainCalculations {


  
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

  
  SAM_SPECIALIZATION(vector<double> getFleetLogFbar(dataSet<double>&, confSet&, vector<int>&, vector<double>&, int, int));
  SAM_SPECIALIZATION(vector<TMBad::ad_aug> getFleetLogFbar(dataSet<TMBad::ad_aug>&, confSet&, vector<int>&, vector<TMBad::ad_aug>&, int, int));
  

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
  
  SAM_SPECIALIZATION(double getFleetCatch(dataSet<double>&, confSet&, vector<int>&, array<double>&, vector<double>&, int, int, int, int));
  SAM_SPECIALIZATION(TMBad::ad_aug getFleetCatch(dataSet<TMBad::ad_aug>&, confSet&, vector<int>&, array<TMBad::ad_aug>&, vector<TMBad::ad_aug>&, int, int, int, int));
  
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

  SAM_SPECIALIZATION(double getFleetLanding(dataSet<double>&, confSet&, vector<int>&, array<double>&, vector<double>&, int, int, int, int));
  SAM_SPECIALIZATION(TMBad::ad_aug getFleetLanding(dataSet<TMBad::ad_aug>&, confSet&, vector<int>&, array<TMBad::ad_aug>&, vector<TMBad::ad_aug>&, int, int, int, int));
  
  // Begining of next year
  template<class Type>
  Type getSSB(dataSet<Type>& dat, confSet& conf, vector<int>& cFleets, Recruitment<Type> &recruit, array<Type>& logN, vector<Type>& logF, int y, int a0, int a1, bool rel = false){
    // Current SSB for recruitment
    int yn = std::min(y+1, dat.propMat.dim[0]-1);
    Type logThisSSB = R_NegInf;
    int ys = std::max(yn-conf.minAge, 0); // Year of birth for predicted logN
    for(int a = 0; a < logN.dim[0]; a++){
      // if(dat.propF(ys,a) > 0) // needs fleet index
      // 	Rf_warning("For next year ssb constraints, prediction of recruitment assumes no fishing before spawning.");
      logThisSSB = logspace_add_SAM(logThisSSB, logN(a,ys) + log(dat.propMat(ys,a)) + log(dat.stockMeanWeight(ys,a) + dat.propM(ys,a) * dat.natMor(ys,a)));
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
	logssb = logspace_add_SAM(logssb, logNp(a-conf.minAge) + log(dat.propMat(yn,a-conf.minAge)) + log(dat.stockMeanWeight(yn,a-conf.minAge)) + dat.propM(ys,a) * dat.natMor(ys,a));
      if(dat.propMat(yn,a-conf.minAge) > 0)
	llssb = logspace_add_SAM(llssb, logN(a-conf.minAge,y) + log(dat.propMat(y,a-conf.minAge)) + log(dat.stockMeanWeight(y,a-conf.minAge)) + dat.propM(ys,a) * dat.natMor(ys,a));
    }
    //
    if(rel)
      return logssb - llssb;
    return logssb;    
  };

  SAM_SPECIALIZATION(double getSSB(dataSet<double>&, confSet&, vector<int>&, Recruitment<double>&, array<double>&, vector<double>&, int, int, int, bool));
  SAM_SPECIALIZATION(TMBad::ad_aug getSSB(dataSet<TMBad::ad_aug>&, confSet&, vector<int>&, Recruitment<TMBad::ad_aug>&, array<TMBad::ad_aug>&, vector<TMBad::ad_aug>&, int, int, int, bool));
  

  // Begining of next year
  template<class Type>
  Type getTSB(dataSet<Type>& dat, confSet& conf, vector<int>& cFleets, Recruitment<Type> &recruit, array<Type>& logN, vector<Type>& logF, int y, int a0, int a1, bool rel = false){
    // Current SSB for recruitment
    int yn = std::min(y+1, dat.propMat.dim[0]-1);
    Type logThisSSB = R_NegInf;
    int ys = std::max(yn-conf.minAge, 0); // Year of birth for predicted logN
    for(int a = 0; a < logN.dim[0]; a++){
      // if(dat.propF(ys,a) > 0) // needs fleet index
      // 	Rf_warning("For next year tsb constraints, prediction of recruitment assumes no fishing before spawning.");
      logThisSSB = logspace_add_SAM(logThisSSB, logN(a,ys) + log(dat.propMat(ys,a)) + log(dat.stockMeanWeight(ys,a)) + dat.propM(ys,a) * dat.natMor(ys,a));
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

  SAM_SPECIALIZATION(double getTSB(dataSet<double>&, confSet&, vector<int>&, Recruitment<double>&, array<double>&, vector<double>&, int, int, int, bool));
  SAM_SPECIALIZATION(TMBad::ad_aug getTSB(dataSet<TMBad::ad_aug>&, confSet&, vector<int>&, Recruitment<TMBad::ad_aug>&, array<TMBad::ad_aug>&, vector<TMBad::ad_aug>&, int, int, int, bool));


  /*

   */
  template<class Type>
  Type getHCRBiomass(dataSet<Type>& dat, confSet& conf, vector<int>& cFleets, array<Type>& logN, vector<Type>& logF, int y, int a0, int a1, int lag, int bioType){
    Type logBiomass = R_NegInf;
    int ys = std::max(y-lag, 0);
    for(int a = a0-conf.minAge; a < a1-conf.minAge; a++){
      Type lb = logN(a,ys) + log(dat.stockMeanWeight(ys,a));
      if(bioType == 0){		// SSB
	if(dat.propMat(ys,a) > 0){
	  Type Fa = 0.0;
	  if(lag > 0){
	    for(int ii = 0; ii < cFleets.size(); ++ii){
	      int f = cFleets(ii);
	      if(dat.propF(ys,a,f) > 0)
		if(conf.keyLogFsta(f,a) > (-1))
		  Fa += dat.propF(ys,a,f) * exp(logF(conf.keyLogFsta(f,a)));
	    }
	  }
	  // else if(dat.propF(ys,a) > 0 && lag >= 0){ // propF needs fleet index
	  //   Rf_error("With lag 0, SSB without fishing before spawning time is used for the HCR.");
	  // }
	  lb += log(dat.propMat(ys,a)) + Fa + dat.propM(ys,a) * dat.natMor(ys,a);
	}
      }else if(bioType == 1){	// TSB
	// Do no more
      }else{
	Rf_error("Unknown biomass type for HCR");
      }
      logBiomass = logspace_add_SAM(logBiomass, lb);
    }
    return exp(logBiomass);
  }

  SAM_SPECIALIZATION(double getHCRBiomass(dataSet<double>&, confSet&, vector<int>&, array<double>&, vector<double>&, int, int, int, int, int));
  SAM_SPECIALIZATION(TMBad::ad_aug getHCRBiomass(dataSet<TMBad::ad_aug>&, confSet&, vector<int>&, array<TMBad::ad_aug>&, vector<TMBad::ad_aug>&, int, int, int, int, int));


  
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
  
  struct ForecastF : NewtonFunctor {
    dataSet<ad> dat;
    confSet conf;
    vector<int> cFleets;
    Recruitment<ad> recruit;
    FConstraintList<ad> cstrs;
    vector<ad> lastLogF;
    array<ad> logN;
    int y;

    ForecastF(dataSet<ad> dat_, confSet conf_, vector<int> cFleets_, Recruitment<ad> recruit_, FConstraintList<ad> cstrs_, vector<ad> lastLogF_, array<ad> logN_, int y_) : dat(dat_), conf(conf_), cFleets(cFleets_), recruit(recruit_), cstrs(cstrs_), lastLogF(lastLogF_), logN(logN_), y(y_) {};

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

	}else if(cstr.cstr == ConstraintType::Constrain_HCR){
	  // Give error if the relative flag makes it this far:
	  RELATIVE_CONSTRAINTS(
			       // Last year
			       Rf_error("A HCR target cannot be relative");
			       ,
			       // Total
			       Rf_error("A HCR target cannot be relative");
			       ,
			       // Fleet
			       Rf_error("A HCR target cannot be relative");
			       );
  
	  /* cstr.settings should be:
	     0: Indicator of biomass type (TSB,SSB,...)
	     1: Indicator of SSB year (-1,0,1)
	     2: Indicator of target Type (F, Catch, Landing,...)
	     3: Amin for TSB/SSB
	     4: Amax for TSB/SSB
	     5: F_origin for HCR
	     6: F_cap for HCR
	     7: B_origin for HCR
	     8: B_cap for HCR
	     9: B_trigger for HCR
	   */
	  int biomassType = CppAD::Integer(cstr.settings(0));
	  int biomassLag = CppAD::Integer(cstr.settings(1));
	  int targetType = CppAD::Integer(cstr.settings(2));
	  int bioA0 = CppAD::Integer(cstr.settings(3));
	  int bioA1 = CppAD::Integer(cstr.settings(4));
	  // Get current biomass
	  ad BforHCR = getHCRBiomass(dat, conf, cFleets, logN, lastLogF, y, bioA0, bioA1, biomassLag, biomassType);
	  // Get optimization target
	  vector<ad> hcrConf = cstr.settings.segment(4,6); // Start to early and overwrite
	  // Insert target
	  hcrConf(0) = exp(cstr.target);
	  ad trgt = hcr(BforHCR, hcrConf);
	  // kappa
	  if(targetType == 0){	  // F target
	    if(cstr.fleet == (-1)){	// Total
	      ad tmp = logFbar - trgt;
	      kappa += tmp * tmp;
	    }else{		// Fleet
	      ad tmp = fleetLogFbar(cstr.fleet) - trgt;
	      kappa += tmp * tmp;
	    }
	  }else if(targetType == 1){ // Catch target
	    ad logC = getFleetCatch(dat, conf, cFleets, logN, newLogF, y, cstr.Amin, cstr.Amax, cstr.fleet);
	    ad tmp = logC - trgt;
	    kappa += tmp * tmp;
	  }else if(targetType == 2){ // Landing target
	    ad logL = getFleetLanding(dat, conf, cFleets, logN, newLogF, y, cstr.Amin, cstr.Amax, cstr.fleet);
	    ad tmp = logL - trgt;
	    kappa += tmp * tmp;
	  }else{
	    Rf_error("Unknown target type for HCR");
	  }
	}else{
	  Rf_error("Constraint type not implemented");
	}
      }
      return kappa;
    }
   
  };

};				// End of namespace
#endif

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
			      newton::newton_config& cfg)SOURCE({
  
  paraSet<TMBad::ad_aug> parad(par);
  Recruitment<TMBad::ad_aug> recruit = makeRecruitmentFunction(conf,parad);

  vector<int> cFleets = getCatchFleets(dat.fleetTypes);

  // Should be deleted by NewtonWrapper
  std::shared_ptr<NewtonFunctor> p_fc(new ConstrainCalculations::ForecastF(dat,conf,cFleets,recruit,cstrs,lastLogF,logN,y));
  
  // vector<double> s0(cFleets.size());
  // s0.setConstant(0);
  // vector<Type> start(tryStart(fc,s0,200));
  vector<Type> start(cFleets.size());
  start.setConstant(0.0);
  
  vector<Type> res = SAM_Newton(p_fc, start, cfg);
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
				})

SAM_SPECIALIZATION(vector<double> calculateNewFVec(dataSet<double>&, confSet&, paraSet<double>&, FConstraintList<double>&, vector<double>&, array<double>&, int, newton::newton_config&));
  SAM_SPECIALIZATION(vector<TMBad::ad_aug> calculateNewFVec(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, FConstraintList<TMBad::ad_aug>&, vector<TMBad::ad_aug>&, array<TMBad::ad_aug>&, int, newton::newton_config&));

