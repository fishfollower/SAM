SAM_DEPENDS(define)
SAM_DEPENDS(logspace)
SAM_DEPENDS(recruitment)
SAM_DEPENDS(convenience)
SAM_DEPENDS(newton)
SAM_DEPENDS(hcr)
SAM_DEPENDS(incidence)
SAM_DEPENDS(derived)
SAM_DEPENDS(predn)
SAM_DEPENDS(extend_array)
SAM_DEPENDS(equilibrium)
SAM_DEPENDS(predf)

HEADER(
enum ConstraintType {
		     Constrain_Fbar = 0,
		     Constrain_Catch = 1,
		     Constrain_SSB = 2,
		     Constrain_TSB = 3,
		     Constrain_Landing = 4,
		     Constrain_KeepRelF = 5,
		     Constrain_HCR = 6,
		     Constrain_ERB = 7,
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
  int useNonLinearityCorrection;

  inline FConstraint() :
    Amin(),
    Amax(),
    fleet(),
    relative(),
    cstr(),
    target(),
    settings(),
    useNonLinearityCorrection()
  {}

  FConstraint(SEXP x);

  template<class T>
  inline FConstraint(const FConstraint<T>& x) : Amin(x.Amin),
						Amax(x.Amax),
						fleet(x.fleet),
						relative(x.relative),
						cstr(x.cstr),
						target(x.target),
						settings(x.settings),
						useNonLinearityCorrection(x.useNonLinearityCorrection) {}
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
	   useNonLinearityCorrection = Rf_asLogical(getListElement(x,"useNonLinearityCorrection",  &Rf_isLogical));
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
  FConstraintList(SEXP x) : vector<FConstraint<Type> >(Rf_length(x)){ 
    //(*this).resize(Rf_length(x));
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
	  fbar(i) += exp(logF(conf.keyLogFsta(f,a-conf.minAge))) * (dat.sampleTimesEnd(f) - dat.sampleTimesStart(f));
      }
      fbar(i) /= Type(a1-a0 + 1.0);
    }
    return log(fbar);
  };

  
  SAM_SPECIALIZATION(vector<double> getFleetLogFbar(dataSet<double>&, confSet&, vector<int>&, vector<double>&, int, int));
  SAM_SPECIALIZATION(vector<TMBad::ad_aug> getFleetLogFbar(dataSet<TMBad::ad_aug>&, confSet&, vector<int>&, vector<TMBad::ad_aug>&, int, int));
  

  template<class Type>
  Type getFleetCatch(dataSet<Type>& dat, confSet& conf, array<Type>& logN, array<Type>& logF, MortalitySet<Type>& mort, int y, int a0, int a1, int fleet){
    Type logCat = R_NegInf;
    // Type cat = 0.0;
    int f0 = 0;
    int f1 = conf.keyLogFsta.dim(0)-1;
    if(fleet > (-1)){
      f0 = fleet;
      f1 = fleet;
    }
    for(int f = f0; f <= f1; ++f){
      if(dat.fleetTypes(f) == 0){
	for(int a=a0; a<=a1; a++){
	  if(conf.keyLogFsta(f,a-conf.minAge) > (-1)){
	    Type lc = logN(a-conf.minAge,y) + mort.logFleetSurvival_before(a-conf.minAge,y,f) + log(mort.fleetCumulativeIncidence(a-conf.minAge,y,f)) + log( dat.catchMeanWeight(y, a-conf.minAge, f));
	    logCat = logspace_add_SAM(logCat, lc);
	    // cat += cc;
	  }
	}
      }
    }
    // return log(cat);
    return logCat;
  }
  
  SAM_SPECIALIZATION(double getFleetCatch(dataSet<double>&, confSet&, array<double>&, array<double>&, MortalitySet<double>&, int, int, int, int));
  SAM_SPECIALIZATION(TMBad::ad_aug getFleetCatch(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>&, int, int, int, int));
  
  template<class Type>
  Type getFleetLanding(dataSet<Type>& dat, confSet& conf, array<Type>& logN, array<Type>& logF, MortalitySet<Type>& mort, int y, int a0, int a1, int fleet){
    Type logCat = R_NegInf;
    // Type cat = 0.0;
    int f0 = 0;
    int f1 = conf.keyLogFsta.dim(0)-1;
    if(fleet > (-1)){
      f0 = fleet;
      f1 = fleet;
    }
    for(int f = f0; f <= f1; ++f){
      if(dat.fleetTypes(f) == 0){
	for(int a=a0; a<=a1; a++){
	  if(conf.keyLogFsta(f,a-conf.minAge) > (-1)){
	    Type LF = dat.landFrac(y,a-conf.minAge,f);
	    Type LW = dat.landMeanWeight(y,a-conf.minAge,f);
	    if(LF > 0 && LW > 0){
	      Type lc = logN(a-conf.minAge,y) + mort.logFleetSurvival_before(a-conf.minAge,y,f) + log(mort.fleetCumulativeIncidence(a-conf.minAge,y,f)) + log(LF) + log(LW);
	      logCat = logspace_add_SAM(logCat, lc);
	    }
	    // cat += cc;
	  }
	}
      }
    }
    // return log(cat);
    return logCat;
  }

  SAM_SPECIALIZATION(double getFleetLanding(dataSet<double>&, confSet&, array<double>&, array<double>&, MortalitySet<double>&, int, int, int, int));
  SAM_SPECIALIZATION(TMBad::ad_aug getFleetLanding(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>&, int, int, int, int));
  
  template<class Type>
  Type getLogSSB(dataSet<Type>& dat, confSet& conf, array<Type>& logN, array<Type>& logF, MortalitySet<Type> mort, int y, int a0, int a1){
    Type logssb = R_NegInf;
    int yu = std::min(y,dat.propMat.dim(0));
    for(int a=a0; a<=a1; ++a){
      if(dat.propMat(yu,a-conf.minAge) > 0){
	Type lssbNew = logN(a-conf.minAge,y) + log(mort.ssbSurvival_before(a-conf.minAge,yu)) + log(dat.propMat(yu,a-conf.minAge)) + log(dat.stockMeanWeight(yu,a-conf.minAge));
	logssb = logspace_add_SAM(logssb, lssbNew);
      }
    }
    return logssb;
  };

  SAM_SPECIALIZATION(double getLogSSB(dataSet<double>&, confSet&, array<double>&, array<double>&, MortalitySet<double>, int, int, int));
  SAM_SPECIALIZATION(TMBad::ad_aug getLogSSB(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>, int, int, int));

 template<class Type>
 Type getLogERB(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logN, array<Type>& logF, MortalitySet<Type> mort, int y, int a0, int a1){
   Type logssb = R_NegInf;
   Type lssb0 = Type(R_NegInf);
   Type lerb0 = Type(R_NegInf); 
    int yu = std::min(y,dat.propMat.dim(0));
    for(int a=a0; a<=a1; ++a){
      if(dat.propMat(yu,a-conf.minAge) > 0){
	Type lssbNew = logN(a-conf.minAge,y) + log(mort.ssbSurvival_before(a-conf.minAge,yu)) + log(dat.propMat(yu,a-conf.minAge)) + exp(par.logFecundityScaling) * log(dat.stockMeanWeight(yu,a-conf.minAge));
	logssb = logspace_add_SAM(logssb, lssbNew);
		// Y=0 ERB
	Type lerb0New = logN(a-conf.minAge,0) + log(mort.ssbSurvival_before(a-conf.minAge,0)) + log(dat.propMat(0,a-conf.minAge)) + exp(par.logFecundityScaling) * log(dat.stockMeanWeight(0,a-conf.minAge));
	lerb0 = logspace_add_SAM(lerb0, lerb0New);
	// Y=0 SSB
	Type lssb0New = logN(a-conf.minAge,0) + log(mort.ssbSurvival_before(a-conf.minAge,0)) + log(dat.propMat(0,a-conf.minAge)) + log(dat.stockMeanWeight(0,a-conf.minAge));
	lssb0 = logspace_add_SAM(lssb0, lssb0New);
      }
    }
    return logssb - lerb0 + lssb0;
  };

  SAM_SPECIALIZATION(double getLogERB(dataSet<double>&, confSet&, paraSet<double>&, array<double>&, array<double>&, MortalitySet<double>, int, int, int));
  SAM_SPECIALIZATION(TMBad::ad_aug getLogERB(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>, int, int, int));
  
  

  template<class Type>
  Type getLogTSB(dataSet<Type>& dat, confSet& conf, array<Type>& logN, array<Type>& logF, MortalitySet<Type> mort, int y, int a0, int a1){
    Type logtsb = R_NegInf;
    int yu = std::min(y,dat.propMat.dim(0));
    for(int a=a0; a<=a1; ++a){
      Type ltsbNew = logN(a-conf.minAge,y) + log(dat.stockMeanWeight(yu,a-conf.minAge));
      logtsb = logspace_add_SAM(logtsb, ltsbNew);
    }
    return logtsb;
  };

SAM_SPECIALIZATION(double getLogTSB(dataSet<double>&, confSet&, array<double>&, array<double>&, MortalitySet<double>, int, int, int));
SAM_SPECIALIZATION(TMBad::ad_aug getLogTSB(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>, int, int, int));

  
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


  struct ForecastF_Fun {
    dataSet<ad> dat;
    confSet conf;
    paraSet<ad> par;
    vector<int> cFleets;
    Recruitment<ad> recruit;
    FConstraint<ad> cstr;
    array<ad> logN;
    array<ad> historicalLogF;
    array<ad> logitFseason;
    int y;

    ForecastF_Fun(dataSet<ad> dat_, confSet conf_, paraSet<ad> par_, vector<int> cFleets_, Recruitment<ad> recruit_, FConstraint<ad> cstr_, array<ad> logN_, array<ad> histLogF_, array<ad> lfs_, int y_) : dat(dat_), conf(conf_), par(par_), cFleets(cFleets_), recruit(recruit_), cstr(cstr_), logN(logN_), historicalLogF(histLogF_), logitFseason(lfs_), y(y_) {};

    ad getTarget(const vector<ad>& nlf){
      vector<ad> newLogF = nlf;
      int compareLag = CppAD::Integer(cstr.settings(0));
      // New F values
      vector<ad> fleetLogFbar = getFleetLogFbar(dat, conf, cFleets, newLogF, cstr.Amin, cstr.Amax);
      ad logFbar = logspace_sum(fleetLogFbar);
	
      // Predict year y and y+1, add one year
      // logisticFseason is already updated at this point;
      array<ad> hLogF2 = historicalLogF;
      hLogF2.col(y) = newLogF;
      hLogF2.col(y+1) = newLogF;
      MortalitySet<ad> mort(dat, conf, par, hLogF2, logitFseason);
      array<ad> logN2 = logN;
      logN2.col(y) = predNFun(dat,conf,par,logN2,hLogF2,recruit,mort,y);
      logN2.col(y+1) = predNFun(dat,conf,par,logN2,hLogF2,recruit,mort,y+1);

      // Add constraint
      if(cstr.cstr == ConstraintType::Constrain_Fbar){
	vector<ad> lastLogF = historicalLogF.col(y-compareLag);   
	// Previous F values
	vector<ad> lastFleetLogFbar = getFleetLogFbar(dat, conf, cFleets, lastLogF, cstr.Amin, cstr.Amax);
	ad lastLogFbar = logspace_sum(lastFleetLogFbar);

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
			       );
	  return trgt;
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
			       );
	  return trgt;
	}

      }else if(cstr.cstr == ConstraintType::Constrain_Catch){
	
	ad trgt = cstr.target;
	RELATIVE_CONSTRAINTS(
			     // Last year
			     ad logCL = getFleetCatch(dat, conf, logN2, hLogF2, mort, y-compareLag, cstr.Amin, cstr.Amax, cstr.fleet);
			     trgt += logCL;
			     ,
			     // Total
			     ad logCL = getFleetCatch(dat, conf, logN2, hLogF2, mort, y, cstr.Amin, cstr.Amax, -1);
			     trgt += logCL;
			     ,
			     // Fleet
			     ad logCL = getFleetCatch(dat, conf, logN2, hLogF2, mort, y, cstr.Amin, cstr.Amax, cstr.relative);
			     trgt += logCL;
			     );	
	return trgt;

      }else if(cstr.cstr == ConstraintType::Constrain_SSB){
	SAM_ASSERT(cstr.relative <= -2, "SSB constraints can only be relative to last year");
	//ad logB = getSSB(dat,conf, cFleets, recruit, logN, newLogF, historicalLogF, y, cstr.Amin, cstr.Amax, cstr.relative == (-2));
	ad trgt = cstr.target;
	if(cstr.relative == (-2))
	  trgt -= getLogSSB(dat, conf, logN2, hLogF2, mort, y-compareLag, cstr.Amin, cstr.Amax);
	return trgt;

      }else if(cstr.cstr == ConstraintType::Constrain_ERB){
	SAM_ASSERT(cstr.relative <= -2, "ERB constraints can only be relative to last year");
	//ad logB = getSSB(dat,conf, cFleets, recruit, logN, newLogF, historicalLogF, y, cstr.Amin, cstr.Amax, cstr.relative == (-2));
	ad trgt = cstr.target;
	if(cstr.relative == (-2))
	  trgt -= getLogERB(dat, conf,par, logN2, hLogF2, mort, y-compareLag, cstr.Amin, cstr.Amax);
	return trgt;
	  
	
      }else if(cstr.cstr == ConstraintType::Constrain_TSB){
	SAM_ASSERT(cstr.relative <= -2, "TSB constraints can only be relative to last year");
	ad trgt = cstr.target;
	if(cstr.relative == (-2))
	  trgt -= getLogTSB(dat, conf, logN2, hLogF2, mort, y-compareLag, cstr.Amin, cstr.Amax);
	return trgt;
	  
      }else if(cstr.cstr == ConstraintType::Constrain_Landing){

	ad trgt = cstr.target;
	RELATIVE_CONSTRAINTS(
			     // Last year
			     ad logLL = getFleetLanding(dat, conf, logN2, hLogF2, mort, y-compareLag, cstr.Amin, cstr.Amax, cstr.fleet);
			     trgt += logLL;
			     ,
			     // Total
			     ad logLL = getFleetLanding(dat, conf, logN2, hLogF2, mort, y, cstr.Amin, cstr.Amax, -1);
			     trgt += logLL;
			     ,
			     // Fleet
			     ad logLL = getFleetLanding(dat, conf, logN2, hLogF2, mort, y, cstr.Amin, cstr.Amax, cstr.relative);
			     trgt += logLL;
			     )	
	  return trgt;

      }else if(cstr.cstr == ConstraintType::Constrain_KeepRelF){
	SAM_ASSERT(cstr.fleet >= 0,"Keep relative Fbar can only be used for fleets, not total");
	SAM_ASSERT(cstr.relative >= 0 && cstr.fleet != cstr.relative, "Keep relative can only be used relative to other fleets in same year");
	vector<ad> lastLogF = historicalLogF.col(y-compareLag);   
	// Previous F values
	vector<ad> lastFleetLogFbar = getFleetLogFbar(dat, conf, cFleets, lastLogF, cstr.Amin, cstr.Amax);	
	return (lastFleetLogFbar(cstr.fleet) - lastFleetLogFbar(cstr.relative));

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
	   0: Indicator of SSB year (-1,0,1)
	   1: Indicator of biomass type (TSB,SSB,...)
	   2: Indicator of target Type (F, Catch, Landing,...)
	   3: Amin for TSB/SSB
	   4: Amax for TSB/SSB
	   5: F_origin for HCR
	   6: F_cap for HCR
	   7: B_origin for HCR
	   8: B_cap for HCR
	   9: B_trigger for HCR
	   10: Soft cap?
	*/
	int biomassType = CppAD::Integer(cstr.settings(1));
	//int targetType = CppAD::Integer(cstr.settings(2));
	int bioA0 = CppAD::Integer(cstr.settings(3));
	int bioA1 = CppAD::Integer(cstr.settings(4));
	// Get current biomass
	ad BforHCR = 0.0;
	ad Btrigger = exp(cstr.target);
	if(biomassType == 0){	// SSB
	  BforHCR = exp(getLogSSB(dat, conf, logN2, hLogF2, mort, y-compareLag, bioA0, bioA1));
	}else if(biomassType ==1){ // TSB
	  BforHCR = exp(getLogTSB(dat, conf, logN2, hLogF2, mort, y-compareLag, bioA0, bioA1));
	}else if(biomassType == 2){ // % of B0
	  BforHCR = exp(getLogSSB(dat, conf, logN2, hLogF2, mort, y-compareLag, bioA0, bioA1));
	  ad b0 = B0_i(dat,conf,par,hLogF2,y-compareLag,0,100);
	  Btrigger *= b0;
	  // }else if(biomassType == 3){ // % of Bmsy
	  //   BforHCR = exp(getLogSSB(dat, conf, logN2, hLogF2, mort, y-compareLag, bioA0, bioA1));
	}else if(biomassType == 3){ // ERB
	  BforHCR = exp(getLogERB(dat, conf, par, logN2, hLogF2, mort, y-compareLag, bioA0, bioA1));
	}else{
	  Rf_error("Wrong biomass type for HCR");
	}
	// Get optimization target
	vector<ad> hcrConf = cstr.settings.segment(4,6); // Start to early and overwrite
	// Insert target
	hcrConf(0) = Btrigger;
	ad trgt = hcr(BforHCR, hcrConf);
	return trgt;
      }else{
	Rf_error("Constraint type not implemented");
      }
      return 0.0;
    }
    
    ad operator()(const vector<ad>& nlf){
      vector<ad> logN1eps = nlf.segment(0,logN.rows());
      vector<ad> logN2eps = nlf.segment(logN.rows(),logN.rows());
      vector<ad> newLogF = nlf.segment(2*logN.rows(), historicalLogF.rows());
      // Return transformation of logF used to set constraints
      	if(cstr.cstr == ConstraintType::Constrain_NONE){
  	  return 0.0;
  	}
  	int compareLag = CppAD::Integer(cstr.settings(0));
  	// New F values
  	vector<ad> fleetLogFbar = getFleetLogFbar(dat, conf, cFleets, newLogF, cstr.Amin, cstr.Amax);
  	ad logFbar = logspace_sum(fleetLogFbar);
	
  	// Predict year y and y+1, add one year
  	// logisticFseason is already updated at this point;
  	array<ad> hLogF2 = historicalLogF;
  	hLogF2.col(y) = newLogF;
  	hLogF2.col(y+1) = newLogF;
  	MortalitySet<ad> mort(dat, conf, par, hLogF2, logitFseason);
  	array<ad> logN2 = logN;
  	logN2.col(y) = predNFun(dat,conf,par,logN2,hLogF2,recruit,mort,y) + logN1eps;
  	logN2.col(y+1) = predNFun(dat,conf,par,logN2,hLogF2,recruit,mort,y+1) + logN2eps;
	
  	// Add constraint
  	if(cstr.cstr == ConstraintType::Constrain_Fbar){
  	  // Previous F values
  	  
  	  // ad lastLogFbar = logspace_sum(lastFleetLogFbar);

  	  if(cstr.fleet == (-1)){	// Total F
	     return logFbar;
  	  }else{
	    vector<ad> lastLogF = historicalLogF.col(y-compareLag);   
  	    vector<ad> lastFleetLogFbar = getFleetLogFbar(dat, conf, cFleets, lastLogF, cstr.Amin, cstr.Amax);
	    return fleetLogFbar(cstr.fleet);
  	  }

  	}else if(cstr.cstr == ConstraintType::Constrain_Catch){

  	  ad logC = getFleetCatch(dat, conf, logN2, hLogF2, mort, y, cstr.Amin, cstr.Amax, cstr.fleet);
  	  return logC;
	  
  	}else if(cstr.cstr == ConstraintType::Constrain_SSB){
  	  SAM_ASSERT(cstr.relative <= -2, "SSB constraints can only be relative to last year");
  	  //ad logB = getSSB(dat,conf, cFleets, recruit, logN, newLogF, historicalLogF, y, cstr.Amin, cstr.Amax, cstr.relative == (-2));
  	  ad logB = getLogSSB(dat, conf, logN2, hLogF2, mort, y+1, cstr.Amin, cstr.Amax);
  	  return logB;
	  
  	}else if(cstr.cstr == ConstraintType::Constrain_TSB){
  	  SAM_ASSERT(cstr.relative <= -2, "TSB constraints can only be relative to last year");
  	  ad logB = getLogTSB(dat, conf, logN2, hLogF2, mort, y+1, cstr.Amin, cstr.Amax);
  	  return logB;
	  
  	}else if(cstr.cstr == ConstraintType::Constrain_Landing){

  	  ad logL = getFleetLanding(dat, conf, logN2, hLogF2, mort, y, cstr.Amin, cstr.Amax, cstr.fleet);
  	  return logL;

  	}else if(cstr.cstr == ConstraintType::Constrain_KeepRelF){

	  SAM_ASSERT(cstr.fleet >= 0,"Keep relative Fbar can only be used for fleets, not total");
	  SAM_ASSERT(cstr.relative >= 0 && cstr.fleet != cstr.relative, "Keep relative can only be used relative to other fleets in same year");	  
  	  return (fleetLogFbar(cstr.fleet) - fleetLogFbar(cstr.relative));

	}else if(cstr.cstr == ConstraintType::Constrain_HCR){
	  int targetType = CppAD::Integer(cstr.settings(2));
  	  if(targetType == 0){	  // F target
  	    if(cstr.fleet == (-1)){	// Total
  	      return logFbar;
  	    }else{		// Fleet
	      return fleetLogFbar(cstr.fleet);
  	    }
  	  }else if(targetType == 1){ // Catch target
  	    ad logC = getFleetCatch(dat, conf, logN2, hLogF2, mort, y, cstr.Amin, cstr.Amax, cstr.fleet);
	    return logC;
  	  }else if(targetType == 2){ // Landing target
  	    ad logL = getFleetLanding(dat, conf, logN2, hLogF2, mort, y, cstr.Amin, cstr.Amax, cstr.fleet);
	    return logL;
  	  }else{
  	    Rf_error("Unknown target type for HCR");
  	  }
  	}else{
  	  Rf_error("Constraint type not implemented");
  	}
  	return 0.0;
    }
    
  };

  
  struct ForecastF : NewtonFunctor {
    dataSet<ad> dat;
    confSet conf;
    paraSet<ad> par;
    vector<int> cFleets;
    Recruitment<ad> recruit;
    FConstraintList<ad> cstrs;
    array<ad> logN;
    array<ad> historicalLogF;
    array<ad> logitFseason;
    int y;
    matrix<ad> fvar;
    matrix<ad> nvar;
    vector<ad> logSel;
    
    ForecastF(dataSet<ad> dat_, confSet conf_, paraSet<ad> par_, vector<int> cFleets_, Recruitment<ad> recruit_, FConstraintList<ad> cstrs_, array<ad> logN_, array<ad> histLogF_, array<ad> lfs_, int y_, matrix<ad> fvar_, matrix<ad> nvar_, vector<ad> logSel_) : dat(dat_), conf(conf_), par(par_), cFleets(cFleets_), recruit(recruit_), cstrs(cstrs_), logN(logN_), historicalLogF(histLogF_), logitFseason(lfs_), y(y_), fvar(fvar_), nvar(nvar_), logSel(logSel_) {};
    

    ad operator()(const vector<ad>& logFs){
      ad kappa = 0.0;
      //
      // matrix<ad> llogF = toFleetMatrix(dat, conf, lastLogF, (vector<ad>)(logFs * ad(0.0)));
      // vector<ad> llfs = (vector<ad>)(lastLogF - lastLFB);
      // matrix<ad> logF = toFleetMatrix(dat, conf, lastLogF, logFs);

      vector<ad> newLogF = logSel;
   
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

	// Functor
	ForecastF_Fun F(dat, conf, par, cFleets, recruit, cstr, logN, historicalLogF, logitFseason, y);
	// Target
	vector<ad> nlf(2*logN.rows() + historicalLogF.rows());
	nlf.setZero();
	nlf.segment(2*logN.rows(), historicalLogF.rows()) = newLogF;
	ad trgt = F.getTarget(newLogF);	  
	// Transformed F
	ad transF = F(nlf);

	// Stochastic correction
	ad v = 0.0;
	if(cstr.useNonLinearityCorrection){
	  matrix<ad> H = autodiff::hessian(F, nlf);
	  matrix<ad> Sigma(2*logN.rows() + historicalLogF.rows(),2*logN.rows() + historicalLogF.rows());
	  Sigma.setZero();
	  Sigma.block(0,0,logN.rows(),logN.rows()) = nvar;
	  Sigma.block(logN.rows(),logN.rows(),logN.rows(),logN.rows()) = nvar;
	  Sigma.block(2*logN.rows(),2*logN.rows(),historicalLogF.rows(),historicalLogF.rows()) = fvar; 
	  //ad cx1 = (newLogF*(vector<ad>(H*newLogF))).sum();
	  ad tr = (H.vec() * Sigma.vec()).sum();
	  v = 0.5 * tr; // (tr + cx1);
	}
	ad tmp = transF + v - trgt;
	kappa += tmp * tmp;
	  
      }
	// int compareLag = CppAD::Integer(cstr.settings(0));
	// // New F values
	// vector<ad> fleetLogFbar = getFleetLogFbar(dat, conf, cFleets, newLogF, cstr.Amin, cstr.Amax);
	// ad logFbar = logspace_sum(fleetLogFbar);
	
	// // Predict year y and y+1, add one year
	// // logisticFseason is already updated at this point;
	// array<ad> hLogF2 = historicalLogF;
	// hLogF2.col(y) = newLogF;
	// hLogF2.col(y+1) = newLogF;
	// MortalitySet<ad> mort(dat, conf, par, hLogF2, logitFseason);
	// array<ad> logN2 = logN;
	// logN2.col(y) = predNFun(dat,conf,par,logN2,hLogF2,recruit,mort,y);
	// logN2.col(y+1) = predNFun(dat,conf,par,logN2,hLogF2,recruit,mort,y+1);

	// // Add constraint
	// if(cstr.cstr == ConstraintType::Constrain_Fbar){
	//   vector<ad> lastLogF = historicalLogF.col(y-compareLag);   
	//   // Previous F values
	//   vector<ad> lastFleetLogFbar = getFleetLogFbar(dat, conf, cFleets, lastLogF, cstr.Amin, cstr.Amax);
	//   ad lastLogFbar = logspace_sum(lastFleetLogFbar);

	//   if(cstr.fleet == (-1)){	// Total F
	//     ad trgt = cstr.target;
	//     RELATIVE_CONSTRAINTS(
	// 			 // Last year
	// 			 trgt += lastLogFbar;
	// 			 ,
	// 			 // Total
	// 			 Rf_error("A total F constraint can only be relative to last year");
	// 			 ,
	// 			 // Fleet
	// 			 Rf_error("A total F constraint can only be relative to last year");
	// 			 )
	//     ad tmp = logFbar - trgt;
	//     kappa += tmp * tmp;
	//   }else{
	//     ad trgt = cstr.target;
	//     RELATIVE_CONSTRAINTS(
	// 			 // Last year
	// 			 trgt += lastFleetLogFbar(cstr.fleet);
	// 			 ,
	// 			 // Total
	// 			 trgt += logFbar;
	// 			 ,
	// 			 // Fleet
	// 			 trgt += fleetLogFbar(cstr.relative);
	// 			 )
	//       ad tmp = fleetLogFbar(cstr.fleet) - trgt;
	//     kappa += tmp * tmp;
	//   }

	// }else if(cstr.cstr == ConstraintType::Constrain_Catch){

	//   ad logC = getFleetCatch(dat, conf, logN2, hLogF2, mort, y, cstr.Amin, cstr.Amax, cstr.fleet);
	//   ad trgt = cstr.target;
	//   RELATIVE_CONSTRAINTS(
	// 		       // Last year
	// 		       ad logCL = getFleetCatch(dat, conf, logN2, hLogF2, mort, y-compareLag, cstr.Amin, cstr.Amax, cstr.fleet);
	// 		       trgt += logCL;
	// 		       ,
	// 		       // Total
	// 		       ad logCL = getFleetCatch(dat, conf, logN2, hLogF2, mort, y, cstr.Amin, cstr.Amax, -1);
	// 		       trgt += logCL;
	// 		       ,
	// 		       // Fleet
	// 		       ad logCL = getFleetCatch(dat, conf, logN2, hLogF2, mort, y, cstr.Amin, cstr.Amax, cstr.relative);
	// 		       trgt += logCL;
	// 		       )	
	//   ad tmp = logC - trgt;
	//   kappa += tmp * tmp;

	// }else if(cstr.cstr == ConstraintType::Constrain_SSB){
	//   SAM_ASSERT(cstr.relative <= -2, "SSB constraints can only be relative to last year");
	//   //ad logB = getSSB(dat,conf, cFleets, recruit, logN, newLogF, historicalLogF, y, cstr.Amin, cstr.Amax, cstr.relative == (-2));
	//   ad logB = getLogSSB(dat, conf, logN2, hLogF2, mort, y+1, cstr.Amin, cstr.Amax);
	//   if(cstr.relative == (-2))
	//     logB -= getLogSSB(dat, conf, logN2, hLogF2, mort, y-compareLag, cstr.Amin, cstr.Amax);
	//   ad tmp = logB - cstr.target;
	//   kappa += tmp * tmp;
	  
	// }else if(cstr.cstr == ConstraintType::Constrain_TSB){
	//   SAM_ASSERT(cstr.relative <= -2, "TSB constraints can only be relative to last year");
	//   ad logB = getLogTSB(dat, conf, logN2, hLogF2, mort, y+1, cstr.Amin, cstr.Amax);
	//   if(cstr.relative == (-2))
	//     logB -= getLogTSB(dat, conf, logN2, hLogF2, mort, y-compareLag, cstr.Amin, cstr.Amax);
	//   ad tmp = logB - cstr.target;
	//   kappa += tmp * tmp;
	  
	// }else if(cstr.cstr == ConstraintType::Constrain_Landing){

	//   ad logL = getFleetLanding(dat, conf, logN2, hLogF2, mort, y, cstr.Amin, cstr.Amax, cstr.fleet);
	//   ad trgt = cstr.target;
	//   RELATIVE_CONSTRAINTS(
	// 		       // Last year
	// 		       ad logLL = getFleetLanding(dat, conf, logN2, hLogF2, mort, y-compareLag, cstr.Amin, cstr.Amax, cstr.fleet);
	// 		       trgt += logLL;
	// 		       ,
	// 		       // Total
	// 		       ad logLL = getFleetLanding(dat, conf, logN2, hLogF2, mort, y, cstr.Amin, cstr.Amax, -1);
	// 		       trgt += logLL;
	// 		       ,
	// 		       // Fleet
	// 		       ad logLL = getFleetLanding(dat, conf, logN2, hLogF2, mort, y, cstr.Amin, cstr.Amax, cstr.relative);
	// 		       trgt += logLL;
	// 		       )	
	//   ad tmp = logL - trgt;
	//   kappa += tmp * tmp;

	// }else if(cstr.cstr == ConstraintType::Constrain_KeepRelF){
	//   vector<ad> lastLogF = historicalLogF.col(y-compareLag);	  
	//   // Previous F values
	//   vector<ad> lastFleetLogFbar = getFleetLogFbar(dat, conf, cFleets, lastLogF, cstr.Amin, cstr.Amax);
	//   SAM_ASSERT(cstr.fleet >= 0,"Keep relative Fbar can only be used for fleets, not total");
	//   SAM_ASSERT(cstr.relative >= 0 && cstr.fleet != cstr.relative, "Keep relative can only be used relative to other fleets in same year");
	//   ad tmp = (lastFleetLogFbar(cstr.fleet) - lastFleetLogFbar(cstr.relative)) - (fleetLogFbar(cstr.fleet) - fleetLogFbar(cstr.relative));
	//   kappa += tmp * tmp;

	// }else if(cstr.cstr == ConstraintType::Constrain_HCR){
	//   // Give error if the relative flag makes it this far:
	//   RELATIVE_CONSTRAINTS(
	// 		       // Last year
	// 		       Rf_error("A HCR target cannot be relative");
	// 		       ,
	// 		       // Total
	// 		       Rf_error("A HCR target cannot be relative");
	// 		       ,
	// 		       // Fleet
	// 		       Rf_error("A HCR target cannot be relative");
	// 		       );
  
	//   /* cstr.settings should be:
	//      0: Indicator of SSB year (-1,0,1)
	//      1: Indicator of biomass type (TSB,SSB,...)
	//      2: Indicator of target Type (F, Catch, Landing,...)
	//      3: Amin for TSB/SSB
	//      4: Amax for TSB/SSB
	//      5: F_origin for HCR
	//      6: F_cap for HCR
	//      7: B_origin for HCR
	//      8: B_cap for HCR
	//      9: B_trigger for HCR
	//      10: Soft cap?
	//    */
	//   int biomassType = CppAD::Integer(cstr.settings(1));
	//   int targetType = CppAD::Integer(cstr.settings(2));
	//   int bioA0 = CppAD::Integer(cstr.settings(3));
	//   int bioA1 = CppAD::Integer(cstr.settings(4));
	//   // Get current biomass
	//   ad BforHCR = 0.0;
	//   ad Btrigger = exp(cstr.target);
	//   if(biomassType == 0){	// SSB
	//     BforHCR = exp(getLogSSB(dat, conf, logN2, hLogF2, mort, y-compareLag, bioA0, bioA1));
	//   }else if(biomassType ==1){ // TSB
	//     BforHCR = exp(getLogTSB(dat, conf, logN2, hLogF2, mort, y-compareLag, bioA0, bioA1));
	//   }else if(biomassType == 2){ // % of B0
	//     BforHCR = exp(getLogSSB(dat, conf, logN2, hLogF2, mort, y-compareLag, bioA0, bioA1));
	//     ad b0 = B0_i(dat,conf,par,hLogF2,y-compareLag,0,100);
	//     Btrigger *= b0;
	//   // }else if(biomassType == 3){ // % of Bmsy
	//   //   BforHCR = exp(getLogSSB(dat, conf, logN2, hLogF2, mort, y-compareLag, bioA0, bioA1));
	//   }else{
	//     Rf_error("Wrong biomass type for HCR");
	//   }
	//   // Get optimization target
	//   vector<ad> hcrConf = cstr.settings.segment(4,6); // Start to early and overwrite
	//   // Insert target
	//   hcrConf(0) = Btrigger;
	//   ad trgt = hcr(BforHCR, hcrConf);
	//   array<ad> logitFseason(0,historicalLogF.dim[1],0);
	  
	//   // kappa
	//   if(targetType == 0){	  // F target
	//     if(cstr.fleet == (-1)){	// Total
	//       ad tmp = logFbar - trgt;
	//       kappa += tmp * tmp;
	//     }else{		// Fleet
	//       ad tmp = fleetLogFbar(cstr.fleet) - trgt;
	//       kappa += tmp * tmp;
	//     }
	//   }else if(targetType == 1){ // Catch target
	//     ad logC = getFleetCatch(dat, conf, logN2, hLogF2, mort, y, cstr.Amin, cstr.Amax, cstr.fleet);
	//     ad tmp = logC - trgt;
	//     kappa += tmp * tmp;
	//   }else if(targetType == 2){ // Landing target
	//     ad logL = getFleetLanding(dat, conf, logN2, hLogF2, mort, y, cstr.Amin, cstr.Amax, cstr.fleet);
	//     ad tmp = logL - trgt;
	//     kappa += tmp * tmp;
	//   }else{
	//     Rf_error("Unknown target type for HCR");
	//   }
	// }else{
	//   Rf_error("Constraint type not implemented");
	// }      
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
			      array<Type>& logN,
			      array<Type>& logF,
			      array<Type>& logitFseason,
			      vector<int>& aveYears,
			      matrix<Type>& fvar,
			      matrix<Type> nvar,
			      vector<Type>& ICESrec,
			      vector<Type>& logSel,
			      int y,
			      newton::newton_config& cfg)SOURCE({
  
  paraSet<TMBad::ad_aug> parad(par);  
  Recruitment<TMBad::ad_aug> recruit = makeRecruitmentFunction(conf,parad);
  if(ICESrec.size() == 2){
    confSet c2(conf);
    recruit = makeICESrecruitment(TMBad::ad_aug(ICESrec(0)),TMBad::ad_aug(ICESrec(1)));
    nvar(0,0) = ICESrec(1) * ICESrec(1);
  }

  vector<int> cFleets = getCatchFleets(dat.fleetTypes);

  // Make data one longer
  dataSet<Type> newDat = dat;
  array<Type> logN2;
  array<Type> logF2;
  array<Type> logitFseason2;
  if(y >= logN.dim(1)-1){
    int nMYears = dat.propMat.dim(0);
    int nYears = 1;
    // propMat
    extendArray(newDat.propMat, nMYears, nYears, aveYears, par.meanLogitMO, conf.keyMatureMean, 1, true);
    // stockMeanWeight
    extendArray(newDat.stockMeanWeight, nMYears, nYears, aveYears, par.meanLogSW, conf.keyStockWeightMean, 0, true);
    // catchMeanWeight
    extendArray(newDat.catchMeanWeight, nMYears, nYears, aveYears, par.meanLogCW, conf.keyCatchWeightMean, 0, true);
    // natMor
    extendArray(newDat.natMor, nMYears, nYears, aveYears, par.meanLogNM, conf.keyMortalityMean, 0, true);
    // landFrac (No biopar process)
    extendArray(newDat.landFrac, nMYears, nYears, aveYears, true);
    // disMeanWeight (No biopar process)
    extendArray(newDat.disMeanWeight, nMYears, nYears, aveYears, true);
    // landMeanWeight (No biopar process)
    extendArray(newDat.landMeanWeight, nMYears, nYears, aveYears, true);
    // propF (No biopar process)
    extendArray(newDat.propF, nMYears, nYears, aveYears, true);
    // propM (No biopar process)
    extendArray(newDat.propM, nMYears, nYears, aveYears, true);
    newDat.noYears = nYears+nMYears;

  
    // Make logN, histLogF, and logitFseason one longer in the year direction
    logN2 = array<Type>(logN.dim(0),logN.dim(1)+nYears);
    for(int i = 0; i < logN2.dim(0); ++i){
      for(int j = 0; j < logN2.dim(1); ++j){
	int j2 = std::min(j,logN.dim(1)-1);
	logN2(i,j) = logN(i,j2);
      }
    }
    logF2 = array<Type>(logF.dim(0),logF.dim(1)+nYears);
    for(int i = 0; i < logF2.dim(0); ++i){
      for(int j = 0; j < logF2.dim(1); ++j){
	int j2 = std::min(j,logF.dim(1)-1);
	logF2(i,j) = logF(i,j2);
      }
      logF2(i,logF2.dim(0)-1) = logF(i,logF.dim(0)-1);
    }
    logitFseason2 = array<Type>(logitFseason.dim(0),logitFseason.dim(1)+nYears,logitFseason.dim(2));
    for(int i = 0; i < logitFseason2.dim(0); ++i){
      for(int k = 0; k < logitFseason2.dim(2); ++k){
	for(int j = 0; j < logitFseason2.dim(1); ++j){
	  int j2 = std::min(j,logitFseason.dim(1)-1);
	  logitFseason2(i,j,k) = logitFseason(i,j2,k);
	}      
      }
    }
  }else{
    logN2 = logN;
    logF2 = logF;
    logitFseason2 = logitFseason;
  }

  // Should be deleted by NewtonWrapper
  std::shared_ptr<ConstrainCalculations::ForecastF> p_fc0 = std::make_shared<ConstrainCalculations::ForecastF>(newDat,conf,par,cFleets,recruit,cstrs,logN2,logF2, logitFseason2, y, fvar, nvar, logSel);
  std::shared_ptr<NewtonFunctor> p_fc = std::dynamic_pointer_cast<NewtonFunctor>(p_fc0);
 
  // vector<double> s0(cFleets.size());
  // s0.setConstant(0);
  // vector<Type> start(tryStart(fc,s0,200));
  vector<Type> start(cFleets.size());
  start.setConstant(log(0.1 / cFleets.size()));
  
  vector<Type> res = SAM_Newton(p_fc, start, cfg);
  vector<Type> newLogF = logSel; //logF.col(y-1);// - lastLogFbar;

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

SAM_SPECIALIZATION(vector<double> calculateNewFVec(dataSet<double>&, confSet&, paraSet<double>&, FConstraintList<double>&, array<double>&,array<double>&,array<double>&, vector<int>&, matrix<double>&, matrix<double>, vector<double>&, vector<double>&, int, newton::newton_config&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> calculateNewFVec(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, FConstraintList<TMBad::ad_aug>&, array<TMBad::ad_aug>&,array<TMBad::ad_aug>&,array<TMBad::ad_aug>&,vector<int>&,matrix<TMBad::ad_aug>&, matrix<TMBad::ad_aug>, vector<TMBad::ad_aug>&, vector<TMBad::ad_aug>&, int, newton::newton_config&));

