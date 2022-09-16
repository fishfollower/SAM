// Forward declarations
template <class Type>
Type hcr(Type ssb, vector<Type> hcrConf);
template<class Type>
Type ssb2F_quick(Type ssbval, vector<Type> logFlast, dataSet<Type> dat, confSet conf, paraSet<Type> par, array<Type> logF, array<Type> logN, int i, Type rec_mean);
template<class Type>
Type landing2F_quick(Type landingval, vector<Type> lastF, vector<Type> M, vector<Type> N, vector<Type> w, vector<Type> frac);
template<class Type>
Type catch2F_quick(Type catchval, vector<Type> lastF, vector<Type> M, vector<Type> N, vector<Type> w);
template<class Type>
Type landing2F(Type landingval, vector<Type> lastF, vector<Type> M, vector<Type> N, vector<Type> w, vector<Type> frac);
template<class Type>
Type catch2F(Type catchval, vector<Type> lastF, vector<Type> M, vector<Type> N, vector<Type> w);
// 

template <class Type>
struct forecastSet {

  enum FModelType {
		   asFModel,
		   // useFscale,
		   // useFval,
		   // useCatchval,
		   // useNextssb,
		   // useLandval,
		   constrainedF,
		   findMSY,
		   HCR,
		   customHCR
  };
  enum recModelType {
		     asRecModel,
		     useIID
  };

  enum FSdTimeScaleModel {
			  rwScale,
			  oneScale,
			  zeroScale		      
  };
  
  int nYears;
  int nCatchAverageYears;
  vector<int> aveYears;
  vector<Type> forecastYear;
  vector<FModelType> FModel;
  // vector<Type> target;
  vector<FConstraintList<Type> > constraints;
  newton::newton_config cfg;
  vector<Type> selectivity;
  vector<recModelType> recModel;
  Type logRecruitmentMedian;
  Type logRecruitmentVar;
  vector<FSdTimeScaleModel> fsdTimeScaleModel;
  vector<int> simFlag;
  vector<Type> hcrConf;
  int hcrCurrentSSB;

  matrix<Type> forecastCalculatedMedian;
  vector<Type> forecastCalculatedLogSdCorrection;
  vector<Type> sel;
  vector<Type> selFull;
  Type initialFbar;
  
  void calculateForecast(array<Type>& logF, array<Type>& logN, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, Recruitment<Type>& recruit, MortalitySet<Type>& mort){
  if(nYears == 0){ // no forecast
    return;
  }

  for(int i = 0; i < nYears; ++i){
    updateForecast(i, logF, logN, dat, conf, par, recruit, mort);
  }
  return;

}

  void updateForecast(int i, array<Type>& logF, array<Type>& logN, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, Recruitment<Type>& recruit, MortalitySet<Type>& mort){
  if(nYears == 0){ // no forecast
    return;
  }

    int indx = forecastYear.size() - nYears + i;    
    Type y = forecastYear(indx);    
    Type lastSSB = ssbi(dat, conf, logN, logF, mort, indx-1);
    // Assuming F before spawning is zero
    Type thisSSB = ssbi(dat, conf, logN, logF, mort, indx);
    Type calcF = 0.0;

    vector<Type> lastShortLogF(logF.rows());
    if(i == 0){
      lastShortLogF = log(initialFbar) + log(sel);
    }else{
      lastShortLogF = forecastCalculatedMedian.col(i-1);
    }
      
    vector<Type> lastFullLogF(logN.rows());
    for(int j = 0; j < logN.rows(); ++j){
      lastFullLogF(j) = R_NegInf;
      for(int f = 0; f < conf.keyLogFsta.dim(0); ++f)
	if(conf.keyLogFsta(f,j)>(-1)){
	  lastFullLogF(j) = logspace_add_SAM(lastFullLogF(j),lastShortLogF(conf.keyLogFsta(f,j)));
	}
    }


    switch(fsdTimeScaleModel(i)) {
    case rwScale:
      forecastCalculatedLogSdCorrection(i) = sqrt(y);
      break;
    case oneScale:
      forecastCalculatedLogSdCorrection(i) = 1.0;
      break;
    case zeroScale:
      forecastCalculatedLogSdCorrection(i) = 1e-6;
      break;
    default:
      Rf_error("F time scale model not implemented");
    }

    // Calculate CV correction
    switch(FModel(i)) { // target is not used. F is a random walk
    case asFModel:
      if(fsdTimeScaleModel(i) != oneScale)
	Rf_warning("F time scale model is ignored when the F model is used for forecasting.");
      forecastCalculatedLogSdCorrection(i) = 1.0;
      forecastCalculatedMedian.col(i) = logF.col(indx - 1);
      break;
    case constrainedF:
      forecastCalculatedMedian.col(i) = calculateNewFVec(dat,conf,par,constraints(i),lastShortLogF,logN,indx, cfg);    
      break;
    // case useFscale: // target is an F scale of previous F
    //   if(i == 0){
    // 	forecastCalculatedMedian.col(i) = log(target(i)) + log(initialFbar) + log(sel);
    //   }else{
    // 	forecastCalculatedMedian.col(i) = log(target(i)) + (vector<Type>)forecastCalculatedMedian.col(i-1);
    //   }
    //   break;
    // case useFval: // target is F value	
    //   forecastCalculatedMedian.col(i) = log(target(i)) + log(sel);
    //   break;
    // case useCatchval: // target is a catch value in weight
    //   calcF = catch2F_quick((Type)target(i),
    // 			      exp(lastFullLogF),
    // 			      (vector<Type>)dat.natMor.matrix().row(indx),
    // 			      exp((vector<Type>)logN.col(indx)),
    // 			      (vector<Type>)dat.catchMeanWeight.matrix().row(indx));
    //   if(i == 0){
    // 	forecastCalculatedMedian.col(i) = log(calcF) + log(initialFbar) + log(sel);
    //   }else{
    // 	forecastCalculatedMedian.col(i) = log(calcF) + (vector<Type>)forecastCalculatedMedian.col(i-1);
    //   }
    //   break;
    // case useNextssb:
    //   calcF = ssb2F_quick((Type)target(i),
    // 			  log(sel), //lastShortLogF,
    // 			  dat,
    // 			  conf,
    // 			  par,
    // 			  logF,
    // 			  logN,
    // 			  mort,
    // 			  indx,
    // 			  logRecruitmentMedian);
    //   forecastCalculatedMedian.col(i) = log(calcF) + log(sel);
    //   break;
    // case useLandval:
    // 	calcF = landing2F_quick((Type)target(i),
    // 				exp(lastFullLogF),
    // 				(vector<Type>)dat.natMor.matrix().row(indx),
    // 				exp((vector<Type>)logN.col(indx)),
    // 				(vector<Type>)dat.landMeanWeight.matrix().row(indx),
    // 				(vector<Type>)dat.landFrac.matrix().row(indx));
    //   if(i == 0){
    // 	forecastCalculatedMedian.col(i) = log(calcF) + log(initialFbar) + log(sel);
    //   }else{
    // 	forecastCalculatedMedian.col(i) = log(calcF) + (vector<Type>)forecastCalculatedMedian.col(i-1);
    //   }
    //   break;
    case findMSY:
      forecastCalculatedMedian.col(i) = par.logFScaleMSY + log(initialFbar) + log(sel); // 
      break;
    case HCR:
      // TODO: Allow hcr forecast based on predicted SSB
      if(hcrCurrentSSB){
	if(sum((vector<Type>)dat.propF.matrix().row(indx)) > 0)
	  Rf_error("currentSSB in HCR can only be used when propF is zero.");
	forecastCalculatedMedian.col(i) = hcr(thisSSB, hcrConf) + log(sel); //      
      }else{
	forecastCalculatedMedian.col(i) = hcr(lastSSB, hcrConf) + log(sel); //
      }
      break;
    case customHCR:
      Rf_error("Forecast type not implemented");
      break;
    default:
      Rf_error("Forecast type not implemented");
    }
    return;
}
  
  forecastSet() : nYears(0),
			nCatchAverageYears(0),
			aveYears(0),
			forecastYear(0),
		  FModel(0),
			// target(0),
		  constraints(0),
		  cfg(),
		  selectivity(0),
			recModel(0),
			logRecruitmentMedian(0),
			logRecruitmentVar(0),
			fsdTimeScaleModel(0),
			simFlag(0),
			hcrConf(0),
			hcrCurrentSSB(0),
			forecastCalculatedMedian(0,0),
			forecastCalculatedLogSdCorrection(0),
			sel(0),
			selFull(0),
			initialFbar(0) {};
  
  forecastSet(SEXP x) // : nYears(Rf_asInteger(getListElement(x,"nYears"))),
		      // 	nCatchAverageYears(0),
		      // 	aveYears(0),
		      // 	forecastYear(0),
		      // 	FModel(static_cast<FModelType>(0)),
		      // 	target(0),
		      // 	selectivity(0),
		      // 	recModel(static_cast<recModelType>(0)),
		      // 	logRecruitmentMedian(0),
		      // 	logRecruitmentVar(0),
		      // 	fsdTimeScaleModel(static_cast<FSdTimeScaleModel>(0)),
		      // 	simFlag(0),
		      // 	hcrConf(0),
		      // 	hcrCurrentSSB(0),
		      // 	forecastCalculatedMedian(0,0),
		      // 	forecastCalculatedLogSdCorrection(0),
		      // 	sel(0),
		      // 	selFull(0),
		      // 	initialFbar(0)
  {
    // If nYears is NULL or 0; only set nYears to 0 -> no forecast
    if(Rf_isNull(getListElement(x,"nYears")) ||
       (int)*REAL(getListElement(x,"nYears")) == 0){
      nYears = 0;
      nCatchAverageYears = 0;
      aveYears = vector<int>(0);
      forecastYear = vector<Type>(0);
      FModel = vector<FModelType>(0);
      // target = vector<Type>(0);
      constraints = vector<FConstraintList<Type> >(0);
      cfg = newton::newton_config();
      selectivity = vector<Type>(0);
      recModel = vector<recModelType>(0);
      logRecruitmentMedian = 0;
      logRecruitmentVar = 0;
      fsdTimeScaleModel = vector<FSdTimeScaleModel>(0);
      simFlag = vector<int>(0);
      hcrConf = vector<Type>(0);
      hcrCurrentSSB = 0;
      forecastCalculatedMedian = matrix<Type>(0,0);
      forecastCalculatedLogSdCorrection = vector<Type>(0);
      sel = vector<Type>(0);
      selFull = vector<Type>(0);
      initialFbar = 0;
    }else{
      using tmbutils::asArray;
      nYears = (int)*REAL(getListElement(x,"nYears"));
      nCatchAverageYears = (int)*REAL(getListElement(x,"nCatchAverageYears"));
      aveYears = asVector<int>(getListElement(x,"aveYears"));
      forecastYear = asVector<Type>(getListElement(x,"forecastYear"));
      vector<int> FModelTmp = asVector<int>(getListElement(x,"FModel"));
      FModel = vector<FModelType>(FModelTmp.size());
      for(int i = 0; i < FModel.size(); ++i)
	FModel(i) = static_cast<FModelType>(FModelTmp(i));
      // Fval = asVector<Type>(getListElement(x,"Fval"));
      // Fscale = asVector<Type>(getListElement(x,"Fscale"));
      // target = asVector<Type>(getListElement(x,"target"));
      SEXP ctmp = PROTECT(getListElement(x,"constraints"));
      constraints = vector<FConstraintList<Type> >(Rf_length(ctmp));
      for(int i = 0; i < Rf_length(ctmp); ++i){
	SEXP tmp = VECTOR_ELT(ctmp, i);
	constraints(i) = FConstraintList<Type>(tmp);
      }
      UNPROTECT(1);
      cfg = newton::newton_config(getListElement(x,"cfg"));
      selectivity = asVector<Type>(getListElement(x,"selectivity"));
      vector<int> recModelTmp = asVector<int>(getListElement(x,"recModel"));
      recModel = vector<recModelType>(recModelTmp.size());
      for(int i = 0; i < recModel.size(); ++i)
	recModel(i) = static_cast<recModelType>(recModelTmp(i));
      logRecruitmentMedian = (Type)*REAL(getListElement(x,"logRecruitmentMedian"));
      logRecruitmentVar = (Type)*REAL(getListElement(x,"logRecruitmentVar"));
      vector<int> fsdTimeScaleModelTmp = asVector<int>(getListElement(x,"fsdTimeScaleModel"));
      fsdTimeScaleModel = vector<FSdTimeScaleModel>(fsdTimeScaleModelTmp.size());
      for(int i = 0; i < fsdTimeScaleModel.size(); ++i)
	fsdTimeScaleModel(i) = static_cast<FSdTimeScaleModel>(fsdTimeScaleModelTmp(i));
      simFlag = asVector<int>(getListElement(x,"simFlag"));
      hcrConf = asVector<Type>(getListElement(x,"hcrConf"));
      hcrCurrentSSB = Rf_asInteger(getListElement(x,"hcrCurrentSSB"));
    }
  };

  template<class T>
  forecastSet(const forecastSet<T>& x) : nYears(x.nYears),
					 nCatchAverageYears(x.nCatchAverageYears),
					 aveYears(x.aveYears),
					 forecastYear(x.forecastYear),
					 FModel(x.FModel),
					 // target(x.target),
					 constraints(x.constraints),
					 cfg(x.cfg),
					 selectivity(x.selectivity),
					 recModel(x.recModel),
					 logRecruitmentMedian(x.logRecruitmentMedian),
					 logRecruitmentVar(x.logRecruitmentVar),
					 fsdTimeScaleModel(x.fsdTimeScaleModel),
					 simFlag(x.simFlag),
					 hcrConf(x.hcrConf),
					 hcrCurrentSSB(x.hcrCurrentSSB),
					 forecastCalculatedMedian(x-forecastCalculatedMedian),
					 forecastCalculatedLogSdCorrection(x.forecastCalculatedLogSdCorrection),
					 sel(x.sel),
					 selFull(x.selFull),
					 initialFbar(x.initialFbar) {}

};



template <class Type>
void extendArray_2D(array<Type>& x, int nModelYears, int nForecastYears, vector<int> aveYears, bool keepModelYears = true){
  vector<int> dim = x.dim;
  array<Type> tmp((int)keepModelYears*nModelYears+nForecastYears,dim(1));
  tmp.setZero();
  vector<Type> ave(dim(1));
  ave.setZero();
  Type nave = aveYears.size();

  // Calculate average
  for(int i = 0; i < aveYears.size(); ++i){
    if(aveYears(i) < dim(0)){    
      for(int j = 0; j < dim(1); ++j){
	ave(j) += x(aveYears(i), j);
      }
    }else{
      nave -= 1.0;
    }
  }

  if(nave == 0)
    Rf_error("ave.years does not cover the data period.");
  
  ave /= nave;

  // Insert values in tmp
  for(int i = 0; i < tmp.dim(0); ++i){
    for(int j = 0; j < tmp.dim(1); ++j){
      if(keepModelYears && i < dim(0)){ // Take value from x
	tmp(i,j) = x(i,j);
      }else{ // Take value from ave
	tmp(i,j) = ave(j);
      }
    }
  }
  // Overwrite x
  // NOTE: x must be resized first, otherwise x=tmp will not work.
  x.initZeroArray(tmp.dim);
  x = tmp;
  return;
}


template <class Type>
void extendArray_3D(array<Type>& x, int nModelYears, int nForecastYears, vector<int> aveYears, bool keepModelYears = true){
  vector<int> dim = x.dim;
  array<Type> tmp((int)keepModelYears*nModelYears+nForecastYears,dim(1),dim(2));
  tmp.setZero();
  array<Type> ave(dim(1),dim(2));
  ave.setZero();
  Type nave = aveYears.size();

  // Calculate average
  for(int i = 0; i < aveYears.size(); ++i){
    if(aveYears(i) < dim(0)){
      for(int k = 0; k < dim(2); ++k){
	for(int j = 0; j < dim(1); ++j){
	  ave(j,k) += x(aveYears(i), j, k);
	}
      }
    }else{
      nave -= 1.0;
    }
  }

  if(nave == 0)
    Rf_error("ave.years does not cover the data period.");
  
  ave /= nave;

  // Insert values in tmp
  for(int i = 0; i < tmp.dim(0); ++i){
    for(int j = 0; j < tmp.dim(1); ++j){
      for(int k = 0; k < tmp.dim(2); ++k){
	if(keepModelYears && i < dim(0)){ // Take value from x
	  tmp(i,j,k) = x(i,j,k);
	}else{ // Take value from ave
	  tmp(i,j,k) = ave(j,k);
	}
      }
    }
  }
  // Overwrite x
  // NOTE: x must be resized first, otherwise x=tmp will not work.
  x.initZeroArray(tmp.dim);
  x = tmp;
  return;
}

template <class Type>
void extendArray(array<Type>& x, int nModelYears, int nForecastYears, vector<int> aveYears, bool keepModelYears = true){
  if(x.dim.size() == 2){
    extendArray_2D(x, nModelYears, nForecastYears, aveYears, keepModelYears);
  }else if(x.dim.size() == 3){
    extendArray_3D(x, nModelYears, nForecastYears, aveYears, keepModelYears);
  }else{
    Rf_error("extendArray is only implemented for arrays of dimension 2 and 3");
  }
  return;
}



template <class Type>
void prepareForForecast(forecastSet<Type>& forecast, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, array<Type>& logN, Recruitment<Type>& recruit){
  int nFYears = forecast.nYears;
  int nMYears = dat.noYears;
  if(nFYears == 0)
    return;
  vector<int> aveYears = forecast.aveYears;
  // propMat
  extendArray(dat.propMat, nMYears, nFYears, aveYears);
  // stockMeanWeight
  extendArray(dat.stockMeanWeight, nMYears, nFYears, aveYears);
  // catchMeanWeight
  extendArray(dat.catchMeanWeight, nMYears, nFYears, aveYears);
  // natMor
  extendArray(dat.natMor, nMYears, nFYears, aveYears);
  // landFrac
  extendArray(dat.landFrac, nMYears, nFYears, aveYears);
  // disMeanWeight
  extendArray(dat.disMeanWeight, nMYears, nFYears, aveYears);
  // landMeanWeight
  extendArray(dat.landMeanWeight, nMYears, nFYears, aveYears);
  // propF
  extendArray(dat.propF, nMYears, nFYears, aveYears);
  // propM
  extendArray(dat.propM, nMYears, nFYears, aveYears);

  // Prepare forecastCalculated...
  forecast.forecastCalculatedMedian = matrix<Type>(logF.rows(), nFYears);
  forecast.forecastCalculatedMedian.setZero();
  forecast.forecastCalculatedLogSdCorrection = vector<Type>(nFYears);
  forecast.forecastCalculatedLogSdCorrection.setZero();

  // Calculate initial Fbar
  int fbarFirst = conf.fbarRange(0) - conf.minAge;
  int fbarLast = conf.fbarRange(1) - conf.minAge;
  forecast.initialFbar = 0.0;
  array<Type> totF = totFFun(conf, logF);
  for(int a = fbarFirst; a <= fbarLast; ++a){  
    //forecast.initialFbar += exp(logF(conf.keyLogFsta(0,a),forecast.forecastYear.size() - nFYears - 1));
    forecast.initialFbar += totF(a,forecast.forecastYear.size() - nFYears - 1);
  }
  forecast.initialFbar /= Type(fbarLast - fbarFirst + 1);

  // Calculate selectivity
  forecast.sel = vector<Type>(logF.rows());
  forecast.sel.setZero();
  forecast.selFull = vector<Type>(logN.rows());
  forecast.selFull.setZero();

  // Correct input selectivity to have Fbar == 1
  Type inputFbar = 0.0;
  if(forecast.selectivity.size() > 0){
    // Rcout << "Using custom selectivity!\n";
    for(int a = fbarFirst; a <= fbarLast; ++a){  
      if(forecast.selectivity.size() == logF.rows()){
	for(int f = 0; f < conf.keyLogFsta.dim(0); ++f)
	  if(conf.keyLogFsta(f,a) > (-1))
	    inputFbar += forecast.selectivity(conf.keyLogFsta(f,a));
      }else if(forecast.selectivity.size() == logN.rows()){
	inputFbar += forecast.selectivity(a);
      }else{
	Rf_error("Wrong size of selectivity. Must match logF or logN array.");
      }
    }
    inputFbar /= Type(fbarLast - fbarFirst + 1);
    if(inputFbar != 1.0){
      // Rf_warning("The input selectivity was re-scaled to have Fbar equal to one.");
      forecast.selectivity /= inputFbar;
    }
  }

  for(int j = 0; j < logN.rows(); ++j){
    for(int f = 0; f < conf.keyLogFsta.dim(0); ++f)
      if(conf.keyLogFsta(f,j)>(-1)){
	if(forecast.selectivity.size() == 0){
	  Type v = exp(logF(conf.keyLogFsta(f,j),forecast.forecastYear.size() - nFYears - 1)) / forecast.initialFbar;
	  forecast.selFull(j) += v;
	  forecast.sel(conf.keyLogFsta(f,j)) = v;
	}else if(forecast.selectivity.size() == logF.rows()){
	  forecast.selFull(j) += forecast.selectivity(conf.keyLogFsta(f,j));
	  forecast.sel(conf.keyLogFsta(f,j)) = forecast.selectivity(conf.keyLogFsta(f,j));
	}else if(forecast.selectivity.size() == logN.rows()){
	  forecast.selFull(j) = forecast.selectivity(j);
	  forecast.sel(conf.keyLogFsta(f,j)) = forecast.selFull(j);
	}else{
	  Rf_error("Wrong size of selectivity. Must match logF or logN array.");
	}
      }
  }
  
  return;  
}

