SAM_DEPENDS(define)
SAM_DEPENDS(toF)
SAM_DEPENDS(recruitment)
SAM_DEPENDS(incidence)
SAM_DEPENDS(derived)
SAM_DEPENDS(logspace)
SAM_DEPENDS(hcr)
SAM_DEPENDS(predf)
SAM_DEPENDS(extend_array)


HEADER(
       template <class Type>
struct forecastSet {

  enum FModelType {
		   asFModel = 0,
		   // useFscale,
		   // useFval,
		   // useCatchval,
		   // useNextssb,
		   // useLandval,
		   constrainedF = 1,
		   fastFixedF = 2,
		   findMSY = 3,
		   HCR = 4,
		   customHCR = 99
  };
  enum recModelType {
		     asRecModel,
		     useIID
  };

  enum FSdTimeScaleModel {
			  rwScale,
			  oneScale,
			  zeroScale,
			  fixedDeviation
  };
  
	 int nYears;
	 int preYears;
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
	 vector<Type> Fdeviation;
	 matrix<Type> FdeviationCov;
	 matrix<Type> FEstCov;
	 int useModelLastN;

  matrix<Type> forecastCalculatedMedian;
  vector<Type> forecastCalculatedLogSdCorrection;
  vector<Type> sel;
  vector<Type> selFull;
  Type initialFbar;
	 matrix<Type> nvar;

	 matrix<Type> cumEpsilon;
  
	 void calculateForecast(array<Type>& logF, array<Type>& logN, array<Type>& logitFseason, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, Recruitment<Type>& recruit, MortalitySet<Type>& mort);
	 void updateForecast(int i, array<Type>& logF, array<Type>& logN, array<Type>& logitFseason, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, Recruitment<Type>& recruit, MortalitySet<Type>& mort, int sim);
  
  forecastSet();
  
  forecastSet(SEXP x);

  template<class T>
  inline forecastSet(const forecastSet<T>& x) : nYears(x.nYears),
						preYears(x.preYears),
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
						Fdeviation(x.Fdeviation),
						FdeviationCov(x.FdeviationCov),
						FEstCov(x.FEstCov),
						useModelLastN(x.useModelLastN),
						forecastCalculatedMedian(x.forecastCalculatedMedian),
						forecastCalculatedLogSdCorrection(x.forecastCalculatedLogSdCorrection),
						sel(x.sel),
						selFull(x.selFull),
						initialFbar(x.initialFbar),
						cumEpsilon(x.cumEpsilon) {}

}
       );

SOURCE(
	 template <class Type>
	 void forecastSet<Type>::calculateForecast(array<Type>& logF, array<Type>& logN, array<Type>& logitFseason, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, Recruitment<Type>& recruit, MortalitySet<Type>& mort){
	   if(nYears == 0){ // no forecast
	     return;
	   }

	   for(int i = 0; i < nYears; ++i){
	     updateForecast(i, logF, logN, logitFseason, dat, conf, par, recruit, mort, false);
	   }
	   return;
	 }
	 )

SOURCE(
       template <class Type>
       void forecastSet<Type>::updateForecast(int i, array<Type>& logF, array<Type>& logN, array<Type>& logitFseason, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, Recruitment<Type>& recruit, MortalitySet<Type>& mort, int sim){
	 if(nYears == 0){ // no forecast
	   return;
	 }
	 int indx = preYears + i;//i starts at 0,forecastYear.size() - nYears + i;    
	 Type y = forecastYear(indx);    
	 Type lastSSB = ssbi(dat, conf, logN, logF, mort, indx-1);
	 // Assuming F before spawning is zero
	 Type thisSSB = ssbi(dat, conf, logN, logF, mort, indx);
	 // Type calcF = 0.0;
	 vector<Type> logSelUse = log(sel);
	 // Type logFbarEpsilon = 0.0;
	 vector<Type> lastFepsilon(logF.rows());
	 lastFepsilon.setZero();
	 if(selectivity.size() == 0){
	   Type lastFbarRealized = fbari(dat, conf, logF, indx-1);
	   logSelUse = logF.col(indx-1) - log(lastFbarRealized);
	 }
	 if(i > 0){
	   cumEpsilon.col(i) = (vector<Type>)cumEpsilon.col(i-1) + (vector<Type>)logF.col(indx-1) - (vector<Type>)forecastCalculatedMedian.col(i-1);
	 }else{
	   lastFepsilon = (vector<Type>)logF.col(indx-1) - (vector<Type>)logF.col(indx-2);
	   cumEpsilon.col(i) = lastFepsilon;
	 }
	 
	   
	 // vector<Type> lastShortLogF(logF.rows());
	 // if(i == 0){
	 //   lastShortLogF = log(initialFbar) + logSelUse;
	 // }else{
	 //   lastShortLogF = forecastCalculatedMedian.col(i-1);
	 // }
     
	 // vector<Type> lastFullLogF(logN.rows());
	 // for(int j = 0; j < logN.rows(); ++j){
	 //   lastFullLogF(j) = R_NegInf;
	 //   for(int f = 0; f < conf.keyLogFsta.dim(0); ++f)
	 //     if(conf.keyLogFsta(f,j)>(-1)){
	 //       lastFullLogF(j) = logspace_add_SAM(lastFullLogF(j),lastShortLogF(conf.keyLogFsta(f,j)));
	 //     }
	 // }

	   
	 vector<Type> muF = get_fmu(dat,conf,par,logF);
	 vector<Type> rhoF = get_frho(dat,conf,par,logF);
	 matrix<Type> fvar;
	 if(fsdTimeScaleModel(i) == fixedDeviation && FdeviationCov.cols() > 0){
	   fvar = FdeviationCov;
	 }else if(FEstCov.cols() > 0){
	   fvar = FEstCov;
	 }else{
	   fvar = get_fvar(dat, conf, par, logF);
	 }
	 if(fsdTimeScaleModel(i) == rwScale){
	   fvar *= (i+2);
	   // for(int j = 0; j < logSelUse.size(); ++j)
	   //   logSelUse(j) += cumEpsilon(j,i);
	 }

	 vector<Type> bnd(logF.dim[0]); bnd.setZero();
	 if(!isNA(conf.boundFbar)){
	   Type lastFbar = fbari(dat, conf, logF, indx-1);
	   for(int q = 0; q < bnd.size(); ++q){
	     // Type v = (logF(q,indx-1)-logF(q,0)) / (34.19951893 * conf.rwBoundLogF);
	     Type v = lastFbar / ( 1.668100537 * conf.boundFbar);
	     bnd(q) = 10.0 * v * v * v * v * v * v * v * v * v;
	   }
	 }
	 vector<Type> predF = muF + rhoF * (logF.col(indx-1) - muF) - bnd;

	 
	   vector<Type> ICESrec;
	   if(recModel(i) == useIID){
	     ICESrec = vector<Type>(2);
	     ICESrec(0) = logRecruitmentMedian;
	     ICESrec(1) = sqrt(logRecruitmentVar);
	   }else{
	     ICESrec = vector<Type>(0);
	   }
    
	   // Calculate CV correction
	   switch(FModel(i)) { // target is not used. F is a random walk
	   case asFModel:
	     if(fsdTimeScaleModel(i) != oneScale)
	       Rf_warning("F time scale model is ignored when the F model is used for forecasting.");
	     forecastCalculatedLogSdCorrection(i) = 1.0;
	     forecastCalculatedMedian.col(i) = predF;// logF.col(indx - 1);
	     break;
	   case constrainedF:
	     forecastCalculatedMedian.col(i) = calculateNewFVec(dat,conf,par,constraints(i),logN,logF,logitFseason, aveYears, fvar, nvar, ICESrec, logSelUse, indx, cfg);    
	     break;
	   case fastFixedF:
	     forecastCalculatedMedian.col(i) = constraints(i)(0).target + logSelUse;
	     break;
	   case findMSY:
	     forecastCalculatedMedian.col(i) = par.logFScaleMSY + log(initialFbar) + logSelUse; // 
	     break;
	   case HCR:
	     // TODO: Allow hcr forecast based on predicted SSB
	     if(hcrCurrentSSB){
	       SAM_ASSERT(sum((vector<Type>)dat.propF.matrix().row(indx)) == 0, "currentSSB in HCR can only be used when propF is zero.");
	       forecastCalculatedMedian.col(i) = hcr(thisSSB, hcrConf) + logSelUse; //      
	     }else{
	       forecastCalculatedMedian.col(i) = hcr(lastSSB, hcrConf) + logSelUse; //
	     }
	     break;
	   case customHCR:
	     Rf_error("Forecast type not implemented");
	     break;
	   default:
	     Rf_error("Forecast type not implemented");
	   }

	   switch(fsdTimeScaleModel(i)) {
	   case rwScale:
	     if(!sim){
	       forecastCalculatedLogSdCorrection(i) = sqrt(y+1);
	     }else{
	       forecastCalculatedLogSdCorrection(i) = 1.0;
	       for(int j = 0; j < forecastCalculatedMedian.rows(); ++j)
		 forecastCalculatedMedian(j,i) += cumEpsilon(j,i);
	     }
	     break;
	   case oneScale:
	     forecastCalculatedLogSdCorrection(i) = 1.0;
	     break;
	   case zeroScale:
	     forecastCalculatedLogSdCorrection(i) = 1e-6;
	     break;
	   case fixedDeviation:
	     forecastCalculatedMedian.col(i) = (vector<Type>)forecastCalculatedMedian.col(i) + Fdeviation;
	     break;
	   default:
	     Rf_error("F time scale model not implemented");
	   }

	   
	   return;
	 }
	 )

SOURCE(
	 template <class Type>
	 forecastSet<Type>::forecastSet() : nYears(),
	 preYears(),
		  nCatchAverageYears(),
		  aveYears(),
		  forecastYear(),
		  FModel(),
		  // target(0),
		  constraints(),
		  cfg(),
		  selectivity(),
		  recModel(),
		  logRecruitmentMedian(),
		  logRecruitmentVar(),
		  fsdTimeScaleModel(),
		  simFlag(),
		  hcrConf(),
		  hcrCurrentSSB(),
	 Fdeviation(),
	 FdeviationCov(),
	 FEstCov(),
	 useModelLastN(),
	 forecastCalculatedMedian(),
		  forecastCalculatedLogSdCorrection(),
		  sel(),
		  selFull(),
	 initialFbar(),
	 cumEpsilon() {}
	 )  

SOURCE(
	 template <class Type>
	 forecastSet<Type>::forecastSet(SEXP x) {
	   // If nYears is NULL or 0; only set nYears to 0 -> no forecast
	   if(Rf_isNull(getListElement(x,"nYears")) ||
	      (int)*REAL(getListElement(x,"nYears")) == 0){
	     preYears = 0;
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
	     Fdeviation = vector<Type>(0);
	     FdeviationCov = matrix<Type>();
	     FEstCov = matrix<Type>();
	     useModelLastN = 1;
	     forecastCalculatedMedian = matrix<Type>(0,0);
	     forecastCalculatedLogSdCorrection = vector<Type>(0);	     
	     sel = vector<Type>(0);
	     selFull = vector<Type>(0);
	     initialFbar = 0;
	     cumEpsilon = matrix<Type>(0,0);
	   }else{
	     using tmbutils::asArray;
	     nYears = (int)*REAL(getListElement(x,"nYears"));
	     preYears = (int)*REAL(getListElement(x,"preYears"));
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
	     Fdeviation = asVector<Type>(getListElement(x,"Fdeviation"));
	     FdeviationCov = asMatrix<Type>(getListElement(x,"FdeviationCov"));
	     FEstCov = asMatrix<Type>(getListElement(x,"FEstCov"));
	     useModelLastN = Rf_asInteger(getListElement(x,"useModelLastN"));
	   }
	 }
	 );

SAM_SPECIALIZATION(struct forecastSet<double>);
SAM_SPECIALIZATION(struct forecastSet<TMBad::ad_aug>);







template <class Type>
void prepareForForecast(forecastSet<Type>& forecast, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, array<Type>& logF, array<Type>& logN, Recruitment<Type>& recruit, objective_function<Type>* of)SOURCE({
    if(forecast.nYears == 0)
      return;
    int nFYears = forecast.nYears - (dat.noYears - forecast.preYears);
    int nMYears = dat.noYears;
    vector<int> aveYears = forecast.aveYears;
  // propMat 
  extendArray(dat.propMat, nMYears, nFYears, aveYears, par.meanLogitMO, conf.keyMatureMean, 1, true);
  REPORT_F(dat.propMat, of);
  // stockMeanWeight
  extendArray(dat.stockMeanWeight, nMYears, nFYears, aveYears, par.meanLogSW, conf.keyStockWeightMean, 0, true);
  REPORT_F(dat.stockMeanWeight, of);
  // catchMeanWeight
  extendArray(dat.catchMeanWeight, nMYears, nFYears, aveYears, par.meanLogCW, conf.keyCatchWeightMean, 0, true);
  REPORT_F(dat.catchMeanWeight, of);
  // natMor
  extendArray(dat.natMor, nMYears, nFYears, aveYears, par.meanLogNM, conf.keyMortalityMean, 0, true);
  REPORT_F(dat.natMor, of);
  // landFrac (No biopar process)
  extendArray(dat.landFrac, nMYears, nFYears, aveYears, true);
  REPORT_F(dat.landFrac, of);
  // disMeanWeight (No biopar process)
  extendArray(dat.disMeanWeight, nMYears, nFYears, aveYears, true);
  REPORT_F(dat.disMeanWeight, of);
  // landMeanWeight (No biopar process)
  extendArray(dat.landMeanWeight, nMYears, nFYears, aveYears, true);
  REPORT_F(dat.landMeanWeight, of);
  // propF (No biopar process)
  extendArray(dat.propF, nMYears, nFYears, aveYears, true);
  REPORT_F(dat.propF, of);
  // propM (No biopar process)
  extendArray(dat.propM, nMYears, nFYears, aveYears, true);
  REPORT_F(dat.propM, of);
  
  // Prepare forecastCalculated...
  forecast.forecastCalculatedMedian = matrix<Type>(logF.rows(), forecast.nYears);
  forecast.forecastCalculatedMedian.setZero();
  forecast.forecastCalculatedLogSdCorrection = vector<Type>(forecast.nYears);
  forecast.forecastCalculatedLogSdCorrection.setZero();

  forecast.cumEpsilon = matrix<Type>(logF.rows(), forecast.nYears);
  forecast.cumEpsilon.setZero();

  forecast.nvar = matrix<Type>(logN.rows(), logN.rows());
  forecast.nvar.setZero();
  for(int i = 0; i < forecast.nvar.rows(); ++i)
    forecast.nvar(i,i) = exp(2.0 * par.logSdLogN(conf.keyVarLogN(i)));

  // Calculate initial Fbar
  int fbarFirst = conf.fbarRange(0) - conf.minAge;
  int fbarLast = conf.fbarRange(1) - conf.minAge;
  // forecast.initialFbar = 0.0;
  // array<Type> totF = totFFun(dat,conf, logF);
  // for(int a = fbarFirst; a <= fbarLast; ++a){  
  //   //forecast.initialFbar += exp(logF(conf.keyLogFsta(0,a),forecast.forecastYear.size() - nFYears - 1));
  //   forecast.initialFbar += totF(a,forecast.forecastYear.size() - nFYears - 1);
  // }
  // forecast.initialFbar /= Type(fbarLast - fbarFirst + 1);
  forecast.initialFbar = fbari(dat,conf,logF,forecast.forecastYear.size() - forecast.nYears - 1,false);
  
  // Calculate selectivity
  forecast.sel = vector<Type>(logF.rows());
  forecast.sel.setZero();
  forecast.selFull = vector<Type>(logN.rows());
  forecast.selFull.setZero();

  // Correct input selectivity to have Fbar == 1
  Type inputFbar = 0.0;
  if(forecast.selectivity.size() > 0){
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
	  Type v = exp(logF(conf.keyLogFsta(f,j),forecast.forecastYear.size() - forecast.nYears - 1)) / forecast.initialFbar;
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
  })

SAM_SPECIALIZATION(void prepareForForecast(forecastSet<double>&, dataSet<double>&, confSet&, paraSet<double>&, array<double>&, array<double>&, Recruitment<double>&, objective_function<double>*));
SAM_SPECIALIZATION(void prepareForForecast(forecastSet<TMBad::ad_aug>&, dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, Recruitment<TMBad::ad_aug>&, objective_function<TMBad::ad_aug>*));

