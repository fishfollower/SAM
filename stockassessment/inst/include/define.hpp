template <class Type>
struct forecastSet;

template <class Type>
struct dataSet;

struct confSet;

template <class Type>
struct paraSet;



#define REPORT_F(name,F)					\
  if(isDouble<Type>::value && F->current_parallel_region<0) {	\
    defineVar(install(#name),					\
	      PROTECT(asSEXP(name)),F->report);			\
    UNPROTECT(1);						\
  }


//This function returns a vector with matrices based on a list of matrices from R
template<class Type>
struct listMatrixFromR : vector<matrix<Type> > {

  listMatrixFromR() : vector<matrix<Type> >() {};
  listMatrixFromR(int n) : vector<matrix<Type> >(n) {};
  listMatrixFromR(SEXP x){ 
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asMatrix<Type>(sm);
    }
  }

  template<class T>
  listMatrixFromR<T> cast() const {
    int n = (*this).size();
    listMatrixFromR<T> d(n);
    for(int i = 0; i < n; ++i){
      matrix<T> tmp = (*this)(i).template cast<T>();
      d(i) = tmp;
    }
    return d;
  }

  listMatrixFromR<Type>& operator=(const listMatrixFromR<Type>& rhs) {
    (*this).resize(rhs.size());
    for(int i = 0; i < rhs.size(); ++i){
      matrix<Type> tmp = rhs(i);
      (*this)(i) = tmp;
    }
    return *this;
  }

  
};


#define ADREPORT_F(name,F) F->reportvector.push(name,#name);

#define SIMULATE_F(F)				\
  if(isDouble<Type>::value && F->do_simulate)

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

bool isNAINT(int x){
  return NA_INTEGER==x;
}


template <class Type>
struct forecastSet {

  enum FModelType {
		   asFModel,
		   useFscale,
		   useFval,
		   useCatchval,
		   useNextssb,
		   useLandval,
		   findMSY
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
  vector<Type> target;
  vector<Type> selectivity;
  matrix<Type> forecastCalculatedMedian;
  vector<Type> forecastCalculatedLogSdCorrection;
  vector<recModelType> recModel;
  Type logRecruitmentMedian;
  Type logRecruitmentVar;
  vector<FSdTimeScaleModel> fsdTimeScaleModel;
  vector<int> simFlag;
  int uniroot;

  void calculateForecast(array<Type>& logF, array<Type>& logN, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par); // Defined after dataSet and confSet
  
  forecastSet() : nYears(0) {};
  
  forecastSet(SEXP x){
    // If nYears is NULL or 0; only set nYears to 0 -> no forecast
    if(Rf_isNull(getListElement(x,"nYears")) ||
       (int)*REAL(getListElement(x,"nYears")) == 0){
      nYears = 0;
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
      target = asVector<Type>(getListElement(x,"target"));
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
      uniroot = (int)*REAL(getListElement(x,"uniroot"));
  
    }
  };

  forecastSet<Type>& operator=(const forecastSet<Type>& rhs) {
    nYears = rhs.nYears;
    if(nYears == 0)
      return *this;
    nCatchAverageYears = rhs.nCatchAverageYears;
    aveYears = rhs.aveYears;
    forecastYear = rhs.forecastYear;
    FModel = rhs.FModel;
    target = rhs.target;
    selectivity = rhs.selectivity;
    forecastCalculatedMedian = rhs.forecastCalculatedMedian;
    forecastCalculatedLogSdCorrection = rhs.forecastCalculatedLogSdCorrection;
    recModel = rhs.recModel;
    logRecruitmentMedian = rhs.logRecruitmentMedian;
    logRecruitmentVar = rhs.logRecruitmentVar;
    fsdTimeScaleModel = rhs.fsdTimeScaleModel;
    simFlag = rhs.simFlag;
    uniroot = rhs.uniroot;
    return *this;
  }

  template<class T>
  forecastSet<T> cast() const {
    forecastSet<T> d;
    d.nYears = nYears;	// int
    if(nYears == 0)
      return d;
    d.nCatchAverageYears = nCatchAverageYears; // int
    d.aveYears = aveYears;	// <int>
    d.forecastYear = forecastYear.template cast<T>();
    // d.FModel = FModel; // <int>
    d.FModel = vector<typename forecastSet<T>::FModelType>(FModel.size());
    for(int i = 0; i < FModel.size(); ++i)
      d.FModel(i) = static_cast<typename forecastSet<T>::FModelType>((int)FModel(i));
    d.target = target.template cast<T>();
    d.selectivity = selectivity.template cast<T>();
    d.forecastCalculatedMedian = forecastCalculatedMedian.template cast<T>();
    d.forecastCalculatedLogSdCorrection = forecastCalculatedLogSdCorrection.template cast<T>();
    // d.recModel = recModel; // <int>
    d.recModel = vector<typename forecastSet<T>::recModelType>(recModel.size());
    for(int i = 0; i < recModel.size(); ++i)
      d.recModel(i) = static_cast<typename forecastSet<T>::recModelType>((int)recModel(i));
    d.logRecruitmentMedian = T(logRecruitmentMedian);
    d.logRecruitmentVar = T(logRecruitmentVar);
    // d.fsdTimeScaleModel = fsdTimeScaleModel; // <int>
    d.fsdTimeScaleModel = vector<typename forecastSet<T>::FSdTimeScaleModel>(fsdTimeScaleModel.size());
    for(int i = 0; i < fsdTimeScaleModel.size(); ++i)
      d.fsdTimeScaleModel(i) = static_cast<typename forecastSet<T>::FSdTimeScaleModel>((int)fsdTimeScaleModel(i));

    d.simFlag = simFlag; // <int>
    d.uniroot = uniroot;
    return d;    
  }

};

template<class Type>
struct referencepointSet {

  enum CatchType {
		  totalCatch,
		  landings,
		  discard
  };
  
  int nYears;
  vector<int> aveYears;
  vector<int> selYears;
  vector<Type> Fsequence;
  vector<Type> xPercent;
  CatchType catchType;

  referencepointSet() : nYears(0) {};
  
  referencepointSet(SEXP x){
    // If nYears is NULL or 0; only set nYears to 0 -> no forecast
    if(Rf_isNull(getListElement(x,"nYears")) ||
       (int)*REAL(getListElement(x,"nYears")) == 0){
      nYears = 0;
    }else{
      using tmbutils::asArray;
      nYears = (int)*REAL(getListElement(x,"nYears"));
      aveYears = asVector<int>(getListElement(x,"aveYears"));
      selYears = asVector<int>(getListElement(x,"selYears"));
      Fsequence = asVector<Type>(getListElement(x,"Fsequence"));
      xPercent = asVector<Type>(getListElement(x,"xPercent"));
      catchType = static_cast<CatchType>((int)*REAL(getListElement(x,"catchType")));
    }
  }
  
  referencepointSet<Type>& operator=(const referencepointSet<Type>& rhs) {
    nYears = rhs.nYears;
    if(nYears == 0)
      return *this;
    aveYears = rhs.aveYears;
    selYears = rhs.selYears;
    Fsequence = rhs.Fsequence;
    xPercent = rhs.xPercent;
    catchType = rhs.catchType;
    return *this;
  }

  template<class T>
  referencepointSet<T> cast() const {
    referencepointSet<T> d;
    d.nYears = nYears;
    if(nYears == 0)
      return d;
    d.aveYears = aveYears;
    d.selYears = selYears;
    d.Fsequence = Fsequence.template cast<T>();
    d.xPercent = xPercent.template cast<T>();
    d.catchType = static_cast<typename referencepointSet<T>::CatchType>((int)catchType);
    return d;    
  }
  

};


template <class Type>
struct dataSet{
  int noFleets;
  vector<int> fleetTypes; 
  vector<Type> sampleTimes;
  int noYears;
  vector<Type> years;
  vector<int> minAgePerFleet;
  vector<int> maxAgePerFleet;
  int nobs;
  array<int> idx1;
  array<int> idx2;
  array<int> idxCor;
  array<int> aux;
  vector<Type> logobs;
  vector<Type> weight;
  // data_indicator<vector<Type>,Type> keep;
  array<Type> propMat;
  array<Type> stockMeanWeight; 
  array<Type> catchMeanWeight;
  array<Type> natMor;
  array<Type> landFrac;
  array<Type> disMeanWeight;
  array<Type> landMeanWeight;
  array<Type> propF;
  array<Type> propM;
  listMatrixFromR<Type> corList;
  forecastSet<Type> forecast;
  referencepointSet<Type> referencepoint;

  dataSet() {};

  dataSet(SEXP x) {
    using tmbutils::asArray;
    noFleets = (int)*REAL(getListElement(x,"noFleets"));
    fleetTypes = asVector<int>(getListElement(x,"fleetTypes"));
    sampleTimes = asVector<Type>(getListElement(x,"sampleTimes"));
    noYears = (int)*REAL(getListElement(x,"noYears"));
    years = asVector<Type>(getListElement(x,"years"));
    minAgePerFleet = asVector<int>(getListElement(x,"minAgePerFleet"));
    maxAgePerFleet = asVector<int>(getListElement(x,"maxAgePerFleet"));
    nobs = (int)*REAL(getListElement(x,"nobs"));
    idx1 = asArray<int>(getListElement(x,"idx1"));
    idx2 = asArray<int>(getListElement(x,"idx2"));
    idxCor = asArray<int>(getListElement(x,"idxCor"));
    aux = asArray<int>(getListElement(x,"aux"));
    logobs = asVector<Type>(getListElement(x,"logobs"));
    weight = asVector<Type>(getListElement(x,"weight"));
    propMat = asArray<Type>(getListElement(x,"propMat"));
    stockMeanWeight = asArray<Type>(getListElement(x,"stockMeanWeight"));
    catchMeanWeight = asArray<Type>(getListElement(x,"catchMeanWeight"));
    natMor = asArray<Type>(getListElement(x,"natMor"));
    landFrac = asArray<Type>(getListElement(x,"landFrac"));
    disMeanWeight = asArray<Type>(getListElement(x,"disMeanWeight"));
    landMeanWeight = asArray<Type>(getListElement(x,"landMeanWeight"));
    propF = asArray<Type>(getListElement(x,"propF"));
    propM = asArray<Type>(getListElement(x,"propM"));
    corList = listMatrixFromR<Type>(getListElement(x,"corList"));
    forecast = forecastSet<Type>(getListElement(x,"forecast"));
    referencepoint = referencepointSet<Type>(getListElement(x,"referencepoint"));
  };

  dataSet<Type>& operator=(const dataSet<Type>& rhs) {
    noFleets = rhs.noFleets;
    fleetTypes = rhs.fleetTypes;
    sampleTimes = rhs.sampleTimes;
    noYears = rhs.noYears;
    years = rhs.years;
    minAgePerFleet = rhs.minAgePerFleet;
    maxAgePerFleet = rhs.maxAgePerFleet;
    nobs = rhs.nobs;
    idx1 = rhs.idx1;
    idx2 = rhs.idx2;
    idxCor = rhs.idxCor;
    aux = rhs.aux;
    logobs = rhs.logobs;
    weight = rhs.weight;
    // keep = rhs.keep;
    propMat = rhs.propMat;
    stockMeanWeight = rhs.stockMeanWeight; 
    catchMeanWeight = rhs.catchMeanWeight;
    natMor = rhs.natMor;
    landFrac = rhs.landFrac;
    disMeanWeight = rhs.disMeanWeight;
    landMeanWeight = rhs.landMeanWeight;
    propF = rhs.propF;
    propM = rhs.propM;
    corList = rhs.corList;
    forecast = rhs.forecast;
    referencepoint = rhs.referencepoint;
    return *this;
  };

  template<class T>
  dataSet<T> cast() const {
    dataSet<T> d;
    d.noFleets = noFleets; // int
    d.fleetTypes = fleetTypes; //<int> 
    d.sampleTimes = sampleTimes.template cast<T>();
    d.noYears = noYears;	// int
    d.years = years.template cast<T>();
    d.minAgePerFleet = minAgePerFleet;//<int> 
    d.maxAgePerFleet = maxAgePerFleet;//<int> 
    d.nobs = nobs;// int
    d.idx1 = idx1;//<int> 
    d.idx2 = idx2;//<int> 
    d.idxCor = idxCor;//<int> 
    d.aux = aux;//<int> 
    d.logobs = logobs.template cast<T>();
    d.weight = weight.template cast<T>();
    // keep = keep;
    // The array must be resized before copying
    d.propMat.initZeroArray(propMat.dim);
    d.propMat = propMat.template cast<T>();
    d.stockMeanWeight.initZeroArray(stockMeanWeight.dim);
    d.stockMeanWeight = stockMeanWeight.template cast<T>(); 
    d.catchMeanWeight.initZeroArray(catchMeanWeight.dim);
    d.catchMeanWeight = catchMeanWeight.template cast<T>();
    d.natMor.initZeroArray(natMor.dim);
    d.natMor = natMor.template cast<T>();
    d.landFrac.initZeroArray(landFrac.dim);
    d.landFrac = landFrac.template cast<T>();
    d.disMeanWeight.initZeroArray(disMeanWeight.dim);
    d.disMeanWeight = disMeanWeight.template cast<T>();
    d.landMeanWeight.initZeroArray(landMeanWeight.dim);
    d.landMeanWeight = landMeanWeight.template cast<T>();
    d.propF.initZeroArray(propF.dim);
    d.propF = propF.template cast<T>();
    d.propM.initZeroArray(propM.dim);
    d.propM = propM.template cast<T>();
    d.corList = corList.template cast<T>();
    d.forecast = forecast.template cast<T>();
    d.referencepoint = referencepoint.template cast<T>();
    return d;    
  }
  
};




template <class Type>
void extendArray(array<Type>& x, int nModelYears, int nForecastYears, vector<int> aveYears, bool keepModelYears = true){
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
void prepareForForecast(dataSet<Type>& dat){
  int nFYears = dat.forecast.nYears;
  int nMYears = dat.noYears;
  if(nFYears == 0)
    return;
  vector<int> aveYears = dat.forecast.aveYears;
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
  return;  
}





struct confSet{
  int minAge;
  int maxAge;
  vector<int> maxAgePlusGroup;
  array<int> keyLogFsta;
  int corFlag;
  array<int> keyLogFpar;
  array<int> keyQpow;
  array<int> keyVarF;
  vector<int> keyVarLogN; 
  array<int> keyVarObs;
  vector<int> obsCorStruct; 
  array<int> keyCorObs;
  int stockRecruitmentModelCode;
  vector<double> constRecBreaks;
  int noScaledYears;
  vector<int> keyScaledYears;
  matrix<int> keyParScaledYA;
  vector<int> fbarRange;
  vector<int> keyBiomassTreat;
  vector<int> simFlag; 
  int resFlag; 
  vector<int> obsLikelihoodFlag;
  int fixVarToWeight;
  double fracMixF;
  double fracMixN;
  vector<double> fracMixObs;
  
  confSet() {};

  confSet(SEXP x){
    using tmbutils::asArray;
    minAge = (int)*REAL(getListElement(x,"minAge"));
    maxAge = (int)*REAL(getListElement(x,"maxAge"));
    maxAgePlusGroup = asVector<int>(getListElement(x,"maxAgePlusGroup"));
    keyLogFsta = asArray<int>(getListElement(x,"keyLogFsta"));
    corFlag = (int)*REAL(getListElement(x,"corFlag"));
    keyLogFpar = asArray<int>(getListElement(x,"keyLogFpar"));
    keyQpow = asArray<int>(getListElement(x,"keyQpow"));
    keyVarF = asArray<int>(getListElement(x,"keyVarF"));
    keyVarLogN = asVector<int>(getListElement(x,"keyVarLogN"));
    keyVarObs = asArray<int>(getListElement(x,"keyVarObs"));
    obsCorStruct = asVector<int>(getListElement(x,"obsCorStruct"));
    keyCorObs = asArray<int>(getListElement(x,"keyCorObs"));
    stockRecruitmentModelCode = (int)*REAL(getListElement(x,"stockRecruitmentModelCode"));
    constRecBreaks = asVector<double>(getListElement(x,"constRecBreaks"));
    noScaledYears = (int)*REAL(getListElement(x,"noScaledYears"));
    keyScaledYears = asVector<int>(getListElement(x,"keyScaledYears"));
    keyParScaledYA = asMatrix<int>(getListElement(x,"keyParScaledYA"));
    fbarRange = asVector<int>(getListElement(x,"fbarRange"));
    keyBiomassTreat = asVector<int>(getListElement(x,"keyBiomassTreat"));
    simFlag = asVector<int>(getListElement(x,"simFlag"));
    resFlag = (int)*REAL(getListElement(x,"resFlag"));
    obsLikelihoodFlag = asVector<int>(getListElement(x,"obsLikelihoodFlag"));
    fixVarToWeight = (int)*REAL(getListElement(x,"fixVarToWeight"));
    fracMixF = (double)*REAL(getListElement(x,"fracMixF"));
    fracMixN = (double)*REAL(getListElement(x,"fracMixN"));
    fracMixObs = asVector<double>(getListElement(x,"fracMixObs"));
  };

  confSet& operator=(const confSet& rhs) {
    minAge = rhs.minAge;
    maxAge = rhs.maxAge; 
    maxAgePlusGroup = rhs.maxAgePlusGroup;
    keyLogFsta = rhs.keyLogFsta;
    keyLogFpar = rhs.keyLogFpar;
    corFlag = rhs.corFlag;
    keyQpow = rhs.keyQpow;
    keyVarF = rhs.keyVarF;
    keyVarLogN = rhs.keyVarLogN;
    keyVarObs = rhs.keyVarObs;
    obsCorStruct = rhs.obsCorStruct;
    keyCorObs = rhs.keyCorObs;
    stockRecruitmentModelCode = rhs.stockRecruitmentModelCode;
    constRecBreaks = rhs.constRecBreaks;
    noScaledYears = rhs.noScaledYears;
    keyScaledYears = rhs.keyScaledYears;
    keyParScaledYA = rhs.keyParScaledYA;
    fbarRange = rhs.fbarRange;
    keyBiomassTreat = rhs.keyBiomassTreat; 
    simFlag = rhs.simFlag;
    resFlag = rhs.resFlag;
    obsLikelihoodFlag = rhs.obsLikelihoodFlag;
    fixVarToWeight = rhs.fixVarToWeight;
    fracMixF = rhs.fracMixF;
    fracMixN = rhs.fracMixN;
    fracMixObs = rhs.fracMixObs;
    return *this;
  };
};

template <class Type>
struct paraSet{
  vector<Type> logFpar; 
  vector<Type> logQpow; 
  vector<Type> logSdLogFsta; 
  vector<Type> logSdLogN; 
  vector<Type> logSdLogObs;
  vector<Type> logSdLogTotalObs;
  vector<Type> transfIRARdist;
  vector<Type> sigmaObsParUS;
  vector<Type> rec_pars; 
  vector<Type> itrans_rho; 
  vector<Type> logScale;
  vector<Type> logitReleaseSurvival;
  vector<Type> logitRecapturePhi; 
  vector<Type> sepFalpha;   
  vector<Type> sepFlogitRho;   
  vector<Type> sepFlogSd;   

  Type logFScaleMSY;
  Type implicitFunctionDelta;
  Type logScaleFmsy;
  Type logScaleFmax;
  Type logScaleF01;
  Type logScaleFcrash;
  Type logScaleFext;
  vector<Type> logScaleFxPercent;
  Type logScaleFlim;

  paraSet() {};
  
  paraSet(SEXP x){
    logFpar = asVector<Type>(getListElement(x,"logFpar"));
    logQpow = asVector<Type>(getListElement(x,"logQpow"));
    logSdLogFsta = asVector<Type>(getListElement(x,"logSdLogFsta"));
    logSdLogN = asVector<Type>(getListElement(x,"logSdLogN"));
    logSdLogObs = asVector<Type>(getListElement(x,"logSdLogObs"));
    logSdLogTotalObs = asVector<Type>(getListElement(x,"logSdLogTotalObs"));
    transfIRARdist = asVector<Type>(getListElement(x,"transfIRARdist"));
    sigmaObsParUS = asVector<Type>(getListElement(x,"sigmaObsParUS"));
    rec_pars = asVector<Type>(getListElement(x,"rec_pars"));
    itrans_rho = asVector<Type>(getListElement(x,"itrans_rho"));
    logScale = asVector<Type>(getListElement(x,"logScale"));
    logitReleaseSurvival = asVector<Type>(getListElement(x,"logitReleaseSurvival"));
    logitRecapturePhi = asVector<Type>(getListElement(x,"logitRecapturePhi"));
    sepFalpha = asVector<Type>(getListElement(x,"sepFalpha"));
    sepFlogitRho = asVector<Type>(getListElement(x,"sepFlogitRho"));
    sepFlogSd = asVector<Type>(getListElement(x,"sepFlogSd"));
    logFScaleMSY = (Type)Rf_asReal(getListElement(x,"logFScaleMSY"));
    implicitFunctionDelta = (Type)Rf_asReal(getListElement(x,"implicitFunctionDelta"));
    logScaleFmsy = (Type)Rf_asReal(getListElement(x,"logScaleFmsy"));
    logScaleFmax = (Type)Rf_asReal(getListElement(x,"logScaleFmax"));
    logScaleF01 = (Type)Rf_asReal(getListElement(x,"logScaleF01"));
    logScaleFcrash = (Type)Rf_asReal(getListElement(x,"logScaleFcrash"));
    logScaleFext = (Type)Rf_asReal(getListElement(x,"logScaleFext"));
    logScaleFxPercent = asVector<Type>(getListElement(x,"logScaleFxPercent"));
    logScaleFlim = (Type)Rf_asReal(getListElement(x,"logScaleFlim"));
  }

  paraSet<Type>& operator=(const paraSet<Type>& rhs) {
    logFpar = rhs.logFpar; 
    logQpow = rhs.logQpow; 
    logSdLogFsta = rhs.logSdLogFsta; 
    logSdLogN = rhs.logSdLogN; 
    logSdLogObs = rhs.logSdLogObs;
    logSdLogTotalObs = rhs.logSdLogTotalObs;
    transfIRARdist = rhs.transfIRARdist;
    sigmaObsParUS = rhs.sigmaObsParUS;
    rec_pars = rhs.rec_pars; 
    itrans_rho = rhs.itrans_rho; 
    logScale = rhs.logScale;
    logitReleaseSurvival = rhs.logitReleaseSurvival;   
    logitRecapturePhi = rhs.logitRecapturePhi;
    sepFalpha = rhs.sepFalpha;
    sepFlogitRho = rhs.sepFlogitRho;
    sepFlogSd = rhs.sepFlogSd;
    logFScaleMSY = rhs.logFScaleMSY;
    implicitFunctionDelta = rhs.implicitFunctionDelta;
    logScaleFmsy = rhs.logScaleFmsy;
    logScaleFmax = rhs.logScaleFmax;
    logScaleF01 = rhs.logScaleF01;
    logScaleFcrash = rhs.logScaleFcrash;
    logScaleFext = rhs.logScaleFext;
    logScaleFxPercent = rhs.logScaleFxPercent;
    logScaleFlim = rhs.logScaleFlim;
    return *this;

  }

  template<class T>
  paraSet<T> cast() const {
    paraSet<T> d;
    d.logFpar = logFpar.template cast<T>(); 
    d.logQpow = logQpow.template cast<T>(); 
    d.logSdLogFsta = logSdLogFsta.template cast<T>(); 
    d.logSdLogN = logSdLogN.template cast<T>(); 
    d.logSdLogObs = logSdLogObs.template cast<T>();
    d.logSdLogTotalObs = logSdLogTotalObs.template cast<T>();
    d.transfIRARdist = transfIRARdist.template cast<T>();
    d.sigmaObsParUS = sigmaObsParUS.template cast<T>();
    d.rec_pars = rec_pars.template cast<T>(); 
    d.itrans_rho = itrans_rho.template cast<T>(); 
    d.logScale = logScale.template cast<T>();
    d.logitReleaseSurvival = logitReleaseSurvival.template cast<T>();   
    d.logitRecapturePhi = logitRecapturePhi.template cast<T>();
    d.sepFalpha = sepFalpha.template cast<T>();
    d.sepFlogitRho = sepFlogitRho.template cast<T>();
    d.sepFlogSd = sepFlogSd.template cast<T>();
    d.logFScaleMSY = T(logFScaleMSY);
    d.implicitFunctionDelta = T(implicitFunctionDelta);
    d.logScaleFmsy = T(logScaleFmsy);
    d.logScaleFmax = T(logScaleFmax);
    d.logScaleF01 = T(logScaleF01);
    d.logScaleFcrash = T(logScaleFcrash);
    d.logScaleFext = T(logScaleFext);
    d.logScaleFxPercent = logScaleFxPercent.template cast<T>();
    d.logScaleFlim = T(logScaleFlim);

    return d;    
  }

};


template<class Type>
void forecastSet<Type>::calculateForecast(array<Type>& logF, array<Type>& logN, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par){
  if(nYears == 0){ // no forecast
    return;
  }
  
  //array<Type>& natMor, array<Type>& catchMeanWeight, ){
  matrix<Type> natMorT = dat.natMor.matrix().transpose();
  matrix<Type> catchMeanWeightT = dat.catchMeanWeight.matrix().transpose();
  matrix<Type> landMeanWeightT = dat.landMeanWeight.matrix().transpose();
  matrix<Type> landFracT = dat.landFrac.matrix().transpose();
  forecastCalculatedMedian = matrix<Type>(logF.rows(), nYears);
  forecastCalculatedMedian.setZero();
  forecastCalculatedLogSdCorrection = vector<Type>(nYears);
  forecastCalculatedLogSdCorrection.setZero();
  int fbarFirst = conf.fbarRange(0) - conf.minAge;
  int fbarLast = conf.fbarRange(1) - conf.minAge;
  Type initialFbar = 0.0;
  for(int a = fbarFirst; a <= fbarLast; ++a){  
    initialFbar += exp(logF(conf.keyLogFsta(0,a),forecastYear.size() - nYears - 1));
  }
  initialFbar /= Type(fbarLast - fbarFirst + 1);

  vector<Type> sel(logF.rows());
  sel.setZero();
  vector<Type> selFull(logN.rows());
  selFull.setZero();

  // Correct input selectivity to have Fbar == 1
  Type inputFbar = 0.0;
  if(selectivity.size() > 0){
    Rcout << "Using custom selectivity!\n";
    for(int a = fbarFirst; a <= fbarLast; ++a){  
      if(selectivity.size() == logF.rows()){
	inputFbar += selectivity(conf.keyLogFsta(0,a));
      }else if(selectivity.size() == logN.rows()){
	inputFbar += selectivity(a);
      }else{
	Rf_error("Wrong size of selectivity. Must match logF or logN array.");
      }
    }
    inputFbar /= Type(fbarLast - fbarFirst + 1);
    if(inputFbar != 1.0){
      Rf_warning("The input selectivity was re-scaled to have Fbar equal to one.");
      selectivity /= inputFbar;
    }
  }
    

    
  for(int j = 0; j < logN.rows(); ++j){
    if(conf.keyLogFsta(0,j)>(-1)){
      if(selectivity.size() == 0){
	selFull(j) = exp(logF(conf.keyLogFsta(0,j),forecastYear.size() - nYears - 1)) / initialFbar;
	sel(conf.keyLogFsta(0,j)) = selFull(j);
      }else if(selectivity.size() == logF.rows()){
	selFull(j) = selectivity(conf.keyLogFsta(0,j));
      }else if(selectivity.size() == logN.rows()){
	selFull(j) = selectivity(j);
      }else{
	Rf_error("Wrong size of selectivity. Must match logF or logN array.");
      }
      sel(conf.keyLogFsta(0,j)) = selFull(j);
    }
  }
        
  for(int i = 0; i < nYears; ++i){
    int indx = forecastYear.size() - nYears + i;
    Type y = forecastYear(indx);
    Type calcF = 0.0;
 
    vector<Type> lastFullLogF(logN.rows());
    for(int j = 0; j < logN.rows(); ++j){
      if(conf.keyLogFsta(0,j)>(-1)){
	if(i == 0){
	  lastFullLogF(j) = log(initialFbar) + log(selFull(j));
	}else{
	  lastFullLogF(j) = forecastCalculatedMedian(conf.keyLogFsta(0,j),i-1);
	}
      }else{
	lastFullLogF(j)= 0.0; 
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
    case useFscale: // target is an F scale of previous F
      if(i == 0){
	forecastCalculatedMedian.col(i) = log(target(i)) + log(initialFbar) + log(sel);
      }else{
	forecastCalculatedMedian.col(i) = log(target(i)) + (vector<Type>)forecastCalculatedMedian.col(i-1);
      }
      break;
    case useFval: // target is F value	
      forecastCalculatedMedian.col(i) = log(target(i)) + log(sel);
      break;
    case useCatchval: // target is a catch value in weight
      if(uniroot){
	calcF = catch2F((Type)target(i), exp(lastFullLogF), (vector<Type>)natMorT.col(indx), exp((vector<Type>)logN.col(indx)), (vector<Type>)catchMeanWeightT.col(indx));
      }else{
	calcF = catch2F_quick((Type)target(i), exp(lastFullLogF), (vector<Type>)natMorT.col(indx), exp((vector<Type>)logN.col(indx)), (vector<Type>)catchMeanWeightT.col(indx));
      }
      if(i == 0){
	forecastCalculatedMedian.col(i) = log(calcF) + log(initialFbar) + log(sel);
      }else{
	forecastCalculatedMedian.col(i) = log(calcF) + (vector<Type>)forecastCalculatedMedian.col(i-1);
      }
      break;
    case useNextssb:	
      Rf_error("Forecast type not implemented");
    case useLandval:
      if(uniroot){
	calcF = landing2F((Type)target(i), exp(lastFullLogF), (vector<Type>)natMorT.col(indx), exp((vector<Type>)logN.col(indx)), (vector<Type>)landMeanWeightT.col(indx), (vector<Type>)landFracT.col(indx));
      }else{
	calcF = landing2F_quick((Type)target(i), exp(lastFullLogF), (vector<Type>)natMorT.col(indx), exp((vector<Type>)logN.col(indx)), (vector<Type>)landMeanWeightT.col(indx), (vector<Type>)landFracT.col(indx));
      }
      if(i == 0){
	forecastCalculatedMedian.col(i) = log(calcF) + log(initialFbar) + log(sel);
      }else{
	forecastCalculatedMedian.col(i) = log(calcF) + (vector<Type>)forecastCalculatedMedian.col(i-1);
      }
      break;
    case findMSY:
      forecastCalculatedMedian.col(i) = par.logFScaleMSY + log(initialFbar) + log(sel); // 
      break;
    default:
      Rf_error("Forecast type not implemented");
    }
  }
    
}









template<class Type>
Type logspace_add_p (Type logx, Type logy, Type p) {
  return log((Type(1)-p)*exp(logy-logx)+p)+logx; // the order of x and y is taylored for this application 
}

template<class Type>
Type logdrobust(Type x, Type p){
  Type ld1=dnorm(x,Type(0.0),Type(1.0),true);
  if(p<Type(1.0e-16)){
    return ld1;
  }else{
    Type ld2=dt(x,Type(3),true);
    Type logres=logspace_add_p(ld2,ld1,p);
    return logres;
  }
}
VECTORIZE2_tt(logdrobust)

template <class Type>
class MVMIX_t{
  Type halfLogDetS;         
  Type p1;                  /*fraction t3*/
  matrix<Type> Sigma;       
  vector<Type> sd;
  matrix<Type> L_Sigma;
  matrix<Type> inv_L_Sigma;
public:
  MVMIX_t(){}
  MVMIX_t(matrix<Type> Sigma_, Type p1_){
    setSigma(Sigma_);
    p1=p1_;
  }
  matrix<Type> cov(){return Sigma;}
  void setSigma(matrix<Type> Sigma_){
    Sigma = Sigma_;
    sd = sqrt(vector<Type>(Sigma.diagonal()));
    Eigen::LLT<Eigen::Matrix<Type,Eigen::Dynamic,Eigen::Dynamic> > llt(Sigma);
    L_Sigma = llt.matrixL();
    vector<Type> D=L_Sigma.diagonal();
    halfLogDetS = sum(log(D));
    inv_L_Sigma = atomic::matinv(L_Sigma);
  }
  void setSigma(matrix<Type> Sigma_, Type p1_){
    setSigma(Sigma_);
    p1=p1_;
  }
  /** \brief Evaluate the negative log density */
  Type operator()(vector<Type> x){
    vector<Type> z = inv_L_Sigma*x;
    return -sum(logdrobust(z,p1))+halfLogDetS;
  }
  Type operator()(vector<Type> x, vector<Type> keep){
    matrix<Type> S = Sigma;
    vector<Type> not_keep = Type(1.0) - keep;
    for(int i = 0; i < S.rows(); i++){
      for(int j = 0; j < S.cols(); j++){
	S(i,j) = S(i,j) * keep(i) * keep(j);
      }
      //S(i,i) += not_keep(i) * pow((Type(1)-p1)*sqrt(Type(0.5)/M_PI)+p1*(Type(1)/M_PI),2); //(t(1))
      S(i,i) += not_keep(i) * pow((Type(1)-p1)*sqrt(Type(0.5)/M_PI)+p1*(Type(2)/(M_PI*sqrt(Type(3)))),2);
    }
    return MVMIX_t<Type>(S,p1)(x * keep);
  }

  vector<Type> simulate() {
    int siz = Sigma.rows();
    vector<Type> x(siz);
    for(int i=0; i<siz; ++i){
      Type u = runif(0.0,1.0);
      if(u<p1){
        x(i) = rt(3.0);
      }else{
        x(i) = rnorm(0.0,1.0);
      }
    }
    x = L_Sigma*x;
    return x;
  }
};

template <class Type>
MVMIX_t<Type> MVMIX(matrix<Type> Sigma, Type p1){
  return MVMIX_t<Type>(Sigma,p1);
}


