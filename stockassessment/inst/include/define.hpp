#pragma once
#ifndef SAM_DEFINE_HPP
#define SAM_DEFINE_HPP

#define SAM_NegInf -20.0
#define SAM_NIZero -10.0
#define SAM_Zero exp(SAM_NIZero)


// #ifndef TMBAD_FRAMEWORK
// #ifndef CPPAD_FRAMEWORK
// #define CPPAD_FRAMEWORK
// #endif
// #endif


// int SEXP2intSAM(SEXP x, int x0 = 0){
//   if(!Rf_isNull(x)){
//     return Rf_asInteger(x);
//   }
//   return x0; 
// }



template <class Type>
vector<Type> SEXP2vecSAM(SEXP x){
  if(!Rf_isNull(x)){
    return asVector<Type>(x);
  }
  return vector<Type>(0); 
}

template <class Type>
Type SEXP2scalarSAM(SEXP x, Type x0){
  if(!Rf_isNull(x)){
    return (Type)Rf_asReal(x);
  }
  return x0; 
}

// template<class Type>
// Type logspace_sum(vector<Type> logx){
//   Type r = R_NegInf;
//   for(int i = 0; i < logx.size(); ++i)
//     r = logspace_add(r, logx(i));
//   return r;
// }


template <class Type>
struct forecastSet;

template <class Type>
struct dataSet;

struct confSet;

template <class Type>
struct paraSet;

template <class Type>
struct Recruitment;

template <class Type>
struct MortalitySet;

template <class Type>
Type hcr(Type ssb, vector<Type> hcrConf);

template <class Type>
array<Type> totFFun(confSet &conf, array<Type> &logF);

#define REPORT_F(name,F)					\
  if(isDouble<Type>::value && F->current_parallel_region<0) {	\
    Rf_defineVar(Rf_install(#name),					\
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
  listMatrixFromR(const listMatrixFromR<T>& other) : vector<matrix<Type> >(other.size()) {
    for(int i = 0; i < other.size(); ++i)
      (*this)(i) = matrix<T>(other(i));
  }
  
  // template<class T>
  // listMatrixFromR<T> cast() const {
  //   int n = (*this).size();
  //   listMatrixFromR<T> d(n);
  //   for(int i = 0; i < n; ++i){
  //     matrix<T> tmp = (*this)(i).template cast<T>();
  //     d(i) = tmp;
  //   }
  //   return d;
  // }

  // listMatrixFromR<Type>& operator=(const listMatrixFromR<Type>& rhs) {
  //   (*this).resize(rhs.size());
  //   for(int i = 0; i < rhs.size(); ++i){
  //     matrix<Type> tmp = rhs(i);
  //     (*this)(i) = tmp;
  //   }
  //   return *this;
  // }

  
};


#define ADREPORT_F(name,F) F->reportvector.push(name,#name);

#define SIMULATE_F(F)				\
  if(isDouble<Type>::value && F->do_simulate)

#define NOT_SIMULATE_F(F)				\
  if(!(isDouble<Type>::value && F->do_simulate))


template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

#ifndef WITH_LIBTMB
bool isNAINT(int x){
  return R_NaInt == x; //NA_INTEGER==x;
}
#endif

template <class Type>
struct forecastSet {

  enum FModelType {
		   asFModel,
		   useFscale,
		   useFval,
		   useCatchval,
		   useNextssb,
		   useLandval,
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
  vector<Type> target;
  vector<Type> selectivity;
  vector<recModelType> recModel;
  Type logRecruitmentMedian;
  Type logRecruitmentVar;
  vector<FSdTimeScaleModel> fsdTimeScaleModel;
  vector<int> simFlag;
  int uniroot;
  vector<Type> hcrConf;
  int hcrCurrentSSB;

  matrix<Type> forecastCalculatedMedian;
  vector<Type> forecastCalculatedLogSdCorrection;
  vector<Type> sel;
  vector<Type> selFull;
  Type initialFbar;
  
  void calculateForecast(array<Type>& logF, array<Type>& logN, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, Recruitment<Type>& recruit, MortalitySet<Type>& mort); // Defined after dataSet and confSet
  void updateForecast(int i, array<Type>& logF, array<Type>& logN, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, Recruitment<Type>& recruit, MortalitySet<Type>& mort); // Defined after dataSet and confSet
  
  forecastSet() : nYears(0),
			nCatchAverageYears(0),
			aveYears(0),
			forecastYear(0),
			FModel(0),
			target(0),
			selectivity(0),
			recModel(0),
			logRecruitmentMedian(0),
			logRecruitmentVar(0),
			fsdTimeScaleModel(0),
			simFlag(0),
			uniroot(0),
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
		      // 	uniroot(0),
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
      target = vector<Type>(0);
      selectivity = vector<Type>(0);
      recModel = vector<recModelType>(0);
      logRecruitmentMedian = 0;
      logRecruitmentVar = 0;
      fsdTimeScaleModel = vector<FSdTimeScaleModel>(0);
      simFlag = vector<int>(0);
      uniroot = 0;
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
					 target(x.target),
					 selectivity(x.selectivity),
					 recModel(x.recModel),
					 logRecruitmentMedian(x.logRecruitmentMedian),
					 logRecruitmentVar(x.logRecruitmentVar),
					 fsdTimeScaleModel(x.fsdTimeScaleModel),
					 simFlag(x.simFlag),
					 uniroot(x.uniroot),
					 hcrConf(x.hcrConf),
					 hcrCurrentSSB(x.hcrCurrentSSB),
					 forecastCalculatedMedian(x-forecastCalculatedMedian),
					 forecastCalculatedLogSdCorrection(x.forecastCalculatedLogSdCorrection),
					 sel(x.sel),
					 selFull(x.selFull),
					 initialFbar(x.initialFbar) {}

  // forecastSet<Type>& operator=(const forecastSet<Type>& rhs) {
  //   nYears = rhs.nYears;
  //   if(nYears == 0)
  //     return *this;
  //   nCatchAverageYears = rhs.nCatchAverageYears;
  //   aveYears = rhs.aveYears;
  //   forecastYear = rhs.forecastYear;
  //   FModel = rhs.FModel;
  //   target = rhs.target;
  //   selectivity = rhs.selectivity;
  //   recModel = rhs.recModel;
  //   logRecruitmentMedian = rhs.logRecruitmentMedian;
  //   logRecruitmentVar = rhs.logRecruitmentVar;
  //   fsdTimeScaleModel = rhs.fsdTimeScaleModel;
  //   simFlag = rhs.simFlag;
  //   uniroot = rhs.uniroot;
  //   hcrConf = rhs.hcrConf;
  //   hcrCurrentSSB = rhs.hcrCurrentSSB;
  //   forecastCalculatedMedian = rhs.forecastCalculatedMedian;
  //   forecastCalculatedLogSdCorrection = rhs.forecastCalculatedLogSdCorrection;
  //   sel = rhs.sel;
  //   selFull = rhs.selFull;
  //   initialFbar = rhs.initialFbar;

  //   return *this;
  // }

  // template<class T>
  // forecastSet<T> cast() const {
  //   forecastSet<T> d;
  //   d.nYears = nYears;	// int
  //   if(nYears == 0)
  //     return d;
  //   d.nCatchAverageYears = nCatchAverageYears; // int
  //   d.aveYears = aveYears;	// <int>
  //   d.forecastYear = forecastYear.template cast<T>();
  //   // d.FModel = FModel; // <int>
  //   d.FModel = vector<typename forecastSet<T>::FModelType>(FModel.size());
  //   for(int i = 0; i < FModel.size(); ++i)
  //     d.FModel(i) = static_cast<typename forecastSet<T>::FModelType>((int)FModel(i));
  //   d.target = target.template cast<T>();
  //   d.selectivity = selectivity.template cast<T>();
  //   // d.recModel = recModel; // <int>
  //   d.recModel = vector<typename forecastSet<T>::recModelType>(recModel.size());
  //   for(int i = 0; i < recModel.size(); ++i)
  //     d.recModel(i) = static_cast<typename forecastSet<T>::recModelType>((int)recModel(i));
  //   d.logRecruitmentMedian = T(logRecruitmentMedian);
  //   d.logRecruitmentVar = T(logRecruitmentVar);
  //   // d.fsdTimeScaleModel = fsdTimeScaleModel; // <int>
  //   d.fsdTimeScaleModel = vector<typename forecastSet<T>::FSdTimeScaleModel>(fsdTimeScaleModel.size());
  //   for(int i = 0; i < fsdTimeScaleModel.size(); ++i)
  //     d.fsdTimeScaleModel(i) = static_cast<typename forecastSet<T>::FSdTimeScaleModel>((int)fsdTimeScaleModel(i));

  //   d.simFlag = simFlag; // <int>
  //   d.uniroot = uniroot;
  //   d.hcrConf = hcrConf.template cast<T>();
  //   d.hcrCurrentSSB = hcrCurrentSSB;
  //   d.forecastCalculatedMedian = forecastCalculatedMedian.template cast<T>();
  //   d.forecastCalculatedLogSdCorrection = forecastCalculatedLogSdCorrection.template cast<T>();
  //   d.sel = sel.template cast<T>();
  //   d.selFull = selFull.template cast<T>();
  //   d.initialFbar = T(initialFbar);
  //   return d;    
  // }

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
  array<int> sumKey;

  dataSet() = default;

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
    sumKey = asArray<int>(getListElement(x,"sumKey"));
  };

  template<class T>
  dataSet(const dataSet<T> &x) :
    noFleets(x.noFleets),
    fleetTypes(x.fleetTypes),
    sampleTimes(x.sampleTimes),
    noYears(x.noYears), 	
    years(x.years),
    minAgePerFleet(x.minAgePerFleet),
    maxAgePerFleet(x.maxAgePerFleet),
    nobs(x.nobs),
    idx1(x.idx1, x.idx1.dim),
    idx2(x.idx2, x.idx2.dim),
    idxCor(x.idxCor, x.idxCor.dim),
    aux(x.aux, x.aux.dim),
    logobs(x.logobs),		
    weight(x.weight),  // Good
    propMat(x.propMat,x.propMat.dim), //(x.propMat),
    stockMeanWeight(x.stockMeanWeight, x.stockMeanWeight.dim),
    catchMeanWeight(x.catchMeanWeight,x.catchMeanWeight.dim),
    natMor(x.natMor, x.natMor.dim),
    landFrac(x.landFrac, x.landFrac.dim),
    disMeanWeight(x.disMeanWeight,x.disMeanWeight.dim), //x.disMeanWeight),
    landMeanWeight(x.landMeanWeight, x.landMeanWeight.dim), //x.landMeanWeight),
    propF(x.propF, x.propF.dim), //x.propF),
    propM(x.propM, x.propM.dim),
    corList(x.corList),
    sumKey(x.sumKey, x.sumKey.dim) {}
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

struct confSet{
  int minAge;
  int maxAge;
  vector<int> maxAgePlusGroup;
  array<int> keyLogFsta;
  vector<int> corFlag;
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
  vector<double> fracMixN;
  vector<double> fracMixObs;
  array<int> predVarObsLink;
  int stockWeightModel;
  vector<int> keyStockWeightMean;
  vector<int> keyStockWeightObsVar;
  int catchWeightModel;
  matrix<int> keyCatchWeightMean;
  matrix<int> keyCatchWeightObsVar;
  int matureModel;
  vector<int> keyMatureMean;
  int mortalityModel;
  vector<int> keyMortalityMean;
  vector<int> keyMortalityObsVar;
  matrix<int> keyXtraSd;
  vector<int> logNMeanCorrection;

  confSet() = default;

  confSet(SEXP x){
    using tmbutils::asArray;
    minAge = (int)*REAL(getListElement(x,"minAge"));
    maxAge = (int)*REAL(getListElement(x,"maxAge"));
    maxAgePlusGroup = asVector<int>(getListElement(x,"maxAgePlusGroup"));
    keyLogFsta = asArray<int>(getListElement(x,"keyLogFsta"));
    corFlag = asVector<int>(getListElement(x,"corFlag"));
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
    fracMixN = asVector<double>(getListElement(x,"fracMixN"));
    fracMixObs = asVector<double>(getListElement(x,"fracMixObs"));
    predVarObsLink = asArray<int>(getListElement(x,"predVarObsLink"));
    stockWeightModel = (int)*REAL(getListElement(x,"stockWeightModel"));
    keyStockWeightMean = asVector<int>(getListElement(x,"keyStockWeightMean"));
    keyStockWeightObsVar = asVector<int>(getListElement(x,"keyStockWeightObsVar"));
    catchWeightModel = (int)*REAL(getListElement(x,"catchWeightModel"));
    keyCatchWeightMean = asVector<int>(getListElement(x,"keyCatchWeightMean"));
    keyCatchWeightObsVar = asVector<int>(getListElement(x,"keyCatchWeightObsVar"));
    matureModel = (int)*REAL(getListElement(x,"matureModel"));
    keyMatureMean = asVector<int>(getListElement(x,"keyMatureMean"));
    mortalityModel = (int)*REAL(getListElement(x,"mortalityModel"));
    keyMortalityMean = asVector<int>(getListElement(x,"keyMortalityMean"));
    keyMortalityObsVar = asVector<int>(getListElement(x,"keyMortalityObsVar"));
    keyXtraSd = asMatrix<int>(getListElement(x,"keyXtraSd"));
    logNMeanCorrection = asVector<int>(getListElement(x,"logNMeanCorrection"));
  }

  confSet(const confSet &other) :
    minAge(other.minAge),
    maxAge(other.maxAge),
    maxAgePlusGroup(other.maxAgePlusGroup),
    keyLogFsta(other.keyLogFsta),
    corFlag(other.corFlag),
    keyLogFpar(other.keyLogFpar),
    keyQpow(other.keyQpow),
    keyVarF(other.keyVarF),
    keyVarLogN(other.keyVarLogN),
    keyVarObs(other.keyVarObs),
    obsCorStruct(other.obsCorStruct),
    keyCorObs(other.keyCorObs),
    stockRecruitmentModelCode(other.stockRecruitmentModelCode),
    constRecBreaks(other.constRecBreaks),
    noScaledYears(other.noScaledYears),
    keyScaledYears(other.keyScaledYears),
    keyParScaledYA(other.keyParScaledYA),
    fbarRange(other.fbarRange),
    keyBiomassTreat(other.keyBiomassTreat),
    simFlag(other.simFlag),
    resFlag(other.resFlag),
    obsLikelihoodFlag(other.obsLikelihoodFlag),
    fixVarToWeight(other.fixVarToWeight),
    fracMixF(other.fracMixF),
    fracMixN(other.fracMixN),
    fracMixObs(other.fracMixObs),
    predVarObsLink(other.predVarObsLink),
    stockWeightModel(other.stockWeightModel),
    keyStockWeightMean(other.keyStockWeightMean),
    keyStockWeightObsVar(other.keyStockWeightObsVar),
    catchWeightModel(other.catchWeightModel),
    keyCatchWeightMean(other.keyCatchWeightMean),
    keyCatchWeightObsVar(other.keyCatchWeightObsVar),
    matureModel(other.matureModel),
    keyMatureMean(other.keyMatureMean),
    mortalityModel(other.mortalityModel),
    keyMortalityMean(other.keyMortalityMean),
    keyMortalityObsVar(other.keyMortalityObsVar),
    keyXtraSd(other.keyXtraSd),
    logNMeanCorrection(other.logNMeanCorrection)
  {}
  
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
  vector<Type> predVarObs;
  Type logFScaleMSY;
  Type implicitFunctionDelta;

  vector<Type> logPhiSW; 
  vector<Type> logSdProcLogSW;
  vector<Type> meanLogSW; 
  vector<Type> logSdLogSW; 
  matrix<Type> logPhiCW; 
  vector<Type> logSdProcLogCW;
  vector<Type> meanLogCW; 
  vector<Type> logSdLogCW; 
  vector<Type> logPhiMO; 
  vector<Type> logSdProcLogitMO;
  vector<Type> meanLogitMO; 
  vector<Type> logSdMO;
  vector<Type> logPhiNM; 
  vector<Type> logSdProcLogNM;
  vector<Type> meanLogNM; 
  vector<Type> logSdLogNM;
  vector<Type> logXtraSd;

  Type splinePenalty;

  paraSet() = default;
  
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

    logPhiSW = asVector<Type>(getListElement(x,"logPhiSW")); 
    logSdProcLogSW = asVector<Type>(getListElement(x,"logSdProcLogSW"));
    meanLogSW  = asVector<Type>(getListElement(x,"meanLogSW")); 
    logSdLogSW = asVector<Type>(getListElement(x,"logSdLogSW"));
    logPhiCW = asMatrix<Type>(getListElement(x,"logPhiCW")); 
    logSdProcLogCW = asVector<Type>(getListElement(x,"logSdProcLogCW"));
    meanLogCW  = asVector<Type>(getListElement(x,"meanLogCW")); 
    logSdLogCW = asVector<Type>(getListElement(x,"logSdLogCW"));
    logPhiMO = asVector<Type>(getListElement(x,"logPhiMO")); 
    logSdProcLogitMO = asVector<Type>(getListElement(x,"logSdProcLogitMO"));
    meanLogitMO  = asVector<Type>(getListElement(x,"meanLogitMO")); 
    logSdMO = asVector<Type>(getListElement(x,"logSdMO"));
    logPhiNM = asVector<Type>(getListElement(x,"logPhiNM")); 
    logSdProcLogNM = asVector<Type>(getListElement(x,"logSdProcLogNM"));
    meanLogNM  = asVector<Type>(getListElement(x,"meanLogNM")); 
    logSdLogNM = asVector<Type>(getListElement(x,"logSdLogNM"));
    logXtraSd = asVector<Type>(getListElement(x,"logXtraSd"));

    splinePenalty = (Type)Rf_asReal(getListElement(x,"splinePenalty"));
  }

  template<class T>
  paraSet(const paraSet<T> &other) :
     logFpar(other.logFpar), 
    logQpow(other.logQpow), 
    logSdLogFsta(other.logSdLogFsta), 
    logSdLogN(other.logSdLogN), 
    logSdLogObs(other.logSdLogObs),
    logSdLogTotalObs(other.logSdLogTotalObs),
    transfIRARdist(other.transfIRARdist),
    sigmaObsParUS(other.sigmaObsParUS),
    rec_pars(other.rec_pars), 
    itrans_rho(other.itrans_rho), 
    logScale(other.logScale),
    logitReleaseSurvival(other.logitReleaseSurvival),   
    logitRecapturePhi(other.logitRecapturePhi),
    sepFalpha(other.sepFalpha),
    sepFlogitRho(other.sepFlogitRho),
    sepFlogSd(other.sepFlogSd),
    logFScaleMSY(other.logFScaleMSY),
    implicitFunctionDelta(other.implicitFunctionDelta),
    logPhiSW(other.logPhiSW), 
    logSdProcLogSW(other.logSdProcLogSW),
    meanLogSW(other.meanLogSW), 
    logSdLogSW(other.logSdLogSW), 
    logPhiCW(other.logPhiCW), 
    logSdProcLogCW(other.logSdProcLogCW),
    meanLogCW(other.meanLogCW), 
    logSdLogCW(other.logSdLogCW),  
    logPhiMO(other.logPhiMO), 
    logSdProcLogitMO(other.logSdProcLogitMO),
    meanLogitMO(other.meanLogitMO), 
    logSdMO(other.logSdMO), 
    logPhiNM(other.logPhiNM),
    logSdProcLogNM(other.logSdProcLogNM),
    meanLogNM(other.meanLogNM),
    logSdLogNM(other.logSdLogNM),
    logXtraSd(other.logXtraSd),
    splinePenalty(other.splinePenalty)  {}

};


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



// Forward declarations
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



template<class Type>
void forecastSet<Type>::calculateForecast(array<Type>& logF, array<Type>& logN, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, Recruitment<Type>& recruit, MortalitySet<Type>& mort){
  if(nYears == 0){ // no forecast
    return;
  }

  for(int i = 0; i < nYears; ++i){
    updateForecast(i, logF, logN, dat, conf, par, recruit, mort);
  }
  return;

}

template<class Type>
void forecastSet<Type>::updateForecast(int i, array<Type>& logF, array<Type>& logN, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, Recruitment<Type>& recruit, MortalitySet<Type>& mort){
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
	  lastFullLogF(j) = logspace_add2(lastFullLogF(j),lastShortLogF(conf.keyLogFsta(f,j)));
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
      calcF = catch2F_quick((Type)target(i),
			      exp(lastFullLogF),
			      (vector<Type>)dat.natMor.matrix().row(indx),
			      exp((vector<Type>)logN.col(indx)),
			      (vector<Type>)dat.catchMeanWeight.matrix().row(indx));
      if(i == 0){
	forecastCalculatedMedian.col(i) = log(calcF) + log(initialFbar) + log(sel);
      }else{
	forecastCalculatedMedian.col(i) = log(calcF) + (vector<Type>)forecastCalculatedMedian.col(i-1);
      }
      break;
    case useNextssb:
      calcF = ssb2F_quick((Type)target(i),
			  log(sel), //lastShortLogF,
			  dat,
			  conf,
			  par,
			  logF,
			  logN,
			  mort,
			  indx,
			  logRecruitmentMedian);
      forecastCalculatedMedian.col(i) = log(calcF) + log(sel);
      break;
    case useLandval:
  	calcF = landing2F_quick((Type)target(i),
				exp(lastFullLogF),
				(vector<Type>)dat.natMor.matrix().row(indx),
				exp((vector<Type>)logN.col(indx)),
				(vector<Type>)dat.landMeanWeight.matrix().row(indx),
				(vector<Type>)dat.landFrac.matrix().row(indx));
      if(i == 0){
	forecastCalculatedMedian.col(i) = log(calcF) + log(initialFbar) + log(sel);
      }else{
	forecastCalculatedMedian.col(i) = log(calcF) + (vector<Type>)forecastCalculatedMedian.col(i-1);
      }
      break;
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
  vector<Type> p1;                  /*fraction t3*/
  matrix<Type> Sigma;       
  vector<Type> sd;
  matrix<Type> L_Sigma;
  matrix<Type> inv_L_Sigma;
public:
  MVMIX_t(){}
  MVMIX_t(matrix<Type> Sigma_, Type p1_){
    setSigma(Sigma_);
    p1=vector<Type>(Sigma_.rows());
    p1.setConstant(p1_);
  }
  MVMIX_t(matrix<Type> Sigma_, vector<Type> p1_){
    setSigma(Sigma_);
    p1=p1_;
  }
  MVMIX_t(matrix<Type> Sigma_, Type p1_, bool useAtomic){
    setSigma(Sigma_, useAtomic);
    p1=p1_;
  }
  matrix<Type> cov(){return Sigma;}
  void setSigma(matrix<Type> Sigma_, bool useAtomic = true){
    Sigma = Sigma_;
    p1 = vector<Type>(Sigma_.rows());
    p1.setZero();
    sd = sqrt(vector<Type>(Sigma.diagonal()));
    Eigen::LLT<Eigen::Matrix<Type,Eigen::Dynamic,Eigen::Dynamic> > llt(Sigma);
    L_Sigma = llt.matrixL();
    vector<Type> D=L_Sigma.diagonal();
    halfLogDetS = sum(log(D));
    // if(useAtomic){
      //inv_L_Sigma = atomic::matinv(L_Sigma);
    // }else{
      inv_L_Sigma = L_Sigma.inverse();
    // }
  }
  void setSigma(matrix<Type> Sigma_, Type p1_){
    setSigma(Sigma_);
    p1.setConstant(p1_);
  }
  /** \brief Evaluate the negative log density */
  Type operator()(vector<Type> x){
    vector<Type> z = inv_L_Sigma*x;
    return -sum(logdrobust(z,p1))+halfLogDetS;
  }
  Type operator()(vector<Type> x, vector<Type> keep){
    matrix<Type> S = Sigma;
    vector<Type> p = p1; 
    vector<Type> not_keep = Type(1.0) - keep;
    for(int i = 0; i < S.rows(); i++){
      for(int j = 0; j < S.cols(); j++){
	S(i,j) = S(i,j) * keep(i) * keep(j);
      }
      S(i,i) += not_keep(i) * pow((Type(1)-p(i))*sqrt(Type(0.5)/M_PI)+p(i)*(Type(2)/(M_PI*sqrt(Type(3)))),2);
    }
    return MVMIX_t<Type>(S,p)(x * keep);
  }

  vector<Type> simulate() {
    int siz = Sigma.rows();
    vector<Type> x(siz);
    for(int i=0; i<siz; ++i){
      Type u = runif(0.0,1.0);
      if(u<p1(i)){
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

template <class Type>
Type findLinkV(Type k, int n=0){
  // small helper function to solve exp(v)-exp(k-v/2)-1=0 for v
  Type v = log(exp(k)+Type(1));
  for(int i=0; i<n; ++i){
    v -= (exp(v)-exp(k-0.5*v)-1.0)/(exp(v)+.5*exp(k-0.5*v));
  }
  return v;
}

template<class Type>
struct referencepointSet {

  enum CatchType {
		  totalCatch,
		  landings,
		  discard
  };

  int nYears;
  int rpType;
  vector<int> aveYears;
  vector<int> selYears;
  vector<Type> logCustomSel;
  vector<Type> xVal;
  CatchType catchType;
  vector<Type> logF0;		// Starting value for optimization
  vector<Type> logSel;

  referencepointSet() = default;

  referencepointSet(int nYears_,
		    int CT,
		    int i,
		    array<Type> logF,
		    confSet conf) :
    nYears(nYears_), rpType(-99), aveYears(1), selYears(1), logCustomSel(0), xVal(0), catchType(static_cast<typename referencepointSet<Type>::CatchType>(CT)), logF0(0), logSel(0) {
    aveYears(0) = i;
    selYears(0) = i;
    setLogSelectivity(logF,conf);
  }
  
  referencepointSet(SEXP x) {
    if(!Rf_isNull(getListElement(x,"nYears"))){
      nYears = Rf_asInteger(getListElement(x,"nYears"));
    }else{
      nYears = 0;
    }
    if(!Rf_isNull(getListElement(x,"rpType"))){
      rpType = Rf_asInteger(getListElement(x,"rpType"));
    }else{
      rpType = -99;
    }
    if(!Rf_isNull(getListElement(x,"aveYears"))){
      aveYears = asVector<int>(getListElement(x,"aveYears"));
    }else{
      aveYears = vector<int>(0);
    }
    if(!Rf_isNull(getListElement(x,"selYears"))){
      selYears = asVector<int>(getListElement(x,"selYears"));
    }else{
      selYears = vector<int>(0);
    }
    if(!Rf_isNull(getListElement(x,"logCustomSel"))){
      logCustomSel = asVector<Type>(getListElement(x,"logCustomSel"));
    }else{
      logCustomSel = vector<Type>(0);
    }
    if(!Rf_isNull(getListElement(x,"xVal"))){
      xVal = asVector<Type>(getListElement(x,"xVal"));
    }else{
      xVal = vector<Type>(0);
    }
    if(!Rf_isNull(getListElement(x,"catchType"))){
      catchType = static_cast<typename referencepointSet<Type>::CatchType>(Rf_asInteger(getListElement(x,"catchType")));
    }else{
      catchType = static_cast<typename referencepointSet<Type>::CatchType>(0);
    }
    if(!Rf_isNull(getListElement(x,"logF0"))){
      logF0 = asVector<Type>(getListElement(x,"logF0"));
    }else{
      logF0 = vector<Type>(0);
    }
    logSel = vector<Type>(0);

  }

  template<class T>
  referencepointSet(const referencepointSet<T>& other) :
    nYears(other.nYears),
    rpType(other.rpType),
    aveYears(other.aveYears),
    selYears(other.selYears),
    logCustomSel(other.logCustomSel),
    xVal(other.xVal),
    catchType(static_cast<typename referencepointSet<Type>::CatchType>((int)other.catchType)),
    logF0(other.logF0),
    logSel(other.logSel)
  {}    
 
  vector<Type> getLogSelectivity(){
    if(logSel.size() == 0)
      Rf_error("logSelectivity not set. Call with logF and conf.");
    return logSel;
  }
  
  void setLogSelectivity(array<Type>& logF, confSet& conf){
    if(logCustomSel.size() == logF.rows()){
      logSel = logCustomSel;
      return;
    }else if(logCustomSel.size() > 0){
      Rf_error("Size of logCustomSel does not match logF");
    }
    vector<Type> logfbartmp(selYears.size());
    logfbartmp.setConstant(R_NegInf);
    Type logfsum = R_NegInf;
    logSel = vector<Type>(logF.rows());
    logSel.setConstant(R_NegInf);
    for(int y = 0; y < selYears.size(); ++y){
      for(int i = 0; i < logSel.size(); ++i){
	logSel(i) = logspace_add2(logSel(i), logF(i,selYears(y)));
      }
      for(int a = conf.fbarRange(0); a <= conf.fbarRange(1); a++){
	for(int f = 0; f < conf.keyLogFsta.dim(0); ++f)
	  if(conf.keyLogFsta(f,a-conf.minAge) > (-1))
	    logfbartmp(y) = logspace_add2(logfbartmp(y), logF(conf.keyLogFsta(f,a-conf.minAge),selYears(y)));
      }
      logfbartmp(y) -= log(Type(conf.fbarRange(1)-conf.fbarRange(0)+1));
      logfsum = logspace_add2(logfsum, logfbartmp(y));
    }
    logSel -= logfsum;
    return;
  }

  Type logFbar(array<Type>& logF, confSet& conf){
    vector<Type> logfbartmp(selYears.size());
    logfbartmp.setConstant(R_NegInf);
    Type logfsum = R_NegInf;
     for(int y = 0; y < selYears.size(); ++y){
       for(int a = conf.fbarRange(0); a <= conf.fbarRange(1); a++){  
	 for(int f = 0; f < conf.keyLogFsta.dim(0); ++f)
	   if(conf.keyLogFsta(f,a-conf.minAge) > (-1))
	     logfbartmp(y) = logspace_add2(logfbartmp(y), logF(conf.keyLogFsta(f,a-conf.minAge),selYears(y)));
       }
       logfbartmp(y) -= log(Type(conf.fbarRange(1)-conf.fbarRange(0)+1));
       logfsum = logspace_add2(logfsum, logfbartmp(y));
     }
     return logfsum;
  }
};



template<class Type>
struct referencepointList : vector<referencepointSet<Type> > {

  newton::newton_config cfg;

  referencepointList() :
    vector<referencepointSet<Type> >(0), cfg(){}
  referencepointList(int i) :
    vector<referencepointSet<Type> >(i), cfg(){}
  referencepointList(SEXP x) :
    vector<referencepointSet<Type> >(Rf_length(x)), cfg(Rf_getAttrib(x,Rf_install("newton_config"))) {
    for(int i = 0; i < Rf_length(x); ++i){
      (*this)(i) = referencepointSet<Type>(VECTOR_ELT(x,i));
    }
  }

  template<class T>
  referencepointList(const listMatrixFromR<T>& other) : vector<referencepointSet<Type> >(other.size()), cfg(other.cfg) {
    for(int i = 0; i < other.size(); ++i)
      (*this)(i) = referencepointSet<T>(other(i));
  }
  
};


#endif
