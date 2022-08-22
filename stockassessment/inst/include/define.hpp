#pragma once
#ifndef SAM_DEFINE_HPP
#define SAM_DEFINE_HPP

#define SAM_NegInf -20.0
#define SAM_NIZero -20.0
#define SAM_Zero exp(SAM_NIZero)

#define SAM_ASSERT(x,msg) if(!(x)) Rf_error(msg);


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



// template <class Type>
// vector<Type> SEXP2vecSAM(SEXP x){
//   if(!Rf_isNull(x)){
//     return asVector<Type>(x);
//   }
//   return vector<Type>(0); 
// }

// template <class Type>
// Type SEXP2scalarSAM(SEXP x, Type x0){
//   if(!Rf_isNull(x)){
//     return (Type)Rf_asReal(x);
//   }
//   return x0; 
// }

// template<class Type>
// Type logspace_sum(vector<Type> logx){
//   Type r = R_NegInf;
//   for(int i = 0; i < logx.size(); ++i)
//     r = logspace_add(r, logx(i));
//   return r;
// }


// template <class Type>
// struct forecastSet;

// template <class Type>
// struct dataSet;

// struct confSet;

// template <class Type>
// struct paraSet;

// template <class Type>
// struct Recruitment;

// template <class Type>
// struct MortalitySet;



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

bool isNAINT(int x){
  return R_NaInt == x; //NA_INTEGER==x;
}

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
    keyCatchWeightMean = asMatrix<int>(getListElement(x,"keyCatchWeightMean"));
    keyCatchWeightObsVar = asMatrix<int>(getListElement(x,"keyCatchWeightObsVar"));
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
	logSel(i) = logspace_add_SAM(logSel(i), logF(i,selYears(y)));
      }
      for(int a = conf.fbarRange(0); a <= conf.fbarRange(1); a++){
	for(int f = 0; f < conf.keyLogFsta.dim(0); ++f)
	  if(conf.keyLogFsta(f,a-conf.minAge) > (-1))
	    logfbartmp(y) = logspace_add_SAM(logfbartmp(y), logF(conf.keyLogFsta(f,a-conf.minAge),selYears(y)));
      }
      logfbartmp(y) -= log(Type(conf.fbarRange(1)-conf.fbarRange(0)+1));
      logfsum = logspace_add_SAM(logfsum, logfbartmp(y));
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
	     logfbartmp(y) = logspace_add_SAM(logfbartmp(y), logF(conf.keyLogFsta(f,a-conf.minAge),selYears(y)));
       }
       logfbartmp(y) -= log(Type(conf.fbarRange(1)-conf.fbarRange(0)+1));
       logfsum = logspace_add_SAM(logfsum, logfbartmp(y));
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
