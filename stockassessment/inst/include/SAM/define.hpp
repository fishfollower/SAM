

//This function returns a vector with matrices based on a list of matrices from R
HEADER(
template<class Type>
struct listMatrixFromR : vector<matrix<Type> > {

  listMatrixFromR();
  listMatrixFromR(int n);
  listMatrixFromR(SEXP x);

  template<class T>
  inline listMatrixFromR(const listMatrixFromR<T>& other) : vector<matrix<Type> >(other.size()) {
    for(int i = 0; i < other.size(); ++i)
      (*this)(i) = matrix<Type>(other(i));
  }
};
       )

SOURCE(
template<class Type>
listMatrixFromR<Type>::listMatrixFromR() : vector<matrix<Type> >() {};
	 )

SOURCE(
template<class Type>
listMatrixFromR<Type>::listMatrixFromR(int n) : vector<matrix<Type> >(n) {};
	 )

SOURCE(
	 template<class Type>
	 listMatrixFromR<Type>::listMatrixFromR(SEXP x){ 
	   (*this).resize(LENGTH(x));
	   for(int i=0; i<LENGTH(x); i++){
	     SEXP sm = VECTOR_ELT(x, i);
	     (*this)(i) = asMatrix<Type>(sm);
	   }
	 }
	 )

SAM_SPECIALIZATION(struct listMatrixFromR<double>);
SAM_SPECIALIZATION(struct listMatrixFromR<TMBad::ad_aug>);



HEADER(
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
  vector<int> minWeek;
  vector<int> maxWeek;
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

  inline dataSet() = default;

  dataSet(SEXP x);

  template<class T>
  inline dataSet(const dataSet<T> &x) :
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
    minWeek(x.minWeek),
    maxWeek(x.maxWeek),
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
});

SOURCE(
    template<class Type>
      dataSet<Type>::dataSet(SEXP x){
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
      minWeek = asVector<int>(getListElement(x,"minWeek"));
      maxWeek = asVector<int>(getListElement(x,"maxWeek"));
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
    )

SAM_SPECIALIZATION(struct dataSet<double>);
SAM_SPECIALIZATION(struct dataSet<TMBad::ad_aug>);

HEADER(
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
  vector<int> keyVarLogP;
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
  vector<int> logNMeanAssumption;
  int initState;

  inline confSet() = default;

  confSet(SEXP x);

  confSet(const confSet &other);
});

SOURCE(
	 confSet::confSet(SEXP x){
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
	   keyVarLogP = asVector<int>(getListElement(x,"keyVarLogP"));
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
	   logNMeanAssumption = asVector<int>(getListElement(x,"logNMeanAssumption"));
	   initState = (int)*REAL(getListElement(x,"initState"));
	 }
	 )

SOURCE(
	 confSet::confSet(const confSet &other) :
	 minAge(other.minAge),
	 maxAge(other.maxAge),
	 maxAgePlusGroup(other.maxAgePlusGroup),
	 keyLogFsta(other.keyLogFsta),
	 corFlag(other.corFlag),
	 keyLogFpar(other.keyLogFpar),
	 keyQpow(other.keyQpow),
	 keyVarF(other.keyVarF),
	 keyVarLogN(other.keyVarLogN),
	 keyVarLogP(other.keyVarLogP),
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
	 logNMeanAssumption(other.logNMeanAssumption),
	 initState(other.initState)
	 {}
	 );

HEADER(
template <class Type>
struct paraSet{
  vector<Type> logFpar; 
  vector<Type> logQpow; 
  vector<Type> logSdLogFsta; 
  vector<Type> logSdLogN; 
  vector<Type> logSdLogP;
  vector<Type> logSdLogObs;
  vector<Type> logSdLogTotalObs;
  vector<Type> transfIRARdist;
  vector<Type> sigmaObsParUS;
  vector<Type> rec_pars; 
  vector<Type> itrans_rho; 
  vector<Type> rhop;
  vector<Type> logScale;
  vector<Type> logitReleaseSurvival;   
  vector<Type> logitRecapturePhi; 
  vector<Type> logAlphaSCB;  
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

  vector<Type> initF;
  vector<Type> initN;

  Type splinePenalty;

  inline paraSet() = default;
  
  paraSet(SEXP x);

  template<class T>
  inline paraSet(const paraSet<T> &other) :
     logFpar(other.logFpar), 
    logQpow(other.logQpow), 
    logSdLogFsta(other.logSdLogFsta), 
    logSdLogN(other.logSdLogN), 
    logSdLogP(other.logSdLogP), 
    logSdLogObs(other.logSdLogObs),
    logSdLogTotalObs(other.logSdLogTotalObs),
    transfIRARdist(other.transfIRARdist),
    sigmaObsParUS(other.sigmaObsParUS),
    rec_pars(other.rec_pars), 
    itrans_rho(other.itrans_rho), 
    rhop(other.rhop),
    logScale(other.logScale),
    logitReleaseSurvival(other.logitReleaseSurvival),   
    logitRecapturePhi(other.logitRecapturePhi),
    logAlphaSCB(other.logAlphaSCB),
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
     initF(other.initF),
     initN(other.initN),
    splinePenalty(other.splinePenalty)  {}

});


SOURCE(
	 template<class Type>
	 paraSet<Type>::paraSet(SEXP x){
	   logFpar = asVector<Type>(getListElement(x,"logFpar"));
	   logQpow = asVector<Type>(getListElement(x,"logQpow"));
	   logSdLogFsta = asVector<Type>(getListElement(x,"logSdLogFsta"));
	   logSdLogN = asVector<Type>(getListElement(x,"logSdLogN"));
	   logSdLogP = asVector<Type>(getListElement(x,"logSdLogP"));
	   logSdLogObs = asVector<Type>(getListElement(x,"logSdLogObs"));
	   logSdLogTotalObs = asVector<Type>(getListElement(x,"logSdLogTotalObs"));
	   transfIRARdist = asVector<Type>(getListElement(x,"transfIRARdist"));
	   sigmaObsParUS = asVector<Type>(getListElement(x,"sigmaObsParUS"));
	   rec_pars = asVector<Type>(getListElement(x,"rec_pars"));
	   itrans_rho = asVector<Type>(getListElement(x,"itrans_rho"));
	   rhop = asVector<Type>(getListElement(x,"rhop"));
	   logScale = asVector<Type>(getListElement(x,"logScale"));
	   logitReleaseSurvival = asVector<Type>(getListElement(x,"logitReleaseSurvival"));
	   logitRecapturePhi = asVector<Type>(getListElement(x,"logitRecapturePhi"));
	   logAlphaSCB = asVector<Type>(getListElement(x,"logAlphaSCB"));
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
	   initF = asVector<Type>(getListElement(x,"initF"));
	   initN = asVector<Type>(getListElement(x,"initN"));

	   splinePenalty = (Type)Rf_asReal(getListElement(x,"splinePenalty"));
	 }
	 );

SAM_SPECIALIZATION(struct paraSet<double>);
SAM_SPECIALIZATION(struct paraSet<TMBad::ad_aug>);

