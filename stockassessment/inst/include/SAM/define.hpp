

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
	 listMatrixFromR<Type>::listMatrixFromR(SEXP x) : vector<matrix<Type> >(Rf_length(x)){ 
	   //(*this).resize(LENGTH(x));
	   for(int i=0; i<Rf_length(x); i++){
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
  vector<Type> sampleTimesStart;
  vector<Type> sampleTimesEnd;
  int noYears;
  vector<Type> years;
  vector<int> minAgePerFleet;
  vector<int> maxAgePerFleet;
  vector<int> fleetCovarianceSize;
  int nobs;
  array<int> idx1;
  array<int> idx2;
  array<int> idxCor;
  vector<int> minWeek;
  vector<int> maxWeek;
  array<int> aux;
  array<Type> auxData;
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

  inline dataSet() :
    noFleets(),
    fleetTypes(),
    sampleTimesStart(),
    sampleTimesEnd(),
    noYears(), 	
    years(),
    minAgePerFleet(),
    maxAgePerFleet(),
    fleetCovarianceSize(),
    nobs(),
    idx1(),
    idx2(),
    idxCor(),
    minWeek(),
    maxWeek(),
    aux(),
    auxData(),
    logobs(),		
    weight(),  // Good
    propMat(), //(x.propMat),
    stockMeanWeight(),
    catchMeanWeight(),
    natMor(),
    landFrac(),
    disMeanWeight(), //x.disMeanWeight),
    landMeanWeight(), //x.landMeanWeight),
    propF(), //x.propF),
    propM(),
    corList(),
    sumKey() {}

  dataSet(SEXP x);

  template<class T>
  inline dataSet(const dataSet<T> &x) :
    noFleets(x.noFleets),
    fleetTypes(x.fleetTypes),
    sampleTimesStart(x.sampleTimesStart),
    sampleTimesEnd(x.sampleTimesEnd),
    noYears(x.noYears), 	
    years(x.years),
    minAgePerFleet(x.minAgePerFleet),
    maxAgePerFleet(x.maxAgePerFleet),
    fleetCovarianceSize(x.fleetCovarianceSize),
    nobs(x.nobs),
    idx1(x.idx1, x.idx1.dim),
    idx2(x.idx2, x.idx2.dim),
    idxCor(x.idxCor, x.idxCor.dim),
    minWeek(x.minWeek),
    maxWeek(x.maxWeek),
    aux(x.aux, x.aux.dim),
    auxData(x.auxData, x.auxData.dim),
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
    sumKey(x.sumKey, x.sumKey.dim)
  {}
});

SOURCE(
       template<class Type>
       dataSet<Type>::dataSet() :
       noFleets(),
       fleetTypes(),
       sampleTimes(),
       noYears(), 	
       years(),
       minAgePerFleet(),
       maxAgePerFleet(),
       nobs(),
       idx1(),
       idx2(),
       idxCor(),
       minWeek(),
       maxWeek(),
       aux(),
       logobs(),		
       weight(),  // Good
       propMat(), //(x.propMat),
       stockMeanWeight(),
       catchMeanWeight(),
       natMor(),
       landFrac(),
       disMeanWeight(), //x.disMeanWeight),
       landMeanWeight(), //x.landMeanWeight),
       propF(), //x.propF),
       propM(),
       corList(),
       sumKey() {};
       )

SOURCE(
    template<class Type>
      dataSet<Type>::dataSet(SEXP x){
      using tmbutils::asArray;
      noFleets = Rf_asInteger(getListElement(x,"noFleets", &isNumericScalar));
      fleetTypes = asVector<int>(getListElement(x,"fleetTypes",  &Rf_isNumeric));
      sampleTimesStart = asVector<Type>(getListElement(x,"sampleTimesStart",  &Rf_isNumeric));
      sampleTimesEnd = asVector<Type>(getListElement(x,"sampleTimesEnd",  &Rf_isNumeric));
      noYears = Rf_asInteger(getListElement(x,"noYears", &isNumericScalar));
      years = asVector<Type>(getListElement(x,"years",  &Rf_isNumeric));
      minAgePerFleet = asVector<int>(getListElement(x,"minAgePerFleet",  &Rf_isNumeric));
      maxAgePerFleet = asVector<int>(getListElement(x,"maxAgePerFleet",  &Rf_isNumeric));
      fleetCovarianceSize = asVector<int>(getListElement(x,"fleetCovarianceSize",  &Rf_isNumeric));
      nobs = Rf_asInteger(getListElement(x,"nobs", &isNumericScalar));
      idx1 = asArray<int>(getListElement(x,"idx1", &Rf_isArray));
      idx2 = asArray<int>(getListElement(x,"idx2", &Rf_isArray));
      idxCor = asArray<int>(getListElement(x,"idxCor", &Rf_isArray));
      minWeek = asVector<int>(getListElement(x,"minWeek",  &Rf_isNumeric));
      maxWeek = asVector<int>(getListElement(x,"maxWeek",  &Rf_isNumeric));
      aux = asArray<int>(getListElement(x,"aux", &Rf_isArray));
      auxData = asArray<Type>(getListElement(x,"auxData", &Rf_isArray));
      logobs = asVector<Type>(getListElement(x,"logobs",  &Rf_isNumeric));
      weight = asVector<Type>(getListElement(x,"weight",  &Rf_isNumeric));
      propMat = asArray<Type>(getListElement(x,"propMat", &Rf_isArray));
      stockMeanWeight = asArray<Type>(getListElement(x,"stockMeanWeight", &Rf_isArray));
      catchMeanWeight = asArray<Type>(getListElement(x,"catchMeanWeight", &Rf_isArray));
      natMor = asArray<Type>(getListElement(x,"natMor", &Rf_isArray));
      landFrac = asArray<Type>(getListElement(x,"landFrac", &Rf_isArray));
      disMeanWeight = asArray<Type>(getListElement(x,"disMeanWeight", &Rf_isArray));
      landMeanWeight = asArray<Type>(getListElement(x,"landMeanWeight", &Rf_isArray));
      propF = asArray<Type>(getListElement(x,"propF", &Rf_isArray));
      propM = asArray<Type>(getListElement(x,"propM", &Rf_isArray));
      corList = listMatrixFromR<Type>(getListElement(x,"corList"));
      sumKey = asArray<int>(getListElement(x,"sumKey", &Rf_isArray));
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
  array<int> keyLogFmu;
  array<int> keyLogFrho;
  array<int> keyLogFbound_kappa;
  array<int> keyLogFbound_alpha;
  array<int> keyLogFbound_tau;
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
  vector<int> fixVarToWeight;
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
  array<int> keyLogFseason;
  vector<double> seasonTimes;
  vector<int> isFishingSeason;
  double seasonFirstYear;
  int seasonFixedEffect;
  
  confSet();

  confSet(SEXP x);

  confSet(const confSet &other);
});

SOURCE(
       confSet::confSet(SEXP x){
	 using tmbutils::asArray;
	 minAge = Rf_asInteger(getListElement(x,"minAge", &isNumericScalar));
	 maxAge = Rf_asInteger(getListElement(x,"maxAge", &isNumericScalar));
	 maxAgePlusGroup = asVector<int>(getListElement(x,"maxAgePlusGroup", &Rf_isNumeric));
	 keyLogFsta = asArray<int>(getListElement(x,"keyLogFsta", &Rf_isArray));
	 keyLogFmu = asArray<int>(getListElement(x,"keyLogFmu", &Rf_isArray));
	 keyLogFrho = asArray<int>(getListElement(x,"keyLogFrho", &Rf_isArray));
	 keyLogFbound_kappa = asArray<int>(getListElement(x,"keyLogFbound_kappa", &Rf_isArray));
	 keyLogFbound_alpha = asArray<int>(getListElement(x,"keyLogFbound_alpha", &Rf_isArray));
	 keyLogFbound_tau = asArray<int>(getListElement(x,"keyLogFbound_tau", &Rf_isArray));
	 corFlag = asVector<int>(getListElement(x,"corFlag", &Rf_isNumeric));
	 keyLogFpar = asArray<int>(getListElement(x,"keyLogFpar", &Rf_isArray));
	 keyQpow = asArray<int>(getListElement(x,"keyQpow", &Rf_isArray));
	 keyVarF = asArray<int>(getListElement(x,"keyVarF", &Rf_isArray));
	 keyVarLogN = asVector<int>(getListElement(x,"keyVarLogN", &Rf_isNumeric));
	 keyVarLogP = asVector<int>(getListElement(x,"keyVarLogP", &Rf_isNumeric));
	 keyVarObs = asArray<int>(getListElement(x,"keyVarObs", &Rf_isArray));
	 obsCorStruct = asVector<int>(getListElement(x,"obsCorStruct", &Rf_isNumeric));
	 keyCorObs = asArray<int>(getListElement(x,"keyCorObs", &Rf_isArray));
	 stockRecruitmentModelCode = Rf_asInteger(getListElement(x,"stockRecruitmentModelCode", &isNumericScalar));
	 constRecBreaks = asVector<double>(getListElement(x,"constRecBreaks", &Rf_isNumeric));
	 noScaledYears = Rf_asInteger(getListElement(x,"noScaledYears", &isNumericScalar));
	 keyScaledYears = asVector<int>(getListElement(x,"keyScaledYears", &Rf_isNumeric));
	 keyParScaledYA = asMatrix<int>(getListElement(x,"keyParScaledYA", &Rf_isMatrix));
	 fbarRange = asVector<int>(getListElement(x,"fbarRange", &Rf_isNumeric));
	 keyBiomassTreat = asVector<int>(getListElement(x,"keyBiomassTreat", &Rf_isNumeric));
	 simFlag = asVector<int>(getListElement(x,"simFlag", &Rf_isNumeric));
	 resFlag = Rf_asInteger(getListElement(x,"resFlag", &isNumericScalar));
	 obsLikelihoodFlag = asVector<int>(getListElement(x,"obsLikelihoodFlag", &Rf_isNumeric));
	 fixVarToWeight = asVector<int>(getListElement(x,"fixVarToWeight", &Rf_isNumeric)); 
	 fracMixF = Rf_asReal(getListElement(x,"fracMixF", &isNumericScalar));
	 fracMixN = asVector<double>(getListElement(x,"fracMixN", &Rf_isNumeric));
	 fracMixObs = asVector<double>(getListElement(x,"fracMixObs", &Rf_isNumeric));
	 predVarObsLink = asArray<int>(getListElement(x,"predVarObsLink", &Rf_isArray));
	 stockWeightModel = Rf_asInteger(getListElement(x,"stockWeightModel", &isNumericScalar));
	 keyStockWeightMean = asVector<int>(getListElement(x,"keyStockWeightMean", &Rf_isNumeric));
	 keyStockWeightObsVar = asVector<int>(getListElement(x,"keyStockWeightObsVar", &Rf_isNumeric));
	 catchWeightModel = Rf_asInteger(getListElement(x,"catchWeightModel", &isNumericScalar));
	 keyCatchWeightMean = asMatrix<int>(getListElement(x,"keyCatchWeightMean", &Rf_isMatrix));
	 keyCatchWeightObsVar = asMatrix<int>(getListElement(x,"keyCatchWeightObsVar", &Rf_isMatrix));
	 matureModel = Rf_asInteger(getListElement(x,"matureModel", &isNumericScalar));
	 keyMatureMean = asVector<int>(getListElement(x,"keyMatureMean", &Rf_isNumeric));
	 mortalityModel = Rf_asInteger(getListElement(x,"mortalityModel", &isNumericScalar));
	 keyMortalityMean = asVector<int>(getListElement(x,"keyMortalityMean", &Rf_isNumeric));
	 keyMortalityObsVar = asVector<int>(getListElement(x,"keyMortalityObsVar", &Rf_isNumeric));
	 keyXtraSd = asMatrix<int>(getListElement(x,"keyXtraSd", &Rf_isMatrix));
	 logNMeanAssumption = asVector<int>(getListElement(x,"logNMeanAssumption", &Rf_isNumeric));
	 initState = Rf_asInteger(getListElement(x,"initState", &isNumericScalar));
	 keyLogFseason = asArray<int>(getListElement(x,"keyLogFseason", &Rf_isArray));
	 seasonTimes = asVector<double>(getListElement(x,"seasonTimes", &Rf_isNumeric));
	 isFishingSeason = asVector<int>(getListElement(x,"isFishingSeason", &Rf_isNumeric));
	 seasonFirstYear = Rf_asReal(getListElement(x,"seasonFirstYear", &isNumericScalar));
	 seasonFixedEffect = Rf_asInteger(getListElement(x,"seasonFixedEffect", &isNumericScalar));
       }
       )

SOURCE(
	 confSet::confSet() :
	 minAge(),
	 maxAge(),
	 maxAgePlusGroup(),
	 keyLogFsta(),
	 keyLogFmu(),
	 keyLogFrho(),
	 keyLogFbound_kappa(),
	 keyLogFbound_alpha(),
	 keyLogFbound_tau(),
	 corFlag(),
	 keyLogFpar(),
	 keyQpow(),
	 keyVarF(),
	 keyVarLogN(),
	 keyVarLogP(),
	 keyVarObs(),
	 obsCorStruct(),
	 keyCorObs(),
	 stockRecruitmentModelCode(),
	 constRecBreaks(),
	 noScaledYears(),
	 keyScaledYears(),
	 keyParScaledYA(),
	 fbarRange(),
	 keyBiomassTreat(),
	 simFlag(),
	 resFlag(),
	 obsLikelihoodFlag(),
	 fixVarToWeight(),
	 fracMixF(),
	 fracMixN(),
	 fracMixObs(),
	 predVarObsLink(),
	 stockWeightModel(),
	 keyStockWeightMean(),
	 keyStockWeightObsVar(),
	 catchWeightModel(),
	 keyCatchWeightMean(),
	 keyCatchWeightObsVar(),
	 matureModel(),
	 keyMatureMean(),
	 mortalityModel(),
	 keyMortalityMean(),
	 keyMortalityObsVar(),
	 keyXtraSd(),
	 logNMeanAssumption(),
	 initState(),
	 keyLogFseason(),
	 seasonTimes(),
	 isFishingSeason(),
	 seasonFirstYear(),
	 seasonFixedEffect()
	 {}
	 );

SOURCE(
	 confSet::confSet(const confSet &other) :
	 minAge(other.minAge),
	 maxAge(other.maxAge),
	 maxAgePlusGroup(other.maxAgePlusGroup),
	 keyLogFsta(other.keyLogFsta),
	 keyLogFmu(other.keyLogFmu),
	 keyLogFrho(other.keyLogFrho),
	 keyLogFbound_kappa(other.keyLogFbound_kappa),
	 keyLogFbound_alpha(other.keyLogFbound_alpha),
	 keyLogFbound_tau(other.keyLogFbound_tau),
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
	 initState(other.initState),
	 keyLogFseason(other.keyLogFseason),
	 seasonTimes(other.seasonTimes),
	 isFishingSeason(other.isFishingSeason),
	 seasonFirstYear(other.seasonFirstYear),
	 seasonFixedEffect(other.seasonFixedEffect)
	 {}
	 );

SOURCE(
	 confSet::confSet() :
	 minAge(),
	 maxAge(),
	 maxAgePlusGroup(),
	 keyLogFsta(),
	 corFlag(),
	 keyLogFpar(),
	 keyQpow(),
	 keyVarF(),
	 keyVarLogN(),
	 keyVarLogP(),
	 keyVarObs(),
	 obsCorStruct(),
	 keyCorObs(),
	 stockRecruitmentModelCode(),
	 constRecBreaks(),
	 noScaledYears(),
	 keyScaledYears(),
	 keyParScaledYA(),
	 fbarRange(),
	 keyBiomassTreat(),
	 simFlag(),
	 resFlag(),
	 obsLikelihoodFlag(),
	 fixVarToWeight(),
	 fracMixF(),
	 fracMixN(),
	 fracMixObs(),
	 predVarObsLink(),
	 stockWeightModel(),
	 keyStockWeightMean(),
	 keyStockWeightObsVar(),
	 catchWeightModel(),
	 keyCatchWeightMean(),
	 keyCatchWeightObsVar(),
	 matureModel(),
	 keyMatureMean(),
	 mortalityModel(),
	 keyMortalityMean(),
	 keyMortalityObsVar(),
	 keyXtraSd(),
	 logNMeanAssumption(),
	 initState()
	 {}
	 );

HEADER(
template <class Type>
struct paraSet{
  vector<Type> logFpar; 
  vector<Type> logQpow; 
  vector<Type> logSdLogFsta; 
  vector<Type> muF; 
  vector<Type> trans_rho_F; 
  vector<Type> boundF_kappa; 
  vector<Type> boundF_alpha; 
  vector<Type> boundF_tau; 
  vector<Type> logSdLogN; 
  vector<Type> logSdLogP;
  vector<Type> logSdLogObs;
  vector<Type> logSdLogTotalObs;
  vector<Type> transfIRARdist;
  vector<Type> sigmaObsParUS;
  vector<Type> rec_pars;
  vector<Type> rec_transphi;
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

  matrix<Type> seasonMu;
  vector<Type> seasonLogitRho;
  vector<Type> seasonLogSd;

  Type splinePenalty;

  inline paraSet() :
    logFpar(), 
    logQpow(), 
    logSdLogFsta(), 
    muF(), 
    trans_rho_F(),
    boundF_kappa(),
    boundF_alpha(),
    boundF_tau(),
    logSdLogN(), 
    logSdLogP(), 
    logSdLogObs(),
    logSdLogTotalObs(),
    transfIRARdist(),
    sigmaObsParUS(),
    rec_pars(),
     rec_transphi(),
    itrans_rho(), 
    rhop(),
    logScale(),
    logitReleaseSurvival(),   
    logitRecapturePhi(),
    logAlphaSCB(),
    sepFalpha(),
    sepFlogitRho(),
    sepFlogSd(),
    logFScaleMSY(),
    implicitFunctionDelta(),
    logPhiSW(), 
    logSdProcLogSW(),
    meanLogSW(), 
    logSdLogSW(), 
    logPhiCW(), 
    logSdProcLogCW(),
    meanLogCW(), 
    logSdLogCW(),  
    logPhiMO(), 
    logSdProcLogitMO(),
    meanLogitMO(), 
    logSdMO(), 
    logPhiNM(),
    logSdProcLogNM(),
    meanLogNM(),
    logSdLogNM(),
    logXtraSd(),
     initF(),
     initN(),
    seasonMu(),
    seasonLogitRho(),
    seasonLogSd(),
    splinePenalty()  {}
  
  paraSet(SEXP x);

  template<class T>
  inline paraSet(const paraSet<T> &other) :
     logFpar(other.logFpar), 
    logQpow(other.logQpow), 
    logSdLogFsta(other.logSdLogFsta), 
    muF(other.muF), 
    trans_rho_F(other.trans_rho_F),
     boundF_kappa(other.boundF_kappa),
     boundF_alpha(other.boundF_alpha),
     boundF_tau(other.boundF_tau),
    logSdLogN(other.logSdLogN), 
    logSdLogP(other.logSdLogP), 
    logSdLogObs(other.logSdLogObs),
    logSdLogTotalObs(other.logSdLogTotalObs),
    transfIRARdist(other.transfIRARdist),
    sigmaObsParUS(other.sigmaObsParUS),
    rec_pars(other.rec_pars),
     rec_transphi(other.rec_transphi),
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
    seasonMu(other.seasonMu),
    seasonLogitRho(other.seasonLogitRho),
    seasonLogSd(other.seasonLogSd),
    splinePenalty(other.splinePenalty)  {}

});

SOURCE(
       template<class Type>
       paraSet<Type>::paraSet() :
       logFpar(), 
       logQpow(), 
       logSdLogFsta(), 
       logSdLogN(), 
       logSdLogP(), 
       logSdLogObs(),
       logSdLogTotalObs(),
       transfIRARdist(),
       sigmaObsParUS(),
       rec_pars(), 
       itrans_rho(), 
       rhop(),
       logScale(),
       logitReleaseSurvival(),   
       logitRecapturePhi(),
       logAlphaSCB(),
       sepFalpha(),
       sepFlogitRho(),
       sepFlogSd(),
       logFScaleMSY(),
       implicitFunctionDelta(),
       logPhiSW(), 
       logSdProcLogSW(),
       meanLogSW(), 
       logSdLogSW(), 
       logPhiCW(), 
       logSdProcLogCW(),
       meanLogCW(), 
       logSdLogCW(),  
       logPhiMO(), 
       logSdProcLogitMO(),
       meanLogitMO(), 
       logSdMO(), 
       logPhiNM(),
       logSdProcLogNM(),
       meanLogNM(),
       logSdLogNM(),
       logXtraSd(),
       initF(),
       initN(),
       splinePenalty()  {};
       )

SOURCE(
	 template<class Type>
	 paraSet<Type>::paraSet(SEXP x){
	   logFpar = asVector<Type>(getListElement(x,"logFpar", &Rf_isNumeric));
	   logQpow = asVector<Type>(getListElement(x,"logQpow", &Rf_isNumeric));
	   logSdLogFsta = asVector<Type>(getListElement(x,"logSdLogFsta", &Rf_isNumeric));
	   muF = asVector<Type>(getListElement(x,"muF", &Rf_isNumeric));
	   trans_rho_F = asVector<Type>(getListElement(x,"trans_rho_F", &Rf_isNumeric));
	   boundF_kappa = asVector<Type>(getListElement(x,"boundF_kappa", &Rf_isNumeric));
	   boundF_alpha = asVector<Type>(getListElement(x,"boundF_alpha", &Rf_isNumeric));
	   boundF_tau = asVector<Type>(getListElement(x,"boundF_tau", &Rf_isNumeric));
	   logSdLogN = asVector<Type>(getListElement(x,"logSdLogN", &Rf_isNumeric));
	   logSdLogP = asVector<Type>(getListElement(x,"logSdLogP", &Rf_isNumeric));
	   logSdLogObs = asVector<Type>(getListElement(x,"logSdLogObs", &Rf_isNumeric));
	   logSdLogTotalObs = asVector<Type>(getListElement(x,"logSdLogTotalObs", &Rf_isNumeric));
	   transfIRARdist = asVector<Type>(getListElement(x,"transfIRARdist", &Rf_isNumeric));
	   sigmaObsParUS = asVector<Type>(getListElement(x,"sigmaObsParUS", &Rf_isNumeric));
	   rec_pars = asVector<Type>(getListElement(x,"rec_pars", &Rf_isNumeric));
	   rec_transphi = asVector<Type>(getListElement(x,"rec_transphi", &Rf_isNumeric));
	   itrans_rho = asVector<Type>(getListElement(x,"itrans_rho", &Rf_isNumeric));
	   rhop = asVector<Type>(getListElement(x,"rhop", &Rf_isNumeric));
	   logScale = asVector<Type>(getListElement(x,"logScale", &Rf_isNumeric));
	   logitReleaseSurvival = asVector<Type>(getListElement(x,"logitReleaseSurvival", &Rf_isNumeric));
	   logitRecapturePhi = asVector<Type>(getListElement(x,"logitRecapturePhi", &Rf_isNumeric));
	   logAlphaSCB = asVector<Type>(getListElement(x,"logAlphaSCB", &Rf_isNumeric));
	   sepFalpha = asVector<Type>(getListElement(x,"sepFalpha", &Rf_isNumeric));
	   sepFlogitRho = asVector<Type>(getListElement(x,"sepFlogitRho", &Rf_isNumeric));
	   sepFlogSd = asVector<Type>(getListElement(x,"sepFlogSd", &Rf_isNumeric));
	   logFScaleMSY = (Type)Rf_asReal(getListElement(x,"logFScaleMSY", &isNumericScalar));
	   implicitFunctionDelta = (Type)Rf_asReal(getListElement(x,"implicitFunctionDelta", &isNumericScalar));

	   logPhiSW = asVector<Type>(getListElement(x,"logPhiSW", &Rf_isNumeric)); 
	   logSdProcLogSW = asVector<Type>(getListElement(x,"logSdProcLogSW", &Rf_isNumeric));
	   meanLogSW  = asVector<Type>(getListElement(x,"meanLogSW", &Rf_isNumeric)); 
	   logSdLogSW = asVector<Type>(getListElement(x,"logSdLogSW", &Rf_isNumeric));
	   logPhiCW = asMatrix<Type>(getListElement(x,"logPhiCW", &Rf_isMatrix)); 
	   logSdProcLogCW = asVector<Type>(getListElement(x,"logSdProcLogCW", &Rf_isNumeric));
	   meanLogCW  = asVector<Type>(getListElement(x,"meanLogCW", &Rf_isNumeric)); 
	   logSdLogCW = asVector<Type>(getListElement(x,"logSdLogCW", &Rf_isNumeric));
	   logPhiMO = asVector<Type>(getListElement(x,"logPhiMO", &Rf_isNumeric)); 
	   logSdProcLogitMO = asVector<Type>(getListElement(x,"logSdProcLogitMO", &Rf_isNumeric));
	   meanLogitMO  = asVector<Type>(getListElement(x,"meanLogitMO", &Rf_isNumeric)); 
	   logSdMO = asVector<Type>(getListElement(x,"logSdMO", &Rf_isNumeric));
	   logPhiNM = asVector<Type>(getListElement(x,"logPhiNM", &Rf_isNumeric)); 
	   logSdProcLogNM = asVector<Type>(getListElement(x,"logSdProcLogNM", &Rf_isNumeric));
	   meanLogNM  = asVector<Type>(getListElement(x,"meanLogNM", &Rf_isNumeric)); 
	   logSdLogNM = asVector<Type>(getListElement(x,"logSdLogNM", &Rf_isNumeric));
	   logXtraSd = asVector<Type>(getListElement(x,"logXtraSd", &Rf_isNumeric));
	   initF = asVector<Type>(getListElement(x,"initF", &Rf_isNumeric));
	   initN = asVector<Type>(getListElement(x,"initN", &Rf_isNumeric));
	   seasonMu = asMatrix<Type>(getListElement(x,"seasonMu", &Rf_isMatrix));
	   seasonLogitRho = asVector<Type>(getListElement(x,"seasonLogitRho", &Rf_isNumeric));
	   seasonLogSd = asVector<Type>(getListElement(x,"seasonLogSd", &Rf_isNumeric));

	   splinePenalty = (Type)Rf_asReal(getListElement(x,"splinePenalty", &isNumericScalar));
	 }
	 );

SAM_SPECIALIZATION(struct paraSet<double>);
SAM_SPECIALIZATION(struct paraSet<TMBad::ad_aug>);

