#define REPORT_F(name,F)					\
if(isDouble<Type>::value && F->current_parallel_region<0) {     \
  defineVar(install(#name),                                     \
            asSEXP_protect(name),F->report);                    \
  UNPROTECT(1);                                                 \
}

#define ADREPORT_F(name,F) F->reportvector.push(name,#name);

#define SIMULATE_F(F)						\
if(isDouble<Type>::value && F->do_simulate)

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

bool isNAINT(int x){
  return NA_INTEGER==x;
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
  array<int> aux;
  vector<Type> logobs;
  vector<Type> weight;
  //data_indicator<vector<Type>,Type> keep;
  array<Type> propMat;
  array<Type> stockMeanWeight; 
  array<Type> catchMeanWeight;
  array<Type> natMor;
  array<Type> landFrac;
  array<Type> disMeanWeight;
  array<Type> landMeanWeight;
  array<Type> propF;
  array<Type> propM;
};

struct confSet{
  int minAge;
  int maxAge;
  int maxAgePlusGroup;
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
  int noScaledYears;
  vector<int> keyScaledYears;
  array<int> keyParScaledYA;
  vector<int> fbarRange;
  vector<int> keyBiomassTreat;
  int simFlag; 
  int resFlag; 
  vector<int> obsLikelihoodFlag;
  int fixVarToWeight;
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
  vector<Type> rec_loga; 
  vector<Type> rec_logb; 
  vector<Type> itrans_rho; 
  vector<Type> logScale;
  vector<Type> logitReleaseSurvival;   
  vector<Type> logitRecapturePhi;   
};
