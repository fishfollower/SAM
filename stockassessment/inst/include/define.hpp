#define REPORT_F(name,F)					\
if(isDouble<Type>::value && F->current_parallel_region<0) {     \
  defineVar(install(#name),                                     \
            PROTECT(asSEXP(name)),F->report);			\
  UNPROTECT(1);                                                 \
}


//This function returns a vector with matrices based on a list of matrices from R
template<class Type>
struct listMatrixFromR : vector<matrix<Type> > {
  listMatrixFromR(SEXP x){ 
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asMatrix<Type>(sm);
    }
  }
};

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
  vector<matrix<Type> > corList;
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
    return *this;
  };
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
  vector<double> constRecBreaks;
  int noScaledYears;
  vector<int> keyScaledYears;
  matrix<int> keyParScaledYA;
  vector<int> fbarRange;
  vector<int> keyBiomassTreat;
  int simFlag; 
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
    maxAgePlusGroup = (int)*REAL(getListElement(x,"maxAgePlusGroup"));
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
    simFlag = (int)*REAL(getListElement(x,"simFlag"));
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
  vector<Type> rec_loga; 
  vector<Type> rec_logb; 
  vector<Type> itrans_rho; 
  vector<Type> logScale;
  vector<Type> logitReleaseSurvival;   
  vector<Type> logitRecapturePhi;   
  vector<Type> sepFalpha;   
  vector<Type> sepFlogitRho;   
  vector<Type> sepFlogSd;   
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
