SAM_DEPENDS(define)
SAM_DEPENDS(recruitment)

#ifndef WITH_SAM_LIB

namespace reproductive_fun {

  template<class Type>
  SparseMatrix<Type> Leslie_i(dataSet<Type> &dat, confSet &conf, paraSet<Type>& par, Recruitment<Type>& recruit, int y){
  
    Type max_dSR = recruit.maxGradient(); // For ricker and BH
    int nAges = conf.maxAge - conf.minAge + 1;

    typedef Triplet<Type> T;
    std::vector<T> tripletList;  
    tripletList.reserve(nAges + (nAges-1));
    // Fecundity in first row
    for(int i = 0; i < nAges; ++i){
      Type f = max_dSR * (dat.stockMeanWeight(y,i)) * (dat.propMat(y,i));
      tripletList.push_back(T(0,i,f));
    }
    // Survival from age class i-1 to i in absence of fishing
    for(int i = 1; i < nAges; ++i){
      Type Mi = dat.natMor(y,i-1);
      tripletList.push_back(T(i,i-1,exp(-Mi)));
    }
    // Plus group
    if(conf.maxAgePlusGroup(0)==1){    
      Type Mi = dat.natMor(y,nAges-1);
      tripletList.push_back(T(nAges-1,nAges-1,exp(-Mi)));
    }
    SparseMatrix<Type> m(nAges,nAges);
    m.setFromTriplets(tripletList.begin(), tripletList.end());
    return m;
  }

  SAM_SPECIALIZATION(SparseMatrix<double> Leslie_i(dataSet<double>&, confSet&, paraSet<double>&, Recruitment<double>&, int));
  SAM_SPECIALIZATION(SparseMatrix<TMBad::ad_aug> Leslie_i(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, Recruitment<TMBad::ad_aug>&, int));

  template<class Type>
  Type norm2(vector<Type> x){
    Type r = 0.0;
    for(int i = 0; i < x.size(); ++i)
      r += x(i) * x(i);
    return sqrt(r);
  }

  SAM_SPECIALIZATION(double norm2(vector<double> x));
  SAM_SPECIALIZATION(TMBad::ad_aug norm2(vector<TMBad::ad_aug> x));

  template<class Type>
  Type maxEigenValue(SparseMatrix<Type> A){
    // Power iteration with fixed number of iterations
    // Assuming A is square
    vector<Type> b(A.cols());
    b.setConstant(1.0);
    for(int i = 0; i < 3 * b.size(); ++i){
      vector<Type> v = A * b;
      b = v / norm2(v);
    }
    vector<Type> v = A * b;
    return (b * v).sum();
  }

  SAM_SPECIALIZATION(double maxEigenValue(SparseMatrix<double> x));
  SAM_SPECIALIZATION(TMBad::ad_aug maxEigenValue(SparseMatrix<TMBad::ad_aug> x));


  template<class Type>
  Type rmax_i(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, Recruitment<Type>& recruit, int y){
    // For Ricker, BH and any CMP model, it is the derivative at 0, but is the largest derivative the correct generalization? What should it be for constant recruitment models? 
    // Implement for other models
    if(conf.stockRecruitmentModelCode == 0 ||
       conf.stockRecruitmentModelCode == 3 ||
       conf.stockRecruitmentModelCode == 62 ||
       conf.stockRecruitmentModelCode == 64 ||
       conf.stockRecruitmentModelCode == 65){
      return R_NaReal;
    }
    Eigen::SparseMatrix<Type> L = Leslie_i(dat,conf,par, recruit,y);
    Type lambda = maxEigenValue(L);
    return log(lambda);
  }

  SAM_SPECIALIZATION(double rmax_i(dataSet<double>&, confSet&, paraSet<double>&, Recruitment<double>&, int));
  SAM_SPECIALIZATION(TMBad::ad_aug rmax_i(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, Recruitment<TMBad::ad_aug>&, int));


  template<class Type>
  Type generationLength_i(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, int y){
    // Following https://www.fao.org/3/a0212e/A0212E13.htm
    Type max_dSR = 1.0; // Cancels out! maxDiffSR(dat,conf,par,y); // For ricker and BH
    Type R0 = 0.0;
    Type Gx = 0.0;
    int nAges = conf.maxAge - conf.minAge + 1;

    Type loglx = 0.0;
    for(int i = 0; i < nAges; ++i){
      loglx -= dat.natMor(y,i);
      Type f = max_dSR * (dat.stockMeanWeight(y,i)) * (dat.propMat(y,i));
      R0 += exp(loglx) * f;
      Gx += exp(loglx) * f * Type(i + conf.minAge);
    }  
    return Gx / R0;
  }


  SAM_SPECIALIZATION(double generationLength_i(dataSet<double>&, confSet&, paraSet<double>&, int));
  SAM_SPECIALIZATION(TMBad::ad_aug generationLength_i(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&,int));

  
}

#endif


template<class Type>
vector<Type> rmax(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, Recruitment<Type>& recruit)SOURCE({
  int timeSteps = dat.natMor.dim(0);
  vector<Type> v(timeSteps);
  v.setZero();
  for(int y = 0; y < timeSteps; ++y)
    v(y) = reproductive_fun::rmax_i(dat, conf, par, recruit, y);
  return v;
  });

SAM_SPECIALIZATION(vector<double> rmax(dataSet<double>&, confSet&, paraSet<double>&, Recruitment<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> rmax(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, Recruitment<TMBad::ad_aug>&));


template<class Type>
vector<Type> generationLength(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par)SOURCE({
  int timeSteps = dat.natMor.dim(0);
  vector<Type> v(timeSteps);
  v.setZero();
  for(int y = 0; y < timeSteps; ++y)
    v(y) = reproductive_fun::generationLength_i(dat, conf, par, y);
  return v;
  });



SAM_SPECIALIZATION(vector<double> generationLength(dataSet<double>&, confSet&, paraSet<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> generationLength(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&));
