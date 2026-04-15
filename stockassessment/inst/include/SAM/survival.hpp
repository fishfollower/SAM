SAM_DEPENDS(convenience)
SAM_DEPENDS(define)
SAM_DEPENDS(incidence)

// Implementation of DOI: 10.4054/DemRes.2013.29.41
// Formulas in appendix
// NOTE: Special cases for constant hazard-at-age!
// Actual calculations are in MortalitySet


////////////////////////////////////////////////////////////////////////////////
////////////////////////// Reporting functions /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Life expectanxy at birth (assuming maximum attainable age is no more than 10 times model maxAge and F=0 before first model age)
template<class Type>
vector<Type> loglifeexpectancy(MortalitySet<Type> &mort)SOURCE({
    int timeSteps = mort.maxYear - mort.minYear + 1;
    vector<Type> v(timeSteps);
    v.setZero();
    for(int y = 0; y < timeSteps; ++y)
      v(y) = mort.logTemporaryLifeExpectancy(y, 0, 10 * mort.maxAge+1);
    return v;
  })


  SAM_SPECIALIZATION(vector<double> loglifeexpectancy(MortalitySet<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> loglifeexpectancy(MortalitySet<TMBad::ad_aug>&));


// Life expectancy at recruitment (assuming maximum attainable age is no more than 10 times model maxAge)
template<class Type>
vector<Type> loglifeexpectancyRec(MortalitySet<Type> &mort)SOURCE({
    int timeSteps = mort.maxYear - mort.minYear + 1;
    vector<Type> v(timeSteps);
    v.setZero();
    for(int y = 0; y < timeSteps; ++y)
      v(y) = logspace_add_SAM(mort.logTemporaryLifeExpectancy(y, mort.minAge, 10 * mort.maxAge+1.0), (Type)log(mort.minAge));
    return v;
  })


  SAM_SPECIALIZATION(vector<double> loglifeexpectancyRec(MortalitySet<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> loglifeexpectancyRec(MortalitySet<TMBad::ad_aug>&));


// Life expectancy at age matrix (assuming maximum attainable age is no more than 10 times model maxAge)
template<class Type>
matrix<Type> loglifeexpectancyAge(MortalitySet<Type>& mort, bool give_log DEFARG(= false), bool predicted DEFARG(= true))SOURCE({
    int timeSteps = mort.maxYear - mort.minYear + 1;
    matrix<Type> v(timeSteps, mort.maxAge+1);
    v.setZero();
    for(int y = 0; y < timeSteps; ++y)
      for(int a = 0; a < v.cols(); ++a)
	v(y,a) = logspace_add_SAM(mort.logTemporaryLifeExpectancy(y, mort.minAge, 10 * mort.maxAge+1.0), (Type)log(a));
    return v;
  })


  SAM_SPECIALIZATION(matrix<double> loglifeexpectancyAge(MortalitySet<double>&, bool, bool));
SAM_SPECIALIZATION(matrix<TMBad::ad_aug> loglifeexpectancyAge(MortalitySet<TMBad::ad_aug>&, bool, bool));


// Years lost to fishing between the start of age minAge and the end of age maxAge
template<class Type>
vector<Type> logYearsLostFishing(MortalitySet<Type> &mort)SOURCE({
    int timeSteps = mort.maxYear - mort.minYear + 1;
    vector<Type> v(timeSteps);
    v.setZero();
    for(int y = 0; y < timeSteps; ++y)
      v(y) = mort.logYearsLostFishing(-1, y, mort.minAge, mort.maxAge+1.0);
    return v;
  })

  SAM_SPECIALIZATION(vector<double> logYearsLostFishing(MortalitySet<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> logYearsLostFishing(MortalitySet<TMBad::ad_aug>&));

template<class Type>
matrix<Type> logYearsLostFishingFleet(MortalitySet<Type> &mort)SOURCE({
    int timeSteps = mort.maxYear - mort.minYear + 1;
    matrix<Type> v(timeSteps, mort.logCIF_F_breakpoints.dim(2));
    v.setZero();
    for(int y = 0; y < timeSteps; ++y)
      for(int f = 0; f < v.cols(); ++f){
	v(y,f) = mort.logYearsLostFishing(f, y, mort.minAge, mort.maxAge+1.0);
      }
    return v;
  })

  SAM_SPECIALIZATION(matrix<double> logYearsLostFishingFleet(MortalitySet<double>&));
SAM_SPECIALIZATION(matrix<TMBad::ad_aug> logYearsLostFishingFleet(MortalitySet<TMBad::ad_aug>&));



// Years lost to other causes between the start of age minAge and the end of age maxAge
template<class Type>
vector<Type> logYearsLostOther(MortalitySet<Type> &mort)SOURCE({
   int timeSteps = mort.maxYear - mort.minYear + 1;
    vector<Type> v(timeSteps);
    v.setZero();
    for(int y = 0; y < timeSteps; ++y)
      v(y) = mort.logYearsLostOther(-1, y, mort.minAge, mort.maxAge+1.0);
    return v;
  })

  SAM_SPECIALIZATION(vector<double> logYearsLostOther(MortalitySet<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> logYearsLostOther(MortalitySet<TMBad::ad_aug>&));




// Years lost to other causes between the start of age minAge and the end of age maxAge
template<class Type>
matrix<Type> logYearsLostOtherRisk(MortalitySet<Type> &mort)SOURCE({
   int timeSteps = mort.maxYear - mort.minYear + 1;
   matrix<Type> v(timeSteps, mort.logCIF_M_breakpoints.dim(2));
    v.setZero();
    for(int y = 0; y < timeSteps; ++y)
      for(int r = 0; r < v.cols(); ++r)
	v(y,r) = mort.logYearsLostOther(r, y, mort.minAge, mort.maxAge+1.0);
    return v;
  })

  SAM_SPECIALIZATION(matrix<double> logYearsLostOtherRisk(MortalitySet<double>&));
  SAM_SPECIALIZATION(matrix<TMBad::ad_aug> logYearsLostOtherRisk(MortalitySet<TMBad::ad_aug>&));


// Temporary life expectancy (for balance) between the start of age minAge and the end of age maxAge
template<class Type>
vector<Type> logTemporaryLifeExpectancy(MortalitySet<Type> &mort)SOURCE({
    int timeSteps = mort.maxYear - mort.minYear + 1;
    vector<Type> v(timeSteps);
    v.setZero();
    for(int y = 0; y < timeSteps; ++y)
      v(y) = mort.logTemporaryLifeExpectancy(y, mort.minAge, mort.maxAge + 1.0);
    return v;
  })


  SAM_SPECIALIZATION(vector<double> logTemporaryLifeExpectancy(MortalitySet<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> logTemporaryLifeExpectancy(MortalitySet<TMBad::ad_aug>&));


