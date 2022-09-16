#pragma once
#ifndef SAM_SURVIVAL_HPP
#define SAM_SURVIVAL_HPP

// Implementation of DOI: 10.4054/DemRes.2013.29.41
// Formulas in appendix
// NOTE: Special cases for constant hazard-at-age!


////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Basic quantities ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Survival until the end of age a1 given survival until the end of age a0-1 assuming hazard rates from year y
template<class Type>
Type survivalFunction_i(dataSet<Type> &dat, confSet &conf, array<Type> &logF,
		      int y, int a0, int a1){
  // Survival until minAge
  Type CH = 0.0;
  for(int a = a0; a <= a1; ++a){
    int j = std::min(std::max(a-conf.minAge,0), dat.natMor.cols()-1);		// Cohort age index
    int i = std::min(y, dat.natMor.rows()-1);	// Cohort year index
    Type M = dat.natMor(i,j);
    Type F = 0.0;
    for(int f = 0; f < conf.keyLogFsta.dim(0); ++f)
      if(conf.keyLogFsta(f,j)>(-1) && a >= conf.minAge)
	F += exp(logF(conf.keyLogFsta(f,j),i));
    Type Z = M + F;
    CH += Z;		// Always one year
  }
  return exp(-CH);
}

// Temporary life expectancy between the start of age a0 and the end of age a1 assuming hazard rates of year y
template<class Type>
Type temporaryLifeExpectancy_i(dataSet<Type> &dat, confSet &conf, array<Type> &logF,
			     int y, int a0, int a1){
  Type le = 0.0;
  for(int a = a0; a <= a1; ++a){
    int j = std::min(std::max(a-conf.minAge,0), dat.natMor.cols()-1);		// Cohort age index
    int i = std::min(y, dat.natMor.rows()-1);	// Cohort year index
    Type M = dat.natMor(i,j);
    Type F = 0.0;
    for(int f = 0; f < conf.keyLogFsta.dim(0); ++f)
      if(conf.keyLogFsta(f,j)>(-1) && a >= conf.minAge)
	F += exp(logF(conf.keyLogFsta(f,j),i));
    Type Z = M + F;
    // Survival until the beginning of current age, i.e., end of last age (should be one for a=a0)
    Type p = survivalFunction_i(dat, conf, logF, y, a0, a-1);
    le += p / Z * (1.0 - exp(-Z));
  }
  return le;
}


// Probability of dying from fishing between the start of age a0 and the end of age a1 assuming hazard rates of year y
template<class Type>
Type cumulativeIncidenceFishing_i(dataSet<Type> &dat, confSet &conf, array<Type> &logF,
			     int y, int a0, int a1){
  Type q = 0.0;
  for(int a = a0; a <= a1; ++a){
    int j = std::min(std::max(a-conf.minAge,0), dat.natMor.cols()-1);		// Cohort age index
    int i = std::min(y, dat.natMor.rows()-1);	// Cohort year index
    Type M = dat.natMor(i,j);
    Type F = 0.0;
    for(int f = 0; f < conf.keyLogFsta.dim(0); ++f)
      if(conf.keyLogFsta(f,j)>(-1) && a >= conf.minAge)
	F += exp(logF(conf.keyLogFsta(f,j),i));
    Type Z = M + F;
    // Survival until the beginning of current age, i.e., end of last age (should be one for a=a0)
    Type p = survivalFunction_i(dat, conf, logF, y, a0, a-1);
    q += p * F / Z * (1.0 - exp(-Z));
  }
  return q;
}


// Probability of dying from other causes between the start of age a0 and the end of age a1 assuming hazard rates of year y
template<class Type>
Type cumulativeIncidenceOther_i(dataSet<Type> &dat, confSet &conf, array<Type> &logF,
			     int y, int a0, int a1){
  Type q = 0.0;
  for(int a = a0; a <= a1; ++a){
    int j = std::min(std::max(a-conf.minAge,0), dat.natMor.cols()-1);		// Cohort age index
    int i = std::min(y, dat.natMor.rows()-1);	// Cohort year index
    Type M = dat.natMor(i,j);
    Type F = 0.0;
    for(int f = 0; f < conf.keyLogFsta.dim(0); ++f)
      if(conf.keyLogFsta(f,j)>(-1) && a >= conf.minAge)
	F += exp(logF(conf.keyLogFsta(f,j),i));
    Type Z = M + F;
    // Survival until the beginning of current age, i.e., end of last age (should be one for a=a0)
    Type p = survivalFunction_i(dat, conf, logF, y, a0, a-1);
    q += p * M / Z * (1.0 - exp(-Z));
  }
  return q;
}

template<class Type>
Type yearsLostFishing_i(dataSet<Type> &dat, confSet &conf, array<Type> &logF,
			     int y, int a0, int a1){
  Type yl = 0.0;
  for(int a = a0; a <= a1; ++a){
    int j = std::min(std::max(a-conf.minAge,0), dat.natMor.cols()-1);		// Cohort age index
    int i = std::min(y, dat.natMor.rows()-1);	// Cohort year index
    Type M = dat.natMor(i,j);
    Type F = 0.0;
    for(int f = 0; f < conf.keyLogFsta.dim(0); ++f)
      if(conf.keyLogFsta(f,j)>(-1) && a >= conf.minAge)
	F += exp(logF(conf.keyLogFsta(f,j),i));
    Type Z = M + F;
    // Survival until the beginning of current age, i.e., end of last age (should be one for a=a0)
    Type p = survivalFunction_i(dat, conf, logF, y, a0, a-1);
    Type q = cumulativeIncidenceFishing_i(dat, conf, logF, y, a0, a-1);
    yl += q + p * F / Z * (1.0 - 1.0 / Z * (1.0 - exp(-Z)));
  }
  return yl;
}


template<class Type>
Type yearsLostOther_i(dataSet<Type> &dat, confSet &conf, array<Type> &logF,
			     int y, int a0, int a1){
  Type yl = 0.0;
  for(int a = a0; a <= a1; ++a){
    int j = std::min(std::max(a-conf.minAge,0), dat.natMor.cols()-1);		// Cohort age index
    int i = std::min(y, dat.natMor.rows()-1);	// Cohort year index
    Type M = dat.natMor(i,j);
    Type F = 0.0;
    for(int f = 0; f < conf.keyLogFsta.dim(0); ++f)
      if(conf.keyLogFsta(f,j)>(-1) && a >= conf.minAge)
	F += exp(logF(conf.keyLogFsta(f,j),i));
    Type Z = M + F;
    // Survival until the beginning of current age, i.e., end of last age (should be one for a=a0)
    Type p = survivalFunction_i(dat, conf, logF, y, a0, a-1);
    Type q = cumulativeIncidenceOther_i(dat, conf, logF, y, a0, a-1);
    yl += q + p * M / Z * (1.0 - 1.0 / Z * (1.0 - exp(-Z)));
  }
  return yl;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////// Reporting functions /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Life expectanxy at birth (assuming maximum attainable age is no more than 10 times model maxAge)
template<class Type>
vector<Type> lifeexpectancy(dataSet<Type> &dat, confSet &conf, array<Type> &logF){
  int timeSteps = dat.natMor.dim(0);
  vector<Type> v(timeSteps);
  v.setZero();
  for(int y = 0; y < timeSteps; ++y)
    v(y) = temporaryLifeExpectancy_i(dat, conf, logF, y, 0, 10 * conf.maxAge);
  return v;
}

// Life expectancy at recruitment (assuming maximum attainable age is no more than 10 times model maxAge)
template<class Type>
vector<Type> lifeexpectancyRec(dataSet<Type> &dat, confSet &conf, array<Type> &logF){
  int timeSteps = dat.natMor.dim(0);
  vector<Type> v(timeSteps);
  v.setZero();
  for(int y = 0; y < timeSteps; ++y)
    v(y) = temporaryLifeExpectancy_i(dat, conf, logF, y, conf.minAge, 10 * conf.maxAge) + conf.minAge;
  return v;
}



// Life expectancy at age matrix (assuming maximum attainable age is no more than 10 times model maxAge)
template<class Type>
matrix<Type> lifeexpectancyAge(dataSet<Type> &dat, confSet &conf, array<Type> &logF, bool give_log = false, bool predicted = true){
  int timeSteps = dat.natMor.dim(0);
  matrix<Type> v(timeSteps, conf.maxAge+1);
  v.setZero();
  for(int y = 0; y < timeSteps; ++y)
    for(int a = 0; a < v.cols(); ++a)
      v(y,a) = temporaryLifeExpectancy_i(dat, conf, logF, y, a, 10 * conf.maxAge) + a;
  return v;
}



// Years lost to fishing between the start of age minAge and the end of age maxAge
template<class Type>
vector<Type> yearsLostFishing(dataSet<Type> &dat, confSet &conf, array<Type> &logF){
  int timeSteps = dat.natMor.dim(0);
  vector<Type> v(timeSteps);
  v.setZero();
  for(int y = 0; y < timeSteps; ++y)
    v(y) = yearsLostFishing_i(dat, conf, logF, y, conf.minAge, conf.maxAge);
  return v;
}


// Years lost to other causes between the start of age minAge and the end of age maxAge
template<class Type>
vector<Type> yearsLostOther(dataSet<Type> &dat, confSet &conf, array<Type> &logF){
  int timeSteps = dat.natMor.dim(0);
  vector<Type> v(timeSteps);
  v.setZero();
  for(int y = 0; y < timeSteps; ++y)
    v(y) = yearsLostOther_i(dat, conf, logF, y, conf.minAge, conf.maxAge);
  return v;
}


// Temporary life expectancy (for balance) between the start of age minAge and the end of age maxAge
template<class Type>
vector<Type> temporaryLifeExpectancy(dataSet<Type> &dat, confSet &conf, array<Type> &logF){
  int timeSteps = dat.natMor.dim(0);
  vector<Type> v(timeSteps);
  v.setZero();
  for(int y = 0; y < timeSteps; ++y)
    v(y) = temporaryLifeExpectancy_i(dat, conf, logF, y, conf.minAge, conf.maxAge);
  return v;
}



#endif	
