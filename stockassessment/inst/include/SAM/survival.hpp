SAM_DEPENDS(convenience)
SAM_DEPENDS(define)


// Implementation of DOI: 10.4054/DemRes.2013.29.41
// Formulas in appendix
// NOTE: Special cases for constant hazard-at-age!


////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Basic quantities ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// Survival until the end of age a1 given survival until the end of age a0-1 assuming hazard rates from year y
template<class Type>
Type survivalFunction_i(dataSet<Type> &dat, confSet &conf, array<Type> &logF,
			int y, int a0, int a1)SOURCE({
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
			  })

SAM_SPECIALIZATION(double survivalFunction_i(dataSet<double>&, confSet&, array<double>&, int, int, int));
SAM_SPECIALIZATION(TMBad::ad_aug survivalFunction_i(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, int, int, int));

// Temporary life expectancy between the start of age a0 and the end of age a1 assuming hazard rates of year y
template<class Type>
Type temporaryLifeExpectancy_i(dataSet<Type> &dat, confSet &conf, array<Type> &logF,
			       int y, int a0, int a1)SOURCE({
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
				 })

SAM_SPECIALIZATION(double temporaryLifeExpectancy_i(dataSet<double>&, confSet&, array<double>&, int, int, int));
SAM_SPECIALIZATION(TMBad::ad_aug temporaryLifeExpectancy_i(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, int, int, int));

// Probability of dying from fishing by fleet f between the start of age a0 and the end of age a1 assuming hazard rates of year y
template<class Type>
Type cumulativeIncidenceFishingFleet_i(dataSet<Type> &dat, confSet &conf, array<Type> &logF,
				       int y, int a0, int a1, int f)SOURCE({
					   Type q = 0.0;
					   for(int a = a0; a <= a1; ++a){
					     int j = std::min(std::max(a-conf.minAge,0), dat.natMor.cols()-1);		// Cohort age index
					     int i = std::min(y, dat.natMor.rows()-1);	// Cohort year index
					     Type M = dat.natMor(i,j);
					     Type F = 0.0;
					     if(conf.keyLogFsta(f,j)>(-1) && a >= conf.minAge)
					       F += exp(logF(conf.keyLogFsta(f,j),i));
					     Type Z = M;
					     for(int ff = 0; ff < conf.keyLogFsta.dim(0); ++ff)
					       if(conf.keyLogFsta(ff,j)>(-1) && a >= conf.minAge)
						 Z += exp(logF(conf.keyLogFsta(ff,j),i));
					     // Survival until the beginning of current age, i.e., end of last age (should be one for a=a0)
					     Type p = survivalFunction_i(dat, conf, logF, y, a0, a-1);
					     q += p * F / Z * (1.0 - exp(-Z));
					   }
					   return q;
					 })


SAM_SPECIALIZATION(double cumulativeIncidenceFishingFleet_i(dataSet<double>&, confSet&, array<double>&, int, int, int, int));
SAM_SPECIALIZATION(TMBad::ad_aug cumulativeIncidenceFishingFleet_i(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, int, int, int, int));


// Probability of dying from fishing between the start of age a0 and the end of age a1 assuming hazard rates of year y
template<class Type>
Type cumulativeIncidenceFishing_i(dataSet<Type> &dat, confSet &conf, array<Type> &logF,
				  int y, int a0, int a1)SOURCE({
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
				    })


SAM_SPECIALIZATION(double cumulativeIncidenceFishing_i(dataSet<double>&, confSet&, array<double>&, int, int, int));
SAM_SPECIALIZATION(TMBad::ad_aug cumulativeIncidenceFishing_i(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, int, int, int));



// Probability of dying from other causes between the start of age a0 and the end of age a1 assuming hazard rates of year y
template<class Type>
Type cumulativeIncidenceOther_i(dataSet<Type> &dat, confSet &conf, array<Type> &logF,
				int y, int a0, int a1)SOURCE({
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
				  })


SAM_SPECIALIZATION(double cumulativeIncidenceOther_i(dataSet<double>&, confSet&, array<double>&, int, int, int));
SAM_SPECIALIZATION(TMBad::ad_aug cumulativeIncidenceOther_i(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, int, int, int));


template<class Type>
Type yearsLostFishingFleet_i(dataSet<Type> &dat, confSet &conf, array<Type> &logF,
			     int y, int a0, int a1, int f)SOURCE({
				 Type yl = 0.0;
				 for(int a = a0; a <= a1; ++a){
				   int j = std::min(std::max(a-conf.minAge,0), dat.natMor.cols()-1);		// Cohort age index
				   int i = std::min(y, dat.natMor.rows()-1);	// Cohort year index
				   Type M = dat.natMor(i,j);
				   Type F = 0.0;
				   if(conf.keyLogFsta(f,j)>(-1) && a >= conf.minAge)
				     F += exp(logF(conf.keyLogFsta(f,j),i));
				   Type Z = M;
				   for(int ff = 0; ff < conf.keyLogFsta.dim(0); ++ff)
				     if(conf.keyLogFsta(ff,j)>(-1) && a >= conf.minAge)
				       Z += exp(logF(conf.keyLogFsta(ff,j),i));
				   // Survival until the beginning of current age, i.e., end of last age (should be one for a=a0)
				   Type p = survivalFunction_i(dat, conf, logF, y, a0, a-1);
				   Type q = cumulativeIncidenceFishing_i(dat, conf, logF, y, a0, a-1);
				   yl += q + p * F / Z * (1.0 - 1.0 / Z * (1.0 - exp(-Z)));
				 }
				 return yl;
			       })


SAM_SPECIALIZATION(double yearsLostFishingFleet_i(dataSet<double>&, confSet&, array<double>&, int, int, int, int));
SAM_SPECIALIZATION(TMBad::ad_aug yearsLostFishingFleet_i(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, int, int, int, int));



template<class Type>
Type yearsLostFishing_i(dataSet<Type> &dat, confSet &conf, array<Type> &logF,
			int y, int a0, int a1)SOURCE({
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
			  })


SAM_SPECIALIZATION(double yearsLostFishing_i(dataSet<double>&, confSet&, array<double>&, int, int, int));
SAM_SPECIALIZATION(TMBad::ad_aug yearsLostFishing_i(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, int, int, int));


template<class Type>
Type yearsLostOther_i(dataSet<Type> &dat, confSet &conf, array<Type> &logF,
		      int y, int a0, int a1)SOURCE({
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
			})


SAM_SPECIALIZATION(double yearsLostOther_i(dataSet<double>&, confSet&, array<double>&, int, int, int));
SAM_SPECIALIZATION(TMBad::ad_aug yearsLostOther_i(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, int, int, int));



////////////////////////////////////////////////////////////////////////////////
////////////////////////// Reporting functions /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Life expectanxy at birth (assuming maximum attainable age is no more than 10 times model maxAge)
template<class Type>
vector<Type> lifeexpectancy(dataSet<Type> &dat, confSet &conf, array<Type> &logF)SOURCE({
    int timeSteps = dat.natMor.dim(0);
    vector<Type> v(timeSteps);
    v.setZero();
    for(int y = 0; y < timeSteps; ++y)
      v(y) = temporaryLifeExpectancy_i(dat, conf, logF, y, 0, 10 * conf.maxAge);
    return v;
  })


  SAM_SPECIALIZATION(vector<double> lifeexpectancy(dataSet<double>&, confSet&, array<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> lifeexpectancy(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&));


// Life expectancy at recruitment (assuming maximum attainable age is no more than 10 times model maxAge)
template<class Type>
vector<Type> lifeexpectancyRec(dataSet<Type> &dat, confSet &conf, array<Type> &logF)SOURCE({
    int timeSteps = dat.natMor.dim(0);
    vector<Type> v(timeSteps);
    v.setZero();
    for(int y = 0; y < timeSteps; ++y)
      v(y) = temporaryLifeExpectancy_i(dat, conf, logF, y, conf.minAge, 10 * conf.maxAge) + conf.minAge;
    return v;
  })


  SAM_SPECIALIZATION(vector<double> lifeexpectancyRec(dataSet<double>&, confSet&, array<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> lifeexpectancyRec(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&));


// Life expectancy at age matrix (assuming maximum attainable age is no more than 10 times model maxAge)
template<class Type>
matrix<Type> lifeexpectancyAge(dataSet<Type> &dat, confSet &conf, array<Type> &logF, bool give_log DEFARG(= false), bool predicted DEFARG(= true))SOURCE({
    int timeSteps = dat.natMor.dim(0);
    matrix<Type> v(timeSteps, conf.maxAge+1);
    v.setZero();
    for(int y = 0; y < timeSteps; ++y)
      for(int a = 0; a < v.cols(); ++a)
	v(y,a) = temporaryLifeExpectancy_i(dat, conf, logF, y, a, 10 * conf.maxAge) + a;
    return v;
  })


  SAM_SPECIALIZATION(matrix<double> lifeexpectancyAge(dataSet<double>&, confSet&, array<double>&, bool, bool));
SAM_SPECIALIZATION(matrix<TMBad::ad_aug> lifeexpectancyAge(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&, bool, bool));


// Years lost to fishing between the start of age minAge and the end of age maxAge
template<class Type>
vector<Type> yearsLostFishing(dataSet<Type> &dat, confSet &conf, array<Type> &logF)SOURCE({
    int timeSteps = dat.natMor.dim(0);
    vector<Type> v(timeSteps);
    v.setZero();
    for(int y = 0; y < timeSteps; ++y)
      v(y) = yearsLostFishing_i(dat, conf, logF, y, conf.minAge, conf.maxAge);
    return v;
  })

  SAM_SPECIALIZATION(vector<double> yearsLostFishing(dataSet<double>&, confSet&, array<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> yearsLostFishing(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&));

template<class Type>
matrix<Type> yearsLostFishingFleet(dataSet<Type> &dat, confSet &conf, array<Type> &logF)SOURCE({
    vector<int> cFleets = getCatchFleets(dat.fleetTypes);
    int timeSteps = dat.natMor.dim(0);
    matrix<Type> v(timeSteps, cFleets.size());
    v.setZero();
    for(int y = 0; y < timeSteps; ++y)
      for(int j = 0; j < cFleets.size(); ++j){
	int f = cFleets[j];
	v(y,j) = yearsLostFishingFleet_i(dat, conf, logF, y, conf.minAge, conf.maxAge,f);
      }
    return v;
  })

  SAM_SPECIALIZATION(matrix<double> yearsLostFishingFleet(dataSet<double>&, confSet&, array<double>&));
SAM_SPECIALIZATION(matrix<TMBad::ad_aug> yearsLostFishingFleet(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&));



// Years lost to other causes between the start of age minAge and the end of age maxAge
template<class Type>
vector<Type> yearsLostOther(dataSet<Type> &dat, confSet &conf, array<Type> &logF)SOURCE({
    int timeSteps = dat.natMor.dim(0);
    vector<Type> v(timeSteps);
    v.setZero();
    for(int y = 0; y < timeSteps; ++y)
      v(y) = yearsLostOther_i(dat, conf, logF, y, conf.minAge, conf.maxAge);
    return v;
  })

  SAM_SPECIALIZATION(vector<double> yearsLostOther(dataSet<double>&, confSet&, array<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> yearsLostOther(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&));


// Temporary life expectancy (for balance) between the start of age minAge and the end of age maxAge
template<class Type>
vector<Type> temporaryLifeExpectancy(dataSet<Type> &dat, confSet &conf, array<Type> &logF)SOURCE({
    int timeSteps = dat.natMor.dim(0);
    vector<Type> v(timeSteps);
    v.setZero();
    for(int y = 0; y < timeSteps; ++y)
      v(y) = temporaryLifeExpectancy_i(dat, conf, logF, y, conf.minAge, conf.maxAge);
    return v;
  })


  SAM_SPECIALIZATION(vector<double> temporaryLifeExpectancy(dataSet<double>&, confSet&, array<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> temporaryLifeExpectancy(dataSet<TMBad::ad_aug>&, confSet&, array<TMBad::ad_aug>&));

