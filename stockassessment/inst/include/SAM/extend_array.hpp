SAM_DEPENDS(logspace)


#ifndef WITH_SAM_LIB
namespace extend_fun {

  template <class Type>
  void extendArray_2D(array<Type>& x, int nModelYears, int nForecastYears, vector<int> aveYears, vector<Type> meanVec, vector<int> keyMeanVec, int meanType = 0, bool keepModelYears = true){
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

    SAM_ASSERT(nave > 0, "ave.years does not cover the data period.");
  
    ave /= nave;

    // Insert values in tmp
    for(int i = 0; i < tmp.dim(0); ++i){
      for(int j = 0; j < tmp.dim(1); ++j){
	if(keepModelYears && i < dim(0)){ // Take value from x
	  tmp(i,j) = x(i,j);
	}else if(meanVec.size() == 0){ // Take value from ave
	  tmp(i,j) = ave(j);
	}else if(meanType == 0){
	  tmp(i,j) = exp(meanVec(keyMeanVec(j)));
	}else if(meanType == 1){
	  tmp(i,j) = invlogit(meanVec(keyMeanVec(j)));
	}else{
	  Rf_error("Wrong mean type in extendArray");
	}
      }
    }
    // Overwrite x
    // NOTE: x must be resized first, otherwise x=tmp will not work.
    x.initZeroArray(tmp.dim);
    x = tmp;
    return;
  }

  SAM_SPECIALIZATION(void extendArray_2D(array<double>&, int, int, vector<int>, vector<double>, vector<int>, int, bool));
  SAM_SPECIALIZATION(void extendArray_2D(array<TMBad::ad_aug>&, int, int, vector<int>, vector<TMBad::ad_aug>, vector<int>, int, bool));

  template <class Type>
  void extendArray_3D(array<Type>& x, int nModelYears, int nForecastYears, vector<int> aveYears, vector<Type> meanVec, matrix<int> keyMeanVec, int meanType = 0,  bool keepModelYears = true){
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

    // if(nave == 0)
    //   Rf_error("ave.years does not cover the data period.");
    SAM_ASSERT(nave > 0, "ave.years does not cover the data period.");
  
    ave /= nave;

    // Insert values in tmp
    for(int i = 0; i < tmp.dim(0); ++i){
      for(int j = 0; j < tmp.dim(1); ++j){
	for(int k = 0; k < tmp.dim(2); ++k){
	  // if(keepModelYears && i < dim(0)){ // Take value from x
	  //   tmp(i,j,k) = x(i,j,k);
	  // }else{ // Take value from ave
	  //   tmp(i,j,k) = ave(j,k);
	  // }
	if(keepModelYears && i < dim(0)){ // Take value from x
	  tmp(i,j,k) = x(i,j,k);
	}else if(meanVec.size() == 0){ // Take value from ave
	  tmp(i,j,k) = ave(j,k);
	}else if(meanType == 0){
	  tmp(i,j,k) = exp(meanVec(keyMeanVec(k,j)));
	}else if(meanType == 1){
	  tmp(i,j,k) = invlogit(meanVec(keyMeanVec(k,j)));
	}else{
	  Rf_error("Wrong mean type in extendArray");
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

  SAM_SPECIALIZATION(void extendArray_3D(array<double>&, int, int, vector<int>, vector<double>, matrix<int>, int, bool));
  SAM_SPECIALIZATION(void extendArray_3D(array<TMBad::ad_aug>&, int, int, vector<int>, vector<TMBad::ad_aug>, matrix<int>, int, bool));


} // end of namespace forecast_fun
#endif



template <class Type>
void extendArray(array<Type>& x, int nModelYears, int nForecastYears, vector<int> aveYears, bool keepModelYears DEFARG(= true))SOURCE({
    vector<Type> a(0);
    matrix<int> b(0,0);
  if(x.dim.size() == 2){
    extend_fun::extendArray_2D(x, nModelYears, nForecastYears, aveYears, a, b.vec(), 0, keepModelYears);
  }else if(x.dim.size() == 3){
    extend_fun::extendArray_3D(x, nModelYears, nForecastYears, aveYears, a, b, 0, keepModelYears);
  }else{
    Rf_error("extendArray is only implemented for arrays of dimension 2 and 3");
  }
  return;
  }
  );

SAM_SPECIALIZATION(void extendArray(array<double>&, int, int, vector<int>, bool));
SAM_SPECIALIZATION(void extendArray(array<TMBad::ad_aug>&, int, int, vector<int>, bool));


template <class Type>
void extendArray(array<Type>& x, int nModelYears, int nForecastYears, vector<int> aveYears, vector<Type> meanVec, vector<int> keyMeanVec, int meanType, bool keepModelYears DEFARG(= true))SOURCE({
  if(x.dim.size() == 2){
    extend_fun::extendArray_2D(x, nModelYears, nForecastYears, aveYears, meanVec, keyMeanVec, meanType, keepModelYears);
  }else if(x.dim.size() == 3){
    Rf_warning("extendArray: dimensions of x and keyMeanVec does not match.");
    matrix<int> kmv(1,keyMeanVec.size());
    kmv.row(0) = keyMeanVec;
    extend_fun::extendArray_3D(x, nModelYears, nForecastYears, aveYears, meanVec, kmv, meanType, keepModelYears);
  }else{
    Rf_error("extendArray is only implemented for arrays of dimension 2 and 3");
  }
  return;
  }
  );

SAM_SPECIALIZATION(void extendArray(array<double>&, int, int, vector<int>, vector<double>, vector<int>, int, bool));
SAM_SPECIALIZATION(void extendArray(array<TMBad::ad_aug>&, int, int, vector<int>, vector<TMBad::ad_aug>, vector<int>, int, bool));



template <class Type>
void extendArray(array<Type>& x, int nModelYears, int nForecastYears, vector<int> aveYears, vector<Type> meanVec, matrix<int> keyMeanVec, int meanType, bool keepModelYears DEFARG(= true))SOURCE({
  if(x.dim.size() == 2){
    Rf_warning("extendArray: dimensions of x and keyMeanVec does not match.");
    extend_fun::extendArray_2D(x, nModelYears, nForecastYears, aveYears, meanVec, keyMeanVec.vec(), meanType, keepModelYears);
  }else if(x.dim.size() == 3){
    extend_fun::extendArray_3D(x, nModelYears, nForecastYears, aveYears, meanVec, keyMeanVec, meanType, keepModelYears);
  }else{
    Rf_error("extendArray is only implemented for arrays of dimension 2 and 3");
  }
  return;
  }
  );

SAM_SPECIALIZATION(void extendArray(array<double>&, int, int, vector<int>,vector<double>, matrix<int>, int, bool));
SAM_SPECIALIZATION(void extendArray(array<TMBad::ad_aug>&, int, int, vector<int>,vector<TMBad::ad_aug>, matrix<int>, int, bool));
