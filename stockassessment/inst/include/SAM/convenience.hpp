SAM_DEPENDS(logspace)

template<class Type>
bool isNA(Type x) SOURCE({
  return R_IsNA(asDouble(x));
  } )

SAM_SPECIALIZATION(bool isNA(double));
SAM_SPECIALIZATION(bool isNA(TMBad::ad_aug));

  
  bool isNAINT(int x) SOURCE({
  return R_NaInt == x; //NA_INTEGER==x;
    } )

  
template<class Type>
Type squash(Type u) SOURCE( {
    Type eps = 1.0e-6;
    u = (1.0 - eps) * (u - .5) + .5;
    return u;
  } );

SAM_SPECIALIZATION(double squash(double));
SAM_SPECIALIZATION(TMBad::ad_aug squash(TMBad::ad_aug));



template<class Type>
Type logspace_add_p (Type logx, Type logy, Type p) SOURCE({
  return log((Type(1)-p)*exp(logy-logx)+p)+logx; // the order of x and y is taylored for this application 
  });

SAM_SPECIALIZATION(double logspace_add_p(double, double, double));
SAM_SPECIALIZATION(TMBad::ad_aug logspace_add_p(TMBad::ad_aug, TMBad::ad_aug, TMBad::ad_aug));


template<class Type>
Type logdrobust(Type x, Type p)SOURCE({
  Type ld1=dnorm(x,Type(0.0),Type(1.0),true);
  if(p<Type(1.0e-16)){
    return ld1;
  }else{
    Type ld2=dt(x,Type(3),true);
    Type logres=logspace_add_p(ld2,ld1,p);
    return logres;
  }
  })

  
template<class Type>
  vector<Type> logdrobust(vector<Type>& x, vector<Type>& p)SOURCE({
      vector<Type> r(x.size());
      r.setZero();
      SAM_ASSERT(x.size() == p.size(), "wrong sizes in logdrobust");
      for(int i = 0; i < r.size(); ++i)
	r(i) = logdrobust(x(i),p(i));
      return r;
  })



SAM_SPECIALIZATION(double logdrobust(double, double));
SAM_SPECIALIZATION(TMBad::ad_aug logdrobust(TMBad::ad_aug, TMBad::ad_aug));
// Vectorized versions
SAM_SPECIALIZATION(vector<double> logdrobust(vector<double>&, vector<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> logdrobust(vector<TMBad::ad_aug>&, vector<TMBad::ad_aug>&));



/*
  Map from real line to the interval (a,b). d determines the steepness.  
 */
template <class Type>
Type toInterval(Type x, Type a, Type b, Type d) SOURCE({
    return (b - a)/(Type(1.0) + exp(-d * x)) + a;
  })

  SAM_SPECIALIZATION(double toInterval(double, double, double, double));
SAM_SPECIALIZATION(TMBad::ad_aug toInterval(TMBad::ad_aug, TMBad::ad_aug, TMBad::ad_aug, TMBad::ad_aug));


vector<int> getCatchFleets(vector<int> fleetTypes)SOURCE({
  std::vector<int> r;
  for(int i = 0; i < (int)fleetTypes.size(); ++i)
    if(fleetTypes(i) == 0)
      r.push_back(i);
  return vector<int>(r);
  });


template<class Type>
Type softmax(Type x, Type y, Type k DEFARG(=1.0))SOURCE({
  return logspace_add_SAM(k * x, k * y) / k;
  });

SAM_SPECIALIZATION(double softmax(double, double, double));
SAM_SPECIALIZATION(TMBad::ad_aug softmax(TMBad::ad_aug, TMBad::ad_aug, TMBad::ad_aug));
