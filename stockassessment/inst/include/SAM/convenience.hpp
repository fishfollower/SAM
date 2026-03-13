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

template<class Type>
Type softabs(Type x, Type k DEFARG(=0.01))SOURCE({
    return sqrt((x) * (x) + k);
  });

SAM_SPECIALIZATION(double softabs(double, double));
SAM_SPECIALIZATION(TMBad::ad_aug softabs(TMBad::ad_aug, TMBad::ad_aug));


template<class Type>
Type softmax2(Type x, Type y, Type k DEFARG(=0.01))SOURCE({
    return 0.5 * (x + y + softabs(x-y,k));
  });

SAM_SPECIALIZATION(double softmax2(double, double, double));
SAM_SPECIALIZATION(TMBad::ad_aug softmax2(TMBad::ad_aug, TMBad::ad_aug, TMBad::ad_aug));



/*
  For the zero mean AR(p) process
X_t = \sum_{i=1}^p  phi_i X_{t-i}
  The characteristic polynomial is 
  p(z) = 1 - phi_1 * z^1 - phi_2 * z^2 - \cdots - phi_p z^p
  The process is stationary if all roots of p(z) = 0 are outside the unit circle (|r_i| > 1)
  or equivalently, if all roots of the reciprocal polynomial
  q(w) = w^p - phi_1 * w^{p-1} - \cdots - phi_{p-1} * z - phi_p
  are inside the unit circle (|r_i| < 1).
 

*/
template<class Type>
vector<Type> logitroots2ARpar(vector<Type> x)SOURCE({
    vector<Type> xInc(x.size());
    xInc.setZero();
    xInc(0) = x(0);
    for(int i = 1; i < x.size(); ++i)
      xInc(i) = xInc(i-1) + exp(x(i));
  int n = x.size();
  vector<Type> p(n);
  p.setZero();
  vector<Type> coef(n+1);
  coef.setZero();
  coef(n) = 1.0;
  for(int i = 1; i <= n; ++i){
    for(int j = n - i - 1; j < n; ++j){
      if(j >= 0)
	coef(j) += Type(-1.0) * toInterval(xInc(i-1),Type(-1.0),Type(1.0),Type(1.0)) * coef(j+1);
    }
  }
  for(int i = 0; i < n; ++i)
    p(i) = -coef((n-1)-i);
  return p;
  })

  SAM_SPECIALIZATION(vector<double> logitroots2ARpar(vector<double>));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> logitroots2ARpar(vector<TMBad::ad_aug>));


template<class Type>
Type matrix_trace(matrix<Type>& x)SOURCE({
    if(x.cols() != x.rows())
      Rf_error("In matrix_trace: x should be square.");
    Type r = 0.0;
    for(int i = 0; i < x.cols(); ++i)
      r += x(i,i);
    return r;
  }
  )

  SAM_SPECIALIZATION(double matrix_trace(matrix<double>&));
SAM_SPECIALIZATION(TMBad::ad_aug matrix_trace(matrix<TMBad::ad_aug>&));


template<class Type>
Type logRiskHazard(Type x, Type m, Type loga, Type logb, int model)SOURCE({
  // Linear spline: a * (x - m) + (b - a) * pos(x-m)
  // Note that we force a<0 and b>0
  //  return -exp(loga) * (x - m) + exp(logspace_add_SAM(loga, logb)) * softmax2(x, Type(0.0), eps);
    Type eps = 0.00001;
    if(model == 0){ // Bathtub
      // Type v1 = -exp(loga) * (x-m) - loga;
      // Type v2 = exp(logb) * (x-m) - logb;
      // return -exp(logk) + logspace_add_SAM(v1, v2);
      Type x1 = -exp(loga) * (x - m) + exp(logspace_add_SAM(loga, logb)) * softmax2(x-m, Type(0.0), eps);
      Type dtp0 = 0.5* sqrt(exp(3*loga+logb)*eps - 2.*exp(2*loga+2*logb)*eps + exp(loga+3*logb)*eps)/exp(loga+logb);
      // if loga > logb, the dtp is dtp0; if loga < logb, dtp is -dtp0
      Type sgn = 2.0 / (1.0 + exp(-100.0 * (loga - logb))) - 1.0;
      Type tp = (dtp0) * sgn + m;
      Type d = -exp(loga) * (tp - m) + exp(logspace_add_SAM(loga, logb)) * softmax2(tp-m, Type(0.0), eps);
      Type r = x1 - d;
      return log(softmax2(r, Type(0.0), Type(0.00001)));
    }else if(model == 1){ // Sigmoid
      return loga - logspace_add_SAM(Type(0.0), -exp(logb) * (x - m));
    }else if(model == 2){ // Increase hockey
      return logb + log( softmax2(x-m, Type(0.0), eps) );      
    }else if(model == 3){ // Decrease hockey
      return log(-exp(loga) * (x - m) + exp(loga) * softmax2(x-m, Type(0.0), eps));
    }else if(model == 4){ // Exp Increase
      return loga + exp(logb) * x;
    }else if(model == 5){ // Exp Decrease
      return loga - exp(logb) * x;
    }else{
      Rf_error("Wrong hazard model!");
    }
  })


SAM_SPECIALIZATION(double logRiskHazard(double,double,double,double,int));
SAM_SPECIALIZATION(TMBad::ad_aug logRiskHazard(TMBad::ad_aug,TMBad::ad_aug,TMBad::ad_aug,TMBad::ad_aug,int));

