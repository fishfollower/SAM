#ifndef SAM_REPRODUCTIVE_HPP
#define SAM_REPRODUCTIVE_HPP


template<class Type>
struct MDSR_NEWT {
  // Negative log gradient of stock-recruitment function. To be implemented for specializations.
  virtual Type fn(Type logs){
    return (Type)R_NaReal;
  }
  Type gr(Type logs){
    Type h = 0.0001 * softmax(fabs(logs),Type(0.0001), Type(100.0));
    Type g = (-fn(logs + 2.0 * h) + 8.0 * fn(logs + h) - 8.0 * fn(logs-h) + fn(logs-2.0*h)) / (12.0 * h);
    return g;
  }
  Type he(Type logs){
    Type h = 0.0001 * softmax(fabs(logs),Type(0.0001), Type(100.0));
    Type gg = (-fn(logs + 2.0 * h) + 16.0 * fn(logs + h) - 30.0 * fn(logs) + 16.0 * fn(logs-h) - fn(logs-2.0*h)) / (12.0 * h * h);
    return gg;
  }
  Type sign0(Type x){
    return x / (fabs(x) + 1e-8);
  }
  Type softmax(Type x, Type y, Type k = 1.0){
    return logspace_add2(k * x, k * y) / k;
  }
  Type numnewt(Type logs){
    // Minimizing
    // Type fnv = fn(logs);
    Type grv = gr(logs);
    Type hev = he(logs);
    // Type g = (f(logs + h) - f(logs - h)) / (2.0 * h);
    Type s = sign0(grv) * sign0(hev);
    Type y = log(fabs(grv)) - log(softmax(fabs(hev), (Type)0.001, (Type)1000.0)); // Damp the gradient
    return logs - 0.9 * s * exp(y); //softmax(exp(y), 0.5 * fabs(logs), Type(1000.0));
  }

  Type minimize(Type logs0){
    Type sv = logs0;
    for(int i = 0; i < 100; ++i){
      sv = numnewt(sv);
    }
    return sv;
  }
  
};

//
  
template<class Type>
struct MDSR_50 : MDSR_NEWT<Type> {
  Type loga;
  Type logb;
  Type logc;
  Type logd;

  MDSR_50(Type la, Type lb, Type lc, Type ld) : MDSR_NEWT<Type>(), loga(la), logb(lb), logc(lc), logd(ld) {};
  
  Type fn (Type logs) override{
    // ln(a) - 2*ln(S + d) - ln(1 + exp((S - b)/(b*g))) + ln(1 + exp(-1/g)) + ln(-b*d*g*(1 + exp((S - b)/(b*g)))*ln(1 + exp((S - b)/(b*g))) + d*(b*g*ln(1 + exp(-1/g)) + S)*exp((S - b)/(b*g)) + ln(1 + exp(-1/g))*b*d*g + S*(S + 2*d))
    Type v0 = (exp(logs)-exp(logb)) * exp(- (logb + logc));
    Type v1 = logspace_add2(Type(0.0), v0);
    Type v2 = logspace_add2(Type(0.0), -exp(-logc));
    Type xx1 = logb + logd + logc + v1 + log(v1);
    Type xx2 = logd + logspace_add2(logb + logc + log(v2), logs) + v0;
    Type xx3 = log(v2) + logb + logd + logc;
    Type xx4 = 2.0 * logs;
    Type xx5 = log(2.0) + logd + logs;
    Type v3 = logspace_add2(logspace_add2(logspace_add2(logspace_add2(xx1,xx2),xx3),xx4),xx5);
    return -(loga - 2.0 * logspace_add2(logs, logd) - v1 + v2 + v3);
  }
};

template<class Type>
Type mdsr_50(Type loga, Type logb, Type logc, Type logd){
  MDSR_50<Type> f(loga, logb, logc, logd);
  Type sv = f.minimize(logd);
  return exp(-f.fn(sv));
}


// 

template<class Type>
struct MDSR_51 : MDSR_NEWT<Type> {
  Type loga;
  Type logb;
  Type logd;

  MDSR_51(Type la, Type lb, Type ld) : MDSR_NEWT<Type>(), loga(la), logb(lb), logd(ld) {};
  
  Type fn (Type logs) override{
    // ln(a) + ln(S) - b*S - 2*ln(S + d) + ln(-S^2*b - S*b*d + S + 2*d)
    Type v1 = logspace_add2(logs,logd);
    Type xA = (2.0*logs+logb);
    Type xB = (logs+logb+logd);
    Type xx = logspace_add2(xA,xB);
    Type v2 = logspace_add2(logs,log(2.0)+logd);
    Type v3 = logspace_sub(v2,xx);
    return -(loga + logs - exp(logb + logs) - 2.0 * v1 + v3);
  }
};

template<class Type>
Type mdsr_51(Type loga, Type logb, Type logd){
  MDSR_51<Type> f(loga, logb, logd);
  Type sv = f.minimize(logd);
  return exp(-f.fn(sv));
}

// 

template<class Type>
struct MDSR_52 : MDSR_NEWT<Type> {
  Type loga;
  Type logb;
  Type logd;

  MDSR_52(Type la, Type lb, Type ld) : MDSR_NEWT<Type>(), loga(la), logb(lb), logd(ld) {};
  
  Type fn (Type logs) override{
    Type v1 = logspace_add2(logs+logb+logd,logs);
    Type v2 = logspace_add2(v1, (Type)log(2.0) + logd);
    Type v3 = logspace_add2(logs+logb,Type(0.0));
    Type v4 = logspace_add2((Type)logs,(Type)logd);
    // a*S*(S*b*d + S + 2*d)/((S*b + 1)^2*(S + d)^2);
    return -(loga + logs + v2 - 2.0 * v3 - 2.0 * v4);
      //loga + logs + logspace_add2(logspace_add2(logs+logb+logd,logs),log((Type)2.0) + logd) - 2.0 * logspace_add2(logs+logb,Type(0.0)) - Type(2.0) * (Type)logspace_add2((Type)logs,(Type)logd);
  }
};

template<class Type>
Type mdsr_52(Type loga, Type logb, Type logd){
  MDSR_52<Type> f(loga, logb, logd);
  Type sv = f.minimize(logd);
  return exp(-f.fn(sv));
}

// 

template<class Type>
struct MDSR_67 : MDSR_NEWT<Type> {
  Type loga;
  Type logb;
  Type logd;

  MDSR_67(Type la, Type lb, Type ld) : MDSR_NEWT<Type>(), loga(la), logb(lb), logd(ld) {};
  
  Type fn (Type logs) override{
    Type v1 = logspace_add2(logs + logb + logd, Type(0.0));
    Type v2 = logspace_sub2(-exp(logd) * v1, -(exp(logd) + 1.0) * v1 + logs + logb + 2.0 * logd);
    // a*S*(S*b*d + S + 2*d)/((S*b + 1)^2*(S + d)^2);
    return -(loga + v2);
      //loga + logs + logspace_add2(logspace_add2(logs+logb+logd,logs),log((Type)2.0) + logd) - 2.0 * logspace_add2(logs+logb,Type(0.0)) - Type(2.0) * (Type)logspace_add2((Type)logs,(Type)logd);
  }
};

template<class Type>
Type mdsr_67(Type loga, Type logb, Type logd){
  MDSR_67<Type> f(loga, logb, logd);
  Type sv = f.minimize(logd);
  return exp(-f.fn(sv));
}



template<class Type>
Type maxDiffSR(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, int y){
  Type mdsr = R_PosInf;
  Type xx = 0.0;
  switch(conf.stockRecruitmentModelCode){
    case 0: // straight RW
       // Not well defined
    break;
    case 1: //ricker
      // Max derivative is alpha
      mdsr = exp(par.rec_pars(0));
    break;
    case 2:  //BH
      // Max derivative is alpha
      mdsr = exp(par.rec_pars(0));
    break;
  case 3: //Constant mean
    // Derivative is Inf: lim_{h->0} (f(0+h)-f(0))/h where f(0)=0, f(h)=some positive value for any h>0
    mdsr = R_PosInf;
    break;
  case 50:
    mdsr = mdsr_50(par.rec_pars(0),par.rec_pars(1),par.rec_pars(2), par.rec_pars(3));
    break;
  case 51:
    mdsr = mdsr_51(par.rec_pars(0),par.rec_pars(1),par.rec_pars(2));
    break;
  case 52:
    mdsr = mdsr_52(par.rec_pars(0),par.rec_pars(1),par.rec_pars(2));
    break;
  case 60:
    mdsr = exp(par.rec_pars(0));
    break;
  case 61: // Hockey stick
    // Max derivative is slope to breakpoint
    // level / blim
    // Type log_level = rec_pars(0);
    // Type log_blim = rec_pars(1);
    mdsr = exp(par.rec_pars(0) - par.rec_pars(1));
    break;
  case 62: // AR1 (on log-scale)
    // Derivative is Inf: lim_{h->0} (f(0+h)-f(0))/h where f(0)=0, f(h)=some positive value for any h>0
    mdsr = R_PosInf;
   break;
  case 63: //Bent hyperbola / Hockey-stick-like
    /*
      Source: e.g., DOI:10.1093/icesjms/fsq055
      rec_pars(0): log-Blim
      rec_pars(1): log of half the slope from 0 to Blim
      rec_pars(2): log-Smoothness parameter
     */
    mdsr = 2.0 * exp(par.rec_pars(1));
    break;
  case 64: // Power CMP
    // Derivative goes to infinity
    mdsr = R_PosInf;
    break;
  case 65: // Power Non-CMP
    // Maximum reproduction rate goes to infinity for large SSB
    mdsr = R_PosInf;
    break;
  case 66: // Shepherd
    // Optimum for S = exp(log((1+d)/(d-1))/d)
    // Maximum derivative d<=1: a, d>1:
    xx = CppAD::CondExpGt(par.rec_pars(2),Type(0.0),exp(log((1.0+exp(par.rec_pars(2))) / (exp(par.rec_pars(2)) - 1.0)) / exp(par.rec_pars(2))), Type(1e-10));
    mdsr = -(-1.0 + (exp(par.rec_pars(2)) - 1.0)*pow(xx/exp(par.rec_pars(1)),exp(par.rec_pars(2))))*exp(par.rec_pars(0))/pow(pow(xx/exp(par.rec_pars(1)),exp(par.rec_pars(2))) + 1.0,2.0);
    break;
  case 67: // Deriso
    mdsr = mdsr_67(par.rec_pars(0),par.rec_pars(1),par.rec_pars(2));
    break;
  case 68: // Saila-Lorda (Iles 1994)
    
    break;
  case 69: // Sigmoidal Beverton-Holt
    
    break;
  case 90: // Non-increasing spline on log(R/S)
    
    break;
  case 91: // Integrated spline on log(R/S)
    
    break;
  case 92: // Spline on log(R/S)
    
    break;
  default:    
      Rf_error("SR model code not recognized");
    break;
  }
  return mdsr;
}



template<class Type>
Eigen::SparseMatrix<Type> Leslie_i(dataSet<Type> &dat, confSet &conf, paraSet<Type>& par, int y){
  
  Type max_dSR = maxDiffSR(dat,conf,par,y); // For ricker and BH
  int nAges = conf.maxAge - conf.minAge + 1;

  typedef Eigen::Triplet<Type> T;
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
  Eigen::SparseMatrix<Type> m(nAges,nAges);
  m.setFromTriplets(tripletList.begin(), tripletList.end());
  return m;
}

template<class Type>
Type norm2(vector<Type> x){
  Type r = 0.0;
  for(int i = 0; i < x.size(); ++i)
    r += x(i) * x(i);
  return sqrt(r);
}

template<class Type>
Type maxEigenValue(Eigen::SparseMatrix<Type> A){
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



template<class Type>
Type rmax_i(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, int y){
  // For Ricker, BH and any CMP model, it is the derivative at 0, but is the largest derivative the correct generalization? What should it be for constant recruitment models? 
  // Implement for other models
  if(conf.stockRecruitmentModelCode == 0 ||
     conf.stockRecruitmentModelCode == 3 ||
     conf.stockRecruitmentModelCode == 62 ||
     conf.stockRecruitmentModelCode == 64 ||
     conf.stockRecruitmentModelCode == 65){
    return R_NaReal;
  }
  Eigen::SparseMatrix<Type> L = Leslie_i(dat,conf,par,y);
  Type lambda = maxEigenValue(L);
  return log(lambda);
}

template<class Type>
vector<Type> rmax(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par){
  int timeSteps = dat.natMor.dim(0);
  vector<Type> v(timeSteps);
  v.setZero();
  for(int y = 0; y < timeSteps; ++y)
    v(y) = rmax_i(dat, conf, par, y);
  return v;
}



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


template<class Type>
vector<Type> generationLength(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par){
  int timeSteps = dat.natMor.dim(0);
  vector<Type> v(timeSteps);
  v.setZero();
  for(int y = 0; y < timeSteps; ++y)
    v(y) = generationLength_i(dat, conf, par, y);
  return v;
}




#endif
