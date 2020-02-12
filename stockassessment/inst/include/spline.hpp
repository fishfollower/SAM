namespace spline_helper {
  
  template<class Type>
  Type softmax(Type x, Type y, Type k = 1.0){
    return logspace_add(k * x, k * y) / k;
  }

  
  // Kumaraswamy-normal (Kw-normal) density function with special choice of a and b
  template<class Type>
  Type dkwnorm(Type x, Type mu, Type sig, Type gam, bool give_log = false){
    Type a = 1.0 + 0.5 * (gam + sqrt(gam * gam));
    Type b = 1.0 - 0.5 * (gam - sqrt(gam * gam));
    Type lpv = pnorm5(x,mu,sig,Type(1.0), Type(1.0));
    Type lGa = a * lpv;
    Type log_res = log(a) + log(b) + dnorm(x,mu,sig, true) + (a-1.0) * lpv + (b-1.0) * logspace_sub(Type(0.0),(Type)lGa);
     if(give_log)
      return log_res;
    return exp(log_res);
  }

  // Kumaraswamy-normal (Kw-normal) distribution function with special choice of a and b
  template<class Type>
  Type pkwnorm(Type x, Type mu, Type sig, Type gam){
    Type a = 1.0 + 0.5 * (gam + sqrt(gam * gam));
    Type b = 1.0 - 0.5 * (gam - sqrt(gam * gam));
    Type lpv = pnorm5(x,mu,sig,Type(1.0), Type(1.0));
    Type lGa = a * lpv;
    // Type lr = logspace_sub(Type(0.0), (Type)(b * logspace_sub(Type(0.0),lGa)));
    return 1.0 - exp(b * logspace_sub(Type(0.0),lGa));
  }
  
  template<class Type>
  matrix<Type> getSigAndGam(vector<Type> knots){
    // if(CppAD::Variable(knots(0)))
    //   Rf_error("Knots can not be parameters");
    if(knots.size() < 3)
      Rf_error("The spline must have at least three knots.");
    matrix<Type> res(knots.size(),2); // Sigma, Gamma
    res.setZero();
    for(int i = 1; i < knots.size() - 1; ++i){
      // Sigma
      res(i,0) = (knots(i+1) - knots(i-1)) / 3.0;
      // Gamma
      if(fabs(knots(i) - knots(i-1)) < 1e-8){
	res(i,1) = res(i,0);
      }else if(fabs(knots(i) - knots(i+1)) < 1e-8){
	res(i,1) = -res(i,0);
      }else{
	res(i,1) = log(knots(i+1) - knots(i)) - log(knots(i) - knots(i-1));
      }
    }
    res(0,0) = (knots(2) - knots(0)) / 3.0;
    res(knots.size()-1,0) = (knots(knots.size()-1) - knots(knots.size()-1-2)) / 3.0;
    if(fabs(knots(1) - knots(0)) < 1e-8){
      res(0,1) = 0.0;
    }else{
      res(0,1) = res(0,0);
    }
    if(fabs(knots(knots.size()-1) - knots(knots.size()-1-1)) < 1e-8){
      res(knots.size()-1,1) = 0.0;
    }else{
      res(knots.size()-1,1) = -res(knots.size()-1,0);
    }
    return res;
  }
}
  
// Spline using Kw-normal density as basis functions
template<class Type>
Type bcspline(Type x, vector<Type> knots, vector<Type> pars){
  // Type x0 = CppAD::CondExpLt(x, knots(0), knots(0),
  // 			     CppAD::CondExpGt(x, knots(knots.size()-1),
  // 					      knots(knots.size()-1),
  // 					      x)
  // 			     );
  if(knots.size() != pars.size())
    Rf_error("Knots and pars must have same length");
  matrix<Type> sg = spline_helper::getSigAndGam(knots);
  Type res = 0.0;
  for(int i = 0; i < knots.size(); ++i){
    Type tmp = spline_helper::dkwnorm(x, knots(i), sg(i,0), sg(i,1), false);
    res += pars(i) * tmp * sg(i,0);
  }
  return res - spline_helper::softmax(x - knots(knots.size() - 1),(Type)0,(Type)100.0);
}

// Integrated spline using Kw-normal density as basis functions
template<class Type>
Type ibcspline(Type x, vector<Type> knots, vector<Type> pars){
  if(knots.size() != pars.size())
    Rf_error("Knots and pars must have same length");
  // Type x0 = CppAD::CondExpLt(x, knots(0), knots(0),
  // 			     CppAD::CondExpGt(x, knots(knots.size()-1),
  // 					      knots(knots.size()-1),
  // 					      x)
  // 			     );
  matrix<Type> sg = spline_helper::getSigAndGam(knots);
  Type res = 0.0;
  Type tmp;
  for(int i = 0; i < knots.size(); ++i){
    // Should be zero at NegInf - not left endpoint of knot interval
    Type v0 = 0.0; //spline_helper::pkwnorm(knots(0), knots(i), sg(i,0), sg(i,1));
    tmp = spline_helper::pkwnorm(x, knots(i), sg(i,0), sg(i,1));
    res += pars(i) * (tmp - v0);
  }
  return res - spline_helper::softmax(tmp * (x - knots(knots.size() - 1)),(Type)0,(Type)100.0);
}

// Monotonically non-increasing spline using Kw-nomal as basis functions
// Has an extra parameter to allow positive value at left bound
template<class Type>
Type ibcdspline(Type x, vector<Type> knots, vector<Type> pars){
  if(knots.size() + 1 != pars.size())
    Rf_error("Pars must have one more element than knots");
  vector<Type> p2(pars.size() - 1);
  p2.setZero();
  for(int i = 0; i < p2.size(); ++i)
    p2(i) = -exp(pars(i));
  Type r = ibcspline(x, knots, p2);
  return r + pars(pars.size()-1);
}

// Monotonically non-decreasing spline using Kw-nomal as basis functions
// Has an extra parameter to allow negative value at left bound
template<class Type>
Type ibcispline(Type x, vector<Type> knots, vector<Type> pars){
  if(knots.size() + 1 != pars.size())
    Rf_error("Pars must have one more element than knots");
  vector<Type> p2(pars.size() - 1);
  p2.setZero();
  for(int i = 0; i < p2.size(); ++i)
    p2(i) = exp(pars(i));
  Type r = ibcspline(x, knots, p2);
  return r + pars(pars.size()-1);
}


extern "C" {
  SEXP bcsplineR(SEXP x, SEXP knots, SEXP pars){
    double r = bcspline(Rf_asReal(x),
			asVector<double>(knots),
			asVector<double>(pars));
    return asSEXP(r);
  }
  SEXP ibcsplineR(SEXP x, SEXP knots, SEXP pars){
    double r = ibcspline(Rf_asReal(x),
			 asVector<double>(knots),
			 asVector<double>(pars));
    return asSEXP(r);
  }
  SEXP ibcdsplineR(SEXP x, SEXP knots, SEXP pars){
    double r = ibcdspline(Rf_asReal(x),
			  asVector<double>(knots),
			  asVector<double>(pars));
    return asSEXP(r);
  }
  SEXP ibcisplineR(SEXP x, SEXP knots, SEXP pars){
    double r = ibcdspline(Rf_asReal(x),
			  asVector<double>(knots),
			  asVector<double>(pars));
    return asSEXP(r);
  }
}
