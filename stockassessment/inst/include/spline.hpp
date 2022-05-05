#pragma once
#ifndef SAM_SPLINE_HPP
#define SAM_SPLINE_HPP


namespace spline_atomic {

  // template<class Type>
  // Type pkwnorm(Type x, Type mu, Type sig, Type gam){
  //   Type a = 1.0 + 0.5 * (gam + sqrt(gam * gam));
  //   Type b = 1.0 - 0.5 * (gam - sqrt(gam * gam));
  //   Type lpv = pnorm5(x,mu,sig,Type(1.0), Type(1.0));
  //   Type lGa = a * lpv;
  //   // Type lr = logspace_sub(Type(0.0), (Type)(b * logspace_sub(Type(0.0),lGa)));
  //   return 1.0 - exp(b * logspace_sub2(Type(0.0),lGa));
  // }

  template<class Float>
  struct pkwnorm1_t {
    typedef Float Scalar; // Required by integrate
    Float gam;         // Parameters
    Float x0;
    // Evaluate joint density of (u, x)
    Float operator() (Float u) {
      Float a = 1.0 + 0.5 * (gam + sqrt(gam * gam));
      Float b = 1.0 - 0.5 * (gam - sqrt(gam * gam));
      Float lpv = pnorm_atomic::pnorm1_1x(u,Float(1.0), Float(1.0));
      Float lGa = a * lpv;
      return 1.0 - exp(b * rec_atomic::logspace_sub2_raw(Float(0.0),lGa));
    }
    // Integrate latent variable (u) out
    Float integrate(Float x) {
      using gauss_kronrod::integrate;
      Float ans = integrate(*this, x0, x);
      return ans;
    }
  };

   template<class Float>
   Float eval_ipkwnorm1(Float x, Float gam, Float x0) {
     pkwnorm1_t<Float> f = {gam, x0};
    return f.integrate(x);
  }

  TMB_BIND_ATOMIC(fun_ipkwnorm1, 110, eval_ipkwnorm1(x[0], x[1], x[2]))

 //  template<class Float>
 //  struct pkwnorm_t {
 //    typedef Float Scalar; // Required by integrate
 //    Float mu, sig, gam;         // Parameters
 //    // Evaluate joint density of (u, x)
 //    Float operator() (Float u) {
 //      Float a = 1.0 + 0.5 * (gam + sqrt(gam * gam));
 //      Float b = 1.0 - 0.5 * (gam - sqrt(gam * gam));
 //      Float lpv = pnorm_atomic::pnorm5_1(u,mu,sig,Float(1.0), Float(1.0));
 //      Float lGa = a * lpv;
 //      return 1.0 - exp(b * rec_atomic::logspace_sub2_raw(Float(0.0),lGa));
 //    }
 //    // Integrate latent variable (u) out
 //    Float integrate(Float x) {
 //      using gauss_kronrod::integrate;
 //      Float ans = integrate(*this, -INFINITY, x);
 //      return ans;
 //    }
 //  };
  
 // template<class Float>
 // Float eval_ipkwnorm(Float x, Float mu, Float sd, Float gam) {
 //   pkwnorm_t<Float> f = {mu, sd, gam};
 //    return f.integrate(x);
 //  }

 //  TMB_BIND_ATOMIC(fun_ipkwnorm, 1111, eval_ipkwnorm(x[0], x[1], x[2], x[3]))
  
}



namespace spline_helper {

 template<class Type>
 Type ipkwnorm(Type x, Type mu, Type sd, Type gam, double x0) {
    vector<Type> args(4); // Last index reserved for derivative order
    args << (x - mu) / sd, gam, (Type)(x0), 0;
    return spline_atomic::fun_ipkwnorm1(CppAD::vector<Type>(args))[0] * sd;
  }
  
  template<class Type>
  Type softmax(Type x, Type y, Type k = 1.0){
    return logspace_add2(k * x, k * y) / k;
  }

  
  // Kumaraswamy-normal (Kw-normal) density function with special choice of a and b
  template<class Type>
  Type dkwnorm(Type x, Type mu, Type sig, Type gam, bool give_log = false){
    Type a = 1.0 + 0.5 * (gam + sqrt(gam * gam));
    Type b = 1.0 - 0.5 * (gam - sqrt(gam * gam));
    Type lpv = pnorm5(x,mu,sig,Type(1.0), Type(1.0));
    Type lGa = a * lpv;
    Type log_res = log(a) + log(b) + dnorm(x,mu,sig, true) + (a-1.0) * lpv + (b-1.0) * logspace_sub2(Type(0.0),(Type)lGa);
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
    return 1.0 - exp(b * logspace_sub2(Type(0.0),lGa));
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
      if(fabs(knots(i) - knots(i-1)) < 0.001){
	res(i,0) = 0.001 / 3.0;
      }else{
	res(i,0) = (knots(i+1) - knots(i-1)) / 3.0;
      }
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
  // if(knots.size() + 1 != pars.size())
  //   Rf_error("Pars must have one more element than knots");
  matrix<Type> sg = spline_helper::getSigAndGam(knots);
  Type res = 0.0;
  for(int i = 0; i < knots.size(); ++i){
    Type tmp = spline_helper::dkwnorm(x, knots(i), sg(i,0), sg(i,1), false);
    res += pars(i) * tmp * sg(i,0);
  }
  Type t1 = 0.0; //2.5 * spline_helper::softmax(x - knots(knots.size() - 1),(Type)0,(Type)100.0);
  // Type t2 = 0.1 * spline_helper::softmax(knots(0) - x,(Type)0,(Type)100.0);
  return res - t1;// + t2 ;
}

// Integrated spline using Kw-normal density as basis functions
template<class Type>
Type ibcspline(Type x, vector<Type> knots, vector<Type> pars){
  // if(knots.size() + 1 != pars.size())
  //   Rf_error("Pars must have one more element than knots");
  // Type x0 = CppAD::CondExpLt(x, knots(0), knots(0),
  // 			     CppAD::CondExpGt(x, knots(knots.size()-1),
  // 					      knots(knots.size()-1),
  // 					      x)
  // 			     );
  matrix<Type> sg = spline_helper::getSigAndGam(knots);
  Type res = 0.0;
  for(int i = 0; i < knots.size(); ++i){
    // Should be zero at NegInf - not left endpoint of knot interval
    Type v0 = spline_helper::pkwnorm(knots(0), knots(i), sg(i,0), sg(i,1));
    Type tmp = spline_helper::pkwnorm(x, knots(i), sg(i,0), sg(i,1));
    res += pars(i) * (tmp - v0);
  }
  // return res - 1.5 * spline_helper::softmax(tmp * (x - knots(knots.size() - 1)),(Type)0,(Type)100.0);
  Type t1 = 0.0; //2.5 * spline_helper::softmax(x - knots(knots.size() - 1),(Type)0,(Type)100.0);
  //Type t2 = 0.1 * spline_helper::softmax(knots(0) - x,(Type)0,(Type)100.0);
  return res - t1; // + t2 ;
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


// Integrated integrated spline using Kw-normal density as basis functions
template<class Type>
Type iibcspline(Type x, vector<Type> knots, vector<Type> pars){
  // Type x0 = CppAD::CondExpLt(x, knots(0), knots(0),
  // 			     CppAD::CondExpGt(x, knots(knots.size()-1),
  // 					      knots(knots.size()-1),
  // 					      x)
  // 			     );
  matrix<Type> sg = spline_helper::getSigAndGam(knots);
  Type res = 0.0;
  for(int i = 0; i < knots.size(); ++i){
    // Should be zero at NegInf - not left endpoint of knot interval
    Type v0 = spline_helper::pkwnorm(knots(0), knots(i), sg(i,0), sg(i,1));
    Type tmp = spline_helper::ipkwnorm(x, knots(i), sg(i,0), sg(i,1), 0.0);
    res += pars(i) * (tmp - v0 * x);
  }
  // return res - 1.5 * spline_helper::softmax(tmp * (x - knots(knots.size() - 1)),(Type)0,(Type)100.0);
  //Type t1 = 2.5 * spline_helper::softmax(0.5*x*x - knots(knots.size() - 1),(Type)0,(Type)100.0);
  //Type t2 = 0.1 * spline_helper::softmax(knots(0) - x,(Type)0,(Type)100.0);
  return res;// - t1; // + t2 ;
}



// Monotonically non-decreasing negative spline using Kw-nomal as basis functions
// Has an extra parameter to allow negative value at left bound
template<class Type>
Type iibcispline(Type x, vector<Type> knots, vector<Type> pars){
  if(knots.size() + 1 != pars.size())
    Rf_error("Pars must have one more element than knots");
  vector<Type> p2(pars.size() - 1);
  p2.setZero();
  Type p2s = 0.0;
  for(int i = 0; i < p2.size(); ++i){
    p2(i) = exp(pars(i));
    p2s += p2(i);
  }
  Type r = iibcspline(x, knots, p2);
  return r - p2s * x + pars(pars.size()-1);
}



#endif
