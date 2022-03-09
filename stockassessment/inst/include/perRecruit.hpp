#ifndef SAM_PERRECRUIT_HPP
#define SAM_PERRECRUIT_HPP

#ifndef TMBAD_FRAMEWORK
#ifndef CPPAD_FRAMEWORK
#define CPPAD_FRAMEWORK
#endif
#endif
  // template<class Float>
  // Float GrNum_sbh_raw(Float s, Float l, Float a, Float b, Float g, Float h){
  //   Float err = 1e5;
  //   Float ans = 0.0;
  //   Float tol = 1e-12;
  //   Float tab[SAM_RECNTAB][SAM_RECNTAB];    
  //   tab[0][0] = (Fn_sbh_raw((Float)(s+h),(Float)l,(Float)a,(Float)b,(Float)g) - Fn_sbh_raw((Float)(s-h),(Float)l,(Float)a,(Float)b,(Float)g)) / (2.0 * h);
  //   for(int i = 1; i < SAM_RECNTAB; ++i){
  //       tab[i][0] = (Fn_sbh_raw((Float)(s+h),(Float)l,(Float)a,(Float)b,(Float)g) - Fn_sbh_raw((Float)(s-h),(Float)l,(Float)a,(Float)b,(Float)g)) / (2.0 * h);
  // 	Float f = 1.0;
  // 	for(int m = 1; m <= i; ++m){
  // 	  f *= 4.0;
  // 	  tab[i][m] = (tab[i][m-1] * f - tab[i-1][m-1]) / (f - 1.0);
  // 	  Float tmp1 = fabs(tab[i][m] - tab[i][m-1]);
  // 	  Float tmp2 = fabs(tab[i][m] - tab[i-1][i-1]);
  // 	  Float errtmp = 0.5 * (tmp1 + tmp2 + fabs(tmp1-tmp2));
  // 	  if(errtmp < err){
  // 	   ans = tab[i][m];
  // 	    err = errtmp;
  // 	    if(err < tol)
  // 	      goto endloop;
  // 	  }
  // 	}
  // 	if(fabs(tab[i][i] - tab[i-1][i-1]) > 2.0 * err)
  // 	  break;
  //     }
  //   endloop:
  //   return ans;  
  // }

#ifdef TMBAD_FRAMEWORK
#define SAM_NegInf -20.0
#define SAM_NIZero -10.0
#else
#define SAM_NegInf -20.0
#define SAM_NIZero -10.0
#endif

#define SAM_Zero exp(SAM_NIZero)
// exp(SAM_NegInf)


template<class Type>
Type logspace_sum(vector<Type> x){
  Type r = R_NegInf;
  for(int i = 0; i < x.size(); ++i)
    r = logspace_add2(r, x(i));
  return r;
}


namespace rec_atomic {
  /*
   * Namespace with calculations for finding equilibrium biomass numerically
   */


  
  /*
   * Spline functions re-implemented for tiny_ad
   */
  
  // Kumaraswamy-normal (Kw-normal) density function with special choice of a and b
  template<class Float>
  Float dkwnorm_raw(Float x, Float mu, Float sig, Float gam, bool give_log = false){
    Float a = 1.0 + 0.5 * (gam + sqrt(gam * gam));
    Float b = 1.0 - 0.5 * (gam - sqrt(gam * gam));
    Float lpv = pnorm_atomic::pnorm5_raw(x,mu,sig,Float(1.0), Float(1.0));
    Float lGa = a * lpv;
    Float log_res = log(a) + log(b) + dnorm(x,mu,sig, true) + (a-1.0) * lpv + (b-1.0) * rec_atomic::logspace_sub2(Float(0.0),(Float)lGa);
    if(give_log)
      return log_res;
    return exp(log_res);
  }

  // Kumaraswamy-normal (Kw-normal) distribution function with special choice of a and b
  template<class Float>
  Float pkwnorm_raw(Float x, Float mu, Float sig, Float gam){
    Float a = 1.0 + 0.5 * (gam + sqrt(gam * gam));
    Float b = 1.0 - 0.5 * (gam - sqrt(gam * gam));
    Float lpv = pnorm_atomic::pnorm5_raw(x,mu,sig,Float(1.0), Float(1.0));
    Float lGa = a * lpv;
    // Float lr = atomic::robust_utils::logspace_sub(Float(0.0), (Float)(b * atomic::robust_utils::logspace_sub(Float(0.0),lGa)));
    return 1.0 - exp(b * rec_atomic::logspace_sub2(Float(0.0),lGa));
  }


  template<class Float>
  matrix<Float> getSigAndGam_raw(vector<Float> knots){
    // if(CppAD::Variable(knots(0)))
    //   Rf_error("Knots can not be parameters");
    if(knots.size() < 3)
      Rf_error("The spline must have at least three knots.");
    matrix<Float> res(knots.size(),2); // Sigma, Gamma
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

  
  // Spline using Kw-normal density as basis functions
  template<class Float>
  Float bcspline_raw(Float x, vector<Float> knots, vector<Float> pars, matrix<Float> sg){
    if(knots.size() != pars.size())
      Rf_error("Knots and pars must have same length");
    // matrix<Float> sg = getSigAndGam(knots);
    Float res = 0.0;
    for(int i = 0; i < knots.size(); ++i){
      Float tmp = dkwnorm_raw(x, knots(i), sg(i,0), sg(i,1), false);
      res += pars(i) * tmp * sg(i,0);
    }
    return res;
  }

  // Integrated spline using Kw-normal density as basis functions
  template<class Float>
  Float ibcspline_raw(Float x, vector<Float> knots, vector<Float> pars, matrix<Float> sg){
    if(knots.size() != pars.size())
      Rf_error("Knots and pars must have same length");
    // matrix<Float> sg = getSigAndGam_raw(knots);
    Float res = 0.0;
    for(int i = 0; i < knots.size(); ++i){
      // Should be zero at NegInf - not left endpoint of knot interval
      Float v0 = 0.0;
      // Float v0 = pkwnorm_raw(knots(0), knots(i), sg(i,0), sg(i,1));
      Float tmp = pkwnorm_raw(x, knots(i), sg(i,0), sg(i,1));
      res += pars(i) * (tmp - v0);
    }
    return res ;
  }

  // Monotonically non-increasing spline using Kw-nomal as basis functions
  // Has an extra parameter to allow positive value at left bound
  template<class Float>
  Float ibcdspline_raw(Float x, vector<Float> knots, vector<Float> pars, matrix<Float> sg){
    if(knots.size() + 1 != pars.size())
      Rf_error("Pars must have one more element than knots");
    vector<Float> p2(pars.size() - 1);
    p2.setZero();
    for(int i = 0; i < p2.size(); ++i)
      p2(i) = -exp(pars(i));
    Float r = ibcspline_raw(x, knots, p2, sg);
    return r + pars(pars.size()-1);
  }

  // Monotonically non-decreasing spline using Kw-nomal as basis functions
  // Has an extra parameter to allow negative value at left bound
  template<class Float>
  Float ibcispline_raw(Float x, vector<Float> knots, vector<Float> pars, matrix<Float> sg){
    if(knots.size() + 1 != pars.size())
      Rf_error("Pars must have one more element than knots");
    vector<Float> p2(pars.size() - 1);
    p2.setZero();
    for(int i = 0; i < p2.size(); ++i)
      p2(i) = exp(pars(i));
    Float r = ibcspline_raw(x, knots, p2, sg);
    return r + pars(pars.size()-1);
  }

  
  /*
   * End of spline functions
   */

  // DOI: 10.1145/361952.361970
  template<class Float>
  Float lambertW_raw(Float x){
    Float wn;
    if(x > 0){
      Float zn, en;
      if(x < 6.46){
	wn = (x + 4.0 / 3.0 * x * x) / (1.0 + 7.0/3.0 * x + 5.0 / 6.0 * x*x);
      }else{
	wn = log(x);
      }

      do{
	zn = log(x) - log(wn) - wn;
	Float tmp1 = 1.0 + wn;
	en = (zn * (2.0 * (tmp1) * (tmp1 + 2.0/3.0 * zn) - zn)) / ((tmp1) * ( 2.0 * (tmp1) * (tmp1 + 2.0/3.0 * zn) - 2.0 * zn));
	wn = wn * (1.0 + en);
      }while(en > 1e-8);   
    }else{
      Rf_warning("lambertW is only implemented for x>0");
      wn = SAM_Zero;
    }
    return wn;
  }

  TMB_BIND_ATOMIC(lambertW0, 1, lambertW_raw(x[0]))

  /*
   * General functions for finding equilibrium numerically
   */
  
  template<class Float, class Functor>
  struct FUN_FOR_NEWT {
    Functor f;
    
    Float logSR(Float logs){
      return f(logs);
    }
    
    Float gr_sr(Float s){
      typedef atomic::tiny_ad::variable<1, 1, Float> F2;
      F2 tmp(s,0);
      F2 ls = log(tmp);
      F2 r0 = f(ls);
      F2 res = exp(r0);
      return res.getDeriv()[0];
    }

    template<class F>
    F Fn(F logs, F l){
      F tmp = log(l) + f(logs) - logs;
      return tmp * tmp;
    }
    Float Gr(Float logs, Float l){
      typedef atomic::tiny_ad::variable<1, 1, Float> F2;
      F2 tmp(logs,0);
      F2 res = Fn(tmp, (F2)l);
      return res.getDeriv()[0];
    }
    Float He(Float logs, Float l){
      typedef atomic::tiny_ad::variable<2, 1, Float> F2;
      F2 tmp(logs,0);
      F2 res = Fn(tmp, (F2)l);
      return res.getDeriv()[0];
    }
  };


  template<class Float>
  struct NEWTON_RESULT {
    Float objective;
    int niter;
    Float par;
    Float gr;
    Float he;
    Float determinant;
    int posdef;
    Float mgc;
    Float logSR;
    Float grSR;
  };
  
  template<class Float, class Functor>
  NEWTON_RESULT<Float> newton(Functor f0, Float l, Float logs0){
    FUN_FOR_NEWT<Float, Functor> f;
    f.f = f0;
    int maxit = 1000;
    Float grad_tol = 1.0e-5;
    Float c1 = 1.0e-4;
    Float c2 = 0.1;
    Float logs = logs0;
    // Float logsOld = logs;
    int it = 0;
    Float mgc = R_PosInf;
    Float gk = f.Gr(logs, l);
    Float Gk = f.He(logs, l);
    Float fCurrent = f.Fn(logs, l);
    bool posdef = Gk > 1e-7;
    Float m = 1.0 / Gk;
    vector<Float> hybridGuess(4);
    hybridGuess << 1.0, 0.5, 0.25, 0.0625;

    while( it < maxit && mgc > grad_tol ){
      // logsOld = logs;
      int useHybridGuessNum = -1;
      Float tk = 0.9;

      Float dk = -gk;

      if(posdef)
      	dk = -m * gk;

      Float gkdk = gk * dk;
     
      if(hybridGuess.size() > 0){	
	for(int i = 0; i < hybridGuess.size(); ++i){
	  Float logs2 = logs + dk * hybridGuess(i);
	  Float f1 = f.Fn(logs2, l);
	  Float g1dk = f.Gr(logs2, l);
	  Float t1 = c1 * hybridGuess(i) * gkdk;
	  Float t2 = c2 * gkdk;	  
	  if((f1 < fCurrent + t1) &&
	     (g1dk > t2)){
	    useHybridGuessNum = i;
	    break;
	  }
	}
      }

      if(useHybridGuessNum == -1){
	// Use step size from Alg.1 in https://arxiv.org/pdf/1612.06965.pdf
	Float rhok = -gkdk;
        Float dkGkdk = dk * Gk * dk;
	Float deltak = sqrt( dkGkdk );
	Float den = rhok + deltak;
	den *= deltak;
        tk = rhok;	
	tk /= (den + 1e-10);
      }else{
	tk = hybridGuess(useHybridGuessNum);
      }
      logs += tk * dk;
      gk = f.Gr(logs, l);
      Gk = f.He(logs, l);
      fCurrent = f.Fn(logs, l);
      m = 1.0 / Gk;
      posdef = Gk > 1e-7;
      mgc = fabs(gk);
      ++it;
    }
    NEWTON_RESULT<Float> res = {fCurrent, it, logs, gk, Gk, Gk, posdef, mgc, f.logSR(logs), f.gr_sr(exp(logs))};

    return res;
  }

  /*
   * Specialization for the Sigmoidal Beverton-Holt recruitment model
   */
  
  template<class Float>
  struct SBH {
    Float a;
    Float b;
    Float g;
    // Float ls0;
    
    template<class F>
    F operator()(F logs){
      // F fls0 = ls0;
      // F ls2 = atomic::robust_utils::logspace_add(logs,fls0);
      return log((F)a) + (F)g * logs - log(1.0+(F)b*exp((F)g * logs));
    }

  };
  
  template<class Float>
  Float Se_sbh_raw(Float l, Float a, Float b, Float g){
    SBH<Float> f;
    f.a = a; f.b = b; f.g = g;// f.ls0 = -10.0;
    NEWTON_RESULT<Float> r = newton<Float, SBH<Float> >(f, l, log(a) - log(b) + log(2.0));
    if(r.par > log(SAM_Zero) &&
       r.objective < 1.0e-6 &&
       fabs(r.grSR) < 1){ // Found stable equilibrium
      return exp(r.par);
    }else if(r.par > log(SAM_Zero) &&
	     r.objective < 1.0e-6 &&
	     fabs(r.grSR) >= 1){ // Found unstable equilibrium
      Rf_warning("Found unstable equilibrium");
      return exp(r.par);
    }else{
      return SAM_Zero * l;
    }
    return SAM_Zero * l;
  }

  TMB_BIND_ATOMIC(Se_sbh0, 1111, Se_sbh_raw(x[0],x[1],x[2],x[3]))


  /*
   * Specialization for the Saila-Lorda recruitment model
   */
  
 
  template<class Float>
  struct SL {
    Float a;
    Float b;
    Float g;
    
    template<class F>
    F operator()(F logs){
      return log((F)a) + (F)g * logs - (F)b * exp(logs);
    }
  };
  
  template<class Float>
  Float Se_sl_raw(Float l, Float a, Float b, Float g){
    if(g < 1.0){
      return (1.0 - g) / b * lambertW_raw( b / (1.0 - g) * pow(a * l, 1 / (1.0 - g) ));
    }
    SL<Float> f;
    f.a = a; f.b = b; f.g = g;
    NEWTON_RESULT<Float> r = newton<Float, SL<Float> >(f, l, log(g) - log(b) + log(2.0));
    if(r.par > log(SAM_Zero) &&
       r.objective < 1.0e-6 &&
       fabs(r.grSR) < 1){ // Found stable equilibrium
      return exp(r.par);
    }else if(r.par > log(SAM_Zero) &&
	     r.objective < 1.0e-6 &&
	     fabs(r.grSR) >= 1){ // Found unstable equilibrium
      Rf_warning("Found unstable equilibrium");
      return exp(r.par);
    }else{
      return SAM_Zero * l;
    }
    return SAM_Zero * l;
  }

  TMB_BIND_ATOMIC(Se_sl0, 1111, Se_sl_raw(x[0],x[1],x[2],x[3]))


  /*
   * Specialization for the Type 2 depensatory Ricker
   */

  template<class Float>
  struct T2DR {
    Float a;
    Float b;
    Float g;
    
    template<class F>
    F operator()(F logs){
      return log((F)a) + logs - (F)b * exp(logs) + logs - log(exp(logs)+(F)g);
    }
  };
  
  template<class Float>
  Float Se_t2dr_raw(Float l, Float a, Float b, Float g){
    T2DR<Float> f;
    f.a = a; f.b = b; f.g = g;
    Float logSe0 = -log(b) + log(a + log(l));
    NEWTON_RESULT<Float> r = newton<Float, T2DR<Float> >(f, l, logSe0);
    if(r.par > log(SAM_Zero) &&
       r.objective < 1.0e-6 &&
       fabs(r.grSR) < 1){ // Found stable equilibrium
      return exp(r.par);
    }else if(r.par > log(SAM_Zero) &&
	     r.objective < 1.0e-6 &&
	     fabs(r.grSR) >= 1){ // Found unstable equilibrium
      Rf_warning("Found unstable equilibrium");
      return exp(r.par);
    }else{
      return SAM_Zero * l;
    }
    return SAM_Zero * l;
  }

  TMB_BIND_ATOMIC(Se_t2dr0, 1111, Se_t2dr_raw(x[0],x[1],x[2],x[3]))


  /*
   * Specialization for the Type 2 depensatory Hockey stick
   */

  template<class Float>
  struct T2DHS {
    Float a;
    Float b;
    Float g;
    
    template<class F>
    F operator()(F logs){
      return log((F)a) - log((F)b) + log(exp(logs) - 0.5*(exp(logs) - (F)b + fabs(exp(logs) - (F)b)))  +
	logs - log(exp(logs)+(F)g);
    }
  };
  
  template<class Float>
  Float Se_t2dhs_raw(Float l, Float a, Float b, Float g){
    T2DHS<Float> f;
    f.a = a; f.b = b; f.g = g;
    Float logSe0 = log(a) + log(l);
    NEWTON_RESULT<Float> r = newton<Float, T2DHS<Float> >(f, l, logSe0);
    if(r.par > log(SAM_Zero) &&
       r.objective < 1.0e-6 &&
       fabs(r.grSR) < 1){ // Found stable equilibrium
      return exp(r.par);
    }else if(r.par > log(SAM_Zero) &&
	     r.objective < 1.0e-6 &&
	     fabs(r.grSR) >= 1){ // Found unstable equilibrium
      Rf_warning("Found unstable equilibrium");
      return exp(r.par);
    }else{
      return SAM_Zero * l;
    }
    return SAM_Zero * l;
  }

  TMB_BIND_ATOMIC(Se_t2dhs0, 1111, Se_t2dhs_raw(x[0],x[1],x[2],x[3]))


  

  /*
   * Specialization for the non-increasing spline on R/S recruitment model
   */
  
 
  template<class Float>
  struct IBCD {
    vector<Float> knots;
    vector<Float> pars;
    matrix<Float> sg;
    
    template<class F>
    F operator()(F logs){
      return logs + ibcdspline_raw(logs, (vector<F>)knots.template cast<F>(), (vector<F>)pars.template cast<F>(), (matrix<F>)sg.template cast<F>());
    }
  };
  
  template<class Float>
  Float Se_ibcd_raw(CppAD::vector<Float> valIn){
    // Float l, vector<Float> knots, vector<Float> pars){
    vector<Float> val(valIn);
    Float l = val(0);
    int n = (int)trunc((val.size() - 3) / 2.0); // One more par than knot
    vector<Float> pars = val.segment(1,n+1);
    vector<Float> knots = val.segment(n+2,n);
    matrix<Float> sg = getSigAndGam_raw(knots);
    Float sv = 0.0;
    for(int i = 0; i < knots.size(); ++i)
      sv += knots(i);
    sv /= (Float)knots.size();
    IBCD<Float> f;
    f.knots = knots; f.pars = pars; f.sg = sg;
    NEWTON_RESULT<Float> r = newton<Float, IBCD<Float> >(f, l, sv);
    if(r.niter > 990){
      Rf_warning("Close to iteration limit");
      r = newton<Float, IBCD<Float> >(f, l, r.par);
    }
    if(r.par > log(SAM_Zero) &&
       r.objective < 1.0e-6 &&
       fabs(r.grSR) < 1){ // Found stable equilibrium
      return exp(r.par);
    }else if(r.par > log(SAM_Zero) &&
	     r.objective < 1.0e-6 &&
	     fabs(r.grSR) >= 1){ // Found unstable equilibrium
      Rf_warning("Found unstable equilibrium");
      return exp(r.par);
    }else{
      Rf_warning("Found no equilibrium");
      return SAM_Zero * l;	// Keep derivative info
    }
    Rf_warning("Found no equilibrium");
    return SAM_Zero * l;	// Keep derivative info
  }

  TMB_BIND_ATOMIC_FLEX_PART(Se_ibcd0, Se_ibcd_raw(x), 1 + (tx.size() - 3) / 2.0 + 1) // Only derivative of lambda and pars





  /*
   * Specialization for the integrated spline on R/S recruitment model
   */
  
 
  template<class Float>
  struct IBC {
    vector<Float> knots;
    vector<Float> pars;
    matrix<Float> sg;
    
    template<class F>
    F operator()(F logs){
      return logs + ibcspline_raw(logs, (vector<F>)knots.template cast<F>(), (vector<F>)pars.template cast<F>(), (matrix<F>)sg.template cast<F>());
    }
  };
  
  template<class Float>
  Float Se_ibc_raw(CppAD::vector<Float> valIn){
    vector<Float> val(valIn);
    Float l = val(0);
    int n = (int)trunc((val.size() - 2) / 2.0); // Same number of pars and knots
    vector<Float> pars = val.segment(1,n);
    vector<Float> knots = val.segment(n+1,n);
    matrix<Float> sg = getSigAndGam_raw(knots);
    Float sv = 0.0;
    for(int i = 0; i < knots.size(); ++i)
      sv += knots(i);
    sv /= (Float)knots.size();
    IBC<Float> f;
    f.knots = knots; f.pars = pars; f.sg = sg;
    NEWTON_RESULT<Float> r = newton<Float, IBC<Float> >(f, l, sv);
    if(r.niter > 990){
      Rf_warning("Close to iteration limit");
      r = newton<Float, IBC<Float> >(f, l, r.par);
    }    if(r.par > log(SAM_Zero) &&
	    r.objective < 1.0e-6 &&
	    fabs(r.grSR) < 1){ // Found stable equilibrium
      return exp(r.par);
    }else if(r.par > log(SAM_Zero) &&
	     r.objective < 1.0e-6 &&
	     fabs(r.grSR) >= 1){ // Found unstable equilibrium
      Rf_warning("Found unstable equilibrium");
      return exp(r.par);
    }else{
      Rf_warning("Found no equilibrium");
      return SAM_Zero * l;	// Keep derivative info
    }
    Rf_warning("Found no equilibrium");
    return SAM_Zero * l;	// Keep derivative info
  }

  TMB_BIND_ATOMIC_FLEX_PART(Se_ibc0, Se_ibc_raw(x), 1 + (tx.size() - 2) / 2.0) // Only derivative of lambda and pars


  
  /*
   * End of namespace
   */

  
}






template<class Type>
Type lambertW(Type x){
  vector<Type> args(2); // Last index reserved for derivative order
  args[0] = x;
  args[1] = 0;
  return rec_atomic::lambertW0(CppAD::vector<Type>(args))[0];
}

template<class Type>
Type Se_sbh(Type l, Type a, Type b, Type g){
  vector<Type> args(5); // Last index reserved for derivative order
  args[0] = l;
  args[1] = a;
  args[2] = b;
  args[3] = g;
  args[4] = 0;
  return rec_atomic::Se_sbh0(CppAD::vector<Type>(args))[0];
}

template<class Type>
Type Se_sl(Type l, Type a, Type b, Type g){
  vector<Type> args(5); // Last index reserved for derivative order
  args[0] = l;
  args[1] = a;
  args[2] = b;
  args[3] = g;
  args[4] = 0;
  return rec_atomic::Se_sl0(CppAD::vector<Type>(args))[0];
}


template<class Type>
Type Se_t2dr(Type l, Type a, Type b, Type g){
  vector<Type> args(5); // Last index reserved for derivative order
  args[0] = l;
  args[1] = a;
  args[2] = b;
  args[3] = g;
  args[4] = 0;
  return rec_atomic::Se_t2dr0(CppAD::vector<Type>(args))[0];
}

  template<class Type>
  Type Se_t2dhs(Type l, Type a, Type b, Type g){
  vector<Type> args(5); // Last index reserved for derivative order
  args[0] = l;
  args[1] = a;
  args[2] = b;
  args[3] = g;
  args[4] = 0;
  return rec_atomic::Se_t2dhs0(CppAD::vector<Type>(args))[0];
}


template<class Type>
Type Se_ibcd(Type l, vector<Type> knots, vector<Type> pars){
  vector<Type> args(1+knots.size()+pars.size()+1); // Last index reserved for derivative order
  args.setZero();
  args(0) = l;
  args.segment(1,pars.size()) = pars;
  args.segment(pars.size()+1,knots.size()) = knots;
  args(args.size()-1) = 0;
  vector<Type> tmp = rec_atomic::Se_ibcd0(CppAD::vector<Type>(args));
  return tmp(0);
}


template<class Type>
Type Se_ibc(Type l, vector<Type> knots, vector<Type> pars){
  vector<Type> args(1+knots.size()+pars.size()+1); // Last index reserved for derivative order
  args.setZero();
  args(0) = l;
  args.segment(1,pars.size()) = pars;
  args.segment(pars.size()+1,knots.size()) = knots;
  args(args.size()-1) = 0;
  vector<Type> tmp = rec_atomic::Se_ibc0(CppAD::vector<Type>(args));
  return tmp(0);
}



/*
 * Quick version for spline recruitment
 */

template<class Type>
Type softmax(Type x, Type y, Type k = 1.0){
  return logspace_add2(k * x, k * y) / k;
}


template<class Type>
Type sign0(Type x){
  return x / (fabs(x) + 1e-8);
}

template<class Type>
struct IBCD_QUICK {
  vector<Type> knots;
  vector<Type> pars;
  Type l;
  
 
  Type sr(Type s){
    return exp(log(s) + ibcdspline(log(s), knots, pars));
  }


  // Type fix(Type s){
  //   return l *  sr(s);
  // }
  Type fix(Type logs){
    return log(l + 1e-16) + logs + ibcdspline(logs, knots, pars);
  }

  Type f(Type logs){
    // l * R(S) = S
    // log(l) + logR(logS) = logS
    return log(l + 1e-16) + logs + ibcdspline(logs, knots, pars) - (logs); 
  }
  
  Type steffensen(Type logs){
    Type a = f(logs);
    Type b = f(logs + 0.99 * a);
    Type g = b / (0.99 * a) - 1.0;
    return logs - a / g;
  }
  
  Type numnewt(Type logs){
    Type h = 0.0001 * softmax(fabs(logs),Type(0.0001), Type(100.0));
    Type a = f(logs);
    Type g = (-f(logs + 2.0 * h) + 8.0 * f(logs + h) - 8.0 * f(logs-h) + f(logs-2.0*h)) / (12 * h);    
    // Type g = (f(logs + h) - f(logs - h)) / (2.0 * h);
    Type s = sign0(a) * sign0(g);
    Type y = log(fabs(a)) - log(softmax(fabs(g), (Type)0.001, (Type)1000.0)); // Damp the gradient
    return logs - 0.9 * s * exp(y); //softmax(exp(y), 0.5 * fabs(logs), Type(1000.0));
  }


  
};

template<class Type>
Type Se_ibcd_quick(Type l, vector<Type> knots, vector<Type> pars, int N){
  Type sv = 0.0;
  // for(int i = 0; i < knots.size(); ++i)
  //   sv += (knots(i));
  // sv /= (Type)knots.size();
  sv = knots((int)floor(knots.size() / 2));
  //sv = 0.75 * (sv - knots(0)) + knots(0);
  //sv = exp(sv);
  IBCD_QUICK<Type> f = {knots, pars, l};
  for(int i = 0; i < N; ++i){
    sv = f.numnewt(sv); //f.numnewt(sv); //f.fix(sv);
    // sv = softmax(tmp,Type(SAM_NIZero),(Type)10.0);
    // sv = tmp; //0.5 * (tmp + l * 100.0 + sqrt((tmp - l * 100.0) * (tmp-l * 100.0) + 1e-3));
  }
  sv = softmax(sv,Type(SAM_NIZero),(Type)100.0);
  return exp(sv);
}




template<class Type>
struct IBC_QUICK {
  vector<Type> knots;
  vector<Type> pars;
  Type l;
  
 
  Type sr(Type s){
    return exp(log(s) + ibcspline(log(s), knots, pars));
  }


  Type fix(Type logs){
    // return l *  sr(s);
    return log(l + 1e-8) + logs + ibcspline(logs, knots, pars);
  }

  Type f(Type logs){
    // l * R(S) = S
    // log(l) + logR(logS) = logS
    return log(l + 1e-8) + logs + ibcspline(logs, knots, pars) - (logs); 
  }

  Type numnewt(Type logs){
    Type h = 0.0001 * softmax(fabs(logs),Type(0.0001), Type(100.0));
    Type a = f(logs);
    Type g = (-f(logs + 2.0 * h) + 8.0 * f(logs + h) - 8.0 * f(logs-h) + f(logs-2.0*h)) / (12 * h);    
    // Type g = (f(logs + h) - f(logs - h)) / (2.0 * h);
    Type s = sign0(a) * sign0(g);
    Type y = log(fabs(a)) - log(softmax(fabs(g), (Type)0.001, (Type)1000.0)); // Damp the gradient
    return logs - 0.9 * s * exp(y); //softmax(exp(y), 0.5 * fabs(logs), Type(1000.0));
  }



};

template<class Type>
Type Se_ibc_quick(Type l, vector<Type> knots, vector<Type> pars, int N){
  Type sv = 0.0;
  for(int i = 0; i < knots.size(); ++i)
    sv += knots(i);
  sv /= (Type)knots.size();
  sv = 0.75 * (sv - knots(0)) + knots(0);
  //sv = exp(sv);
  IBC_QUICK<Type> f = {knots, pars, l};
  for(int i = 0; i < N; ++i){
    sv = f.numnewt(sv);
  }
  sv = softmax(sv,Type(SAM_NIZero),(Type)100.0);
  return exp(sv);
}


/*
 * Functions to get derivatives of stock recruitment curve
 */

// Functional form

#ifdef TMBAD_FRAMEWORK
template<class Type>
struct F_dFunctionalSR {
  vector<Type> rp;
  int srmc;
  template<class T>
  T operator()(vector<T> x){  // Evaluate function
    vector<T> rp2 = rp.template cast<T>();
    return exp(functionalStockRecruitment(x[0], rp2, srmc));
  }    
};
#endif

template<class Type>
Type dFunctionalSR(Type ssb, vector<Type> rp, int srmc){
  
#ifdef CPPAD_FRAMEWORK
  vector<AD<Type> > rp2(rp.size());
  rp2 = rp.template cast<AD<Type> >();
  CppAD::vector<AD<Type> > x( 1 );
  x[0] = ssb;
  CppAD::Independent(x);
  CppAD::vector<AD<Type> > y( 1 );
  y[0] = exp(functionalStockRecruitment(x[0], rp2, srmc));
  CppAD::ADFun<Type> F(x, y);
  CppAD::vector<Type> x_eval( 1 );
  x_eval[0] = ssb;
  vector<Type> r = F.Jacobian(x_eval);
  return r[0];
#endif
#ifdef TMBAD_FRAMEWORK
  F_dFunctionalSR<Type> Fd = {rp,srmc};
  vector<Type> x(1);
  x(0) = ssb;
  TMBad::ADFun<> G(TMBad::StdWrap<F_dFunctionalSR<Type> ,vector<TMBad::ad_aug> >(Fd), x);
  G = G.JacFun();
  return G(x)[0];
#endif
}

// Spline

#ifdef TMBAD_FRAMEWORK
template<class Type>
struct F_dSplineSR {
  vector<Type> knots;
  vector<Type> rp;
  int srmc;
  template<class T>
  T operator()(vector<T> x){  // Evaluate function
    T logx = log(x[0]);
    vector<T> rp2 = rp.template cast<T>();
    vector<T> kn2 = knots.template cast<T>();
    if(srmc == 90){
      return exp(logx + ibcdspline(logx, kn2, rp2)); 
    }else if(srmc == 91){
      return exp(logx + ibcspline(logx, kn2, rp2)); 
    }else if(srmc == 92){
      return exp(logx + bcspline(logx, kn2, rp2)); 
    }else{
      Rf_error("Not a spline recruitment model");
    }
    return 0.0;
  }    
};
#endif

template<class Type>
Type dSplineSR(Type ssb, vector<Type> knots, vector<Type> rp, int srmc){
  
#ifdef CPPAD_FRAMEWORK
  vector<AD<Type> > rp2(rp.size());
  rp2 = rp.template cast<AD<Type> >();
  vector<AD<Type> > kn2(knots.size());
  kn2 = knots.template cast<AD<Type> >();
  CppAD::vector<AD<Type> > x( 1 );
  x[0] = ssb;
  CppAD::Independent(x);
  CppAD::vector<AD<Type> > y( 1 );
  if(srmc == 90){
    y[0] = exp(log(x[0]) + ibcdspline(log(x[0]), kn2, rp2)); 
  }else if(srmc == 91){
    y[0] = exp(log(x[0]) + ibcspline(log(x[0]), kn2, rp2)); 
  }else if(srmc == 92){
    y[0] = exp(log(x[0]) + bcspline(log(x[0]), kn2, rp2)); 
  }else{
    Rf_error("Not a spline recruitment model");
  }
  CppAD::ADFun<Type> F(x, y);
  CppAD::vector<Type> x_eval( 1 );
  x_eval[0] = ssb;
  vector<Type> r = F.Jacobian(x_eval);
  return r[0];
#endif
#ifdef TMBAD_FRAMEWORK
  F_dSplineSR<Type> Fd = {knots,rp,srmc};
  vector<Type> x(1);
  x(0) = ssb;
  TMBad::ADFun<> G(TMBad::StdWrap<F_dSplineSR<Type> ,vector<TMBad::ad_aug> >(Fd), x);
  G = G.JacFun();
  return G(x)[0];
#endif
}


/*
 * Per-recruit calculations
 */



template<class Type>
struct PERREC_t {
  Type logFbar;
  Type logYPR;
  Type logSPR;
  Type logSe;
  Type logRe;
  Type logYe;
  Type dSR0;
  Type logLifeExpectancy;
  Type logYearsLost;
  Type logDiscYPR;
  Type logDiscYe;
};


template<class T>
T Se_SR52(T SPR, T a, T b, T d){
  T v0 = sqrt(SPR*SPR*a*a - 2.0 * SPR * a * b * d + b*b*d*d - 2 * a * SPR - 2.0 * b * d + 1.0);
  T v1 = (a * SPR - b*d - 1.0 + v0);
  // T v2 = (a * SPR - b*d - 1.0 - v0);
  // v1 must be larger
  return v1 / (2.0 * b);
}

template<class Type, class T>
PERREC_t<T> perRecruit(T logFbar, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, vector<Type>& logSel, vector<int> aveYears, int nYears = 300, int CT = 0){


#ifdef CPPAD_FRAMEWORK
  using CppAD::fabs;
#endif
  
  // Prepare data
  dataSet<T> newDat = dat.template cast<T>();
  int nMYears = dat.noYears;
  // propMat
  extendArray(newDat.propMat, nMYears, nYears, aveYears, false);
   // stockMeanWeight
  extendArray(newDat.stockMeanWeight, nMYears, nYears, aveYears, false);
  // catchMeanWeight
  extendArray(newDat.catchMeanWeight, nMYears, nYears, aveYears, false);
  // natMor
  extendArray(newDat.natMor, nMYears, nYears, aveYears, false);
  // landFrac
  extendArray(newDat.landFrac, nMYears, nYears, aveYears, false);
  // disMeanWeight
  extendArray(newDat.disMeanWeight, nMYears, nYears, aveYears, false);
  // landMeanWeight
  extendArray(newDat.landMeanWeight, nMYears, nYears, aveYears, false);
  // propF
  extendArray(newDat.propF, nMYears, nYears, aveYears, false);
  // propM
  extendArray(newDat.propM, nMYears, nYears, aveYears, false);
  newDat.noYears = nYears;

  // Prepare conf
  confSet newConf = conf;
  //// Set random walk recruitment
  newConf.stockRecruitmentModelCode = 0;

  // Prepare parameters
  paraSet<T> newPar = par.template cast<T>();

  vector<T> logSelT = logSel.template cast<T>();
  // Make logF array
  array<T> logF(logSelT.size(), nYears);
  logF.setZero();
  for(int i = 0; i < nYears; ++i)
    logF.col(i) = logSelT + logFbar;

  // Make logN array - start with one recruit
  int nAge = conf.maxAge - conf.minAge + 1;
  array<T> logN(nAge, nYears);
  logN.setConstant(SAM_NegInf);
  logN(0,0) = 0.0;

  // Run loop over years
  // typename referencepointSet<T>::RecCorrectionType RecCorrection = static_cast<typename referencepointSet<T>::RecCorrectionType>(RC);
  for(int i = 1; i < nYears; ++i){
    // predN
    logN.col(i) = predNFun(newDat, newConf, newPar, logN, logF, i);
    // for(int j = 1; j < logN.rows(); ++j){
      // if(RecCorrection == referencepointSet<T>::RecMean){
      // 	logN(j,i) += 0.5 * exp(2.0 * newPar.logSdLogN(conf.keyVarLogN(j)));
      // }else if(RecCorrection == referencepointSet<T>::RecMedian){
	
      // }else if(RecCorrection == referencepointSet<T>::RecMode){
      // 	logN(j,i) += -exp(newPar.logSdLogN(conf.keyVarLogN(j)));
      // }else{
      // 	Rf_error("Recruitment correction type not implemented");
      // }
    // }
    // remove recruitment
    logN(0,i) = SAM_NegInf;
  }

 
  // Calculate yield
//  #ifdef TMBAD_FRAMEWORK
//  vector<T> logcat(nYears);
//   logcat.setConstant(R_NegInf);
//   switch(newDat.referencepoint.catchType){
//   case referencepointSet<T>::totalCatch:
//     logcat = catchFun(newDat, newConf, logN, logF, true);
//     break;
//   case referencepointSet<T>::landings:
//     logcat = log(landFun(newDat, newConf, logN, logF));
//     break;
//   case referencepointSet<T>::discard:
//     logcat = log(disFun(newDat, newConf, logN, logF));
//     break;
//   default:
//     Rf_error("Unknown reference point catch type.");
//     break;
//   }
//   T logYPR = logspace_sum(logcat); //log(sum(cat) + (T)exp(-12.0));//
// #else
  vector<T> cat(nYears);
  cat.setZero();
  typename referencepointSet<T>::CatchType catchType = static_cast<typename referencepointSet<T>::CatchType>(CT);
  switch(catchType){
  case referencepointSet<T>::totalCatch:
    cat = catchFun(newDat, newConf, logN, logF);
    break;
  case referencepointSet<T>::landings:
    cat = landFun(newDat, newConf, logN, logF);
    break;
  case referencepointSet<T>::discard:
    cat = disFun(newDat, newConf, logN, logF);
    break;
  default:
    Rf_error("Unknown reference point catch type.");
    break;
  }
  T logYPR = log(sum(cat) + (T)exp(SAM_NegInf));//

  T logYLTF = log(yearsLostFishing_i(newDat, newConf, logF, newDat.natMor.dim(0)-1, conf.minAge, conf.maxAge));
  T logLifeExpectancy = log(temporaryLifeExpectancy_i(newDat, newConf, logF, newDat.natMor.dim(0)-1, conf.minAge, 10 * conf.maxAge) + (T)conf.minAge);
    
  // Calculate spawners
  //vector<T> logssb = ssbFun(newDat, newConf, logN, logF, true);
  vector<T> ssb = ssbFun(newDat, newConf, logN, logF);
  T logSPR = log(softmax(sum(ssb),(T)SAM_Zero,(T)1000.0)); //log(sum(ssb)); log(sum(ssb) + (T)exp(-12.0));

 ////////////////////////////////////////////////////////////////////////////////
  // Survival calculations                                                      //
  ////////////////////////////////////////////////////////////////////////////////
  T discYPR = 0.0;
  // T discYe = 0.0;
  // Custom recursive calculation inspired by survival.hpp
  T yl_logp = 0.0;
  T yl_q = 0.0;
  T yl = 0.0;
  // vector<T> catYe = cat * exp(logRe);
  for(int aa = 0; aa < cat.size(); ++aa){
    int j = std::min(std::max(aa-newConf.minAge,0), newDat.natMor.cols()-1);		// Cohort age index
    int i = std::min(aa, newDat.natMor.rows()-1);
    T M = newDat.natMor(i,j);
    T F = 0.0;    
    if(newConf.keyLogFsta(0,j)>(-1) && aa >= newConf.minAge)
      F += exp(logF(newConf.keyLogFsta(0,j),i));
    T Z = M + F;
    yl += yl_q + exp(yl_logp) * F / Z * (1.0 - 1.0 / Z * (1.0 - exp(-Z)));
    discYPR += cat(aa) * exp(-yl);
    // discYe += catYe(aa) * exp(-yl);
    yl_q += exp(yl_logp) * F / Z * (1.0 - exp(-Z));
    yl_logp += -Z;
  }
  T logDiscYPR = log(discYPR + (T)exp(SAM_NegInf));
  // T logDiscYe = log(discYe + (T)exp(SAM_NegInf));
  ////////////////////////////////////////////////////////////////////////////////
  // End survival calculations                                                  //
  ////////////////////////////////////////////////////////////////////////////////

  

  //T logSPR = logspace_sum(logssb);
  T lambda = exp(logSPR); // sum(ssb);

  if(conf.stockRecruitmentModelCode == 0){//  ||
     // conf.stockRecruitmentModelCode == 3){
    PERREC_t<T> res = {logFbar, // logFbar
		       logYPR,	   // logYPR
		       logSPR,	   // logSPR
		       R_NaReal,	   // logSe
		       R_NaReal,		 // logRe
		       R_NaReal,// logYe
		       R_NaReal,// dSR0
		       logLifeExpectancy, // logLifeExpectancy
		       logYLTF, // logYearsLost
		       logDiscYPR, // logDiscYPR
		       R_NaReal	 // logDiscYe
    }; 
    return res;
  }
  
  // Calculate Se
  T Se = SAM_Zero; //R_NegInf;
  
  T dsr0 = 10000.0;
  if(conf.stockRecruitmentModelCode != 0 &&
     conf.stockRecruitmentModelCode != 3 &&
     conf.stockRecruitmentModelCode != 62 &&
     conf.stockRecruitmentModelCode != 65 &&
     (conf.stockRecruitmentModelCode != 68) && // || newPar.rec_pars[2] < 0) &&
     (conf.stockRecruitmentModelCode != 69) &&
     conf.stockRecruitmentModelCode != 90 &&
     conf.stockRecruitmentModelCode != 91 &&
     conf.stockRecruitmentModelCode != 92){ // || newPar.rec_pars[2] < 0 ))
    dsr0 = dFunctionalSR(T(SAM_Zero), newPar.rec_pars, conf.stockRecruitmentModelCode);
  }
  if(conf.stockRecruitmentModelCode == 68 || conf.stockRecruitmentModelCode == 69){
    dsr0 = CppAD::CondExpLt(newPar.rec_pars[2],
			    (T)0.0,
			    dFunctionalSR(T(SAM_Zero),
					  newPar.rec_pars,
					  conf.stockRecruitmentModelCode),
			    dsr0);
  }

  if(conf.stockRecruitmentModelCode == 90 ||
     conf.stockRecruitmentModelCode == 91 ||
     conf.stockRecruitmentModelCode == 92){
    // dSplineSR uses logssb
    dsr0 = dSplineSR((T)SAM_Zero,
		     (vector<T>)newConf.constRecBreaks.template cast<T>(),
		     newPar.rec_pars,
		     conf.stockRecruitmentModelCode);
  }

  T logRecCorrection = 0.0;
  // if(RecCorrection == referencepointSet<T>::RecMean){
  //   logRecCorrection = 0.5 * exp(2.0 * newPar.logSdLogN(conf.keyVarLogN(0)));
  // }else if(RecCorrection == referencepointSet<T>::RecMedian){
  //   logRecCorrection = 0.0;
  // }else if(RecCorrection == referencepointSet<T>::RecMode){
  //   logRecCorrection = -exp(newPar.logSdLogN(conf.keyVarLogN(0)));
  // }else{
  //   Rf_error("Recruitment correction type not implemented");
  // }
    
  switch(conf.stockRecruitmentModelCode){
  case 0: // straight RW 
    Rf_error("Equilibrium SSB not implemented");
    // Projecting forward, the mean on log-scale will always be last years recruitment
    // Not sure it is a good idea to implement
    //Se = lambda * newPar.logN(0,newPar.logN.cols());
    break;
  case 1: //ricker
    Se = exp(-newPar.rec_pars(1)) * (newPar.rec_pars(0) + logRecCorrection + log(lambda)); //log((exp(newPar.rec_pars(0)) * lambda));
    break;
  case 2:  //BH
    // Handle negative values below!
    Se = (exp(newPar.rec_pars(0) + logRecCorrection) * lambda - 1.0) * exp(-newPar.rec_pars(1));
      // CppAD::CondExpGt((newPar.rec_pars(0) + logSPR), T(SAM_Zero),
      // 			   fabs(exp(newPar.rec_pars(0)) * lambda - 1.0) * exp(-newPar.rec_pars(1)),
      // 			   T(SAM_Zero));
    break;
  case 3: //Constant mean
    // Constant recruitment - last year of assessment
    Se = lambda * exp(newPar.rec_pars(newPar.rec_pars.size() - 1) + logRecCorrection);
    break;
  case 50:			// Type 2 depensatory logistic hockey stick
    Rf_error("Equilibrium not implemented yet");
    break;
  case 51:			// Type 2 depensatory Ricker
    Se = Se_t2dr(lambda, exp(newPar.rec_pars(0) + logRecCorrection), exp(newPar.rec_pars(1)), exp(newPar.rec_pars(2)));
    break;
  case 52:			// Type 2 depensatory BH
    Se = Se_SR52(lambda, exp(newPar.rec_pars(0) + logRecCorrection), exp(newPar.rec_pars(1)), exp(newPar.rec_pars(2)));
    break;
  case 53:			// Type 2 depensatory BH
    Se = Se_t2dhs(lambda, exp(newPar.rec_pars(0) + logRecCorrection), exp(newPar.rec_pars(1)), exp(newPar.rec_pars(2)));
    break;
  case 60:			// logistic hockey stick
    Rf_error("Equilibrium not implemented yet");
    break;    
  case 61: // Hockey stick
    Se = lambda * exp(newPar.rec_pars(0) + logRecCorrection);
    break;
  case 62: // AR1 (on log-scale)
    Se = lambda * exp(newPar.rec_pars(0) + logRecCorrection);
    break;
  case 63: //Bent hyperbola / Hockey-stick-like
    Se = (2.0 * sqrt(exp(2.0 * newPar.rec_pars(0)) + exp(2.0 * newPar.rec_pars(2)) / 4.0) / (lambda * exp(newPar.rec_pars(1) + logRecCorrection)) - 2.0 * exp(newPar.rec_pars(0)) - 2.0 * sqrt(exp(2.0 * newPar.rec_pars(0)) + exp(2.0 * newPar.rec_pars(2)) / 4.0)) / ( 1.0 / ((lambda * lambda * exp(2.0 * (newPar.rec_pars(1) + logRecCorrection)))) - 2.0 / (lambda * exp(newPar.rec_pars(1) + logRecCorrection))  );  
    break;
  case 64: // Power CMP
    Se = exp(1.0 / (1.0 - invlogit(newPar.rec_pars(1))) * (newPar.rec_pars(0) + logRecCorrection + log(lambda)));
    break;
  case 65: // Power Non-CMP
    Se = exp(1.0 / (1.0 - (exp(newPar.rec_pars(1)) + 1.0001)) * (newPar.rec_pars(0) + logRecCorrection + log(lambda)));
    break;
  case 66: // Shepherd
    // Se = CppAD::CondExpGt((newPar.rec_pars(0) + logSPR), T(SAM_Zero),
    // 			  exp( newPar.rec_pars(1) + 1.0 / exp(newPar.rec_pars(2)) * log(fabs(exp(newPar.rec_pars(0)) * lambda - 1.0))),
    // 			  T(SAM_Zero));
    Se = exp( newPar.rec_pars(1) + 1.0 / exp(newPar.rec_pars(2)) * log(softmax(exp(newPar.rec_pars(0) + logRecCorrection) * lambda - 1.0,(T)SAM_Zero, (T)100.0)) );
    break;
  case 67: // Deriso
    // Handle negative values below!
    // Se = CppAD::CondExpGt((newPar.rec_pars(0) + logSPR), T(SAM_Zero),
    // 			  fabs(exp(exp(-newPar.rec_pars(2)) * (newPar.rec_pars(0) + logSPR)) - 1.0) * exp(-newPar.rec_pars(1)),
    // 			   T(SAM_Zero));
    Se = (exp(exp(-newPar.rec_pars(2)) * (newPar.rec_pars(0) + logRecCorrection + logSPR)) - 1.0) * exp(-newPar.rec_pars(1) - newPar.rec_pars(2));
    break;
  case 68: // Saila-Lorda (cases: gamma > 1; gamma = 1; gamma < 1)
    Se = Se_sl(lambda, exp(newPar.rec_pars(0) + logRecCorrection), exp(newPar.rec_pars(1)), exp(newPar.rec_pars(2)));
    break;
  case 69: // Sigmoidal Beverton-Holt (cases: gamma < 1; gamma >= 1)
    Se = Se_sbh(lambda, exp(newPar.rec_pars(0) + logRecCorrection), exp(newPar.rec_pars(1)), exp(newPar.rec_pars(2)));
    break;
  case 90: // Non-increasing spline on log R/S
    Se = Se_ibcd_quick(lambda,(vector<T>)newConf.constRecBreaks.template cast<T>(), newPar.rec_pars, dat.referencepoint.optN);
    break;
  case 91: // integrated spline on log R/S
    Se = Se_ibc_quick(lambda,(vector<T>)newConf.constRecBreaks.template cast<T>(), newPar.rec_pars, dat.referencepoint.optN);
    break;
  case 92: // spline on log R/S
    Rf_error("Not implemented yet");
    break;
  case 93:
    Rf_error("Not implemented yet");
    break;
    default:
      Rf_error("SR model code not recognized");
    break;   
  }
  //T logSe = log(Se);
  if(conf.stockRecruitmentModelCode == 63 ||
     conf.stockRecruitmentModelCode == 66 ||
     conf.stockRecruitmentModelCode == 90)
    Se = CppAD::CondExpGt(-logSPR, log(dsr0), (T)SAM_NegInf, Se);
  
  // T logSe = CppAD::CondExpGt(exp(-logSPR), dsr0 - (T)1e-3,
  // 			     log(fabs(Se)) - 3.0 * (exp(-logSPR) - dsr0),
  // 			     log(fabs(Se)));
  T logSe = log(softmax(Se, (T)exp(SAM_NegInf), (T)1000.0));
  
  // T logYe = CppAD::CondExpGt(exp(-logSPR), dsr0,
  // 			     logSe - logSPR + logYPR - 3.0 * (exp(-logSPR) - dsr0),
  // 			     logSe - logSPR + logYPR);
  T logYe = logSe - logSPR + logYPR;
  //log(softmax(exp(logSe - logSPR + logYPR), (T)SAM_Zero, (T)1.0));
  T logRe = logSe - logSPR;
  //log(softmax(exp(logSe - logSPR), (T)SAM_NegInf, (T)1.0));
  T logDiscYe = logDiscYPR + logRe;
 
  
  // Return
  PERREC_t<T> res = {logFbar, // logFbar
		     logYPR,	// logYPR
		     logSPR,	// logSPR
		     logSe,	// logSe
		     logRe,	// logRe
		     logYe,	// logYe
		     dsr0,
		     logLifeExpectancy, // logLifeExpectancy
		     logYLTF, // logYearsLost
		     logDiscYPR,
		     logDiscYe,
  };	// DSR0

  return res;
}


/*
 * struct for reference points
 */

#ifdef TMBAD_FRAMEWORK
template<class Type>
struct REFERENCE_POINTS;
// Structs for REFERENCE_POINTS AD
// YPR
template<class Type>
struct F_FYPR {
  REFERENCE_POINTS<Type>& parent;
  template<class T>
  T operator()(vector<T> x){  // Evaluate function
    return parent.YPR(x[0]);
  }
};

// SPR
template<class Type>
struct F_FSPR {
  REFERENCE_POINTS<Type>& parent;
  template<class T>
  T operator()(vector<T> x){  // Evaluate function
    return parent.SPR(x[0]);
  }
};

// FSR
template<class Type>
struct F_FSR {
  REFERENCE_POINTS<Type>& parent;
  template<class T>
  T operator()(vector<T> x){  // Evaluate function
    return parent.SR(x[0]);
  }
};

#endif




template<class Type>
struct REFERENCE_POINTS {


  

  // Input data
  dataSet<Type>& dat;
  confSet& conf;
  paraSet<Type>& par;

  int nYears;
  vector<int> aveYears;
  vector<int> selYears;

  array<Type>& logN;
  array<Type>& logF;
  
  // Input F

  Type logFsq;			// Status quo
  Type logF0;			// "No" fishing
  Type logFmsy; 		// Maximizes yield
  Type logFmypyl;		// Maximizes yield per year lost
  Type logFmdy;		        // Maximizes discounted yield per year lost
  Type logFmax;			// Maximizes yield per recruit
  vector<Type> logFxdYPR;			// F such that YPR'(0) = 0.1 * YPR'(F)
  Type logFcrash;		// F such that 1/SPR(f) = SR'(0) (i.e. stock crashes [with compensatory recruitment] if slope of spawner-per-recruit in origin is less than slope of stock-recruitment model in origin)
  Type logFext;			// F such that stock dies out - SSB(F) < epsilon and F smallest possible
  vector<Type> logFxPercent;	// F such that SSB is reduced to x% of unfished stock (Se(F) = x/100 * Se(0) )
  //Type logFmed;		// Fishing  mortality  rate  F  corresponding  to  a  SSB/R  equal  to  the  inverse  of  the  50th  percentile of the observed R/SSB
  Type logFlim;			// F such that Se(F) = Blim (for hockey-stick-like stock recruitment only)
  //Type logFpa;		// F that corresponds to logBpa
  matrix<Type> logFmsyRange;		// F that corresponds to yield of x% of MSY
  vector<Type> logFxB0;
  
  // Corresponding SSB??
  Type logBsq;
  Type logB0;
  Type logBmsy;
  Type logBmypyl;
  Type logBmdy;
  Type logBmax;
  vector<Type> logBxdYPR;
  Type logBcrash;
  Type logBext;
  vector<Type> logBxPercent;
  Type logBlim;			// Known from model parameters (for hockey-stick-like stock recruitment only)
  //Type logBpa;			// Ba = Blim * exp(1.645 * sigma) where sigma is the standard deviation of log(SSB) at the start of the year following the terminal year of the assessment if sigma is unknown, 0.2 can be used as default.
  matrix<Type> logBmsyRange;
  vector<Type> logBxB0;

  // Corresponding Recruitment
  Type logRsq;
  Type logR0;
  Type logRmsy;
  Type logRmypyl;
  Type logRmdy;
  Type logRmax;
  vector<Type> logRxdYPR;
  Type logRcrash;
  Type logRext;
  vector<Type> logRxPercent;
  Type logRlim;			// Known from model parameters (for hockey-stick-like stock recruitment only)
  //Type logBpa;			// Ba = Blim * exp(1.645 * sigma) where sigma is the standard deviation of log(SSB) at the start of the year following the terminal year of the assessment if sigma is unknown, 0.2 can be used as default.
  matrix<Type> logRmsyRange;
  vector<Type> logRxB0;

  // Corresponding Yield??
  Type logYsq;
  Type logY0;
  Type logYmsy;
  Type logYmypyl;
  Type logYmdy;
  Type logYmax;
  vector<Type> logYxdYPR;
  Type logYcrash;
  Type logYext;
  vector<Type> logYxPercent;
  Type logYlim;			// Known from model parameters (for hockey-stick-like stock recruit
  matrix<Type> logYmsyRange;
  vector<Type> logYxB0;

  Type logYPRsq;
  Type logYPR0;
  Type logYPRmsy;
  Type logYPRmypyl;
  Type logYPRmdy;
  Type logYPRmax;
  vector<Type> logYPRxdYPR;
  Type logYPRcrash;
  Type logYPRext;
  vector<Type> logYPRxPercent;
  Type logYPRlim;
  matrix<Type> logYPRmsyRange;
  vector<Type> logYPRxB0;

  Type logSPRsq;
  Type logSPR0;
  Type logSPRmsy;
  Type logSPRmypyl;
  Type logSPRmdy;
  Type logSPRmax;
  vector<Type> logSPRxdYPR;
  Type logSPRcrash;
  Type logSPRext;
  vector<Type> logSPRxPercent;
  Type logSPRlim;
  matrix<Type> logSPRmsyRange;
  vector<Type> logSPRxB0;

  
  // Derived values
  vector<Type> logSel;

#ifdef CPPAD_FRAMEWORK
  CppAD::ADFun<Type> FSR;
  CppAD::ADFun<Type> FYPR;
  CppAD::ADFun<Type> FSPR;
#endif
#ifdef TMBAD_FRAMEWORK
  // ADFuns
  TMBad::ADFun<> FSR;
  TMBad::ADFun<> FYPR;
  TMBad::ADFun<> FSPR;
#endif
  
  REFERENCE_POINTS(){}
  REFERENCE_POINTS(dataSet<Type>& dat_,
		   confSet& conf_,
		   paraSet<Type>& par_,
		   array<Type>& logN_,
		   array<Type>& logF_,
		   int fromSQ
		   ): dat(dat_), conf(conf_), par(par_),
		      logN(logN_),
		      logF(logF_)
  {
    nYears = dat.referencepoint.nYears;
    aveYears = dat.referencepoint.aveYears;
    selYears = dat.referencepoint.selYears;

    // Calculate current selectivity and status quo
    vector<Type> logfbartmp(selYears.size());
    logfbartmp.setConstant(R_NegInf);
    Type logfsum = R_NegInf;
    logSel = vector<Type>(logF.rows());
    logSel.setConstant(R_NegInf);

    for(int y = 0; y < selYears.size(); ++y){
      for(int i = 0; i < logSel.size(); ++i)
	logSel(i) = logspace_add2(logSel(i), logF(i,selYears(y)));
      for(int a = conf.fbarRange(0); a <= conf.fbarRange(1); a++){  
	logfbartmp(y) = logspace_add2(logfbartmp(y), logF(conf.keyLogFsta(0,a-conf.minAge),selYears(y)));
      }
      logfbartmp(y) -= log(Type(conf.fbarRange(1)-conf.fbarRange(0)+1));
      logfsum = logspace_add2(logfsum, logfbartmp(y));
    }

    logFsq = logfsum - log(logfbartmp.size());//logfbartmp(logfbartmp.size() - 1); //log(sum(fbartmp)) - log(fbartmp.size());
    logSel -= logfsum;
    
    logBsq = log(Se(exp(logFsq)));
    logRsq = log(Re(exp(logFsq)));
    logYsq = log(yield(exp(logFsq)));
    logYPRsq = log(YPR(exp(logFsq)));
    logSPRsq = log(SPR(exp(logFsq)));

    logF0 = SAM_NegInf;
    logB0 = log(Se(exp(logF0)));
    logR0 = log(Re(exp(logF0)));
    logY0 = log(yield(exp(logF0)));
    logYPR0 = log(YPR(exp(logF0)));
    logSPR0 = log(SPR(exp(logF0)));

    // Calculate actual F values
    if(CppAD::Variable(par.logScaleFmypyl)){
      logFmypyl = par.logScaleFmypyl + logFsq; // logFsq + par.logScaleFmsy;
      logBmypyl = log(Se(exp(logFmypyl)));
      logRmypyl = log(Re(exp(logFmypyl)));
      logYmypyl = log(yield(exp(logFmypyl)));
      logYPRmypyl = log(YPR(exp(logFmypyl)));
      logSPRmypyl = log(SPR(exp(logFmypyl)));
    }else{
      logFmypyl = R_NaReal;// R_NaReal;
      logBmypyl = R_NaReal;
      logRmypyl = R_NaReal;
      logYmypyl = R_NaReal;
      logYPRmypyl = R_NaReal;
      logSPRmypyl = R_NaReal;
    }

        // Calculate actual F values
    if(CppAD::Variable(par.logScaleFmdy)){
      logFmdy = par.logScaleFmdy + logFsq * Type(fromSQ); // logFsq + par.logScaleFmsy;
      logBmdy = log(Se(exp(logFmdy)));
      logRmdy = log(Re(exp(logFmdy)));
      logYmdy = log(yield(exp(logFmdy)));
      logYPRmdy = log(YPR(exp(logFmdy)));
      logSPRmdy = log(SPR(exp(logFmdy)));
    }else{
      logFmdy = R_NaReal;// R_NaReal;
      logBmdy = R_NaReal;
      logRmdy = R_NaReal;
      logYmdy = R_NaReal;
      logYPRmdy = R_NaReal;
      logSPRmdy = R_NaReal;
    }

    
    // Calculate actual F values
    if(CppAD::Variable(par.logScaleFmsy)){
      logFmsy = par.logScaleFmsy + logFsq * Type(fromSQ); // logFsq + par.logScaleFmsy;
      logBmsy = log(Se(exp(logFmsy)));
      logRmsy = log(Re(exp(logFmsy)));
      logYmsy = log(yield(exp(logFmsy)));
      logYPRmsy = log(YPR(exp(logFmsy)));
      logSPRmsy = log(SPR(exp(logFmsy)));
    }else{
      logFmsy = R_NaReal;// R_NaReal;
      logBmsy = R_NaReal;
      logRmsy = R_NaReal;
      logYmsy = R_NaReal;
      logYPRmsy = R_NaReal;
      logSPRmsy = R_NaReal;
    }
    
    // MSY range lower
    logFmsyRange = matrix<Type>(par.logScaleFmsyRange.rows(),par.logScaleFmsyRange.cols());
    logBmsyRange = matrix<Type>(par.logScaleFmsyRange.rows(),par.logScaleFmsyRange.cols());
    logRmsyRange = matrix<Type>(par.logScaleFmsyRange.rows(),par.logScaleFmsyRange.cols());
    logYmsyRange = matrix<Type>(par.logScaleFmsyRange.rows(),par.logScaleFmsyRange.cols());
    logYPRmsyRange = matrix<Type>(par.logScaleFmsyRange.rows(),par.logScaleFmsyRange.cols());
    logSPRmsyRange = matrix<Type>(par.logScaleFmsyRange.rows(),par.logScaleFmsyRange.cols());

    for(int i = 0; i < par.logScaleFmsyRange.cols(); ++i){
      // Lower end of range
      if(CppAD::Variable(par.logScaleFmsy) && CppAD::Variable(par.logScaleFmsyRange(0,i))){
	logFmsyRange(0,i) = logFmsy - exp(-par.logScaleFmsyRange(0,i)); //logFsq + par.logScaleF35;
	logBmsyRange(0,i) = log(Se(exp(logFmsyRange(0,i))));
	logRmsyRange(0,i) = log(Re(exp(logFmsyRange(0,i))));
	logYmsyRange(0,i) = log(yield(exp(logFmsyRange(0,i))));
	logYPRmsyRange(0,i) = log(YPR(exp(logFmsyRange(0,i))));
	logSPRmsyRange(0,i) = log(SPR(exp(logFmsyRange(0,i))));
      }else{
	logFmsyRange(0,i) = R_NaReal;//R_NaReal;
	logBmsyRange(0,i) = R_NaReal;
	logRmsyRange(0,i) = R_NaReal;
	logYmsyRange(0,i) = R_NaReal;
	logYPRmsyRange(0,i) = R_NaReal;
	logSPRmsyRange(0,i) = R_NaReal;
      }
      // Upper end of range
      if(CppAD::Variable(par.logScaleFmsy) && CppAD::Variable(par.logScaleFmsyRange(1,i))){
	logFmsyRange(1,i) = logFmsy + exp(par.logScaleFmsyRange(1,i)); //logFsq + par.logScaleF35;
	logBmsyRange(1,i) = log(Se(exp(logFmsyRange(1,i))));
	logRmsyRange(1,i) = log(Re(exp(logFmsyRange(1,i))));
	logYmsyRange(1,i) = log(yield(exp(logFmsyRange(1,i))));
	logYPRmsyRange(1,i) = log(YPR(exp(logFmsyRange(1,i))));
	logSPRmsyRange(1,i) = log(SPR(exp(logFmsyRange(1,i))));
      }else{
	logFmsyRange(1,i) = R_NaReal;//R_NaReal;
	logBmsyRange(1,i) = R_NaReal;
	logRmsyRange(1,i) = R_NaReal;
	logYmsyRange(1,i) = R_NaReal;
	logYPRmsyRange(1,i) = R_NaReal;
	logSPRmsyRange(1,i) = R_NaReal;
      }
    }

    
    if(CppAD::Variable(par.logScaleFmax)){
      logFmax = par.logScaleFmax + logFsq * Type(fromSQ); //logFsq + par.logScaleFmax;
      logBmax = log(Se(exp(logFmax)));
      logRmax = log(Re(exp(logFmax)));
      logYmax = log(yield(exp(logFmax)));
      logYPRmax = log(YPR(exp(logFmax)));
      logSPRmax = log(SPR(exp(logFmax)));
    }else{
      logFmax = R_NaReal;//R_NaReal;
      logBmax = R_NaReal;
      logRmax = R_NaReal;
      logYmax = R_NaReal;
      logYPRmax = R_NaReal;
      logSPRmax = R_NaReal;
    }


    logFxdYPR = vector<Type>(par.logScaleFxdYPR.size());
    logBxdYPR = vector<Type>(par.logScaleFxdYPR.size());
    logRxdYPR = vector<Type>(par.logScaleFxdYPR.size());
    logYxdYPR = vector<Type>(par.logScaleFxdYPR.size());
    logYPRxdYPR = vector<Type>(par.logScaleFxdYPR.size());
    logSPRxdYPR = vector<Type>(par.logScaleFxdYPR.size());
    for(int i = 0; i < par.logScaleFxdYPR.size(); ++i){
      if(CppAD::Variable(par.logScaleFxdYPR(i))){
	logFxdYPR(i) = par.logScaleFxdYPR(i) + logFsq * Type(fromSQ); //logFsq + par.logScaleF01;
	logBxdYPR(i) = log(Se(exp(logFxdYPR(i))));
	logRxdYPR(i) = log(Re(exp(logFxdYPR(i))));
	logYxdYPR(i) = log(yield(exp(logFxdYPR(i))));
	logYPRxdYPR(i) = log(YPR(exp(logFxdYPR(i))));
	logSPRxdYPR(i) = log(SPR(exp(logFxdYPR(i))));
      }else{
	logFxdYPR(i) = R_NaReal;//R_NaReal;
	logBxdYPR(i) = R_NaReal;
	logRxdYPR(i) = R_NaReal;
	logYxdYPR(i) = R_NaReal;
	logYPRxdYPR(i) = R_NaReal;
	logSPRxdYPR(i) = R_NaReal;
      }
    }

    if(CppAD::Variable(par.logScaleFcrash)){
      logFcrash = par.logScaleFcrash + logFsq * Type(fromSQ); //logFsq + par.logScaleFcrash;
      logBcrash = log(Se(exp(logFcrash)));
      logRcrash = log(Re(exp(logFcrash)));
      logYcrash = log(yield(exp(logFcrash)));
      logYPRcrash = log(YPR(exp(logFcrash)));
      logSPRcrash = log(SPR(exp(logFcrash)));
    }else{
      logFcrash = R_NaReal;//R_NaReal;
      logBcrash = R_NaReal;
      logRcrash = R_NaReal;
      logYcrash = R_NaReal;
      logYPRcrash = R_NaReal;
      logSPRcrash = R_NaReal;
    }

    
    if(CppAD::Variable(par.logScaleFext)){
      logFext = par.logScaleFext + logFsq * Type(fromSQ); //logFsq + par.logScaleFcrash;
      logBext = log(Se(exp(logFext)));
      logRext = log(Re(exp(logFext)));
      logYext = log(yield(exp(logFext)));
      logYPRext = log(YPR(exp(logFext)));
      logSPRext = log(SPR(exp(logFext)));
    }else{
      logFext = R_NaReal;//R_NaReal;
      logBext = R_NaReal;
      logRext = R_NaReal;
      logYext = R_NaReal;
      logYPRext = R_NaReal;
      logSPRext = R_NaReal;
    }

    logFxPercent = vector<Type>(par.logScaleFxPercent.size());
    logBxPercent = vector<Type>(par.logScaleFxPercent.size());
    logRxPercent = vector<Type>(par.logScaleFxPercent.size());
    logYxPercent = vector<Type>(par.logScaleFxPercent.size());
    logYPRxPercent = vector<Type>(par.logScaleFxPercent.size());
    logSPRxPercent = vector<Type>(par.logScaleFxPercent.size());

    for(int i = 0; i < par.logScaleFxPercent.size(); ++i){
      if(CppAD::Variable(par.logScaleFxPercent(i))){
	logFxPercent(i) = par.logScaleFxPercent(i) + logFsq; //logFsq + par.logScaleF35;
	logBxPercent(i) = log(Se(exp(logFxPercent(i))));
	logRxPercent(i) = log(Re(exp(logFxPercent(i))));
	logYxPercent(i) = log(yield(exp(logFxPercent(i))));
	logYPRxPercent(i) = log(YPR(exp(logFxPercent(i))));
	logSPRxPercent(i) = log(SPR(exp(logFxPercent(i))));
      }else{
	logFxPercent(i) = R_NaReal;//R_NaReal;
	logBxPercent(i) = R_NaReal;
	logRxPercent(i) = R_NaReal;
	logYxPercent(i) = R_NaReal;
	logYPRxPercent(i) = R_NaReal;
	logSPRxPercent(i) = R_NaReal;
      }
    }
    
    logFxB0 = vector<Type>(par.logScaleFxB0.size());
    logBxB0 = vector<Type>(par.logScaleFxB0.size());
    logRxB0 = vector<Type>(par.logScaleFxB0.size());
    logYxB0 = vector<Type>(par.logScaleFxB0.size());
    logYPRxB0 = vector<Type>(par.logScaleFxB0.size());
    logSPRxB0 = vector<Type>(par.logScaleFxB0.size());
    for(int i = 0; i < par.logScaleFxB0.size(); ++i){
      if(CppAD::Variable(par.logScaleFxB0(i))){
	logFxB0(i) = par.logScaleFxB0(i) + logFsq * Type(fromSQ); //logFsq + par.logScaleF01;
	logBxB0(i) = log(Se(exp(logFxB0(i))));
	logRxB0(i) = log(Re(exp(logFxB0(i))));
	logYxB0(i) = log(yield(exp(logFxB0(i))));
	logYPRxB0(i) = log(YPR(exp(logFxB0(i))));
	logSPRxB0(i) = log(SPR(exp(logFxB0(i))));
      }else{
	logFxB0(i) = R_NaReal;//R_NaReal;
	logBxB0(i) = R_NaReal;
	logRxB0(i) = R_NaReal;
	logYxB0(i) = R_NaReal;
	logYPRxB0(i) = R_NaReal;
	logSPRxB0(i) = R_NaReal;
      }
    }


    
    if(CppAD::Variable(par.logScaleFlim) &&
       (conf.stockRecruitmentModelCode == 61 ||
	conf.stockRecruitmentModelCode == 63)){    
      logFlim = par.logScaleFlim + logFsq; //logFsq + par.logScaleFlim;
      logYlim = log(yield(exp(logFlim)));
      if(conf.stockRecruitmentModelCode == 61){
	logBlim = par.rec_pars(1);
      }else if(conf.stockRecruitmentModelCode == 63){
	logBlim = par.rec_pars(0);      
      }
      logRlim = log(Re(exp(logFlim)));
      logYPRlim = log(YPR(exp(logFlim)));
      logSPRlim = log(SPR(exp(logFlim)));
    }else if(CppAD::Variable(par.logScaleFcrash) &&
			     !(conf.stockRecruitmentModelCode == 61 ||
			       conf.stockRecruitmentModelCode == 63)){
      logFlim = par.logScaleFcrash + logFsq * Type(fromSQ) + log(0.95); //logFsq + par.logScaleFcrash;
      logBlim = log(Se(exp(logFlim)));
      logRcrash = log(Re(exp(logFlim)));
      logYcrash = log(yield(exp(logFlim)));
      logYPRcrash = log(YPR(exp(logFlim)));
      logSPRcrash = log(SPR(exp(logFlim)));
    }else if(CppAD::Variable(par.logScaleFext) &&
	     (conf.stockRecruitmentModelCode == 61 ||
	      conf.stockRecruitmentModelCode == 63)){
      logFlim = par.logScaleFext + logFsq * Type(fromSQ) + log(0.95); //logFsq + par.logScaleFcrash;
      logBlim = log(Se(exp(logFlim)));
      logRcrash = log(Re(exp(logFlim)));
      logYcrash = log(yield(exp(logFlim)));
      logYPRcrash = log(YPR(exp(logFlim)));
      logSPRcrash = log(SPR(exp(logFlim)));
    }else{
      logFlim = R_NaReal;//R_NaReal;
      if(conf.stockRecruitmentModelCode == 61){
	logBlim = par.rec_pars(1);
      }else if(conf.stockRecruitmentModelCode == 63){
	logBlim = par.rec_pars(0);      
      }else{
	logBlim = R_NaReal;
      }
      logRlim = R_NaReal;
      logYlim = R_NaReal;
      logYPRlim = R_NaReal;
      logSPRlim = R_NaReal;
    }
 
  
#ifdef CPPAD_FRAMEWORK
    // Prepare AD
    vector<Type> Fsqvec(1);
    Fsqvec(0) = exp(-10.0);
    CppAD::vector<AD<Type> > x1( Fsqvec );
    CppAD::vector<AD<Type> > y1( 1 );
    CppAD::Independent(x1);
    y1[0] = YPR(x1);
    FYPR = CppAD::ADFun<Type>(x1, y1);
    
    CppAD::vector<AD<Type> > x2( Fsqvec );
    CppAD::vector<AD<Type> > y2( 1 );
    CppAD::Independent(x2);
    y2[0] = SPR(x2);
    FSPR = CppAD::ADFun<Type>(x2, y2);

    vector<Type> Bsqvec(1);
    Bsqvec(0) = exp(12.0);
    CppAD::vector<AD<Type> > x3( Bsqvec );
    CppAD::vector<AD<Type> > y3( 1 );
    CppAD::Independent(x3);
    y3[0] = SR(x3);
    FSR = CppAD::ADFun<Type>(x3, y3);

#endif
#ifdef TMBAD_FRAMEWORK
    // Prepare AD
    vector<Type> x0(1);
    x0(0) = exp(-10.0);
    // YPR
    F_FYPR<Type> FyprInst = {*this};
    FYPR = TMBad::ADFun<>(TMBad::StdWrap<F_FYPR<Type>,vector<TMBad::ad_aug> >(FyprInst), x0);
    FYPR = FYPR.JacFun();

    // SPR
    F_FSPR<Type> FsprInst = {*this};
    FSPR = TMBad::ADFun<>(TMBad::StdWrap<F_FSPR<Type>,vector<TMBad::ad_aug> >(FsprInst), x0);
    FSPR = FSPR.JacFun();

    // FSR
    F_FSR<Type> FsrInst = {*this};
    FSR = TMBad::ADFun<>(TMBad::StdWrap<F_FSR<Type>,vector<TMBad::ad_aug> >(FsrInst), x0);
    FSR = FSR.JacFun();

#endif

    
  }

  template<class T>
  T logDiscYield(T logFbar){
    PERREC_t<T> r = perRecruit<Type, T>(logFbar, dat, conf, par, logSel, aveYears, nYears,
					//dat.referencepoint.RecCorrection,
					dat.referencepoint.catchType);
    return r.logDiscYe;
  }
  
  template<class T>
  T logYPR(T logFbar){
    PERREC_t<T> r = perRecruit<Type, T>(logFbar, dat, conf, par, logSel, aveYears, nYears,
					//dat.referencepoint.RecCorrection,
					dat.referencepoint.catchType);
    return r.logYPR;
  }

  template<class T>
  T YPR(T Fbar){
    return exp(logYPR(log(Fbar)));
  }

  AD<Type> YPR(CppAD::vector<AD<Type> > Fbar){
    PERREC_t<AD<Type> > r = perRecruit<Type, AD<Type> >(log(Fbar[0]), dat, conf, par, logSel, aveYears, nYears,
							//dat.referencepoint.RecCorrection,
							dat.referencepoint.catchType);
    return exp(r.logYPR);
  }
  
  Type dYPR(Type Fbar){
#ifdef CPPAD_FRAMEWORK
    vector<Type> Fv(1);
    Fv(0) = Fbar;
    CppAD::vector<Type> x_eval( Fv );
    return FYPR.Jacobian(x_eval)[0];
#endif
#ifdef TMBAD_FRAMEWORK
    vector<Type> xx(1);
    xx(0) = Fbar;
    vector<Type> tmp = FYPR(xx);
    return tmp(0);
#endif
     }

  template<class T>
  T SPR(T Fbar){
    PERREC_t<T> r = perRecruit<Type, T>(log(Fbar), dat, conf, par, logSel, aveYears, nYears,
					//dat.referencepoint.RecCorrection,
					dat.referencepoint.catchType);
    return exp(r.logSPR);
  }

  AD<Type> SPR(CppAD::vector<AD<Type> > Fbar){
    PERREC_t<AD<Type> > r = perRecruit<Type, AD<Type> >(log(Fbar[0]), dat, conf, par, logSel, aveYears, nYears,
							//dat.referencepoint.RecCorrection,
							dat.referencepoint.catchType);
    return exp(r.logSPR);
  }

  Type dSPR(Type Fbar){
#ifdef CPPAD_FRAMEWORK
    vector<Type> Fv(1);
    Fv(0) = Fbar;
    CppAD::vector<Type> x_eval( Fv );
    return FSPR.Jacobian(x_eval)[0];
#endif
#ifdef TMBAD_FRAMEWORK
    vector<Type> xx(1);
    xx(0) = Fbar;
    vector<Type> tmp = FSPR(xx);
    return tmp(0);
#endif
  }

  Type Se(Type Fbar){
    PERREC_t<Type> r = perRecruit<Type, Type>(log(Fbar), dat, conf, par, logSel, aveYears, nYears,
					      //dat.referencepoint.RecCorrection,
					      dat.referencepoint.catchType);
    return exp(r.logSe);
  }

  Type Re(Type Fbar){
    PERREC_t<Type> r = perRecruit<Type, Type>(log(Fbar), dat, conf, par, logSel, aveYears, nYears,
					      //dat.referencepoint.RecCorrection,
					      dat.referencepoint.catchType);
    return exp(r.logRe);
  }

  Type yield(Type Fbar){
    PERREC_t<Type> r = perRecruit<Type, Type>(log(Fbar), dat, conf, par, logSel, aveYears, nYears,
					      //dat.referencepoint.RecCorrection,
					      dat.referencepoint.catchType);
    return exp(r.logYe);
  }
  
  Type yearsLost(Type Fbar){
    PERREC_t<Type> r = perRecruit<Type, Type>(log(Fbar), dat, conf, par, logSel, aveYears, nYears,
					      //dat.referencepoint.RecCorrection,
					      dat.referencepoint.catchType);
    return exp(r.logYearsLost);
  }

  Type logYearsLost(Type Fbar){
    PERREC_t<Type> r = perRecruit<Type, Type>(log(Fbar), dat, conf, par, logSel, aveYears, nYears,
					      //dat.referencepoint.RecCorrection,
					      dat.referencepoint.catchType);
    return r.logYearsLost;
  }

  Type lifeexpectancy(Type Fbar){
    PERREC_t<Type> r = perRecruit<Type, Type>(log(Fbar), dat, conf, par, logSel, aveYears, nYears,
					      //dat.referencepoint.RecCorrection,
					      dat.referencepoint.catchType);
    return exp(r.logLifeExpectancy);
  }


  template<class T>
  T SR(T ssb){
    vector<T> rp2 = par.rec_pars.template cast<T>();
    if(conf.stockRecruitmentModelCode == 0){
      return 0.0;
    }else if(conf.stockRecruitmentModelCode == 3){
      return exp(par.rec_pars(par.rec_pars.size() - 1));
    }else if(conf.stockRecruitmentModelCode == 62){
      return exp(par.rec_pars(0));
    }else if(conf.stockRecruitmentModelCode == 90){
      return exp(log(ssb) + ibcdspline(log(ssb),
				       (vector<T>)(conf.constRecBreaks.template cast<T>()),
				       rp2));
    }else if(conf.stockRecruitmentModelCode == 91){
      return exp(log(ssb) + ibcspline(log(ssb),
				      (vector<T>)(conf.constRecBreaks.template cast<T>()),
				      rp2));
    }else if(conf.stockRecruitmentModelCode == 92){
      return exp(log(ssb) + bcspline(log(ssb),
				     (vector<T>)(conf.constRecBreaks.template cast<T>()),
				     rp2));
    }else{
      return exp(functionalStockRecruitment(ssb, rp2, conf.stockRecruitmentModelCode));
    }
  }

  AD<Type> SR(CppAD::vector<AD<Type> > ssb){
    AD<Type> s0 = ssb[0];
    vector<AD<Type> > rp(par.rec_pars.size());
    rp = par.rec_pars.template cast<AD<Type> >();
    if(conf.stockRecruitmentModelCode == 0){
      return 0.0;
    }else if(conf.stockRecruitmentModelCode == 3){
      return exp(rp(rp.size() - 1));
    }else if(conf.stockRecruitmentModelCode == 62){
      return exp(rp(0));
    }else if(conf.stockRecruitmentModelCode == 90){
      return exp(log(s0) + ibcdspline(log(s0),
				      (vector<AD<Type> >)(conf.constRecBreaks.template cast<AD<Type> >()),
				      rp));
    }else if(conf.stockRecruitmentModelCode == 91){
      return exp(log(s0) + ibcspline(log(s0),
				  (vector<AD<Type> >)(conf.constRecBreaks.template cast<AD<Type> >()),
				     rp));
    }else if(conf.stockRecruitmentModelCode == 92){
      return exp(log(s0) + bcspline(log(s0),
				 (vector<AD<Type> >)(conf.constRecBreaks.template cast<AD<Type> >()),
				    rp));
    }else{
      return exp(functionalStockRecruitment(s0, rp, conf.stockRecruitmentModelCode));
    }
  }


  Type dSR(Type ssb){
#ifdef CPPAD_FRAMEWORK
    vector<Type> Fv(1);
    Fv(0) = ssb;
    CppAD::vector<Type> x_eval( Fv );
    return FSR.Jacobian(x_eval)[0];
#endif
#ifdef TMBAD_FRAMEWORK
    vector<Type> xx(1);
    xx(0) = ssb;
    vector<Type> tmp = FSR(xx);
    return tmp(0);
#endif
  }

  
  
  // Calculate "likelihood" contribution to estimate reference points
  Type operator()() {
    Type nll = 0.0;
    if(CppAD::Variable(par.logScaleFmsy)){
      Type tmpYmsy = yield(exp(logFmsy));
      nll -= log(tmpYmsy);
      for(int i = 0; i < par.logScaleFmsyRange.cols(); ++i){
	// Lower end of range
	if(CppAD::Variable(par.logScaleFmsyRange(0,i))){
	  Type tmp = log((yield(exp(logFmsyRange(0,i))))) - (log(dat.referencepoint.MSYRange(i)) + log(tmpYmsy));
	  nll += tmp * tmp;
	}
	      
	// Upper end of range
	if(CppAD::Variable(par.logScaleFmsyRange(1,i))){
	  Type tmp = log((yield(exp(logFmsyRange(1,i))))) - (log(dat.referencepoint.MSYRange(i)) + log(tmpYmsy));
	  nll += tmp * tmp;	  
	}
		
      }
    }

    if(CppAD::Variable(par.logScaleFmypyl)){
      Type tmp = log(yield(exp(logFmypyl))) - logspace_add((Type)0.0,logYearsLost(exp(logFmypyl)) - log((Type)conf.maxAge - (Type)conf.minAge + (Type)1.0));
      nll -= log(tmp);
    }

    if(CppAD::Variable(par.logScaleFmdy)){
      nll -= logDiscYield(logFmdy);
    }

    
    if(CppAD::Variable(par.logScaleFmax)){
      nll -= logYPR(logFmax);
    }
 
    for(int i = 0; i < par.logScaleFxdYPR.size(); ++i){
      if(CppAD::Variable(par.logScaleFxdYPR(i))){
	Type tmp = dat.referencepoint.xdYPR(i) * dYPR(Type(SAM_Zero)) - dYPR(exp(logFxdYPR(i)));
	nll += tmp * tmp;
      }
    }
 
    if(CppAD::Variable(par.logScaleFcrash)){
      Type tmp = dSR(Type(SAM_Zero)) - (1.0 / SPR(exp(logFcrash)));
      nll += tmp * tmp;
    }
 
    if(CppAD::Variable(par.logScaleFext)){
      Type tmp1 = log(Se(exp(logFext)));
      Type tmp2 = exp(logFext);
      nll += tmp1 * tmp1 + 10.0 * tmp2;
      // Type tmp = Se(exp(logFext)) - Se(exp(logFext - 0.01));
      // nll += tmp;
    }
    
    for(int i = 0; i < par.logScaleFxPercent.size(); ++i){
      if(CppAD::Variable(par.logScaleFxPercent(i))){
	Type tmp = dat.referencepoint.xPercent(i) * SPR(Type(SAM_Zero)) - SPR(exp(logFxPercent(i)));
	nll += tmp * tmp;
      }
    }
 

    for(int i = 0; i < par.logScaleFxB0.size(); ++i){
      if(CppAD::Variable(par.logScaleFxB0(i))){
	Type tmp = log(dat.referencepoint.xB0(i)) + logB0 - log(Se(exp(logFxB0(i))));
	nll += tmp * tmp;
      }
    }
 
    if(!isNA(logBlim) && CppAD::Variable(par.logScaleFlim)){
      Type tmp = logBlim - log(Se(exp(logFlim)));
      nll += tmp * tmp;
    }
 
    return par.implicitFunctionDelta * nll;
  }

};

template<class Type>
Type nllReferencepoints(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, array<Type> &logF, objective_function<Type> *of){

  if(dat.referencepoint.nYears == 0)
    return 0.0;

  REFERENCE_POINTS<Type> referencepoint(dat, conf, par, logN, logF, dat.referencepoint.fromSQ);



  if(dat.referencepoint.Fsequence.size() > 0){
    vector<Type> Fseq = dat.referencepoint.Fsequence;
    vector<Type> refpointseq_logYPR(Fseq.size());
    refpointseq_logYPR.setZero();
    vector<Type> refpointseq_logSPR(Fseq.size());
    refpointseq_logSPR.setZero();
    vector<Type> refpointseq_logSe(Fseq.size());
    refpointseq_logSe.setZero();
    vector<Type> refpointseq_logYe(Fseq.size());
    refpointseq_logYe.setZero();
    vector<Type> refpointseq_logDiscYe(Fseq.size());
    refpointseq_logDiscYe.setZero();
    vector<Type> refpointseq_logDiscYPR(Fseq.size());
    refpointseq_logDiscYPR.setZero();
    vector<Type> refpointseq_logRe(Fseq.size());
    refpointseq_logRe.setZero();
    Type refpointseq_dSR0 = 0.0;

    vector<Type> refpointseq_logYearsLost(Fseq.size());
    refpointseq_logYearsLost.setZero();
    vector<Type> refpointseq_logLifeExpectancy(Fseq.size());
    refpointseq_logLifeExpectancy.setZero();

    
    for(int i = 0; i < Fseq.size(); ++i){
      PERREC_t<Type> v = perRecruit<Type, Type>(log(Fseq(i)),
						referencepoint.dat,
						referencepoint.conf,
						referencepoint.par,
						referencepoint.logSel,
						referencepoint.aveYears,
						referencepoint.nYears,
						//dat.referencepoint.RecCorrection,
						dat.referencepoint.catchType);
      refpointseq_logYPR(i) = v.logYPR;
      refpointseq_logSPR(i) = v.logSPR;
      refpointseq_logSe(i) = v.logSe;
      refpointseq_logYe(i) = v.logYe;
      refpointseq_logRe(i) = v.logRe;
      refpointseq_logDiscYe(i) = v.logDiscYe;
      refpointseq_logDiscYPR(i) = v.logDiscYPR;
      refpointseq_logYearsLost(i) = v.logYearsLost;
      refpointseq_logLifeExpectancy(i) = v.logLifeExpectancy;
      if(i == 0)
	refpointseq_dSR0 = v.dSR0;
	
    }

    REPORT_F(refpointseq_logYPR, of);
    REPORT_F(refpointseq_logSPR, of);
    REPORT_F(refpointseq_logSe, of);
    REPORT_F(refpointseq_logYe, of);
    REPORT_F(refpointseq_logRe, of);
    REPORT_F(refpointseq_dSR0, of);
    REPORT_F(refpointseq_logDiscYe, of);
    REPORT_F(refpointseq_logDiscYPR, of);
    REPORT_F(refpointseq_logYearsLost, of);
    REPORT_F(refpointseq_logLifeExpectancy, of);
    ADREPORT_F(refpointseq_logYPR, of);
    ADREPORT_F(refpointseq_logSPR, of);
    ADREPORT_F(refpointseq_logSe, of);
    ADREPORT_F(refpointseq_logYe, of);
    ADREPORT_F(refpointseq_logRe, of);
    ADREPORT_F(refpointseq_dSR0, of);
    ADREPORT_F(refpointseq_logYearsLost, of);
    ADREPORT_F(refpointseq_logLifeExpectancy, of);
  }

  
  ADREPORT_F(referencepoint.logFsq,of);
  ADREPORT_F(referencepoint.logBsq,of);
  ADREPORT_F(referencepoint.logRsq,of);
  ADREPORT_F(referencepoint.logYsq,of);
  ADREPORT_F(referencepoint.logYPRsq,of);
  ADREPORT_F(referencepoint.logSPRsq,of);

  ADREPORT_F(referencepoint.logF0,of);
  ADREPORT_F(referencepoint.logB0,of);
  ADREPORT_F(referencepoint.logR0,of);
  ADREPORT_F(referencepoint.logY0,of);
  ADREPORT_F(referencepoint.logYPR0,of);
  ADREPORT_F(referencepoint.logSPR0,of);

  ADREPORT_F(referencepoint.logFmsy,of);
  ADREPORT_F(referencepoint.logBmsy,of);
  ADREPORT_F(referencepoint.logRmsy,of);
  ADREPORT_F(referencepoint.logYmsy,of);
  ADREPORT_F(referencepoint.logYPRmsy,of);
  ADREPORT_F(referencepoint.logSPRmsy,of);

  
  ADREPORT_F(referencepoint.logFmsyRange,of);
  ADREPORT_F(referencepoint.logBmsyRange,of);
  ADREPORT_F(referencepoint.logRmsyRange,of);
  ADREPORT_F(referencepoint.logYmsyRange,of);
  ADREPORT_F(referencepoint.logYPRmsyRange,of);
  ADREPORT_F(referencepoint.logSPRmsyRange,of);

  ADREPORT_F(referencepoint.logFmax,of);
  ADREPORT_F(referencepoint.logBmax,of);
  ADREPORT_F(referencepoint.logRmax,of);
  ADREPORT_F(referencepoint.logYmax,of);
  ADREPORT_F(referencepoint.logYPRmax,of);
  ADREPORT_F(referencepoint.logSPRmax,of);
 
  ADREPORT_F(referencepoint.logFxdYPR,of);
  ADREPORT_F(referencepoint.logBxdYPR,of);
  ADREPORT_F(referencepoint.logRxdYPR,of);
  ADREPORT_F(referencepoint.logYxdYPR,of);
  ADREPORT_F(referencepoint.logYPRxdYPR,of);
  ADREPORT_F(referencepoint.logSPRxdYPR,of);

  ADREPORT_F(referencepoint.logFcrash,of);
  ADREPORT_F(referencepoint.logBcrash,of);
  ADREPORT_F(referencepoint.logRcrash,of);
  ADREPORT_F(referencepoint.logYcrash,of);
  ADREPORT_F(referencepoint.logYPRcrash,of);
  ADREPORT_F(referencepoint.logSPRcrash,of);

  ADREPORT_F(referencepoint.logFext,of);
  ADREPORT_F(referencepoint.logBext,of);
  ADREPORT_F(referencepoint.logRext,of);
  ADREPORT_F(referencepoint.logYext,of);
  ADREPORT_F(referencepoint.logYPRext,of);
  ADREPORT_F(referencepoint.logSPRext,of);

  
  ADREPORT_F(referencepoint.logFxPercent,of);
  ADREPORT_F(referencepoint.logBxPercent,of);
  ADREPORT_F(referencepoint.logRxPercent,of);
  ADREPORT_F(referencepoint.logYxPercent,of);
  ADREPORT_F(referencepoint.logYPRxPercent,of);
  ADREPORT_F(referencepoint.logSPRxPercent,of);

  ADREPORT_F(referencepoint.logFxB0,of);
  ADREPORT_F(referencepoint.logBxB0,of);
  ADREPORT_F(referencepoint.logRxB0,of);
  ADREPORT_F(referencepoint.logYxB0,of);
  ADREPORT_F(referencepoint.logYPRxB0,of);
  ADREPORT_F(referencepoint.logSPRxB0,of);

  
  ADREPORT_F(referencepoint.logFlim,of);
  ADREPORT_F(referencepoint.logBlim,of);
  ADREPORT_F(referencepoint.logRlim,of);
  ADREPORT_F(referencepoint.logYlim,of);
  ADREPORT_F(referencepoint.logYPRlim,of);
  ADREPORT_F(referencepoint.logSPRlim,of);

  ADREPORT_F(referencepoint.logFmypyl,of);
  ADREPORT_F(referencepoint.logBmypyl,of);
  ADREPORT_F(referencepoint.logRmypyl,of);
  ADREPORT_F(referencepoint.logYmypyl,of);
  ADREPORT_F(referencepoint.logYPRmypyl,of);
  ADREPORT_F(referencepoint.logSPRmypyl,of);

  ADREPORT_F(referencepoint.logFmdy,of);
  ADREPORT_F(referencepoint.logBmdy,of);
  ADREPORT_F(referencepoint.logRmdy,of);
  ADREPORT_F(referencepoint.logYmdy,of);
  ADREPORT_F(referencepoint.logYPRmdy,of);
  ADREPORT_F(referencepoint.logSPRmdy,of);


  // Fbar relative to (last year) reference points
  // vector<Type> logfbar = fbarFun(conf, logF, true);
  // vector<Type> fbar = exp(logfbar);
  vector<Type> fbar = fbarFun(conf, logF);
  vector<Type> logfbar = log(fbar);

  vector<Type> relref_logfbar_fmsy = logfbar - referencepoint.logFmsy;
  ADREPORT_F(relref_logfbar_fmsy,of);
  // vector<Type> relref_logfbar_fmax = logfbar - referencepoint.logFmax;
  // ADREPORT_F(relref_logfbar_fmax,of);
  // vector<Type> relref_logfbar_f01 = logfbar - referencepoint.logF01;
  // ADREPORT_F(relref_logfbar_f01,of);
  // vector<Type> relref_logfbar_f35 = logfbar - referencepoint.logF35;
  // ADREPORT_F(relref_logfbar_f35,of);

  // SSB relative to (last year) reference points
  // vector<Type> logssb = ssbFun(dat, conf, logN, logF, true);
  // vector<Type> ssb = exp(logssb);
  vector<Type> ssb = ssbFun(dat, conf, logN, logF);
  vector<Type> logssb = log(ssb);

  vector<Type> relref_logssb_bmsy = logssb - referencepoint.logBmsy;
  ADREPORT_F(relref_logssb_bmsy,of);
  // vector<Type> relref_logssb_bmax = logssb - referencepoint.logBmax;
  // ADREPORT_F(relref_logssb_bmax,of);
  // vector<Type> relref_logssb_b01 = logssb - referencepoint.logB01;
  // ADREPORT_F(relref_logssb_b01,of);
  // vector<Type> relref_logssb_b35 = logssb - referencepoint.logB35;
  // ADREPORT_F(relref_logssb_b35,of);

  Type ans = referencepoint();
  return ans;
}




// R functions
#ifdef TMBAD_FRAMEWORK
struct F_dFunctionalSR2 {
  long int nrp;
  int srmc;
  template<class T>
  T operator()(vector<T> x){  // Evaluate function
    vector<T> rp2(nrp);
    for(int i = 0; i < nrp; ++i)
      rp2(i) = x(i);
    T ssb = x(nrp);
    return exp(functionalStockRecruitment(ssb, rp2, srmc));
  }
};
#endif
  

#endif
