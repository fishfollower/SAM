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



#define SAM_NegInf -20.0
#define SAM_NIZero -10.0
#define SAM_Zero exp(SAM_NIZero)
#define SAM_RECNTAB 20
// exp(SAM_NegInf)


namespace rec_atomic {

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
      Rf_error("lambertW is only implemented for x>0");
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
      return SAM_Zero;
    }
    return SAM_Zero;
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
      return SAM_Zero;
    }
    return SAM_Zero;
  }

  TMB_BIND_ATOMIC(Se_sl0, 1111, Se_sl_raw(x[0],x[1],x[2],x[3]))


  
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


template<class Type>
struct PERREC_t {
  Type logFbar;
  Type logYPR;
  Type logSPR;
  Type logSe;
  Type logRe;
  Type logYe;
  Type dSR0;
};


template<class Type, class T>
PERREC_t<T> perRecruit(T Fbar, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, vector<Type>& sel, vector<int> aveYears, int nYears = 300){
  
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

  vector<T> selT = sel.template cast<T>();
  // Make logF array
  array<T> logF(selT.size(), nYears);
  logF.setZero();
  for(int i = 0; i < nYears; ++i)
    logF.col(i) += log(selT) + log(Fbar);

  // Make logN array - start with one recruit
  int nAge = conf.maxAge - conf.minAge + 1;
  array<T> logN(nAge, nYears);
  logN.setZero();
  logN += SAM_NegInf;
  logN(0,0) = 0.0;

  // Run loop over years
  for(int i = 1; i < nYears; ++i){
    // predN
    logN.col(i) = predNFun(newDat, newConf, newPar, logN, logF, i);
    // remove recruitment
    logN(0,i) = SAM_NegInf;
  }

 
  // Calculate yield
  vector<T> cat(nYears);
  switch(newDat.referencepoint.catchType){
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
  T logYPR = log(sum(cat));
  // Calculate spawners
  vector<T> ssb = ssbFun(newDat, newConf, logN, logF);
  T logSPR = log(sum(ssb));
  T lambda = sum(ssb); //exp(logSPR);

  if(conf.stockRecruitmentModelCode == 0 ||
     conf.stockRecruitmentModelCode == 3){
    PERREC_t<T> res = {log(Fbar), // logFbar
		       logYPR,	   // logYPR
		       logSPR,	   // logSPR
		       R_NaReal,	   // logSe
		       R_NaReal,		 // logRe
		       R_NaReal,// logYe
		       R_NaReal}; // dSR0
    return res;
  }
  



  // Calculate Se
  T Se = 0.0; //R_NegInf;
  
  T dsr0 = 10000.0;
  if(conf.stockRecruitmentModelCode != 62 &&
     conf.stockRecruitmentModelCode != 65 &&
     conf.stockRecruitmentModelCode != 68 &&
     conf.stockRecruitmentModelCode != 69)
    dsr0 = dFunctionalSR(T(SAM_Zero), newPar.rec_pars, conf.stockRecruitmentModelCode);
  if(conf.stockRecruitmentModelCode == 68 || conf.stockRecruitmentModelCode == 69)
    dsr0 = CppAD::CondExpLt((T)newPar.rec_pars[2], (T)0, (T)dFunctionalSR(T(SAM_Zero), newPar.rec_pars, conf.stockRecruitmentModelCode), (T)dsr0);

  switch(conf.stockRecruitmentModelCode){
  case 0: // straight RW 
    Rf_error("Equilibrium SSB not implemented");
    break;
  case 1: //ricker
    Se = exp(-newPar.rec_pars(1)) * (newPar.rec_pars(0) + log(lambda)); //log((exp(newPar.rec_pars(0)) * lambda));
    break;
  case 2:  //BH
    Se =  CppAD::CondExpGt((newPar.rec_pars(0) + logSPR), T(SAM_Zero),
			   fabs(exp(newPar.rec_pars(0)) * lambda - 1.0) * exp(-newPar.rec_pars(1)),
			   T(SAM_Zero));
    break;
  case 3: //Constant mean
    Rf_error("Equilibrium SSB not implemented");
    break;
  case 61: // Hockey stick
    Se = lambda * exp(newPar.rec_pars(0));
    break;
  case 62: // AR1 (on log-scale)
    Se = lambda * exp(newPar.rec_pars(0));
    break;
  case 63: //Bent hyperbola / Hockey-stick-like
    Se = (2.0 * sqrt(exp(2.0 * newPar.rec_pars(0)) + exp(2.0 * newPar.rec_pars(2)) / 4.0) / (lambda * exp(newPar.rec_pars(1))) - 2.0 * exp(newPar.rec_pars(0)) - 2.0 * sqrt(exp(2.0 * newPar.rec_pars(0)) + exp(2.0 * newPar.rec_pars(2)) / 4.0)) / ( 1.0 / ((lambda * lambda * exp(2.0 * newPar.rec_pars(1)))) - 2.0 / (lambda * exp(newPar.rec_pars(1)))  );  
    break;
  case 64: // Power CMP
    Se = exp(1.0 / (1.0 - invlogit(newPar.rec_pars(1))) * (newPar.rec_pars(0) + log(lambda)));
    break;
  case 65: // Power Non-CMP
    Se = exp(1.0 / (1.0 - (exp(newPar.rec_pars(1)) + 1.0001)) * (newPar.rec_pars(0) + log(lambda)));
    break;
  case 66: // Shepherd
    Se = CppAD::CondExpGt((newPar.rec_pars(0) + logSPR), T(SAM_Zero),
    			  exp( newPar.rec_pars(1) + 1.0 / exp(newPar.rec_pars(2)) * log(fabs(exp(newPar.rec_pars(0)) * lambda - 1.0))),
    			  T(SAM_Zero));
    break;
  case 67: // Deriso
    Se = CppAD::CondExpGt((newPar.rec_pars(0) + logSPR), T(SAM_Zero),
			  fabs(exp(exp(-newPar.rec_pars(2)) * (newPar.rec_pars(0) + logSPR)) - 1.0) * exp(-newPar.rec_pars(1)),
			   T(SAM_Zero));
    break;
  case 68: // Saila-Lorda (cases: gamma > 1; gamma = 1; gamma < 1
    Se = Se_sl(lambda, exp(newPar.rec_pars(0)), exp(newPar.rec_pars(1)), exp(newPar.rec_pars(2)));
    break;
  case 69: // Sigmoidal Beverton-Holt
    Se = Se_sbh(lambda, exp(newPar.rec_pars(0)), exp(newPar.rec_pars(1)), exp(newPar.rec_pars(2)));
    break;
    default:
      error("SR model code not recognized");
    break;   
  }

  //T logSe = log(Se);
  //Type logSe = CppAD::CondExpGt(exp(-logSPR), dsr0, Type(SAM_NegInf), log(Se));
  T logSe = CppAD::CondExpGt(exp(-logSPR), dsr0,
  			     log(fabs(Se)) - 10.0 * (exp(-logSPR) - dsr0),
  			     log(fabs(Se)));

  T logYe = CppAD::CondExpGt(exp(-logSPR), dsr0,
  			     logSe - logSPR + logYPR - 10.0 * (exp(-logSPR) - dsr0),
  			     logSe - logSPR + logYPR);

  // Return
  PERREC_t<T> res = {log(Fbar), // logFbar
		     logYPR,	   // logYPR
		     logSPR,	   // logSPR
		     logSe,	   // logSe
		     logSe - logSPR,		 // logRe
		     logYe,   // logYe
		     dsr0};			// DSR0

  return res;
}


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
  Type logFmax;			// Maximizes yield per recruit
  Type logF01;			// F such that YPR'(0) = 0.1 * YPR'(F)
  Type logFcrash;		// F such that 1/SPR(f) = SR'(0) (i.e. stock crashes [with compensatory recruitment] if slope of spawner-per-recruit in origin is less than slope of stock-recruitment model in origin)
  Type logFext;			// F such that stock dies out - SSB(F) < epsilon and F smallest possible
  vector<Type> logFxPercent;			// F such that SSB is reduced to x% of unfished stock (Se(F) = x/100 * Se(0) )
  //Type logFmed;			// Fishing  mortality  rate  F  corresponding  to  a  SSB/R  equal  to  the  inverse  of  the  50th  percentile of the observed R/SSB
  Type logFlim;			// F such that Se(F) = Blim (for hockey-stick-like stock recruitment only)
  //Type logFpa;			// F that corresponds to logBpa
   

  // Corresponding SSB??
  Type logBsq;
  Type logB0;
  Type logBmsy;
  Type logBmax;
  Type logB01;
  Type logBcrash;
  Type logBext;
  vector<Type> logBxPercent;
  Type logBlim;			// Known from model parameters (for hockey-stick-like stock recruitment only)
  //Type logBpa;			// Ba = Blim * exp(1.645 * sigma) where sigma is the standard deviation of log(SSB) at the start of the year following the terminal year of the assessment if sigma is unknown, 0.2 can be used as default.

  // Corresponding Recruitment
  Type logRsq;
  Type logR0;
  Type logRmsy;
  Type logRmax;
  Type logR01;
  Type logRcrash;
  Type logRext;
  vector<Type> logRxPercent;
  Type logRlim;			// Known from model parameters (for hockey-stick-like stock recruitment only)
  //Type logBpa;			// Ba = Blim * exp(1.645 * sigma) where sigma is the standard deviation of log(SSB) at the start of the year following the terminal year of the assessment if sigma is unknown, 0.2 can be used as default.

  
  // Corresponding Yield??
  Type logYsq;
  Type logY0;
  Type logYmsy;
  Type logYmax;
  Type logY01;
  Type logYcrash;
  Type logYext;
  vector<Type> logYxPercent;
  Type logYlim;			// Known from model parameters (for hockey-stick-like stock recruit

  Type logYPRsq;
  Type logYPR0;
  Type logYPRmsy;
  Type logYPRmax;
  Type logYPR01;
  Type logYPRcrash;
  Type logYPRext;
  vector<Type> logYPRxPercent;
  Type logYPRlim;

  Type logSPRsq;
  Type logSPR0;
  Type logSPRmsy;
  Type logSPRmax;
  Type logSPR01;
  Type logSPRcrash;
  Type logSPRext;
  vector<Type> logSPRxPercent;
  Type logSPRlim;

  
  // Derived values
  vector<Type> sel;

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
		   array<Type>& logF_	        
		   ): dat(dat_), conf(conf_), par(par_),
		      logN(logN_),
		      logF(logF_)
  {
    nYears = dat.referencepoint.nYears;
    aveYears = dat.referencepoint.aveYears;
    selYears = dat.referencepoint.selYears;

    // Calculate current selectivity and status quo
    vector<Type> fbartmp(selYears.size());
    fbartmp.setZero();
    sel = vector<Type>(logF.rows());
    sel.setZero();

    for(int y = 0; y < selYears.size(); ++y){
      sel += exp(logF.col(selYears(y)));
      for(int a = conf.fbarRange(0); a <= conf.fbarRange(1); a++){  
	fbartmp(y) += exp(logF(conf.keyLogFsta(0,a-conf.minAge),selYears(y)));
      }
      fbartmp(y) /= Type(conf.fbarRange(1)-conf.fbarRange(0)+1);
    }

    logFsq = log(sum(fbartmp)) - log(fbartmp.size());
    sel *= 1.0 / sum(fbartmp);
    
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
    if(CppAD::Variable(par.logScaleFmsy)){
      logFmsy = par.logScaleFmsy; // logFsq + par.logScaleFmsy;
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

    
    if(CppAD::Variable(par.logScaleFmax)){
      logFmax = par.logScaleFmax; //logFsq + par.logScaleFmax;
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


    if(CppAD::Variable(par.logScaleF01)){
      logF01 = par.logScaleF01; //logFsq + par.logScaleF01;
      logB01 = log(Se(exp(logF01)));
      logR01 = log(Re(exp(logF01)));
      logY01 = log(yield(exp(logF01)));
      logYPR01 = log(YPR(exp(logF01)));
      logSPR01 = log(SPR(exp(logF01)));
    }else{
      logF01 = R_NaReal;//R_NaReal;
      logB01 = R_NaReal;
      logR01 = R_NaReal;
      logY01 = R_NaReal;
      logYPR01 = R_NaReal;
      logSPR01 = R_NaReal;
    }


    if(CppAD::Variable(par.logScaleFcrash)){
      logFcrash = par.logScaleFcrash; //logFsq + par.logScaleFcrash;
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
      logFext = par.logScaleFext; //logFsq + par.logScaleFcrash;
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
	logFxPercent(i) = par.logScaleFxPercent(i); //logFsq + par.logScaleF35;
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
    
    if(CppAD::Variable(par.logScaleFlim) &&
       (conf.stockRecruitmentModelCode == 61 ||
	conf.stockRecruitmentModelCode == 63)){    
      logFlim = par.logScaleFlim; //logFsq + par.logScaleFlim;
      logYlim = log(yield(exp(logFlim)));
      if(conf.stockRecruitmentModelCode == 61){
	logBlim = par.rec_pars(1);
      }else if(conf.stockRecruitmentModelCode == 63){
	logBlim = par.rec_pars(0);      
      }
      logRlim = log(Re(exp(logFlim)));
      logYPRlim = log(YPR(exp(logFlim)));
      logSPRlim = log(SPR(exp(logFlim)));
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
  T YPR(T Fbar){
    PERREC_t<T> r = perRecruit<Type, T>(Fbar, dat, conf, par, sel, aveYears, nYears);
    return exp(r.logYPR);
  }

  AD<Type> YPR(CppAD::vector<AD<Type> > Fbar){
    PERREC_t<AD<Type> > r = perRecruit<Type, AD<Type> >(Fbar[0], dat, conf, par, sel, aveYears, nYears);
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
    PERREC_t<T> r = perRecruit<Type, T>(Fbar, dat, conf, par, sel, aveYears, nYears);
    return exp(r.logSPR);
  }

  AD<Type> SPR(CppAD::vector<AD<Type> > Fbar){
    PERREC_t<AD<Type> > r = perRecruit<Type, AD<Type> >(Fbar[0], dat, conf, par, sel, aveYears, nYears);
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
    PERREC_t<Type> r = perRecruit<Type, Type>(Fbar, dat, conf, par, sel, aveYears, nYears);
    return exp(r.logSe);
  }

  Type Re(Type Fbar){
    PERREC_t<Type> r = perRecruit<Type, Type>(Fbar, dat, conf, par, sel, aveYears, nYears);
    return exp(r.logRe);
  }

  Type yield(Type Fbar){
    PERREC_t<Type> r = perRecruit<Type, Type>(Fbar, dat, conf, par, sel, aveYears, nYears);
    return exp(r.logYe);
  }

    template<class T>
    T SR(T ssb){
    if(conf.stockRecruitmentModelCode == 0)
      return 0.0;
    if(conf.stockRecruitmentModelCode == 62)
      return exp(par.rec_pars(0));
    vector<T> rp2 = par.rec_pars.template cast<T>();
    return exp(functionalStockRecruitment(ssb, rp2, conf.stockRecruitmentModelCode));
  }

  AD<Type> SR(CppAD::vector<AD<Type> > ssb){
    AD<Type> s0 = ssb[0];
    vector<AD<Type> > rp(par.rec_pars.size());
    rp = par.rec_pars.template cast<AD<Type> >();
    if(conf.stockRecruitmentModelCode == 0)
      return 0.0;
    if(conf.stockRecruitmentModelCode == 62)
      return exp(rp(0));
    return exp(functionalStockRecruitment(s0, rp, conf.stockRecruitmentModelCode));
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
      nll -= log(yield(exp(logFmsy)));
    }

    if(CppAD::Variable(par.logScaleFmax)){
      nll -= log(YPR(exp(logFmax)));
    }

    if(CppAD::Variable(par.logScaleF01)){
     Type tmp = 0.1 * dYPR(Type(SAM_Zero)) - dYPR(exp(logF01));
      nll += tmp * tmp;
    }

    if(CppAD::Variable(par.logScaleFcrash)){
      Type tmp = dSR(Type(SAM_Zero)) - (1.0 / SPR(exp(logFcrash)));
      nll += tmp * tmp;
    }
    
    if(CppAD::Variable(par.logScaleFext)){
      Type tmp1 = Se(exp(logFext));
      Type tmp2 = exp(logFext);
      nll += tmp1 * tmp1 + tmp2 * tmp2;
    }
    
    for(int i = 0; i < par.logScaleFxPercent.size(); ++i){
      if(CppAD::Variable(par.logScaleFxPercent(i))){
	Type tmp = dat.referencepoint.xPercent(i) * SPR(Type(SAM_Zero)) - SPR(exp(logFxPercent(i)));
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

  REFERENCE_POINTS<Type> referencepoint(dat, conf, par, logN, logF);


  if(dat.referencepoint.Fsequence.size() > 0){
    vector<Type> Fseq = dat.referencepoint.Fsequence;
    vector<Type> logYPR(Fseq.size());
    logYPR.setZero();
    vector<Type> logSPR(Fseq.size());
    logSPR.setZero();
    vector<Type> logSe(Fseq.size());
    logSe.setZero();
    vector<Type> logYe(Fseq.size());
    logYe.setZero();
    vector<Type> logRe(Fseq.size());
    logRe.setZero();
    Type dSR0 = 0.0;

    for(int i = 0; i < Fseq.size(); ++i){
      PERREC_t<Type> v = perRecruit<Type, Type>(Fseq(i),
						referencepoint.dat,
						referencepoint.conf,
						referencepoint.par,
						referencepoint.sel,
						referencepoint.aveYears,
						referencepoint.nYears);
      logYPR(i) = v.logYPR;
      logSPR(i) = v.logSPR;
      logSe(i) = v.logSe;
      logYe(i) = v.logYe;
      logRe(i) = v.logRe;
      if(i == 0)
	dSR0 = v.dSR0;
	
    }

    REPORT_F(logYPR, of);
    REPORT_F(logSPR, of);
    REPORT_F(logSe, of);
    REPORT_F(logYe, of);
    REPORT_F(logRe, of);
    REPORT_F(dSR0, of);

    ADREPORT_F(logYPR, of);
    ADREPORT_F(logSPR, of);
    ADREPORT_F(logSe, of);
    ADREPORT_F(logYe, of);
    ADREPORT_F(logRe, of);
    ADREPORT_F(dSR0, of);
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

  ADREPORT_F(referencepoint.logFmax,of);
  ADREPORT_F(referencepoint.logBmax,of);
  ADREPORT_F(referencepoint.logRmax,of);
  ADREPORT_F(referencepoint.logYmax,of);
  ADREPORT_F(referencepoint.logYPRmax,of);
  ADREPORT_F(referencepoint.logSPRmax,of);
 
  ADREPORT_F(referencepoint.logF01,of);
  ADREPORT_F(referencepoint.logB01,of);
  ADREPORT_F(referencepoint.logR01,of);
  ADREPORT_F(referencepoint.logY01,of);
  ADREPORT_F(referencepoint.logYPR01,of);
  ADREPORT_F(referencepoint.logSPR01,of);

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

  ADREPORT_F(referencepoint.logFlim,of);
  ADREPORT_F(referencepoint.logBlim,of);
  ADREPORT_F(referencepoint.logRlim,of);
  ADREPORT_F(referencepoint.logYlim,of);
  ADREPORT_F(referencepoint.logYPRlim,of);
  ADREPORT_F(referencepoint.logSPRlim,of);

  // Fbar relative to (last year) reference points
  vector<Type> fbar = fbarFun(conf, logF);
  vector<Type> logfbar = log(fbar);

  vector<Type> relref_logfbar_fmsy = logfbar - referencepoint.logFmsy;
  ADREPORT_F(relref_logfbar_fmsy,of);
  vector<Type> relref_logfbar_fmax = logfbar - referencepoint.logFmax;
  ADREPORT_F(relref_logfbar_fmax,of);
  vector<Type> relref_logfbar_f01 = logfbar - referencepoint.logF01;
  ADREPORT_F(relref_logfbar_f01,of);
  // vector<Type> relref_logfbar_f35 = logfbar - referencepoint.logF35;
  // ADREPORT_F(relref_logfbar_f35,of);

  // SSB relative to (last year) reference points
    vector<Type> ssb = ssbFun(dat, conf, logN, logF);
  vector<Type> logssb = log(ssb);

  vector<Type> relref_logssb_bmsy = logssb - referencepoint.logBmsy;
  ADREPORT_F(relref_logssb_bmsy,of);
  vector<Type> relref_logssb_bmax = logssb - referencepoint.logBmax;
  ADREPORT_F(relref_logssb_bmsy,of);
  vector<Type> relref_logssb_b01 = logssb - referencepoint.logB01;
  ADREPORT_F(relref_logssb_b01,of);
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
  
extern "C" {

  SEXP perRecruitR(SEXP Fbar, SEXP dat, SEXP conf, SEXP pl, SEXP sel, SEXP aveYears, SEXP nYears){
    dataSet<double> d0(dat);
    confSet c0(conf);
    paraSet<double> p0(pl);
    vector<double> s0 = asVector<double>(sel);
    vector<int> a0 = asVector<int>(aveYears);
    double Fbar0 = Rf_asReal(Fbar);
    int nY0 = Rf_asInteger(nYears);
 
    PERREC_t<double> y = perRecruit<double, double>(Fbar0, d0, c0, p0, s0, a0, nY0);
    const char *resNms[] = {"logF", "logYPR", "logSPR", "logSe", "logRe", "logYe", ""}; // Must end with ""
    SEXP res;
    PROTECT(res = Rf_mkNamed(VECSXP, resNms));
    SET_VECTOR_ELT(res, 0, asSEXP(y.logFbar));
    SET_VECTOR_ELT(res, 1, asSEXP(y.logYPR));
    SET_VECTOR_ELT(res, 2, asSEXP(y.logSPR));
    SET_VECTOR_ELT(res, 3, asSEXP(y.logSe));
    SET_VECTOR_ELT(res, 4, asSEXP(y.logRe));
    SET_VECTOR_ELT(res, 5, asSEXP(y.logYe));

    UNPROTECT(1);    
    return res;

  }


  SEXP stockRecruitmentModelR(SEXP ssb, SEXP rec_pars, SEXP code){
    double b = Rf_asReal(ssb);
    vector<double> rp = asVector<double>(rec_pars);
    int srmc = Rf_asInteger(code);
	
    double v = exp(functionalStockRecruitment(b, rp, srmc));

#ifdef CPPAD_FRAMEWORK
    vector<AD<double> > rp2(rp.size()+1);
    rp2 = rp.template cast<AD<double> >();
    CppAD::vector<AD<double> > x( 1 );
    x[0] = b;
    CppAD::Independent(x);
    CppAD::vector<AD<double> > y( 1 );
    y[0] = exp(functionalStockRecruitment(x[0], rp2, srmc));
    CppAD::ADFun<double> F(x, y);
    CppAD::vector<double> x_eval( 1 );
    x_eval[0] = b;
    vector<double> r = F.Jacobian(x_eval);
#endif
#ifdef TMBAD_FRAMEWORK
   
    F_dFunctionalSR2 Fd = {rp.size(),srmc};
    vector<double> x(rp.size() + 1);
    for(int i = 0; i < rp.size(); ++i)
      x(i) = rp(i);
    x(rp.size()) = b;
    TMBad::ADFun<> G(TMBad::StdWrap<F_dFunctionalSR2,vector<TMBad::ad_aug> >(Fd), x);
    // TMBad::ADFun<> G(Fd,x);
    G = G.JacFun();
    vector<double> r = G(x);
#endif

    const char *resNms[] = {"Recruits", "Gradient", ""}; // Must end with ""
    SEXP res;
    PROTECT(res = Rf_mkNamed(VECSXP, resNms));
    SET_VECTOR_ELT(res, 0, asSEXP(v));
    SET_VECTOR_ELT(res, 1, asSEXP(r));
 
    UNPROTECT(1);    
    return res;
      
  }

  SEXP Se_sbhR(SEXP lambda, SEXP a, SEXP b, SEXP g){
    double r = Se_sbh(Rf_asReal(lambda), Rf_asReal(a), Rf_asReal(b), Rf_asReal(g));
    return asSEXP(r);
  }

  SEXP Se_slR(SEXP lambda, SEXP a, SEXP b, SEXP g){
    double r = Se_sl(Rf_asReal(lambda), Rf_asReal(a), Rf_asReal(b), Rf_asReal(g));
    return asSEXP(r);
  }
  
}
