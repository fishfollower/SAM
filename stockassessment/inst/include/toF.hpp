#pragma once
#ifndef SAM_TOF_HPP
#define SAM_TOF_HPP


enum ConstraintType {
		     Constrain_Fbar,
		     Constrain_Catch,
		     Constrain_SSB,
		     Constrain_TSB,
		     Constrain_Landing,
		     Constrain_KeepRelF,
		     Constrain_NONE = 99
};


template<class Type>
struct FConstraint {
  int Amin;
  int Amax;
  int fleet;
  bool relative;
  ConstraintType cstr;
  Type target;

  FConstraint() = default;

  FConstraint(SEXP x){
    Amin = (int)*REAL(getListElement(x,"Amin"));
    Amax = (int)*REAL(getListElement(x,"Amax"));
    fleet = (int)*REAL(getListElement(x,"fleet"));
    relative = (bool)*REAL(getListElement(x,"relative"));
    cstr = static_cast<ConstraintType>((int)*REAL(getListElement(x,"cstr")));
    target = (Type)*REAL(getListElement(x,"target"));
  };

  template<class T>
  FConstraint(const FConstraint<T>& x) : Amin(x.Amin),
					 Amax(x.Amax),
					 fleet(x.fleet),
					 relative(x.relative),
					 cstr(x.cstr),
					 target(x.target) {};  
};

template<class Type>
struct FConstraintList : vector<FConstraint<Type> > {
  FConstraintList() : vector<FConstraint<Type> >() {};
  FConstraintList(int n) : vector<FConstraint<Type> >(n) {};
  FConstraintList(SEXP x){ 
    (*this).resize(Rf_length(x));
    for(int i=0; i<Rf_length(x); i++){
      (*this)(i) = FConstraint<Type>(VECTOR_ELT(x, i));
    }
  }
  template<class T>
  FConstraintList(const FConstraintList<T>& other) : vector<FConstraint<Type> >(other.size()) {
    for(int i = 0; i < other.size(); ++i)
      (*this)(i) = FConstraint<T>(other(i));
  }
  
};




namespace ConstrainCalculations {



 
  // Newton approach


  std::vector<int> getCatchFleets(vector<int> fleetTypes){
    std::vector<int> r;
    for(int i = 0; i < fleetTypes.size(); ++i)
      if(fleetTypes(i) == 0)
	r.push_back(i);
    return r;
  };
  
  template<class Type>
  vector<Type> getFleetLogFbar(dataSet<Type>& dat, confSet& conf, vector<Type>& logF, int a0, int a1){
    std::vector<int> cFleets = getCatchFleets(dat.fleetTypes);
    vector<Type> fbar(cFleets.size());
    fbar.setConstant(0.0);
    for(int i = 0; i < cFleets.size(); ++i){
      int f = cFleets[i];
      for(int a=a0; a<=a1; a++){
	if(conf.keyLogFsta(f,a-conf.minAge) > (-1))
	  fbar(i) += exp(logF(conf.keyLogFsta(f,a-conf.minAge)));
      }
      fbar(i) /= Type(a1-a0 + 1.0);
    }
    return log(fbar);
  };

  template<class Type>
  Type getFleetCatch(dataSet<Type>& dat, confSet& conf, array<Type>& logN, vector<Type>& logF, int y, int a0, int a1, int fleet){
    Type logCat = R_NegInf;
    std::vector<int> cFleets = getCatchFleets(dat.fleetTypes);
    for(int a=a0; a<=a1; a++){
      Type logZa = log(dat.natMor(y, a-conf.minAge));
      for(int ii = 0; ii < cFleets.size(); ++ii){
	int f = cFleets[ii];
	if(conf.keyLogFsta(f,a-conf.minAge) > (-1))
	  logZa = logspace_add_SAM(logZa, logF(conf.keyLogFsta(f,a-conf.minAge)));
      }
      Type logv = logspace_sub_SAM(Type(0.0), -exp(logZa)) - logZa;
      int f0 = 0;
      int f1 = cFleets.size()-1;
      if(fleet > (-1)){
	f0 = fleet;
	f1 = fleet;
      }
    for(int ii = f0; ii <= f1; ++ii){
	int f = cFleets[ii];
	if(conf.keyLogFsta(f,a-conf.minAge) > (-1)){
	  Type logFI = logv + logF(conf.keyLogFsta(f,a-conf.minAge));
	  Type lc =  logFI + logN(a-conf.minAge,y) + log(dat.catchMeanWeight(y, a-conf.minAge, f));
	  logCat = logspace_add_SAM(logCat, lc);
	}
      }
    }
    return logCat;
  }
  
  template<class Type>
  Type getFleetLanding(dataSet<Type>& dat, confSet& conf, array<Type>& logN, vector<Type>& logF, int y, int a0, int a1, int fleet){
        Type logCat = R_NegInf;
    std::vector<int> cFleets = getCatchFleets(dat.fleetTypes);
    for(int a=a0; a<=a1; a++){
      Type logZa = log(dat.natMor(y, a-conf.minAge));
      for(int ii = 0; ii < cFleets.size(); ++ii){
	int f = cFleets[ii];
	if(conf.keyLogFsta(f,a-conf.minAge) > (-1))
	  logZa = logspace_add_SAM(logZa, logF(conf.keyLogFsta(f,a-conf.minAge)));
      }
      Type logv = logspace_sub_SAM(Type(0.0), -exp(logZa)) - logZa;
      int f0 = 0;
      int f1 = cFleets.size()-1;
      if(fleet > (-1)){
	f0 = fleet;
	f1 = fleet;
      }
    for(int ii = f0; ii <= f1; ++ii){
	int f = cFleets[ii];	
	if(conf.keyLogFsta(f,a-conf.minAge) > (-1)){
	  Type LW = dat.landMeanWeight(y, a-conf.minAge, f);
	  Type LF = dat.landFrac(y,a-conf.minAge,f);
	  if(LW > 0 && LF > 0){
	    Type logFI = logv + logF(conf.keyLogFsta(f,a-conf.minAge));
	    Type lc =  logFI + logN(a-conf.minAge,y) + log(dat.catchMeanWeight(y, a-conf.minAge, f));
	    logCat = logspace_add_SAM(logCat, lc);
	  }
	}
      }
    }
    // Type logCat = R_NegInf;
    // std::vector<int> cFleets = getCatchFleets(dat.fleetTypes);
    // for(int a=a0; a<=a1; a++){
    //   Type logZa = log(dat.natMor(y, a-conf.minAge));
    //   Type logv = logspace_sub_SAM(Type(0.0), -exp(logZa)) - logZa;
    //   for(int ii = 0; ii < cFleets.size(); ++ii){
    // 	int f = cFleets[ii];
    // 	if(conf.keyLogFsta(f,a-conf.minAge) > (-1))
    // 	  logZa = logspace_add_SAM(logZa, logF(conf.keyLogFsta(f,a-conf.minAge)));
    //   }
    //   int f0 = 0;
    //   int f1 = cFleets.size();
    //   if(fleet >= 0){
    // 	f0 = fleet;
    // 	f1 = fleet;
    //   }
    //   for(int ii = f0; ii <= f1; ++ii){
    // 	int f = cFleets[ii];
    // 	Type LW = dat.landMeanWeight(y, a-conf.minAge, f);
    // 	Type LF = dat.landFrac(y,a-conf.minAge,f);
    // 	if(LW > 0 && LF > 0 && conf.keyLogFsta(f,a-conf.minAge) > (-1)){
    // 	  Type logFI = logv + logF(conf.keyLogFsta(f,a-conf.minAge));
    // 	  Type lc =  log(LF) + logFI + logN(a-conf.minAge,y) + log(LW);
    // 	  logCat = logspace_add_SAM(logCat, lc);
    // 	}
    //   }
    // }
    return logCat;
  }

  // Begining of next year
  template<class Type>
  Type getSSB(dataSet<Type>& dat, confSet& conf, Recruitment<Type> &recruit, array<Type>& logN, vector<Type>& logF, int y, int a0, int a1, bool rel = false){
    // Current SSB for recruitment
    std::vector<int> cFleets = getCatchFleets(dat.fleetTypes);
    int yn = std::min(y+1, dat.propMat.dim[0]-1);
    Type logThisSSB = R_NegInf;
    int ys = std::max(yn-conf.minAge, 0); // Year of birth for predicted logN
    for(int a = 0; a < logN.dim[0]; a++){
      logThisSSB = logspace_add_SAM(logThisSSB, logN(a,ys) + log(dat.propMat(ys,a)) + log(dat.stockMeanWeight(ys,a)));
    }
    // Next N
    vector<Type> logNp(logN.dim[0]);
    logNp.setConstant(R_NegInf);
    //// Recruit
    logNp(0) = recruit(logThisSSB,logN(0,y), dat.years(0)+y+1);
    //// Middle
    for(int a = 1; a < logNp.size(); ++a){
      Type logZa = log(dat.natMor(y, a-1));
      Type logv = logspace_sub_SAM(Type(0.0), -exp(logZa)) - logZa;
      for(int ii = 0; ii < cFleets.size(); ++ii){
	int f = cFleets[ii];
	if(conf.keyLogFsta(f,a-1) > (-1))
	  logZa = logspace_add_SAM(logZa, logF(conf.keyLogFsta(f,a-1)));
      }     
      logNp(a) = logN(a-1,y) - exp(logZa);
    }
    //// Plus group
    if(conf.maxAgePlusGroup(0)==1){
      int a = logNp.size()-1;
      Type v1 = logNp(a);
   Type logZa = log(dat.natMor(y, a));
      Type logv = logspace_sub_SAM(Type(0.0), -exp(logZa)) - logZa;
      for(int ii = 0; ii < cFleets.size(); ++ii){
	int f = cFleets[ii];
	if(conf.keyLogFsta(f,a) > (-1))
	  logZa = logspace_add_SAM(logZa, logF(conf.keyLogFsta(f,a)));
      }
      Type v2 = logN(a, y) - exp(logZa);
      logNp(a) = logspace_add_SAM(v1,v2);
    }    
   // Next SSB
    Type logssb = R_NegInf;
    Type llssb = R_NegInf;
    for(int a = a0; a <= a1; a++){
      if(dat.propMat(yn,a-conf.minAge) > 0)
	logssb = logspace_add_SAM(logssb, logNp(a-conf.minAge) + log(dat.propMat(yn,a-conf.minAge)) + log(dat.stockMeanWeight(yn,a-conf.minAge)));
      if(dat.propMat(yn,a-conf.minAge) > 0)
	llssb = logspace_add_SAM(llssb, logN(a-conf.minAge,y) + log(dat.propMat(y,a-conf.minAge)) + log(dat.stockMeanWeight(y,a-conf.minAge)));
    }
    //
    if(rel)
      return logssb - llssb;
    return logssb;    
  };

  // Begining of next year
  template<class Type>
  Type getTSB(dataSet<Type>& dat, confSet& conf, Recruitment<Type> &recruit, array<Type>& logN, vector<Type>& logF, int y, int a0, int a1, bool rel = false){
    // Current SSB for recruitment
    std::vector<int> cFleets = getCatchFleets(dat.fleetTypes);
    int yn = std::min(y+1, dat.propMat.dim[0]-1);
    Type logThisSSB = R_NegInf;
    int ys = std::max(yn-conf.minAge, 0); // Year of birth for predicted logN
    for(int a = 0; a < logN.dim[0]; a++){
      logThisSSB = logspace_add_SAM(logThisSSB, logN(a,ys) + log(dat.propMat(ys,a)) + log(dat.stockMeanWeight(ys,a)));
    }
    // Next N
    vector<Type> logNp(logN.dim[0]);
    logNp.setConstant(R_NegInf);
    //// Recruit
    logNp(0) = recruit(logThisSSB,logN(0,y), dat.years(0)+y+1);
    //// Middle
    for(int a = 1; a < logNp.size(); ++a){
      Type logZa = log(dat.natMor(y, a-1));
      Type logv = logspace_sub_SAM(Type(0.0), -exp(logZa)) - logZa;
      for(int ii = 0; ii < cFleets.size(); ++ii){
	int f = cFleets[ii];
	if(conf.keyLogFsta(f,a-1) > (-1))
	  logZa = logspace_add_SAM(logZa, logF(conf.keyLogFsta(f,a-1)));
      }     
      logNp(a) = logN(a-1,y) - exp(logZa);
    }
    //// Plus group
    if(conf.maxAgePlusGroup(0)==1){
      int a = logNp.size()-1;
      Type v1 = logNp(a);
   Type logZa = log(dat.natMor(y, a));
      Type logv = logspace_sub_SAM(Type(0.0), -exp(logZa)) - logZa;
      for(int ii = 0; ii < cFleets.size(); ++ii){
	int f = cFleets[ii];
	if(conf.keyLogFsta(f,a) > (-1))
	  logZa = logspace_add_SAM(logZa, logF(conf.keyLogFsta(f,a)));
      }
      Type v2 = logN(a, y) - exp(logZa);
      logNp(a) = logspace_add_SAM(v1,v2);
    }    
    // Next TSB
    Type logtsb = R_NegInf;
    Type lltsb = R_NegInf;
    for(int a = a0; a <= a1; a++){
      logtsb = logspace_add_SAM(logtsb, logNp(a-conf.minAge) + log(dat.stockMeanWeight(yn,a-conf.minAge)));
      lltsb = logspace_add_SAM(lltsb, logN(a-conf.minAge,y) + log(dat.stockMeanWeight(y,a-conf.minAge)));
    }
    //
    if(rel)
      return logtsb - lltsb;
    return logtsb;
  };
  
  typedef TMBad::ad_aug ad;

  struct ForecastF {
    dataSet<ad> dat;
    confSet conf;
    Recruitment<ad> recruit;
    FConstraintList<ad> cstrs;
    vector<ad> lastLogF;
    array<ad> logN;
    int y;

    ad operator()(const vector<ad>& logFs){
      ad kappa = 0.0;
      //
      // matrix<ad> llogF = toFleetMatrix(dat, conf, lastLogF, (vector<ad>)(logFs * ad(0.0)));
      // vector<ad> llfs = (vector<ad>)(lastLogF - lastLFB);
      // matrix<ad> logF = toFleetMatrix(dat, conf, lastLogF, logFs);

      vector<ad> newLogF = lastLogF;// - lastLogFbar;
      std::vector<int> cFleets = getCatchFleets(dat.fleetTypes);
 
      vector<bool> done(newLogF.size());
      done.setConstant(false);
      for(int i = 0; i < cFleets.size(); ++i){
	int f = cFleets[i];
	for(int a = 0; a < conf.keyLogFsta.dim(1); ++a){
	  int indx = conf.keyLogFsta(f,a);
	  if(indx > (-1) && !done(indx)){
	    newLogF(indx) += logFs(i);
	    done(indx) = true;
	  }
	}
      }
      
      for(int i = 0; i < cstrs.size(); ++i){
	FConstraint<ad> cstr = cstrs(i);
	if(cstr.cstr == ConstraintType::Constrain_NONE){
	  continue;
	}
	// Previous F values
	vector<ad> lastFleetLogFbar = getFleetLogFbar(dat, conf, lastLogF, cstr.Amin, cstr.Amax);
	ad lastLogFbar = logspace_sum(lastFleetLogFbar);
	// New F values
	vector<ad> fleetLogFbar = getFleetLogFbar(dat, conf, newLogF, cstr.Amin, cstr.Amax);
	ad logFbar = logspace_sum(fleetLogFbar);
	// Add constraint

	if(cstr.cstr == ConstraintType::Constrain_Fbar){

	  if(cstr.fleet == (-1)){	// Total F
	    ad trgt = cstr.target;
	    if(cstr.relative)
	      trgt += lastLogFbar;
	    ad tmp = logFbar - trgt;
	    kappa += tmp * tmp;
	  }else{
	    ad trgt = cstr.target;
	    if(cstr.relative)
	      trgt += lastFleetLogFbar(cstr.fleet);
	    ad tmp = fleetLogFbar(cstr.fleet) - trgt;
	    kappa += tmp * tmp;
	  }

	}else if(cstr.cstr == ConstraintType::Constrain_Catch){

	  ad logC = getFleetCatch(dat, conf, logN, newLogF, y, cstr.Amin, cstr.Amax, cstr.fleet);
	  ad trgt = cstr.target;
	  if(cstr.relative){
	    ad logCL = getFleetCatch(dat, conf, logN, lastLogF, y-1, cstr.Amin, cstr.Amax, cstr.fleet);
	    trgt += logCL;
	  }
	  ad tmp = logC - trgt;
	  kappa += tmp * tmp;

	}else if(cstr.cstr == ConstraintType::Constrain_SSB){
	  
	  ad logB = getSSB(dat,conf, recruit, logN, newLogF, y, cstr.Amin, cstr.Amax, cstr.relative);
	  ad tmp = logB - cstr.target;
	  kappa += tmp * tmp;
	  
	}else if(cstr.cstr == ConstraintType::Constrain_TSB){
	  
	  ad logB = getTSB(dat,conf, recruit, logN, newLogF, y, cstr.Amin, cstr.Amax, cstr.relative);
	  ad tmp = logB - cstr.target;
	  kappa += tmp * tmp;
	  
	}else if(cstr.cstr == ConstraintType::Constrain_Landing){

	  ad logL = getFleetLanding(dat, conf, logN, newLogF, y, cstr.Amin, cstr.Amax, cstr.fleet);
	  ad trgt = cstr.target;
	  if(cstr.relative){
	    ad logLL = getFleetLanding(dat, conf, logN, lastLogF, y-1, cstr.Amin, cstr.Amax, cstr.fleet);
	    trgt += logLL;
	  }
	  ad tmp = logL - trgt;
	  kappa += tmp * tmp;

	}else if(cstr.cstr == ConstraintType::Constrain_KeepRelF){

	  SAM_ASSERT(cstr.fleet >= 0,"Keep relative Fbar can only be used for fleets, not total")
	  ad tmp = (lastFleetLogFbar(cstr.fleet) - lastFleetLogFbar(CppAD::Integer(cstr.target))) - (fleetLogFbar(cstr.fleet) - fleetLogFbar(CppAD::Integer(cstr.target)));
	  kappa += tmp * tmp;	  
	}else{
	  Rf_error("Constraint type not implemented");
	}
      }
      return kappa;
    }
   
  };

};				// End of namespace


template<class Type>
vector<Type> calculateNewFVec(dataSet<Type>& dat,
			      confSet& conf,
			      paraSet<Type>& par,
			      FConstraintList<Type>& cstrs,
			      vector<Type>& lastLogF,
			      array<Type>& logN,
			      int y,
			      newton::newton_config& cfg){

  paraSet<TMBad::ad_aug> parad(par);
  Recruitment<TMBad::ad_aug> recruit = makeRecruitmentFunction(conf,parad);

  std::vector<int> cFleets = ConstrainCalculations::getCatchFleets(dat.fleetTypes);

  // vector<Type> Zero(cFleets.size()); Zero.setZero();
  // matrix<Type> llogF = ConstrainCalculations::toFleetMatrix(dat, conf, lastLogF, Zero);    
  // vector<Type> lastFleetLogFbar = ConstrainCalculations::getFleetLogFbar(conf, llogF, conf.fbarRange[0], conf.fbarRange[1]);
  // Type lastLogFbar = logspace_sum(lastFleetLogFbar);

  ConstrainCalculations::ForecastF fc = {dat, conf, recruit, cstrs, lastLogF, logN, y};

  
  vector<Type> start(cFleets.size());
  start.setConstant(0);
  
  cfg.simplify = false;		// Needed for logspace_add in older versions of TMB
  vector<Type> res = newton::Newton(fc, start, cfg);
  vector<Type> newLogF = lastLogF;// - lastLogFbar;

  vector<bool> done(newLogF.size());
  done.setConstant(false);
  for(int i = 0; i < cFleets.size(); ++i){
    int f = cFleets[i];
    for(int a = 0; a < conf.keyLogFsta.dim(1); ++a){
      int indx = conf.keyLogFsta(f,a);
      if(indx > (-1) && !done(indx)){
	newLogF(indx) += res(i);
	done(indx) = true;
      }
    }
  }
  return newLogF;
}




// Faster versions

// template<class Type>
// struct CATCH2F_QUICK {
//   vector<Type> Flast;
//   vector<Type> M;
//   vector<Type> N;
//   vector<Type> w;
//   vector<Type> frac;
//   Type catchval;
    
//   Type operator()(Type logFScale){
//     vector<Type> Fa = exp(logFScale) * Flast;
//     vector<Type> Z = Fa + M;
//     vector<Type> C = Fa * (Type(1.0) - exp(-Z)) * N / Z * w * frac;
//     return log(catchval) - log(sum(C) + 1e-5) + logFScale; //(sum(C) / catchval) * FScale;
//   }    
// };



// template<class Type>
// Type catch2F_quick(Type catchval, vector<Type> lastF, vector<Type> M, vector<Type> N, vector<Type> w) {
//   Type sv = 0.0;
//   int maxAge = lastF.size();
//   vector<Type> frac(maxAge);
//   frac.setZero();
//   frac += 1.0;
//   CATCH2F_QUICK<Type> f = {lastF, M, N, w, frac, catchval};
//   for(int i = 0; i < 100; ++i){
//     Type tmp = f(sv);
//     sv = tmp;
//   }
//   return exp(sv);
// }


// template<class Type>
// Type landing2F_quick(Type landingval, vector<Type> lastF, vector<Type> M, vector<Type> N, vector<Type> w, vector<Type> frac) {
//   Type sv = 0.0;
//   CATCH2F_QUICK<Type> f = {lastF, M, N, w, frac, landingval};
//   for(int i = 0; i < 100; ++i){
//     Type tmp = f(sv);
//     sv = tmp;
//   }
//   return exp(sv);
// }

// template<class Type>
// struct SSB2F_QUICK {
//   vector<Type> logFlast;
    
//   array<Type> logN;
//   array<Type> logF;

//   MortalitySet<Type> mort;
  
//   confSet cf;
//   dataSet<Type> ds;
//   paraSet<Type> ps;
//   int i;

//   Type ssbval;
//   Type rec_mean;

//   Type f(Type logFScale){
//     Recruitment<Type> recruit = Recruitment<Type>("zero",new Rec_None<Type>());
//     array<Type> lN = logN;
//     array<Type> lF = logF;
//     lF += logFScale;
//     vector<Type> nextN = predNFun(ds, cf, ps, lN, lF, recruit, mort, i);
//     if(!isNA(rec_mean))
//       nextN(0) = rec_mean;
//     int jj = i;
//     if(i + 1 < lN.cols())
//       jj = i+1;
//     Type newSSB = 0.0;
//     for(int q = 0; q < nextN.size(); ++q)
//       newSSB += exp(nextN(q)) * ds.propMat(jj,q) * ds.stockMeanWeight(jj,q);    
//     // if(i + 1 < logN.cols()){      
//     //   lN.col(i+1) = nextN;
//     //   lF.col(i+1) = logFScale + logFlast;
//     //   newSSB = ssbi(ds, cf, lN, lF, i+1);
//     // }else{
//     //   lN.col(i) = nextN;
//     //   lF.col(i) = logFScale + logFlast;
//     //   newSSB = ssbi(ds, cf, lN, lF, i);
//     // }
//     return log(ssbval) - log(newSSB + 1e-5);
//   }

//   Type softmax(Type x, Type y, Type k = 1.0){
//     return logspace_add2(k * x, k * y) / k;
//   }
//   Type sign0(Type x){
//     return x / (fabs(x) + 1e-8);
//   }
//   Type numnewt(Type logs){
//     Type h = 0.0001 * softmax(fabs(logs),Type(0.01), Type(100.0));
//     Type a = f(logs);
//     Type g = (-f(logs + 2.0 * h) + 8.0 * f(logs + h) - 8.0 * f(logs-h) + f(logs-2.0*h)) / (12 * h);
//     // Type g = (f(logs + h) - f(logs - h)) / (2.0 * h);
//     Type s = sign0(a) * sign0(g);
//     Type y = log(fabs(a)) - log(softmax(fabs(g), (Type)0.001, (Type)1000.0)); // Damp the gradient
//     return logs - 0.9 * s * exp(y); //softmax(exp(y), 0.5 * fabs(logs), Type(1000.0));
//   }
  
// };



// template<class Type>
// Type ssb2F_quick(Type ssbval, vector<Type> logFlast, dataSet<Type> dat, confSet conf, paraSet<Type> par, array<Type> logF, array<Type> logN, MortalitySet<Type> mort, int i, Type rec_mean) {
//   Type sv = 0;
//   SSB2F_QUICK<Type> f = {logFlast, logN, logF, mort, conf, dat, par, i, ssbval, rec_mean};
//   for(int j = 0; j < 30; ++j){
//     sv = f.numnewt(sv);
//   }
//   return exp(sv);
// }


#endif
