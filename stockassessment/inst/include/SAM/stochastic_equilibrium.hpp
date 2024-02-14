SAM_DEPENDS(convenience)
SAM_DEPENDS(define)
SAM_DEPENDS(logspace)
SAM_DEPENDS(recruitment)
SAM_DEPENDS(n)
SAM_DEPENDS(extend_array)
SAM_DEPENDS(equilibrium_recycling)

HEADER(
template<class Type>
struct stochasticCalculator {

  dataSet<Type> dat;
  confSet conf;
  paraSet<Type> par;
  array<Type> logFSel;
  array<Type> logFseason;
  int CT;
  MortalitySet<Type> mort;
  // Type logFbar;
  Recruitment<Type> recruit;  
  // vector<Type> logSel;

  stochasticCalculator() : dat(), conf(), par(),logFSel(), logFseason(), CT(), mort(), recruit() {};
  
  stochasticCalculator( dataSet<Type> dat,
			confSet conf,
			paraSet<Type> par,
			array<Type> logFSel_,
			array<Type> logFseason_,
			int CT_
			// Type logFbar,
			// vector<Type>& logSel
			) :
    dat(dat), conf(conf), par(par), logFSel(logFSel_, logFSel_.dim), logFseason(logFseason_,logFseason_.dim), CT(CT_), mort(dat,conf,par,logFSel,logFseason) {
    recruit = makeRecruitmentFunction(conf,par);
  }
  
  template<class T>
  inline stochasticCalculator(const stochasticCalculator<T> &x) :
    dat(x.dat),
    conf(x.conf),
    par(x.par),
    logFSel(x.logFSel,x.logFSel.dim),
    logFseason(x.logFseason,x.logFseason.dim),
    CT(x.CT),
    mort(dat,conf,par,logFSel,logFseason)
  {
    recruit = makeRecruitmentFunction(conf,par);
  }

  Type dSR(Type logssb){
    return recruit.dSR(logssb);
  }

  Type getLogSSB(vector<Type>& logX){
    int y = dat.propM.rows()-1;
    Type logFbar = logX(0);
    vector<Type> logN = logX.segment(1,logX.size()-1);
    array<Type> logF = logFSel;
    logF.col(y) += logFbar;
    mort.updateYear(dat,conf,par,logF,logFseason,y);
    Type logThisSSB=Type(R_NegInf);
    // Calculate log-SSB
    for(int a = conf.minAge; a <= conf.maxAge; ++a){      
      int ax = a - conf.minAge;
      if(dat.propMat(y,ax) > 0){
	Type lssbNew = logN(a) + log(mort.ssbSurvival_before(ax,y)) + log(dat.propMat(y,ax)) + log(dat.stockMeanWeight(y,ax));
	logThisSSB = logspace_add_SAM(logThisSSB, lssbNew);
      }
    }

    return logThisSSB;
  }

  vector<Type> getCumulativeHazard(vector<Type>& logX){
    int y = dat.propM.rows()-1;
   Type logFbar = logX(0);
    vector<Type> logN = logX.segment(1,logX.size()-1);
    array<Type> logF = logFSel;
    logF.col(y) += logFbar;
    mort.updateYear(dat,conf,par,logF,logFseason,y);
    vector<Type> res(logN.size());
    res.setZero();
    res.segment(conf.minAge, conf.maxAge - conf.minAge + 1) = (vector<Type>)mort.cumulativeHazard.col(y);
    return res;
  }
  
  Type getLogRecruitment(vector<Type>& logX){
    return logX(conf.minAge+1);
  }

  vector<Type> predN(vector<Type>& logX){
    int y = dat.propM.rows()-1;
    Type logFbar = logX(0);
    vector<Type> logN = logX.segment(1,logX.size()-1);
    array<Type> logF = logFSel;
    logF.col(y) += logFbar;
    mort.updateYear(dat,conf,par,logF,logFseason,y);

   // From age 0 to maxAge
    vector<Type> predN(logN.size());
    predN.setConstant(R_NegInf);
    
    Type logThisSSB = getLogSSB(logX);
    Type lastLogR = R_NaReal;
    lastLogR = logN(conf.minAge);    
    predN(0) = recruit(logThisSSB, lastLogR, y); // dat.years(0) + i is needed for forecast
    if(par.rec_transphi.size() > 0){
      Rf_warning("Stochastic reference points does not currently account for autocorrelation in recruitment");
    }
  
    switch(conf.logNMeanAssumption(0)){
    case 0:			// Median on natural scale
      predN(0) += 0.0;
      break;
    case 1:			// Mean on natural scale
      predN(0) -= 0.5 * exp(2.0 * par.logSdLogN(conf.keyVarLogN(0)));
      break;
    case 2:			// Mode on natural scale
      predN(0) += exp(2.0 * par.logSdLogN(conf.keyVarLogN(0)));
      break;
    default:
      Rf_error("logNMeanCorrection not implemented.");
      break;    
    }

    vector<Type> cumHaz = getCumulativeHazard(logX);
    for(int j=1; j<logN.size(); ++j){
      // NOTE: cumHaz is zero before conf.minAge
      predN(j)=logN(j-1) - cumHaz(j-1); //totF(j-1,i-1)-dat.natMor(i-1,j-1); 
    }
    if(conf.maxAgePlusGroup(0)==1){// plusgroup adjustment if catches need them 
      Type v1 = predN(logN.size()-1); // Already updated above
      Type v2 = logN(logN.size()-1) - cumHaz(logN.size()-1); //totF(stateDimN-1,i-1) - dat.natMor(i-1,stateDimN-1); // Remaining in plus group from last year
      predN(logN.size()-1) = logspace_add_SAM(v1,v2);
    }

    for(int j=1; j<logN.size(); ++j){
      switch(conf.logNMeanAssumption(1)){
      case 0:			// Median on natural scale
	predN(j) += 0.0;
	break;
      case 1:			// Mean on natural scale
	predN(j) -= 0.5 * exp(2.0 * par.logSdLogN(conf.keyVarLogN(j)));
	break;
      case 2:			// Mode on natural scale
	predN(j) += exp(2.0 * par.logSdLogN(conf.keyVarLogN(j)));
	break;
      default:
	Rf_error("logNMeanCorrection not implemented.");
	break;    
      }
    }
    vector<Type> res(logX.size());
    res(0) = logFbar;
    res.segment(1,logN.size()) = predN;
    return res;
  }

  Type getLogCatch(vector<Type>& logX){
    typename referencepointSet<Type>::CatchType catchType = static_cast<typename referencepointSet<Type>::CatchType>(CT);
    int y = dat.propM.rows()-1;
    Type logFbar = logX(0);
    vector<Type> logN = logX.segment(1,logX.size()-1);
    array<Type> logF = logFSel;
    logF.col(y) += logFbar;
    mort.updateYear(dat,conf,par,logF,logFseason,y);
      //vector<Type> cif = getTotalCIF_F(logN);
    Type logC = R_NegInf;
    for(int a = conf.minAge; a <= conf.maxAge; ++a){
      int ax = a - conf.minAge;
      for(int f = 0; f < conf.keyLogFsta.dim(0); ++f){
	if(dat.fleetTypes(f) == 0){ // Only catch fleets
	  Type tmp;
	  Type tmp0 = logN(a) + mort.logFleetSurvival_before(ax,y,f) + log(mort.fleetCumulativeIncidence(ax,y,f));
	  switch(catchType){
	   case referencepointSet<Type>::totalCatch:
	     tmp = tmp0 + log(dat.catchMeanWeight(y, ax, f));
	     break;
	  case referencepointSet<Type>::landings:
	    tmp = tmp0 + log(dat.landMeanWeight(y,ax,f)) + log(dat.landFrac(y,ax,f));
	    break;
	  case referencepointSet<Type>::discard:
	    tmp = tmp0 + log(dat.disMeanWeight(y,ax,f)) + log(1.0 - (dat.landFrac(y,ax,f)-1e-8));
	    break;
	  default:
	     Rf_error("Unknown reference point catch type.");
	     break;
	  }
	  logC = logspace_add_SAM(logC, tmp);	 
	}
      }
    }
    return logC;
  }
  
  
};
       )

HEADER(
       template<class Type>
       struct S_Functor_N {
	 stochasticCalculator<Type> calc;

	 template<class T>
	 vector<Type> operator()(vector<Type>& logX){
	   return calc.predN(logX);
	 }
	 
	 template<class T>
	 vector<T> operator()(vector<T>& logX){
	   stochasticCalculator<T> c2(calc); 
	   return c2.predN(logX);
	 }
       };
       )


// SOURCE(
//        template<class Type>
//        struct S_Functor_N_i {
// 	 stochasticCalculator<Type> calc;
// 	 int i;
// 	 template<class T>
// 	 T operator()(vector<T>& logX){
// 	   stochasticCalculator<T> c2(calc); 
// 	   return c2.predN(logX)(i);
// 	 }
//        };
//        )



HEADER(
       template<class Type>
       struct S_Functor_S {
	 stochasticCalculator<Type> calc;

	 vector<Type> operator()(vector<Type>& logX){
	   vector<Type> res(1);
	   res(0) = calc.getLogSSB(logX);
	   return res;
	 }
	 
	 template<class T>
	 vector<T> operator()(vector<T>& logX){
	   stochasticCalculator<T> c2(calc);
	   vector<T> res(1);
	   res(0) = c2.getLogSSB(logX);
	   return res;
	 }
       };
       )

HEADER(
       template<class Type>
       struct S_Functor_C {
	 stochasticCalculator<Type> calc;

	 vector<Type> operator()(vector<Type>& logX){
	   vector<Type> res(1);
	   res(0) = calc.getLogCatch(logX);
	   return res;
	 }
	 
	 template<class T>
	 vector<T> operator()(vector<T>& logX){
	   stochasticCalculator<T> c2(calc); 
	   vector<T> res(1);
	   res(0) = c2.getLogCatch(logX);
	   return res;
	 }
       };
       )

HEADER(
       template<class Type>
       struct S_Functor_SPR {
	 stochasticCalculator<Type> calc;

	 vector<Type> operator()(vector<Type>& logX){
	   vector<Type> res(1);
	   res(0) =  calc.getLogSSB(logX) - calc.getLogRecruitment(logX);
	   return res;
	 }
	 
	 template<class T>
	 vector<T> operator()(vector<T>& logX){
	   stochasticCalculator<T> c2(calc); 
	   vector<T> res(1);
	   res(0) =  c2.getLogSSB(logX) - c2.getLogRecruitment(logX);
	   return res;
	 }
       };
       )

HEADER(
       template<class Type>
       struct S_Functor_CPR {
  	 stochasticCalculator<Type> calc;

	 vector<Type> operator()(vector<Type>& logX){
	   vector<Type> res(1);
	   res(0) = calc.getLogCatch(logX) - calc.getLogRecruitment(logX);
	   return res;
	 }
	 
	 template<class T>
	 vector<T> operator()(vector<T>& logX){
	   stochasticCalculator<T> c2(calc); 
	   vector<T> res(1);
	   res(0) = c2.getLogCatch(logX) - c2.getLogRecruitment(logX);
	   return res;
	 }
       };
       )

HEADER(
template<class Fun, class Type>
Type Deriv(Fun F, vector<Type> x0, int i, int outputIndex){
  Type h = 0.0001 * sqrt(x0(i) * x0(i) + (Type)SAM_Zero);
  vector<Type> x1 = x0;
  x1(i) = x0(i) + 2.0 * h;
  Type v1 = F(x1)(outputIndex);
  vector<Type> x2 = x0;
  x2(i) = x0(i) + h;
  Type v2 = F(x2)(outputIndex);
  vector<Type> x3 = x0;
  x3(i) = x0(i) - h;
  Type v3 = F(x3)(outputIndex);
  vector<Type> x4 = x0;
  x4(i) = x0(i) - 2.0 * h;
  Type v4 = F(x4)(outputIndex);
  Type g = (-v1 + 8.0 * v2 - 8.0 * v3 + v4) / (12.0 * h);
  return g;
}
       )

HEADER(
template<class Fun, class Type>
vector<Type> Gradient(Fun F, vector<Type> x0, int outputIndex){
  vector<Type> res(x0.size());
  for(int i = 0; i < res.size(); ++i)
    res(i) = Deriv(F,x0,i,outputIndex);
  return res;
}
       )

HEADER(
template<class Fun, class Type>
vector<vector<Type> > GradientList(Fun F, vector<Type> x0){
  int n = F(x0).size();
  vector<vector<Type> > res(n);
  for(int i = 0; i < n; ++i)
    res(i) = Gradient(F, x0, i);
  return res;
}
       )

HEADER(
template<class Fun, class Type>
Type Deriv2(Fun F, vector<Type> x0, int i, int j, int outputIndex){
  Type h1 = 0.0001 * sqrt(x0(i) * x0(i) + (Type)SAM_Zero);
  Type v0 = F(x0)(outputIndex);
  // with respect to i i
  vector<Type> x1 = x0;
  x1(i) = x0(i) + h1;
  Type v1 = F(x1)(outputIndex);
  vector<Type> x2 = x0;
  x2(i) = x0(i) - h1;
  Type v2 = F(x2)(outputIndex);
  Type vii = (v1 + v2 - 2.0 * v0) / (h1 * h1);
  if(i == j)
    return vii;
  // with respect to jj
  Type h2 = 0.0001 * sqrt(x0(j) * x0(j) + (Type)SAM_Zero);
  x1 = x0;
  x1(j) = x0(j) + h1;
  v1 = F(x1)(outputIndex);
  x2 = x0;
  x2(j) = x0(j) - h1;
  v2 = F(x2)(outputIndex);
  Type vjj = (v1 + v2 - 2.0 * v0) / (h1 * h1);
  // with respect to ij
  x1 = x0;
  x1(i) = x0(i) + h1;
  x1(j) = x0(j) + h2;
  v1 = F(x1)(outputIndex);
  x2 = x0;
  x2(i) = x0(i) - h1;
  x2(j) = x0(j) - h2;
  v2 = F(x2)(outputIndex);
  Type hdi = vii * h1 * h1;
  Type hdj = vjj * h2 * h2;
  return (v1 + v2 - 2.0 * v0 - hdi - hdj) / (2.0 * h1 * h2);
}
       )

HEADER(
template<class Fun, class Type>
matrix<Type> Hessian(Fun F, vector<Type>& x0, int outputIndex){
  // Can reuse calculations if copying from Deriv2 instead of calling it
  matrix<Type> res(x0.size(), x0.size());
  vector<Type> h(x0.size());
  Type v0 = F(x0)(outputIndex);
  // Diagonal first
  for(int i = 0; i < x0.size(); ++i){
    h(i) = 0.0001 * sqrt(x0(i) * x0(i) + (Type)SAM_Zero);
    vector<Type> x1 = x0;
    x1(i) = x0(i) + h(i);
    Type v1 = F(x1)(outputIndex);
    vector<Type> x2 = x0;
    x2(i) = x0(i) - h(i);
    Type v2 = F(x2)(outputIndex);
    Type vii = (v1 + v2 - 2.0 * v0) / (h(i) * h(i));
    res(i,i) = vii;
  }
  // Lower triangle
  for(int i = 0; i < x0.size(); ++i){
    for(int j = 0; j < i; ++j){
      vector<Type> x1 = x0;
      x1(i) = x0(i) + h(i);
      x1(j) = x0(j) + h(j);
      Type v1 = F(x1)(outputIndex);
      vector<Type> x2 = x0;
      x2(i) = x0(i) - h(i);
      x2(j) = x0(j) - h(j);
      Type v2 = F(x2)(outputIndex);
      Type hdi = res(i,i) * h(i) * h(i);
      Type hdj = res(j,j) * h(j) * h(j);
      Type vij = (v1 + v2 - 2.0 * v0 - hdi - hdj) / (2.0 * h(i) * h(j));
      res(i,j) = vij;
      res(j,i) = vij;
    }
  }
  return res;
}
       )

HEADER(
template<class Fun, class Type>
vector<matrix<Type> > HessianList(Fun F, vector<Type> x0){
  int n = F(x0).size();
  vector<matrix<Type> > res(n);
  for(int i = 0; i < n; ++i)
    res(i) = Hessian(F, x0, i);
  return res;
}
       )

HEADER(
template<class Type>
struct DerivOutput {
  vector<Type> value;
  vector<vector<Type> > gradientList;
  vector<matrix<Type> > hessianList;
};
       )

HEADER(
template<class Fun, class Type>
DerivOutput<Type> getDerivOutput(Fun F, vector<Type> x0){
  // Value
  vector<Type> v0 = F(x0);
  vector<vector<Type> > gradientList(v0.size());
  vector<matrix<Type> > hessianList(v0.size());
  matrix<Type> v_p(v0.size(), x0.size());
  v_p.setZero();
  matrix<Type> v_m(v0.size(), x0.size());
  v_m.setZero();
  vector<matrix<Type> > v_pp(x0.size());
  vector<matrix<Type> > v_mm(x0.size());  
  vector<Type> h(x0.size());
  for(int i = 0; i < x0.size(); ++i){
    h(i) = 0.0001 + 0.0001 * sqrt(x0(i) * x0(i) + (Type)SAM_Zero);
    vector<Type> xp = x0;
    xp(i) += h(i);
    vector<Type> xm = x0;
    xm(i) -= h(i);
    v_p.col(i) = F(xp);
    v_m.col(i) = F(xm);
    if(i > 0){
      matrix<Type> v_ppi(v0.size(), i);
      matrix<Type> v_mmi(v0.size(), i);
      for(int j = 0; j < i; ++j){
	vector<Type> xp2 = xp;
	xp2(j) += h(j);
	vector<Type> xm2 = xm;
	xm2(j) -= h(j);
	v_ppi.col(j) = F(xp2);
	v_mmi.col(j) = F(xm2);
      }
      v_pp(i) = v_ppi;
      v_mm(i) = v_mmi;
    }
  }
  // Calculate derivatives per outputIndex
  for(int outputIndex = 0; outputIndex < v0.size(); ++outputIndex){
    vector<Type> G(x0.size());
    matrix<Type> H(x0.size(), x0.size());
    for(int i = 0; i < x0.size(); ++i){
      G(i) = (v_p(outputIndex,i) - v_m(outputIndex,i)) / (2.0 * h(i));
      H(i,i) = (v_p(outputIndex,i) - 2.0 * v0(outputIndex) + v_m(outputIndex,i)) / (h(i) * h(i));
      for(int j = 0; j < i; ++j){
	H(i,j) = ( v_pp(i).col(j)(outputIndex) - v_p(outputIndex,i) - v_p(outputIndex,j) + 2.0 * v0(outputIndex) - v_m(outputIndex,i) - v_m(outputIndex,j) + v_mm(i).col(j)(outputIndex)) / (2.0 * h(i) * h(j));
	H(j,i) = H(i,j);
      }
    }
    gradientList(outputIndex) = G;
    hessianList(outputIndex) = H;
  }
  DerivOutput<Type> r = {v0, gradientList, hessianList};
  return r;
}
       )

HEADER(
template<class Type>
struct StochasticWorkerYearOutput {
  vector<Type> Ex;
  matrix<Type> Sx;
}
       )

HEADER(
template<class Type>
struct EquilibriumRecycler_Stochastic_Worker {

  Type logFbar0;
  confSet conf;
  paraSet<Type> par;
  vector<Type> logSel;
  vector<Type> logN0;
  int nYears;
  int CT;
  int DT; //0: AD, 1: Numeric, 2: simplified numeric
  
  dataSet<Type> newDat;
  array<Type> logFSel;
  array<Type> logitFseason;
  stochasticCalculator<Type> calc;

  vector<Type> Ex0;
  matrix<Type> Sx0;
  matrix<Type> nvar;

  S_Functor_N<Type> F;
  TMBad::ADFun<> G;
  TMBad::ADFun<> H;

  S_Functor_S<Type> F_S;
  TMBad::ADFun<> G_S;
  TMBad::ADFun<> H_S;

  S_Functor_C<Type> F_C;
  TMBad::ADFun<> G_C;
  TMBad::ADFun<> H_C;

  S_Functor_SPR<Type> F_SPR;
  TMBad::ADFun<> G_SPR;
  TMBad::ADFun<> H_SPR;

  S_Functor_CPR<Type> F_CPR;
  TMBad::ADFun<> G_CPR;
  TMBad::ADFun<> H_CPR;

  EquilibriumRecycler_Stochastic_Worker() : logFbar0(),
					    conf(),
					    par(),
					    logSel(),
					    logN0(),
					    nYears(),
					    CT(),
					    DT(),
					    newDat(),
					    logitFseason(),
					    calc(),
					    Ex0(),
					    Sx0(),
					    nvar()
  {}
    
  template<class T>
  EquilibriumRecycler_Stochastic_Worker(T logFbar0_,
					dataSet<T> dat,
					confSet conf_,
					paraSet<T> par_,
					vector<T> logSel_,
					vector<int> aveYears,
					vector<T> logN0_,
					int nYears_,
					int CT_,
					int DT_) :
    logFbar0(logFbar0_),
    conf(conf_),
    par(par_),
    logSel(logSel_),
    logN0(logN0_),
    nYears(nYears_),
    CT(CT_),
    DT(DT_),
    newDat(dat),
    logitFseason(),
    calc(),
    Ex0(),
    Sx0(),
    nvar() {
    // Setup
    if(nYears <= 0)
      Rf_error("nYears must be greater than 0.");
    if(aveYears.size() == 0)
      Rf_error("aveYears must be given.");
    if(logSel.size() != conf.keyLogFsta.maxCoeff()+1)
      Rf_error("Wrong size of selectivity vector");
    // Prepare data
    dataSet<Type> newDat = dat;
    int nMYears = dat.noYears;
    // propMat
    if(par.meanLogitMO.size() > 0)
      Rf_warning("Biopar not implemented for stochastic reference points yet!");
    //extendArray(newDat.propMat, nMYears, nYears, aveYears, par.meanLogitMO, conf.keyMatureMean, 1, false);
    extendArray(newDat.propMat, nMYears, nYears, aveYears, false);
    // stockMeanWeight
    if(par.meanLogSW.size() > 0)
      Rf_warning("Biopar not implemented for stochastic reference points yet!");
    //extendArray(newDat.stockMeanWeight, nMYears, nYears, aveYears, par.meanLogSW, conf.keyStockWeightMean, 0, false);
    extendArray(newDat.stockMeanWeight, nMYears, nYears, aveYears, false);
    // catchMeanWeight
    if(par.meanLogCW.size() > 0)
      Rf_warning("Biopar not implemented for stochastic reference points yet!");
    //extendArray(newDat.catchMeanWeight, nMYears, nYears, aveYears, par.meanLogCW, conf.keyCatchWeightMean, 0, false);
    extendArray(newDat.catchMeanWeight, nMYears, nYears, aveYears, false);
    // natMor
    if(par.meanLogNM.size() > 0)
      Rf_warning("Biopar not implemented for stochastic reference points yet!");
    //extendArray(newDat.natMor, nMYears, nYears, aveYears, par.meanLogNM, conf.keyMortalityMean, 0, false);
    extendArray(newDat.natMor, nMYears, nYears, aveYears, false);
    // landFrac (No biopar process)
    extendArray(newDat.landFrac, nMYears, nYears, aveYears, false);
    // disMeanWeight (No biopar process)
    extendArray(newDat.disMeanWeight, nMYears, nYears, aveYears, false);
    // landMeanWeight (No biopar process)
    extendArray(newDat.landMeanWeight, nMYears, nYears, aveYears, false);
    // propF (No biopar process)
    extendArray(newDat.propF, nMYears, nYears, aveYears, false);
    // propM (No biopar process)
    extendArray(newDat.propM, nMYears, nYears, aveYears, false);
    newDat.noYears = nYears;

    // Recruitment<Type> recruit = makeRecruitmentFunction(conf, par);
    // Rescale to ensure Fbar for logSel is 1
    vector<Type> totF(conf.maxAge-conf.minAge+1);
    totF.setZero();
    Type fbarSel = 0.0;
    for(int a = conf.fbarRange(0); a <= conf.fbarRange(1); ++a){
      int ax = a - conf.minAge;
      Type totF = 0.0;
      for(int f = 0; f < conf.keyLogFsta.rows(); ++f){
	if(conf.keyLogFsta(f,ax) > (-1)){
	  Type Fval = exp(logSel(conf.keyLogFsta(f,ax))) * (dat.sampleTimesEnd(f) - dat.sampleTimesStart(f));
	  totF += Fval;
	}
      }
      fbarSel += totF / Type(conf.fbarRange(1)-conf.fbarRange(0)+1.0);
    }
    logSel -= log(fbarSel);

    // Make logF array
    logFSel = array<Type>(logSel.size(), nMYears + nYears);
    logFSel.setZero();
    for(int i = 0; i < nYears; ++i)
      logFSel.col(i) = logSel;

    // Make logitF season array
    logitFseason = array<Type>(par.seasonMu.rows(), nMYears + nYears, par.seasonMu.cols());
    logitFseason.setZero();
    for(int i = 0; i < nYears; ++i)
      for(int j = 0; j < par.seasonMu.cols(); ++j)
	for(int k = 0; k < par.seasonMu.rows(); ++k)
	  logitFseason(k,i,j) = par.seasonMu(k,j);

    //MortalitySet<Type> mort(newDat, conf, par, logF, logitFseason);

    calc = stochasticCalculator<Type>(newDat, conf, par, logFSel, logitFseason, CT);

    // Initialize
    // Going from age 0 to age conf.maxAge with Fbar as the first element
    Ex0 = vector<Type>(conf.maxAge+2);
    Ex0.setZero();
    Ex0(0) = logFbar0;
    Ex0.segment(conf.minAge+1,conf.maxAge-conf.minAge+1) = logN0;
    for(int i = 1; i < conf.minAge+1; ++i)
      Ex0(i) = logN0(0);

    // Run to deterministic equilibrium
    // for(int y = 0; y < nYears; ++y){
    //   Ex = calc.predN(Ex);
    // }
    // Rcout << "Ex:\n\t" << Ex << "\n";
    int n = Ex0.size();
    Sx0 = matrix<Type>(n, n);
    Sx0.setZero();

    // Sx0.setConstant(0.9);
    // for(int i = 0; i < n; ++i)
    //   Sx0(i,i) = 1.0;

    matrix<Type> nvar0 = get_nvar(conf,par);
    nvar = matrix<Type>(n,n);
    nvar.setZero();
    nvar.block(conf.minAge+1, conf.minAge+1, nvar0.rows(), nvar0.cols()) = nvar0;

    F = {calc};
    if(DT == 0){
    G = TMBad::ADFun<>(TMBad::StdWrap<S_Functor_N<Type>,vector<TMBad::ad_aug> >(F), Ex0);
    // G.optimize();
    G = G.JacFun();
    // G.optimize();
    H = G.JacFun();
    // H.optimize();
    // G = G.atomic();
    // H = H.atomic();
    }
    
    F_S = {calc};
    if(DT == 0){
    G_S = TMBad::ADFun<> (TMBad::StdWrap<S_Functor_S<Type>,vector<TMBad::ad_aug> >(F_S), Ex0);
    // G_S.optimize();
    G_S = G_S.JacFun();
    // G_S.optimize();
    H_S = G_S.JacFun();
    // H_S.optimize();
    // G_S = G_S.atomic();
    // H_S = H_S.atomic();
    }
    
    F_C = {calc};
    if(DT == 0){
      G_C = TMBad::ADFun<>(TMBad::StdWrap<S_Functor_C<Type>,vector<TMBad::ad_aug> >(F_C), Ex0);
      // G_C.optimize();
      G_C = G_C.JacFun();
      // G_C.optimize();
      H_C = G_C.JacFun();
      // H_C.optimize();
      // G_C = G_C.atomic();
      // H_C = H_C.atomic();
    }
    
    F_CPR = {calc};
    if(DT == 0){
      G_CPR = TMBad::ADFun<>(TMBad::StdWrap<S_Functor_CPR<Type>,vector<TMBad::ad_aug> >(F_CPR), Ex0);
    // G_CPR.optimize();
      G_CPR = G_CPR.JacFun();
    // G_CPR.optimize();
      H_CPR = G_CPR.JacFun();
    // H_CPR.optimize();
    // G_CPR = G_CPR.atomic();
    // H_CPR = H_CPR.atomic();
    }
    
    F_SPR = {calc};
    if(DT == 0){
      G_SPR = TMBad::ADFun<>(TMBad::StdWrap<S_Functor_SPR<Type>,vector<TMBad::ad_aug> >(F_SPR), Ex0);
      // G_SPR.optimize();
      G_SPR = G_SPR.JacFun();
      // G_SPR.optimize();
      H_SPR = G_SPR.JacFun();
      // H_SPR.optimize();
      // G_SPR = G_SPR.atomic();
      // H_SPR = H_SPR.atomic();
    }
    return;
  }

  StochasticWorkerYearOutput<Type> updateNext(Type logFbar, vector<Type> Ex, matrix<Type> Sx, int y){
    vector<Type> pEx;
    vector<vector<Type> > gx;
    vector<matrix<Type> > hx;
    if(DT == 0){
      int n = Ex.size();
      pEx = F(Ex);
      // vector<Type> pEx = calc.predN(Ex);
	
      matrix<Type> gx0 = asMatrix((vector<Type>)G(Ex),n,n); // Transposed Jacobian -> each column is a gradient	
      matrix<Type> hx0 = asMatrix((vector<Type>)H(Ex),n,n*n); // Each column is a gradient of an entry in the Jacobian -> first n columns is the first hessian (transposed), next n are the second hessian ...
      vector<vector<Type> > gtmp(n);
      vector<matrix<Type> > htmp(n);
      for(int i = 0; i < n; ++i){
	gtmp(i) = gx0.col(i);
	htmp(i) = hx0.block(0,n*i,n,n);
      }
      gx = gtmp;
      hx = htmp;      
    }else if(DT == 1){
      pEx = F(Ex);
      gx = GradientList(F, Ex);
      hx = HessianList(F, Ex);
    }else if(DT == 2){
      DerivOutput<Type> dx = getDerivOutput(F, Ex);
      pEx = dx.value;
      gx = dx.gradientList;
      hx = dx.hessianList;
    }
    matrix<Type> Vx(Ex.size(), Ex.size());
    Vx.setZero();
    for(int i = 0; i < Ex.size(); ++i){
      vector<Type> gi = gx(i);	
      matrix<Type> Hi = hx(i);
      matrix<Type> m1 = Hi * Sx;
      pEx(i) += 0.5 * matrix_trace(m1);
      for(int j = 0; j <= i; ++j){
	vector<Type> gj = gx(j);
	matrix<Type> Hj = hx(j);
	matrix<Type> m2 = Hj * Sx;
	matrix<Type> m3 = m1 * m2;
	Type tt = matrix_trace(m3);
	Vx(i,j) = (Type)(gi * (Sx * gj)).sum() + 0.5 * tt;
	Vx(j,i) = Vx(i,j);
      }
    }
    vector<Type> ExOut = pEx;
    matrix<Type> SxOut = Vx + nvar;
    StochasticWorkerYearOutput<Type> res = {ExOut, SxOut};
    return res;
  }
 
  
  PERREC_t<Type> summary_along_years(Type logFbar, int outType, int Ntail, Type q){
    vector<Type> Ex = Ex0;
    Ex(0) = logFbar;
    matrix<Type> Sx = Sx0;
    Type logSe = R_NegInf;
    Type logYe = R_NegInf;
    Type logRe = R_NegInf;
    Type logSPR = R_NegInf;
    Type logYPR = R_NegInf;
    
    for(int y = 0; y < nYears; ++y){
      StochasticWorkerYearOutput<Type> r0 = updateNext(logFbar, Ex, Sx, y);
      Ex = r0.Ex;
      Sx = r0.Sx;
      // Add to reporting
      if(y + 1 > nYears - Ntail){
	// Derived values
	//// SSB
	Type val_s = 0;
	vector<Type> g_s;
	matrix<Type> h_s;
	if(DT==0){
	  g_s = (vector<Type>)G_S(Ex);
	  h_s = asMatrix((vector<Type>)H_S(Ex),Ex.size(),Ex.size());
	  val_s = F_S(Ex)(0);
	}else if(DT==1){
	  g_s = Gradient(F_S, Ex, 0);
	  h_s = Hessian(F_S, Ex, 0);
	  val_s = F_S(Ex)(0);
	}else if(DT==2){
	  DerivOutput<Type> dx_S = getDerivOutput(F_S, Ex);
	  g_s = dx_S.gradientList(0);
	  h_s = dx_S.hessianList(0);
	  val_s = dx_S.value(0);
	}
	matrix<Type> m_s = h_s * Sx;
	matrix<Type> m_s2 = m_s * m_s;
	Type E_logSe = val_s + 0.5 * matrix_trace(m_s);
	Type V_logSe = (g_s * (Sx * g_s)).sum() + 0.5 * matrix_trace(m_s2);
	//// Catch
	Type val_c = 0;
	vector<Type> g_c;
	matrix<Type> h_c;
	if(DT==0){
	  g_c = (vector<Type>)G_C(Ex);
	  h_c = asMatrix((vector<Type>)H_C(Ex),Ex.size(),Ex.size());
	  val_c = F_C(Ex)(0);
	}else if(DT==1){
	  g_c = Gradient(F_C, Ex, 0);
	  h_c = Hessian(F_C, Ex, 0);
	  val_c = F_C(Ex)(0);
	}else if(DT==2){
	  DerivOutput<Type> dx_C = getDerivOutput(F_C, Ex);
	  g_c = dx_C.gradientList(0);
	  h_c = dx_C.hessianList(0);
	  val_c = dx_C.value(0);
	}
	matrix<Type> m_c = h_c * Sx;
	matrix<Type> m_c2 = m_c * m_c;
	Type E_logYe = val_c + 0.5 * matrix_trace(m_c);
	Type V_logYe = (g_c * (Sx * g_c)).sum() + 0.5 * matrix_trace(m_c2);

	//// Recruitment
	Type E_logRe = Ex(conf.minAge+1);
	Type V_logRe = Sx(conf.minAge+1,conf.minAge+1);

	//// Catch / Recruitment
	Type val_cpr = 0;
	vector<Type> g_cpr;
	matrix<Type> h_cpr;
	if(DT==0){
	  g_cpr = (vector<Type>)G_CPR(Ex);
	  h_cpr = asMatrix((vector<Type>)H_CPR(Ex),Ex.size(),Ex.size());
	  val_cpr = F_CPR(Ex)(0);
	}else if(DT==1){
	  g_cpr = Gradient(F_CPR, Ex, 0);
	  h_cpr = Hessian(F_CPR, Ex, 0);
	  val_cpr = F_CPR(Ex)(0);
	}else if(DT==2){
	  DerivOutput<Type> dx_CPR = getDerivOutput(F_CPR, Ex);
	  g_cpr = dx_CPR.gradientList(0);
	  h_cpr = dx_CPR.hessianList(0);
	  val_cpr = dx_CPR.value(0);
	}
	matrix<Type> m_cpr = h_cpr * Sx;
	matrix<Type> m_cpr2 = m_cpr * m_cpr;
	Type E_logYPR = val_cpr + 0.5 * matrix_trace(m_cpr);
	Type V_logYPR = (g_cpr * (Sx * g_cpr)).sum() + 0.5 * matrix_trace(m_cpr2);

	//// SSB / Recruitment
	Type val_spr = 0;
	vector<Type> g_spr;
	matrix<Type> h_spr;
	if(DT==0){
	  g_spr = (vector<Type>)G_SPR(Ex);
	  h_spr = asMatrix((vector<Type>)H_SPR(Ex),Ex.size(),Ex.size());
	  val_spr = F_SPR(Ex)(0);
	}else if(DT==1){
	  g_spr = Gradient(F_SPR, Ex, 0);
	  h_spr = Hessian(F_SPR, Ex, 0);
	  val_spr = F_SPR(Ex)(0);
	}else if(DT==2){
	  DerivOutput<Type> dx_SPR = getDerivOutput(F_SPR, Ex);
	  g_spr = dx_SPR.gradientList(0);
	  h_spr = dx_SPR.hessianList(0);
	  val_spr = dx_SPR.value(0);
	}
	matrix<Type> m_spr = h_spr * Sx;
	matrix<Type> m_spr2 = m_spr * m_spr;
	Type E_logSPR = val_spr + 0.5 * matrix_trace(m_spr);
	Type V_logSPR = (g_spr * (Sx * g_spr)).sum() + 0.5 * matrix_trace(m_spr2);
	if(outType == 0){ // Median
	  logSe = logspace_add(logSe, E_logSe);
	  logYe = logspace_add(logYe, E_logYe);
	  logRe = logspace_add(logRe, E_logRe);
	  logSPR = logspace_add(logSPR, E_logSPR);
	  logYPR = logspace_add(logYPR, E_logYPR);
	}else if(outType == 1){ // Mean
	  logSe = logspace_add(logSe, E_logSe + 0.5 * V_logSe);
	  logYe = logspace_add(logYe, E_logYe + 0.5 * V_logYe);
	  logRe = logspace_add(logRe, E_logRe + 0.5 * V_logRe);
	  logSPR = logspace_add(logSPR, E_logSPR + 0.5 * V_logSPR);
	  logYPR = logspace_add(logYPR, E_logYPR + 0.5 * V_logYPR);
	}else if(outType == 2){ // Mode
	  logSe = logspace_add(logSe, E_logSe - V_logSe);
	  logYe = logspace_add(logYe, E_logYe - V_logYe);
	  logRe = logspace_add(logRe, E_logRe - V_logRe);
	  logSPR = logspace_add(logSPR, E_logSPR - V_logSPR);
	  logYPR = logspace_add(logYPR, E_logYPR - V_logYPR);
	}else if(outType == 3){ // Quantile
	  logSe = logspace_add(logSe, qnorm(q, E_logSe, sqrt(V_logSe)));
	  logYe = logspace_add(logYe, qnorm(q, E_logYe, sqrt(V_logYe)));
	  logRe = logspace_add(logRe, qnorm(q, E_logRe, sqrt(V_logRe)));
	  logSPR = logspace_add(logSPR, qnorm(q, E_logSPR, sqrt(V_logSPR)));
	  logYPR = logspace_add(logYPR, qnorm(q, E_logYPR, sqrt(V_logYPR)));
	}
      }
    }

    //// Life years lost
    array<Type> logF = logFSel;
    logF.col(newDat.natMor.dim(0)-1) += logFbar;
    Type logLifeExpectancy = log(temporaryLifeExpectancy_i(newDat, conf, logF, newDat.natMor.dim(0)-1, conf.minAge, 10 * conf.maxAge) + (Type)conf.minAge + SAM_Zero);

    //// Life expectancy
    Type logYLTF = log(yearsLostFishing_i(newDat, conf, logF, newDat.natMor.dim(0)-1, conf.minAge, conf.maxAge) + SAM_Zero);
      
    //// dSR0
    Type dSR0 = calc.dSR((Type)SAM_NegInf);
   
    PERREC_t<Type> r2  = {logFbar,// - log(Ntail),
     logYPR - log(Ntail),
     logSPR - log(Ntail),
     logSe - log(Ntail),
     logRe - log(Ntail),
     logYe - log(Ntail),
     dSR0,  // Only calculated once for now
     logLifeExpectancy, // Only calculated once for now
     logYLTF, // Only calculated once for now
     R_NaReal,
     R_NaReal
   };
   return r2;
  }

  STOCHASTIC_PERREC_t<Type> operator()(Type logFbar){
    // Prepare
    vector<Type> Ex = Ex0;
    Ex(0) = logFbar;
    matrix<Type> Sx = Sx0;
    // Calculate
    // Loop to the end of time
    // Approximations of mean and covariance from Rego et al (2021) DOI: 10.1002/cnm.3535
    vector<Type> lastEx = Ex;
    // int n = Ex.size();
    for(int y = 0; y < nYears; ++y){
      StochasticWorkerYearOutput<Type> r0 = updateNext(logFbar, Ex, Sx, y);
      lastEx = Ex;
      Ex = r0.Ex;
      Sx = r0.Sx;
    }
    // Derived values
    //// SSB
    Type val_s = 0;
    vector<Type> g_s;
    matrix<Type> h_s;
    if(DT==0){
      g_s = (vector<Type>)G_S(Ex);
       h_s = asMatrix((vector<Type>)H_S(Ex),Ex.size(),Ex.size());
       val_s = F_S(Ex)(0);
    }else if(DT==1){
      g_s = Gradient(F_S, Ex, 0);
      h_s = Hessian(F_S, Ex, 0);
      val_s = F_S(Ex)(0);
    }else if(DT==2){
      DerivOutput<Type> dx_S = getDerivOutput(F_S, Ex);
      g_s = dx_S.gradientList(0);
      h_s = dx_S.hessianList(0);
      val_s = dx_S.value(0);
    }
    matrix<Type> m_s = h_s * Sx;
    matrix<Type> m_s2 = m_s * m_s;
    Type E_logSSB = val_s + 0.5 * matrix_trace(m_s);
    Type V_logSSB = (g_s * (Sx * g_s)).sum() + 0.5 * matrix_trace(m_s2);

    //// Catch
    Type val_c = 0;
    vector<Type> g_c;
    matrix<Type> h_c;
    if(DT==0){
      g_c = (vector<Type>)G_C(Ex);
      h_c = asMatrix((vector<Type>)H_C(Ex),Ex.size(),Ex.size());
      val_c = F_C(Ex)(0);
    }else if(DT==1){
      g_c = Gradient(F_C, Ex, 0);
      h_c = Hessian(F_C, Ex, 0);
      val_c = F_C(Ex)(0);
    }else if(DT==2){
      DerivOutput<Type> dx_C = getDerivOutput(F_C, Ex);
      g_c = dx_C.gradientList(0);
      h_c = dx_C.hessianList(0);
      val_c = dx_C.value(0);
    }
    matrix<Type> m_c = h_c * Sx;
    matrix<Type> m_c2 = m_c * m_c;
    Type E_logYield = val_c + 0.5 * matrix_trace(m_c);
    Type V_logYield = (g_c * (Sx * g_c)).sum() + 0.5 * matrix_trace(m_c2);

    //// Recruitment
    Type E_logRecruit = Ex(conf.minAge+1);
    Type V_logRecruit = Sx(conf.minAge+1,conf.minAge+1);

    //// Catch / Recruitment
    Type val_cpr = 0;
    vector<Type> g_cpr;
    matrix<Type> h_cpr;
    if(DT==0){
      g_cpr = (vector<Type>)G_CPR(Ex);
      h_cpr = asMatrix((vector<Type>)H_CPR(Ex),Ex.size(),Ex.size());
      val_cpr = F_CPR(Ex)(0);
    }else if(DT==1){
      g_cpr = Gradient(F_CPR, Ex, 0);
      h_cpr = Hessian(F_CPR, Ex, 0);
      val_cpr = F_CPR(Ex)(0);
    }else if(DT==2){
      DerivOutput<Type> dx_CPR = getDerivOutput(F_CPR, Ex);
      g_cpr = dx_CPR.gradientList(0);
      h_cpr = dx_CPR.hessianList(0);
      val_cpr = dx_CPR.value(0);
    }
    matrix<Type> m_cpr = h_cpr * Sx;
    matrix<Type> m_cpr2 = m_cpr * m_cpr;
    Type E_logYPR = val_cpr + 0.5 * matrix_trace(m_cpr);
    Type V_logYPR = (g_cpr * (Sx * g_cpr)).sum() + 0.5 * matrix_trace(m_cpr2);

    //// SSB / Recruitment
    Type val_spr = 0;
    vector<Type> g_spr;
    matrix<Type> h_spr;
    if(DT==0){
      g_spr = (vector<Type>)G_SPR(Ex);
      h_spr = asMatrix((vector<Type>)H_SPR(Ex),Ex.size(),Ex.size());
      val_spr = F_SPR(Ex)(0);
    }else if(DT==1){
      g_spr = Gradient(F_SPR, Ex, 0);
      h_spr = Hessian(F_SPR, Ex, 0);
      val_spr = F_SPR(Ex)(0);
    }else if(DT==2){
      DerivOutput<Type> dx_SPR = getDerivOutput(F_SPR, Ex);
      g_spr = dx_SPR.gradientList(0);
      h_spr = dx_SPR.hessianList(0);
      val_spr = dx_SPR.value(0);
    }
    matrix<Type> m_spr = h_spr * Sx;
    matrix<Type> m_spr2 = m_spr * m_spr;
    Type E_logSPR = val_spr + 0.5 * matrix_trace(m_spr);
    Type V_logSPR = (g_spr * (Sx * g_spr)).sum() + 0.5 * matrix_trace(m_spr2);

    //// Life years lost
    array<Type> logF = logFSel;
    logF.col(newDat.natMor.dim(0)-1) += logFbar;
    Type E_logLifeExpectancy = log(temporaryLifeExpectancy_i(newDat, conf, logF, newDat.natMor.dim(0)-1, conf.minAge, 10 * conf.maxAge) + (Type)conf.minAge + SAM_Zero);

    //// Life expectancy
    Type E_logYLTF = log(yearsLostFishing_i(newDat, conf, logF, newDat.natMor.dim(0)-1, conf.minAge, conf.maxAge) + SAM_Zero);

    //// dSR0
    Type dSR0 = calc.dSR((Type)SAM_NegInf);
    
  // Collect result
  int nAge = conf.maxAge - conf.minAge + 1;
  STOCHASTIC_PERREC_t<Type> res = {
    logFbar, // Type logFbar;

    E_logYPR,// Type E_logYPR;
    E_logSPR,// Type E_logSPR;
    E_logSSB,// Type E_logSe;
    E_logRecruit,// Type E_logRe;
    E_logYield,// Type E_logYe;
    E_logLifeExpectancy,// Type E_logLifeExpectancy;
    E_logYLTF,// Type E_logYearsLost;

    (Type)Sx(0,0),// Type V_logFbar;
    V_logYPR,// Type V_logYPR;
    V_logSPR,// Type V_logSPR;
    V_logSSB,// Type V_logSe;
    V_logRecruit,// Type V_logRe;
    V_logYield,// Type V_logYe;
    (Type)0.0,// Type V_logLifeExpectancy;
    (Type)0.0,// Type V_logYearsLost;

    (vector<Type>)Ex.segment(conf.minAge+1, nAge),// vector<Type> E_logN;
    (matrix<Type>)Sx.block(conf.minAge+1,conf.minAge+1,nAge,nAge),// matrix<Type> V_logN;

    Ex - lastEx,
    dSR0
  };

  //Return result
  return res;
    
  }
  
}

       )


//////////////////////////////////////////////////////////////////////////////////////////
// For optimizing/reporting natural scale MEDIAN quantities

HEADER(
template<class Type>
struct EquilibriumRecycler_Stochastic_Median : EquilibriumRecycler<Type> {

  EquilibriumRecycler_Stochastic_Worker<Type> wrk;
  
  EquilibriumRecycler_Stochastic_Median();
  EquilibriumRecycler_Stochastic_Median(Type logFbar0,
					dataSet<Type> dat,
					confSet conf,
					paraSet<Type> par,
					vector<Type> logSel,
					vector<int> aveYears,
					vector<Type> logN0,
					int nYears,
					int CT,
					int DT);
  PERREC_t<Type> operator()(Type logFbar);
};
       )

SOURCE(
template<class Type>
EquilibriumRecycler_Stochastic_Median<Type>::EquilibriumRecycler_Stochastic_Median() : EquilibriumRecycler<Type>(), wrk() {};
       )

SOURCE(
       template<class Type>
       EquilibriumRecycler_Stochastic_Median<Type>::EquilibriumRecycler_Stochastic_Median(Type logFbar0,
											  dataSet<Type> dat,
											  confSet conf,
											  paraSet<Type> par,
											  vector<Type> logSel,
											  vector<int> aveYears,
											  vector<Type> logN0,
											  int nYears,
											  int CT,
											  int DT) : EquilibriumRecycler<Type>(), wrk(logFbar0, dat,conf,par,logSel,aveYears,logN0,nYears,CT,DT) {};
       )


SOURCE(
template<class Type>
PERREC_t<Type> EquilibriumRecycler_Stochastic_Median<Type>::operator()(Type logFbar){
  STOCHASTIC_PERREC_t<Type> r = wrk(logFbar);
  PERREC_t<Type> r2  = {r.E_logFbar,
	    r.E_logYPR,
	    r.E_logSPR,
	    r.E_logSe,
	    r.E_logRe,
	    r.E_logYe,
	    r.dSR0,
	    r.E_logLifeExpectancy,
	    r.E_logYearsLost,
	    R_NaReal,
	    R_NaReal
	  };
  return r2;
  }
       )


SAM_SPECIALIZATION(struct EquilibriumRecycler_Stochastic_Median<double>);
SAM_SPECIALIZATION(struct EquilibriumRecycler_Stochastic_Median<TMBad::ad_aug>);


//////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////
// For optimizing/reporting natural scale MEAN quantities

HEADER(
template<class Type>
struct EquilibriumRecycler_Stochastic_Mean : EquilibriumRecycler<Type> {

  EquilibriumRecycler_Stochastic_Worker<Type> wrk;
  
  EquilibriumRecycler_Stochastic_Mean();
  EquilibriumRecycler_Stochastic_Mean(Type logFbar0,
					dataSet<Type> dat,
					confSet conf,
					paraSet<Type> par,
					vector<Type> logSel,
					vector<int> aveYears,
					vector<Type> logN0,
					int nYears,
				      int CT,
				      int DT);
  PERREC_t<Type> operator()(Type logFbar);
};
       )

SOURCE(
template<class Type>
EquilibriumRecycler_Stochastic_Mean<Type>::EquilibriumRecycler_Stochastic_Mean() : EquilibriumRecycler<Type>(), wrk() {};
       )

SOURCE(
       template<class Type>
       EquilibriumRecycler_Stochastic_Mean<Type>::EquilibriumRecycler_Stochastic_Mean(Type logFbar0,
											  dataSet<Type> dat,
											  confSet conf,
											  paraSet<Type> par,
											  vector<Type> logSel,
											  vector<int> aveYears,
											  vector<Type> logN0,
											  int nYears,
										      int CT,
										      int DT) : EquilibriumRecycler<Type>(), wrk(logFbar0, dat,conf,par,logSel,aveYears,logN0,nYears,CT,DT) {};
       )


SOURCE(
template<class Type>
PERREC_t<Type> EquilibriumRecycler_Stochastic_Mean<Type>::operator()(Type logFbar){
  STOCHASTIC_PERREC_t<Type> r = wrk(logFbar);
  PERREC_t<Type> r2  = {r.E_logFbar + 0.5 * r.V_logFbar,
	    r.E_logYPR + 0.5 * r.V_logYPR,
	    r.E_logSPR + 0.5 * r.V_logSPR,
	    r.E_logSe + 0.5 * r.V_logSe,
	    r.E_logRe + 0.5 * r.V_logRe,
	    r.E_logYe + 0.5 * r.V_logYe,
	    r.dSR0,
	    r.E_logLifeExpectancy + 0.5 * r.V_logLifeExpectancy,
	    r.E_logYearsLost + 0.5 * r.V_logYearsLost,
	    R_NaReal,
	    R_NaReal
	  };
  return r2;
  }
       )


SAM_SPECIALIZATION(struct EquilibriumRecycler_Stochastic_Mean<double>);
SAM_SPECIALIZATION(struct EquilibriumRecycler_Stochastic_Mean<TMBad::ad_aug>);


//////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////
// For optimizing/reporting natural scale MODE quantities

HEADER(
template<class Type>
struct EquilibriumRecycler_Stochastic_Mode : EquilibriumRecycler<Type> {

  EquilibriumRecycler_Stochastic_Worker<Type> wrk;
  
  EquilibriumRecycler_Stochastic_Mode();
  EquilibriumRecycler_Stochastic_Mode(Type logFbar0,
					dataSet<Type> dat,
					confSet conf,
					paraSet<Type> par,
					vector<Type> logSel,
					vector<int> aveYears,
					vector<Type> logN0,
					int nYears,
				      int CT,
				      int DT);
  PERREC_t<Type> operator()(Type logFbar);
};
       )

SOURCE(
template<class Type>
EquilibriumRecycler_Stochastic_Mode<Type>::EquilibriumRecycler_Stochastic_Mode() : EquilibriumRecycler<Type>(), wrk() {};
       )

SOURCE(
       template<class Type>
       EquilibriumRecycler_Stochastic_Mode<Type>::EquilibriumRecycler_Stochastic_Mode(Type logFbar0,
											  dataSet<Type> dat,
											  confSet conf,
											  paraSet<Type> par,
											  vector<Type> logSel,
											  vector<int> aveYears,
											  vector<Type> logN0,
											  int nYears,
										      int CT,
										      int DT) : EquilibriumRecycler<Type>(), wrk(logFbar0, dat,conf,par,logSel,aveYears,logN0,nYears,CT,DT) {};
       )


SOURCE(
template<class Type>
PERREC_t<Type> EquilibriumRecycler_Stochastic_Mode<Type>::operator()(Type logFbar){
  STOCHASTIC_PERREC_t<Type> r = wrk(logFbar);
  PERREC_t<Type> r2  = {r.E_logFbar - r.V_logFbar,
	    r.E_logYPR - r.V_logYPR,
	    r.E_logSPR - r.V_logSPR,
	    r.E_logSe - r.V_logSe,
	    r.E_logRe - r.V_logRe,
	    r.E_logYe - r.V_logYe,
	    r.dSR0,
	    r.E_logLifeExpectancy - r.V_logLifeExpectancy,
	    r.E_logYearsLost - r.V_logYearsLost,
	    R_NaReal,
	    R_NaReal
	  };
  return r2;
  }
       )


SAM_SPECIALIZATION(struct EquilibriumRecycler_Stochastic_Mode<double>);
SAM_SPECIALIZATION(struct EquilibriumRecycler_Stochastic_Mode<TMBad::ad_aug>);


//////////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////////////
// For optimizing/reporting natural scale QUANTILE quantities

HEADER(
template<class Type>
struct EquilibriumRecycler_Stochastic_Quantile : EquilibriumRecycler<Type> {

  EquilibriumRecycler_Stochastic_Worker<Type> wrk;
  Type q;
  
  EquilibriumRecycler_Stochastic_Quantile();
  EquilibriumRecycler_Stochastic_Quantile(Type logFbar0,
					dataSet<Type> dat,
					confSet conf,
					paraSet<Type> par,
					vector<Type> logSel,
					vector<int> aveYears,
					vector<Type> logN0,
					int nYears,
					  int CT,
					  int DT,
					  Type q);
  PERREC_t<Type> operator()(Type logFbar);
};
       )

SOURCE(
template<class Type>
EquilibriumRecycler_Stochastic_Quantile<Type>::EquilibriumRecycler_Stochastic_Quantile() : EquilibriumRecycler<Type>(), wrk(), q() {};
       )

SOURCE(
       template<class Type>
       EquilibriumRecycler_Stochastic_Quantile<Type>::EquilibriumRecycler_Stochastic_Quantile(Type logFbar0,
											  dataSet<Type> dat,
											  confSet conf,
											  paraSet<Type> par,
											  vector<Type> logSel,
											  vector<int> aveYears,
											  vector<Type> logN0,
											  int nYears,
											      int CT,
											      int DT,
											      Type q) : EquilibriumRecycler<Type>(), wrk(logFbar0, dat,conf,par,logSel,aveYears,logN0,nYears,CT,DT), q(q) {};
       )


SOURCE(
template<class Type>
PERREC_t<Type> EquilibriumRecycler_Stochastic_Quantile<Type>::operator()(Type logFbar){
  STOCHASTIC_PERREC_t<Type> r = wrk(logFbar);
  PERREC_t<Type> r2  = {qnorm(q,r.E_logFbar,sqrt(r.V_logFbar)),
	    qnorm(q,r.E_logYPR, sqrt(r.V_logYPR)),
	    qnorm(q,r.E_logSPR, sqrt(r.V_logSPR)),
	    qnorm(q,r.E_logSe, sqrt(r.V_logSe)),
	    qnorm(q,r.E_logRe, sqrt(r.V_logRe)),
	    qnorm(q,r.E_logYe, sqrt(r.V_logYe)),
	    r.dSR0,
	    qnorm(q,r.E_logLifeExpectancy, sqrt(r.V_logLifeExpectancy)),
	    qnorm(q,r.E_logYearsLost, sqrt(r.V_logYearsLost)),
	    R_NaReal,
	    R_NaReal
	  };
  return r2;
  }
       )


SAM_SPECIALIZATION(struct EquilibriumRecycler_Stochastic_Quantile<double>);
SAM_SPECIALIZATION(struct EquilibriumRecycler_Stochastic_Quantile<TMBad::ad_aug>);


//////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Stochastic //////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

template<class Type>
STOCHASTIC_PERREC_t<Type> perRecruit_S(const Type& logFbar, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, vector<Type>& logSel, vector<int>& aveYears, vector<Type>& logN0, int nYears DEFARG(= 300), int CT DEFARG(= 0), int DT DEFARG(= 0))SOURCE({

    EquilibriumRecycler_Stochastic_Worker<Type> wrk(logFbar,dat,conf,par,logSel,aveYears,logN0,nYears,CT,DT);
    return wrk(logFbar);

  }
  )


  SAM_SPECIALIZATION(STOCHASTIC_PERREC_t<double> perRecruit_S(const double&, dataSet<double>&, confSet&, paraSet<double>&, vector<double>&, vector<int>&, vector<double>&, int, int, int));
SAM_SPECIALIZATION(STOCHASTIC_PERREC_t<TMBad::ad_aug> perRecruit_S(const TMBad::ad_aug&, dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, vector<TMBad::ad_aug>&, vector<int>&, vector<TMBad::ad_aug>&, int, int, int));
