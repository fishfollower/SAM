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
  MortalitySet<Type> mort;
  // Type logFbar;
  Recruitment<Type> recruit;  
  // vector<Type> logSel;

  stochasticCalculator() : dat(), conf(), par(),logFSel(), logFseason(), mort(), recruit() {};
  
  stochasticCalculator( dataSet<Type> dat,
			confSet conf,
			paraSet<Type> par,
			array<Type> logFSel_,
			array<Type> logFseason_
			// Type logFbar,
			// vector<Type>& logSel
			) :
    dat(dat), conf(conf), par(par), logFSel(logFSel_, logFSel_.dim), logFseason(logFseason_,logFseason_.dim), mort(dat,conf,par,logFSel,logFseason) {
    recruit = makeRecruitmentFunction(conf,par);
  }
  
  template<class T>
  inline stochasticCalculator(const stochasticCalculator<T> &x) :
    dat(x.dat),
    conf(x.conf),
    par(x.par),
    logFSel(x.logFSel,x.logFSel.dim),
    logFseason(x.logFseason,x.logFseason.dim),
    mort(dat,conf,par,logFSel,logFseason)
  {
    recruit = makeRecruitmentFunction(conf,par);
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
	  Type tmp = logN(a) + mort.logFleetSurvival_before(ax,y,f) + log(mort.fleetCumulativeIncidence(ax,y,f)) + log(dat.catchMeanWeight(y, ax, f));
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
template<class Type>
struct EquilibriumRecycler_Stochastic_Worker {

  Type logFbar0;
  confSet conf;
  paraSet<Type> par;
  vector<Type> logSel;
  vector<Type> logN0;
  int nYears;
  int CT;

  
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
    newDat(),
    logitFseason(),
    calc(),
    Ex0(),
    Sx0(),
    nvar() {}
    
  template<class T>
  EquilibriumRecycler_Stochastic_Worker(T logFbar0_,
					dataSet<T> dat,
					confSet conf_,
					paraSet<T> par_,
					vector<T> logSel_,
					vector<int> aveYears,
					vector<T> logN0_,
					int nYears_,
					int CT_) :
    logFbar0(logFbar0_),
    conf(conf_),
    par(par_),
    logSel(logSel_),
    logN0(logN0_),
    nYears(nYears_),
    CT(CT_),
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

   calc = stochasticCalculator<Type>(newDat, conf, par, logFSel, logitFseason);

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
    G = TMBad::ADFun<>(TMBad::StdWrap<S_Functor_N<Type>,vector<TMBad::ad_aug> >(F), Ex0);
    // G.optimize();
    G = G.JacFun();
    // G.optimize();
    H = G.JacFun();
    // H.optimize();
    // G = G.atomic();
    // H = H.atomic();
    
    F_S = {calc};
    G_S = TMBad::ADFun<> (TMBad::StdWrap<S_Functor_S<Type>,vector<TMBad::ad_aug> >(F_S), Ex0);
    // G_S.optimize();
    G_S = G_S.JacFun();
    // G_S.optimize();
    H_S = G_S.JacFun();
    // H_S.optimize();
    // G_S = G_S.atomic();
    // H_S = H_S.atomic();
  
    F_C = {calc};
    G_C = TMBad::ADFun<>(TMBad::StdWrap<S_Functor_C<Type>,vector<TMBad::ad_aug> >(F_C), Ex0);
    // G_C.optimize();
    G_C = G_C.JacFun();
    // G_C.optimize();
    H_C = G_C.JacFun();
    // H_C.optimize();
    // G_C = G_C.atomic();
    // H_C = H_C.atomic();
    
    F_CPR = {calc};
    G_CPR = TMBad::ADFun<>(TMBad::StdWrap<S_Functor_CPR<Type>,vector<TMBad::ad_aug> >(F_CPR), Ex0);
    // G_CPR.optimize();
    G_CPR = G_CPR.JacFun();
    // G_CPR.optimize();
    H_CPR = G_CPR.JacFun();
    // H_CPR.optimize();
    // G_CPR = G_CPR.atomic();
    // H_CPR = H_CPR.atomic();

    F_SPR = {calc};
    G_SPR = TMBad::ADFun<>(TMBad::StdWrap<S_Functor_SPR<Type>,vector<TMBad::ad_aug> >(F_SPR), Ex0);
    // G_SPR.optimize();
    G_SPR = G_SPR.JacFun();
    // G_SPR.optimize();
    H_SPR = G_SPR.JacFun();
    // H_SPR.optimize();
    // G_SPR = G_SPR.atomic();
    // H_SPR = H_SPR.atomic();

    return;
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
    int n = Ex.size();
    for(int y = 0; y < nYears; ++y){
      vector<Type> pEx = F(Ex);
      // vector<Type> pEx = calc.predN(Ex);
      matrix<Type> gx = asMatrix((vector<Type>)G(Ex),n,n); // Transposed Jacobian -> each column is a gradient
      matrix<Type> hx = asMatrix((vector<Type>)H(Ex),n,n*n); // Each column is a gradient of an entry in the Jacobian -> first n columns is the first hessian (transposed), next n are the second hessian ...
      matrix<Type> Vx(Ex.size(), Ex.size());
      Vx.setZero();
      for(int i = 0; i < Ex.size(); ++i){
	// S_Functor_N_i<Type> Fi = {calc, i};
	vector<Type> gi = gx.col(i);
	// vector<Type> gi0 = autodiff::gradient(Fi,Ex);    
	matrix<Type> Hi = hx.block(0,n*i,n,n);
	// matrix<Type> Hi0 = autodiff::hessian(Fi,Ex);
	// Rcout << y << " - " << i << ": " << (gi0-gi).abs().sum() << "  "  << (Hi0-Hi).array().abs().sum() << "\n";
	matrix<Type> m1 = Hi * Sx;
	pEx(i) += 0.5 * matrix_trace(m1);
	for(int j = 0; j <= i; ++j){
	  // S_Functor_N_i<Type> Fj = {calc, j};
	  vector<Type> gj = gx.col(j);
	  // vector<Type> gj0 = autodiff::gradient(Fj,Ex);
	  matrix<Type> Hj = hx.block(0,n*j,n,n);
	  // matrix<Type> Hj0 = autodiff::hessian(Fj,Ex);
	  // Rcout << "\t\t" << "\t" << j << ": " << (gj0-gj).abs().sum() << "  "  << (Hj0-Hj).array().abs().sum() << "\n";
	  matrix<Type> m2 = Hj * Sx;
	  matrix<Type> m3 = m1 * m2;
	  Type tt = matrix_trace(m3);
	  Vx(i,j) = (Type)(gi * (Sx * gj)).sum() + 0.5 * tt;
	  Vx(j,i) = Vx(i,j);
	}
      }
      lastEx = Ex;
      Ex = pEx;
      Sx = Vx + nvar;
    }
  // Derived values
  //// SSB
  vector<Type> g_s = G_S(Ex);
  matrix<Type> h_s = asMatrix((vector<Type>)H_S(Ex),Ex.size(),Ex.size());
  matrix<Type> m_s = h_s * Sx;
  matrix<Type> m_s2 = m_s * m_s;
  Type E_logSSB = F_S(Ex)(0) + 0.5 * matrix_trace(m_s);
  Type V_logSSB = (g_s * (Sx * g_s)).sum() + 0.5 * matrix_trace(m_s2);

  //// Catch
  vector<Type> g_c = G_C(Ex);
  matrix<Type> h_c = asMatrix((vector<Type>)H_C(Ex),n,n);
  matrix<Type> m_c = h_c * Sx;
  matrix<Type> m_c2 = m_c * m_c;
  Type E_logYield = F_C(Ex)(0) + 0.5 * matrix_trace(m_c);
  Type V_logYield = (g_c * (Sx * g_c)).sum() + 0.5 * matrix_trace(m_c2);

  //// Recruitment
  Type E_logRecruit = Ex(conf.minAge);
  Type V_logRecruit = Sx(conf.minAge,conf.minAge);

  //// Catch / Recruitment
  vector<Type> g_cpr = G_CPR(Ex);
  matrix<Type> h_cpr = asMatrix((vector<Type>)H_CPR(Ex),Ex.size(),Ex.size());
  matrix<Type> m_cpr = h_cpr * Sx;
  matrix<Type> m_cpr2 = m_cpr * m_cpr;
  Type E_logYPR = F_CPR(Ex)(0) + 0.5 * matrix_trace(m_cpr);
  Type V_logYPR = (g_cpr * (Sx * g_cpr)).sum() + 0.5 * matrix_trace(m_cpr2);

  //// SSB / Recruitment
  vector<Type> g_spr = G_SPR(Ex);
  matrix<Type> h_spr = asMatrix((vector<Type>)H_SPR(Ex),Ex.size(),Ex.size());
  matrix<Type> m_spr = h_spr * Sx;
  matrix<Type> m_spr2 = m_spr * m_spr;
  Type E_logSPR = F_SPR(Ex)(0) + 0.5 * matrix_trace(m_spr);
  Type V_logSPR = (g_spr * (Sx * g_spr)).sum() + 0.5 * matrix_trace(m_spr2);

  //// Life years lost
  array<Type> logF = logFSel;
  logF.col(newDat.natMor.dim(0)-1) += logFbar;
  Type E_logLifeExpectancy = log(temporaryLifeExpectancy_i(newDat, conf, logF, newDat.natMor.dim(0)-1, conf.minAge, 10 * conf.maxAge) + (Type)conf.minAge + SAM_Zero);

  //// Life expectancy
  Type E_logYLTF = log(yearsLostFishing_i(newDat, conf, logF, newDat.natMor.dim(0)-1, conf.minAge, conf.maxAge) + SAM_Zero);

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

    Ex - lastEx
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
					int CT);
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
											  int CT) : EquilibriumRecycler<Type>(), wrk(logFbar0, dat,conf,par,logSel,aveYears,logN0,nYears,CT) {};
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
	    R_NaReal,
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
					int CT);
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
											  int CT) : EquilibriumRecycler<Type>(), wrk(logFbar0, dat,conf,par,logSel,aveYears,logN0,nYears,CT) {};
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
	    R_NaReal,
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
					int CT);
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
											  int CT) : EquilibriumRecycler<Type>(), wrk(logFbar0, dat,conf,par,logSel,aveYears,logN0,nYears,CT) {};
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
	    R_NaReal,
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
											      Type q) : EquilibriumRecycler<Type>(), wrk(logFbar0, dat,conf,par,logSel,aveYears,logN0,nYears,CT), q(q) {};
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
	    R_NaReal,
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
STOCHASTIC_PERREC_t<Type> perRecruit_S(const Type& logFbar, dataSet<Type>& dat, confSet& conf, paraSet<Type>& par, vector<Type>& logSel, vector<int>& aveYears, vector<Type>& logN0, int nYears DEFARG(= 300), int CT DEFARG(= 0))SOURCE({

    EquilibriumRecycler_Stochastic_Worker<Type> wrk(logFbar,dat,conf,par,logSel,aveYears,logN0,nYears,CT);
    return wrk(logFbar);

  }
  )


SAM_SPECIALIZATION(STOCHASTIC_PERREC_t<double> perRecruit_S(const double&, dataSet<double>&, confSet&, paraSet<double>&, vector<double>&, vector<int>&, vector<double>&, int, int));
SAM_SPECIALIZATION(STOCHASTIC_PERREC_t<TMBad::ad_aug> perRecruit_S(const TMBad::ad_aug&, dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, vector<TMBad::ad_aug>&, vector<int>&, vector<TMBad::ad_aug>&, int, int));
