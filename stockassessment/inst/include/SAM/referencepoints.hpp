SAM_DEPENDS(define)
SAM_DEPENDS(recruitment)
SAM_DEPENDS(refpointset)
SAM_DEPENDS(derived)
SAM_DEPENDS(equilibrium)
SAM_DEPENDS(newton)

#include <memory>


///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Deterministic /////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef WITH_SAM_LIB
namespace referencepoints_helper {

  typedef TMBad::ad_aug ad;

// Deterministic Reference point functors
// template<class Type>
struct RPD_Base : NewtonFunctor {

  dataSet<ad> dat;
  confSet conf;
  paraSet<ad> par;
  referencepointSet<ad> rp;

  RPD_Base() : dat(), conf(), par(), rp() {};
  
  RPD_Base(const dataSet<ad>& dat_,
	   const confSet& conf_,
	   const paraSet<ad>& par_,
	   const referencepointSet<ad>& rp_) :
    dat(dat_),
    conf(conf_),
    par(par_),
    rp(rp_){}

  template<class T>
  RPD_Base(const dataSet<T>& dat_,
	   const confSet& conf_,
	   const paraSet<T>& par_,
	   const referencepointSet<T>& rp_) :
    dat(dat_),
    conf(conf_),
    par(par_),
    rp(rp_){}

  PERREC_t<ad> getPerRec(const ad& logFbar){
    vector<ad> ls = rp.getLogSelectivity();
    PERREC_t<ad> r =  perRecruit_D(logFbar, dat, conf, par, ls, rp.aveYears, rp.nYears, rp.catchType);
    return r;
  }

  template<class T>
  vector<T> par2logF(const vector<T>& x){
    return x;
  }

  virtual ad operator()(const vector<ad>& logFbar) = 0;
  
};


template<class Type>
struct RefPointD_Base {
  
  virtual PERREC_t<Type> getPerRecruit(Type logFbar) = 0;
  virtual vector<Type> optimize(vector<Type> logF0) = 0;
  virtual vector<Type> par2logF(const vector<Type>& x) = 0;
  
};


template<class Type>
class Referencepoint_D {
  std::shared_ptr<RefPointD_Base<Type> > ptr;
  const char* name;
  int id;
  vector<Type> logF;
public:
  Referencepoint_D() : ptr(), name(), id(), logF() {}

  Referencepoint_D(const char* name_,
		   int id_,
		   vector<Type> logF0,
		   RefPointD_Base<Type>* p
		   ) :
    ptr(p), name(name_), id(id_), logF(logF0.size()) {
    logF.setConstant(R_NegInf);
    if(ptr != nullptr){
      logF = ptr->optimize(logF0);
    }else{
      Rf_error("Referencepoint_D: Pointer not initialized.");
    }
  }
  
  ~Referencepoint_D() { ptr.reset(); }

  const char* makeName(const char* rpnm){
    std::string nam("");
    nam.append("referencepoint_");
    std::ostringstream s;
    s << id;
    nam.append(s.str());
    nam.append("_");
    nam.append(name);
    nam.append("_");
    nam.append(rpnm);
    return strdup(nam.data());
  }
  
  void adreport(objective_function<Type> *of){
    vector<Type> logFout = ptr->par2logF(logF);
    int n = logFout.size();
    vector<Type> logYPR(n);
    vector<Type> logSPR(n);
    vector<Type> logSe(n);
    vector<Type> logRe(n);
    vector<Type> logYe(n);
    vector<Type> logLifeExpectancy(n);
    vector<Type> logYearsLost(n);
    vector<Type> logDiscYPR(n);
    vector<Type> logDiscYe(n);
    for(int i = 0; i < n; ++i){
      PERREC_t<Type> pr = ptr->getPerRecruit(logFout(i));
      logYPR(i) = pr.logYPR;
      logSPR(i) = pr.logSPR;
      logSe(i) = pr.logSe;
      logRe(i) = pr.logRe;
      logYe(i) = pr.logYe;
      logLifeExpectancy(i) = pr.logLifeExpectancy;
      logYearsLost(i) = pr.logYearsLost;
      logDiscYPR(i) = pr.logDiscYPR;
      logDiscYe(i) = pr.logDiscYe;
    }    
    of->reportvector.push(logFout, makeName("logF"));
    of->reportvector.push(logYPR, makeName("logYPR"));
    of->reportvector.push(logSPR, makeName("logSPR"));
    of->reportvector.push(logSe, makeName("logSe"));
    of->reportvector.push(logRe, makeName("logRe"));
    of->reportvector.push(logYe, makeName("logYe"));
    of->reportvector.push(logLifeExpectancy, makeName("logLifeExpectancy"));
    of->reportvector.push(logYearsLost, makeName("logYearsLost"));
    of->reportvector.push(logDiscYPR, makeName("logDiscYPR"));
    of->reportvector.push(logDiscYe, makeName("logDiscYe"));
    return;
  }

  void report(vector<Type> x, const char* name, objective_function<Type> *of){
    if(! (isDouble<Type>::value && of -> current_parallel_region<0 ) )
      return;
    SEXP Xval;
    PROTECT( Xval = asSEXP(x) );
    Rf_defineVar(Rf_install(name), Xval, of -> report);
    UNPROTECT(1);      
    return;
  }
  void report(objective_function<Type> *of){
    if(! (isDouble<Type>::value && of -> current_parallel_region<0 ) )
      return;
    vector<Type> logFout = ptr->par2logF(logF);
    int n = logFout.size();
    vector<Type> logYPR(n);
    vector<Type> logSPR(n);
    vector<Type> logSe(n);
    vector<Type> logRe(n);
    vector<Type> logYe(n);
    vector<Type> logLifeExpectancy(n);
    vector<Type> logYearsLost(n);
    vector<Type> logDiscYPR(n);
    vector<Type> logDiscYe(n);
    for(int i = 0; i < n; ++i){
      PERREC_t<Type> pr = ptr->getPerRecruit(logFout(i));
      logYPR(i) = pr.logYPR;
      logSPR(i) = pr.logSPR;
      logSe(i) = pr.logSe;
      logRe(i) = pr.logRe;
      logYe(i) = pr.logYe;
      logLifeExpectancy(i) = pr.logLifeExpectancy;
      logYearsLost(i) = pr.logYearsLost;
      logDiscYPR(i) = pr.logDiscYPR;
      logDiscYe(i) = pr.logDiscYe;
      }
    report(logFout, makeName("logF"), of);
    report(logYPR, makeName("logYPR"), of);
    report(logSPR, makeName("logSPR"), of);
    report(logSe, makeName("logSe"), of);
    report(logRe, makeName("logRe"), of);
    report(logYe, makeName("logYe"), of);
    report(logLifeExpectancy, makeName("logLifeExpectancy"), of);
    report(logYearsLost, makeName("logYearsLost"), of);
    report(logDiscYPR, makeName("logDiscYPR"), of);
    report(logDiscYe, makeName("logDiscYe"), of);
    return;
  }

    
};



///////////////////////////////////////////////////////////////////////////////


  // Enum of implemented reference points
  enum ReferencePointDeterministic {
				    None = -99,
				    FixedF = -1,
				    StatusQuo = 0,
				    MSY = 1,
				    MSYRange = 2,
				    Max = 3,
				    xdYPR = 4,
				    xSPR = 5,
				    xB0 = 6,
				    MYPYLdiv = 7,
				    MYPYLprod = 8,
				    MDY = 9,
				    Crash = 10,
				    Ext = 11,
				    Lim = 12
  };

  
  template<class Type>
  struct RefPointD_Known : RefPointD_Base<Type> {
    dataSet<Type> dat;
    confSet conf;
    paraSet<Type> par;
    referencepointSet<Type> rp;

    RefPointD_Known() : RefPointD_Base<Type>(), dat(), conf(), par(), rp() {};
  
    RefPointD_Known(const dataSet<Type>& dat_,
		    const confSet& conf_,
		    const paraSet<Type>& par_,
		    const referencepointSet<Type>& rp_) :
      RefPointD_Base<Type>(),
      dat(dat_),
      conf(conf_),
      par(par_),
      rp(rp_){}

    PERREC_t<Type> getPerRecruit(Type logFbar){
      vector<Type> ls = rp.getLogSelectivity();
      PERREC_t<Type> r =  perRecruit_D(logFbar, dat, conf, par, ls, rp.aveYears, rp.nYears, rp.catchType);
      return r;
    }

    vector<Type> par2logF(const vector<Type>& x){
      return x;
    }
    virtual vector<Type> optimize(vector<Type> logF0) = 0;
  };

  template<class Type>
  struct RefPointD_Numeric : RefPointD_Base<Type> {

    dataSet<Type> dat;
    confSet conf;
    paraSet<Type> par;
    referencepointSet<Type> rp;

    // Functor f;			// i.e., RPD_MSY. Should be derived from RPD_Base
    std::shared_ptr<RPD_Base> ptr;
    newton::newton_config cfg;
  
    RefPointD_Numeric(const dataSet<Type>& dat_,
		      const confSet& conf_,
		      const paraSet<Type>& par_,
		      const referencepointSet<Type>& rp_,
		      std::shared_ptr<RPD_Base> p_,
		      newton::newton_config cfg_ = newton::newton_config()) :
      RefPointD_Base<Type>(),
      dat(dat_),
      conf(conf_),
      par(par_),
      rp(rp_),
      ptr(p_), cfg(cfg_) {};

  
    PERREC_t<Type> getPerRecruit(Type logFbar){
       vector<Type> ls = rp.getLogSelectivity();
       PERREC_t<Type> r =  perRecruit_D(logFbar, dat, conf, par, ls, rp.aveYears, rp.nYears, rp.catchType);
      return r;
    }

    vector<Type> par2logF(const vector<Type>& x){
      return ptr->par2logF(x);
    };
    vector<Type> optimize(vector<Type> logF0){
      // cfg.simplify = false;	// Needed in older versions of TMB for logspace_add
      std::shared_ptr<NewtonFunctor> pUse = std::static_pointer_cast<NewtonFunctor>(ptr);
      return SAM_Newton(pUse,logF0, cfg);
    }  
  };

  //   using RPD_Base::RPD_Base;	       
#define USING_RPD_BASE_0(NAME)			\
  RPD_##NAME() : RPD_Base() {}			\
  RPD_##NAME(const dataSet<ad>& dat,		\
	     const confSet& conf,		\
	     const paraSet<ad>& par,		\
	     const referencepointSet<ad>& rp) :	\
  RPD_Base(dat,conf,par,rp) {}			\
  template<class T>				\
  RPD_##NAME(const dataSet<T>& dat,		\
	     const confSet& conf,		\
	     const paraSet<T>& par,		\
	     const referencepointSet<T>& rp) :	\
       RPD_Base(dat,conf,par,rp) {}		\
  using RPD_Base::getPerRec

#define USING_RPD_BASE(NAME)			\
  USING_RPD_BASE_0(NAME);			\
  using RPD_Base::par2logF

  //  using RefPointD_Known<Type>::RefPointD_Known;    
#define USING_REFPOINT_KNOWN(NAME)					\
  RefPointD_##NAME(dataSet<Type>& dat,					\
		   confSet& conf,					\
		   paraSet<Type>& par,					\
		   referencepointSet<Type>& rp) :			\
  RefPointD_Known<Type>(dat,conf,par,rp) {};				\
  using RefPointD_Known<Type>::getPerRecruit;				\
  using RefPointD_Known<Type>::par2logF

 

  // Define macro for repeat work
#define MAKE_REFPOINT_D(NAME)						\
  template<class Type>							\
  struct RefPointD_##NAME : RefPointD_Numeric<Type> { \
    RefPointD_##NAME(dataSet<Type>& dat,				\
		     confSet& conf,					\
		     paraSet<Type>& par,				\
		     referencepointSet<Type>& rp,			\
		     newton::newton_config cfg = newton::newton_config()) : \
    RefPointD_Numeric<Type>(dat,conf,par,rp,std::make_shared<RPD_##NAME>(dat,conf,par,rp),cfg) {}; \
  }
  

  //////////////////////////////// Specializations ////////////////////////////////

  //////// FixedF Reference point (Combination of old B0 and Fsequence) ////////

  template<class Type>
  struct RefPointD_FixedF : RefPointD_Known<Type> {

    USING_REFPOINT_KNOWN(FixedF);

    vector<Type> optimize(vector<Type> logF0){
      return logF0;
    }

  };
  
  //////// MSY Reference point ////////


  struct RPD_MSY : RPD_Base {

    USING_RPD_BASE(MSY);
  
    ad operator()(const vector<ad> &x) {
      ad logFbar = x(0);
      PERREC_t<ad> r = getPerRec(logFbar);
      return -r.logYe;
    }
  };

  MAKE_REFPOINT_D(MSY);

  //////// MSYRange Reference point ////////

  // template<class Type>
  struct RPD_MSYrange : RPD_Base {

    USING_RPD_BASE_0(MSYrange);

    // Does not report MSY
    template<class T>
    vector<T> par2logF(const vector<T>& x){
      SAM_ASSERT((x.size()%2)==1,"In reference point MSYrange, length of F must be odd.");
      vector<T>r(x.size()-1);
      r.setConstant(x(0));
      for(int i = 1; i < x.size(); ++i){
	if((i%2)==1){		// lower
	  r(i-1) -= exp(-x(i));
	}else{			// upper
	  r(i-1) += exp(x(i));	  
	}
      }
      return r;
    }
  
    ad operator()(const vector<ad> &x) {
      vector<ad> logFs = par2logF(x);
      SAM_ASSERT( logFs.size() == (2 * this->rp.xVal.size()),"In reference point MSYrange, length of F does not match length of fractions.");
      ad kappa = 0.0;
      // MSY
      PERREC_t<ad> r = getPerRec(x(0));
      kappa += -r.logYe;
      // Ranges
      for(int i = 0; i < logFs.size(); ++i){
	int xvi = (i) / 2;
	PERREC_t<ad> r2 = getPerRec(logFs(i));
	ad tmp = r2.logYe - (log(this->rp.xVal(xvi)) + r.logYe); 
	kappa += tmp * tmp;	
      }
      return kappa;
    }
  };

  MAKE_REFPOINT_D(MSYrange);


  //////// Max Reference point ////////


  // template<class Type>
  struct RPD_Max : RPD_Base {

    USING_RPD_BASE(Max);
  
    ad operator()(const vector<ad> &x) {
      ad logFbar = x(0);
      PERREC_t<ad> r = getPerRec(logFbar);
      return -r.logYPR;
    }
  };

  MAKE_REFPOINT_D(Max);


  //////// %dYPR(0) Reference point (e.g., F~0.1~) ////////


  // template<class Type>
  struct RPD_xdYPR : RPD_Base {

    USING_RPD_BASE(xdYPR);
   
    ad operator()(const vector<ad> &x) {
      SAM_ASSERT(x.size() == this->rp.xVal.size(),"In reference point xdYPR, length of F does not match length of fractions.");
      ad dYPR0 = dYPR(ad(SAM_NegInf), this->dat, this->conf, this->par, this->rp);
      ad kappa = 0.0;
      for(int i = 0; i < x.size(); ++i){
	ad logFbar = x(i);
	ad v = dYPR(logFbar, this->dat,this->conf,this->par,this->rp);
	ad tmp = v - this->rp.xVal(i) * dYPR0;
	kappa += tmp * tmp;
      }
      return kappa;
    }
  };

  MAKE_REFPOINT_D(xdYPR);



  //////// %SPR(0) Reference point (e.g., F~35%SPR~) ////////


  // template<class Type>
  struct RPD_xSPR : RPD_Base {

    USING_RPD_BASE(xSPR);
 
    ad operator()(const vector<ad> &x) {
      SAM_ASSERT(x.size() == this->rp.xVal.size(),"In reference point xSPR, length of F does not match length of fractions.");
      ad SPR0 = spawnersPerRecruit_i(ad(SAM_NegInf), this->dat, this->conf, this->par, this->rp);
      ad kappa = 0.0;
      for(int i = 0; i < x.size(); ++i){
	ad logFbar = x(i);
	ad v = spawnersPerRecruit_i(logFbar, this->dat,this->conf,this->par,this->rp);
	ad tmp = v - this->rp.xVal(i) * SPR0;
	kappa += tmp * tmp;
      }
      return kappa;
    }
  };

  MAKE_REFPOINT_D(xSPR);



  //////// %B(0) Reference point (e.g., F~20%B0~) ////////


  struct RPD_xB0 : RPD_Base {

    USING_RPD_BASE(xB0);
    
    ad operator()(const vector<ad> &x) {
      SAM_ASSERT(x.size() == this->rp.xVal.size(),"In reference point xB0, length of F does not match length of fractions.");
      ad B0 = B0_i(this->dat, this->conf, this->par, this->rp);
      ad kappa = 0.0;
      for(int i = 0; i < x.size(); ++i){
	ad logFbar = x(i);      
	ad v = equilibriumBiomass_i(logFbar, this->dat,this->conf,this->par,this->rp);
	ad tmp = v - this->rp.xVal(i) * B0;
	kappa += tmp * tmp;
      }
      return kappa;
    }
  };

  MAKE_REFPOINT_D(xB0);


  //////// Maximum yield per year lost (v1) ////////


  // template<class Type>
  struct RPD_MYPYLdiv : RPD_Base {

    USING_RPD_BASE(MYPYLdiv);
 
    ad operator()(const vector<ad> &x) {    
      ad logFbar = x(0);
      ad logAgeRange = log((ad)this->conf.maxAge - (ad)this->conf.minAge + (ad)1.0);
      PERREC_t<ad> r = getPerRec(logFbar);
      // Yield / (1 + YearsLost/AgeRange)
      ad tmp = r.logYe - logspace_add_SAM(ad(0.0), r.logYearsLost - logAgeRange);
      return -tmp;
    }
  };

  MAKE_REFPOINT_D(MYPYLdiv);

  //////// Maximum yield per year lost (v2) ////////


  struct RPD_MYPYLprod : RPD_Base {

    USING_RPD_BASE(MYPYLprod);
 
    ad operator()(const vector<ad> &x) {    
      SAM_ASSERT(x.size() == this->rp.xVal.size(),"In reference point MYPYLprod, length of F does not match length of fractions.");
      ad logAgeRange = log((ad)this->conf.maxAge - (ad)this->conf.minAge + (ad)1.0);
	// Yield * (1 - (YearsLost/AgeRange)^d)
	ad kappa = 0.0;    
      for(int i = 0; i < x.size(); ++i){
	ad logFbar = x(i);
	PERREC_t<ad> r = getPerRec(logFbar);
	kappa -= r.logYe + logspace_sub_SAM(ad(0.0), this->rp.xVal(i) * (r.logYearsLost - logAgeRange));
      }
      return kappa;
    }
  };

  MAKE_REFPOINT_D(MYPYLprod);

  //////// Maximum (life year) Discounted Yield ////////


  struct RPD_MDY : RPD_Base {

    USING_RPD_BASE(MDY);
 
    ad operator()(const vector<ad> &x) {    
      ad logFbar = x(0);
      PERREC_t<ad> r = getPerRec(logFbar);
      return -r.logDiscYe;
    }
  };

  MAKE_REFPOINT_D(MDY);


  //////// Crash ////////


  struct RPD_Crash : RPD_Base {

    USING_RPD_BASE(Crash);

    ad operator()(const vector<ad> &x) {    
      ad logFbar = x(0);
      ad logdSR0 = log(getPerRec(SAM_NegInf).dSR0); //log(rec.dSR(Type(SAM_Zero)))
      PERREC_t<ad> r = getPerRec(logFbar);
      ad tmp = logdSR0 - (-r.logSPR);
      return tmp * tmp;
    }
  };

  MAKE_REFPOINT_D(Crash);



  //////// Ext ////////


  struct RPD_Ext : RPD_Base {

    USING_RPD_BASE(Ext);
  
    ad operator()(const vector<ad> &x) {    
      ad logFbar = x(0);
      PERREC_t<ad> r = getPerRec(logFbar);
      return r.logSe * r.logSe;
      // Type h = 0.2;
      // int N = 5;
      // Type v = 0.0;
      // Type d = 0.0;
      // for(int i = 0; i < N; ++i){
      // 	Type xx = - h + 2.0 * h / Type(N-1) * Type(i);
      // 	PERREC_t<Type> r2 = getPerRec(logFbar + xx);
      // 	Type W = 0.75 * (1.0 - (xx / h) * (xx / h));
      // 	// Type ls = TMBad::CondExpLt(r2.logSe, Type(-5.0), Type(-5.0), r2.logSe);
      // 	v += r2.logSe * W;
      // 	d += W;
      // }
      // // PERREC_t<Type> rm = getPerRec(logFbar - 0.01);
      // // PERREC_t<Type> rp = getPerRec(logFbar + 0.01);
      // // Type tmp1 = r.logSe;
      // // Type tmp2 = rm.logSe - rp.logSe;
      // // return tmp1 * tmp1 + tmp2 * tmp2;
      // Type tmp = v / d;
      // return tmp * tmp;
    }
  };

  MAKE_REFPOINT_D(Ext);


  //////// Lim?? ////////

}

#endif

/////////////////////////////// Function to call from main program ////////////////////////////////


template<class Type>
void reportDeterministicReferencePoints(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, array<Type> &logF, Recruitment<Type> &recruit, referencepointList<Type> &referencepoints, objective_function<Type> *of)SOURCE({
  
  if(referencepoints.size() == 0)
    return;

  using namespace referencepoints_helper;
  newton::newton_config cfg = referencepoints.cfg;
  // cfg.simplify = false; 	// Needed for logspace_add
  
  // Report reference points
  for(int i = 0; i < referencepoints.size(); ++i){
    referencepointSet<Type> rps = referencepoints(i);
    rps.setLogSelectivity(logF,conf);
    ReferencePointDeterministic rpt = static_cast<ReferencePointDeterministic>(rps.rpType);
    Referencepoint_D<Type> rp;
  
    if(rpt == ReferencePointDeterministic::None){
      continue;
    }else if(rpt  == ReferencePointDeterministic::FixedF){
      rp = Referencepoint_D<Type>("FixedF",i,rps.logF0, new RefPointD_FixedF<Type>(dat,conf,par,rps));
    }else if(rpt  == ReferencePointDeterministic::StatusQuo){
      if(rps.xVal.size() == 0)
	Rf_error("Referencepoint StatusQuo must have at least one xVal");
      vector<Type> logFsq(rps.xVal.size());
      for(int xi = 0; xi < rps.xVal.size(); ++xi)
	logFsq(xi) = fbari(dat,conf, logF, logF.cols()-1 - CppAD::Integer(rps.xVal(xi)), true);
      rp = Referencepoint_D<Type>("StatusQuo",i,logFsq, new RefPointD_FixedF<Type>(dat,conf,par,rps));
    }else if(rpt  == ReferencePointDeterministic::MSY){
      rp = Referencepoint_D<Type>("MSY",i,rps.logF0, new RefPointD_MSY<Type>(dat,conf,par,rps, cfg));
    }else if(rpt  == ReferencePointDeterministic::MSYRange){
      rp = Referencepoint_D<Type>("MSYRange",i,rps.logF0, new RefPointD_MSYrange<Type>(dat,conf,par,rps, cfg));
    }else if(rpt  == ReferencePointDeterministic::Max){
      rp = Referencepoint_D<Type>("Max",i,rps.logF0, new RefPointD_Max<Type>(dat,conf,par,rps, cfg));
    }else if(rpt  == ReferencePointDeterministic::xdYPR){
      rp = Referencepoint_D<Type>("xdYPR",i,rps.logF0, new RefPointD_xdYPR<Type>(dat,conf,par,rps, cfg));
    }else if(rpt  == ReferencePointDeterministic::xSPR){
      rp = Referencepoint_D<Type>("xSPR",i,rps.logF0, new RefPointD_xSPR<Type>(dat,conf,par,rps, cfg));
    }else if(rpt  == ReferencePointDeterministic::xB0){
      rp = Referencepoint_D<Type>("xB0",i,rps.logF0, new RefPointD_xB0<Type>(dat,conf,par,rps, cfg));
    }else if(rpt  == ReferencePointDeterministic::MYPYLdiv){
      rp = Referencepoint_D<Type>("MYPYLdiv",i,rps.logF0, new RefPointD_MYPYLdiv<Type>(dat,conf,par,rps, cfg));
    }else if(rpt  == ReferencePointDeterministic::MYPYLprod){
      rp = Referencepoint_D<Type>("MYPYLprod",i,rps.logF0, new RefPointD_MYPYLprod<Type>(dat,conf,par,rps, cfg));
    }else if(rpt  == ReferencePointDeterministic::MDY){
      rp = Referencepoint_D<Type>("MDY",i,rps.logF0, new RefPointD_MDY<Type>(dat,conf,par,rps, cfg));
    }else if(rpt  == ReferencePointDeterministic::Crash){
      rp = Referencepoint_D<Type>("Crash",i,rps.logF0, new RefPointD_Crash<Type>(dat,conf,par,rps, cfg));
    }else if(rpt  == ReferencePointDeterministic::Ext){
      rp = Referencepoint_D<Type>("Ext",i,rps.logF0, new RefPointD_Ext<Type>(dat,conf,par,rps, cfg));
    }else if(rpt  == ReferencePointDeterministic::Lim){
      Rf_error("Lim not implemented yet");
    }else{
      Rf_error("Referencepoint type not implemented");
    }
    rp.report(of);
    rp.adreport(of);
  }

  // Report relative reference points
 
    return;
})

SAM_SPECIALIZATION(void reportDeterministicReferencePoints(dataSet<double>&, confSet&, paraSet<double>&, array<double>&, array<double>&, Recruitment<double>&, referencepointList<double>&, objective_function<double>*));
SAM_SPECIALIZATION(void reportDeterministicReferencePoints(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, Recruitment<TMBad::ad_aug>&, referencepointList<TMBad::ad_aug>&, objective_function<TMBad::ad_aug>*));
