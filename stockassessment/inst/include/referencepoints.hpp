#pragma once
#ifndef SAM_REFERENCEPOINT_HPP
#define SAM_REFERENCEPOINT_HPP

#include <memory>

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

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Deterministic /////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

// Deterministic Reference point functors
template<class Type>
struct RPD_Base {

  dataSet<Type> dat;
  confSet conf;
  paraSet<Type> par;
  referencepointSet<Type> rp;

  RPD_Base() = default;
  
  RPD_Base(const dataSet<Type>& dat_,
	   const confSet& conf_,
	   const paraSet<Type>& par_,
	   const referencepointSet<Type>& rp_) :
      dat(dat_),
      conf(conf_),
      par(par_),
      rp(rp_){};

  PERREC_t<Type> getPerRec(const Type& logFbar){
    vector<Type> ls = rp.getLogSelectivity();
    PERREC_t<Type> r =  perRecruit_D(logFbar, dat, conf, par, ls, rp.aveYears, rp.nYears, rp.catchType);
    return r;
  }

  virtual Type operator()(const vector<Type>& logFbar) = 0;
  
};


template<class Type>
struct RefPointD_Base {
  
  virtual PERREC_t<Type> getPerRecruit(Type logFbar) = 0;
  virtual vector<Type> optimize(vector<Type> logF0) = 0;
  
};

template<class Type>
struct RefPointD_Known : RefPointD_Base<Type> {
  dataSet<Type> dat;
  confSet conf;
  paraSet<Type> par;
  referencepointSet<Type> rp;

  RefPointD_Known() = default;
  
  RefPointD_Known(const dataSet<Type>& dat_,
	    const confSet& conf_,
	    const paraSet<Type>& par_,
	    const referencepointSet<Type>& rp_) :
    dat(dat_),
    conf(conf_),
    par(par_),
    rp(rp_){}

  PERREC_t<Type> getPerRecruit(Type logFbar){
    vector<Type> ls = rp.getLogSelectivity();
    PERREC_t<Type> r =  perRecruit_D(logFbar, dat, conf, par, ls, rp.aveYears, rp.nYears, rp.catchType);
    return r;
  }

  virtual vector<Type> optimize(vector<Type> logF0) = 0;
};

template<class Type, class Functor>
struct RefPointD_Numeric : RefPointD_Base<Type> {

  Functor f;			// i.e., RPD_MSY. Should be derived from RPD_Base
  newton::newton_config cfg;
  
  RefPointD_Numeric(Functor f_,
		    newton::newton_config cfg_ = newton::newton_config()) :
    f(f_), cfg(cfg_) {};
  
  PERREC_t<Type> getPerRecruit(Type logFbar){
    PERREC_t<TMBad::ad_aug> r0 = f.getPerRec(logFbar);
    PERREC_t<Type> r = {
			newton::unsafe_cast<Type>(r0.logFbar),
			newton::unsafe_cast<Type>(r0.logYPR),
			newton::unsafe_cast<Type>(r0.logSPR),
			newton::unsafe_cast<Type>(r0.logSe),
			newton::unsafe_cast<Type>(r0.logRe),
			newton::unsafe_cast<Type>(r0.logYe),
			newton::unsafe_cast<Type>(r0.dSR0),
			newton::unsafe_cast<Type>(r0.logLifeExpectancy),
			newton::unsafe_cast<Type>(r0.logYearsLost),
			newton::unsafe_cast<Type>(r0.logDiscYPR),
			newton::unsafe_cast<Type>(r0.logDiscYe)
    };
    return r;
  }

  vector<Type> optimize(vector<Type> logF0){
    cfg.trace = 1;
    return newton::Newton(f,logF0, cfg);
  }  
};



template<class Type>
class Referencepoint_D {
  std::shared_ptr<RefPointD_Base<Type> > ptr;
  const char* name;
  int id;
  vector<Type> logF;
public:
  Referencepoint_D() = default;

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
    int n = logF.size();
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
      PERREC_t<Type> pr = ptr->getPerRecruit(logF(i));
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
    of->reportvector.push(logF, makeName("logF"));
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
    int n = logF.size();
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
      PERREC_t<Type> pr = ptr->getPerRecruit(logF(i));
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
    report(logF, makeName("logF"), of);
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

#define USING_RPD_BASE				\
  using RPD_Base<Type>::RPD_Base;		\
  using RPD_Base<Type>::getPerRec

#define USING_REFPOINT_KNOWN				\
  using RefPointD_Known<Type>::RefPointD_Known;		\
  using RefPointD_Known<Type>::getPerRecruit


// Define macro for repeat work
#define MAKE_REFPOINT_D(NAME)						\
  template<class Type>							\
  struct RefPointD_##NAME : RefPointD_Numeric<Type, RPD_##NAME<TMBad::ad_aug> > { \
    RefPointD_##NAME(dataSet<Type>& dat,				\
		     confSet& conf,					\
		     paraSet<Type>& par,				\
		     referencepointSet<Type>& rp,			\
		     newton::newton_config cfg = newton::newton_config()) : \
    RefPointD_Numeric<Type, RPD_##NAME<TMBad::ad_aug> >(RPD_##NAME<TMBad::ad_aug>(dat,conf,par,rp),cfg) {}; \
  }
  

//////////////////////////////// Specializations ////////////////////////////////

//////// FixedF Reference point (Combination of old B0 and Fsequence) ////////

template<class Type>
struct RefPointD_FixedF : RefPointD_Known<Type> {

  USING_REFPOINT_KNOWN;

  vector<Type> optimize(vector<Type> logF0){
    return logF0;
  }

};
  
//////// MSY Reference point ////////


template<class Type>
struct RPD_MSY : RPD_Base<Type> {

  USING_RPD_BASE;
  
  Type operator()(const vector<Type> &x) {
    Type logFbar = x(0);
    PERREC_t<Type> r = getPerRec(logFbar);
    return -r.logYe;
  }
};

MAKE_REFPOINT_D(MSY);


//////// MSY Reference point ////////


template<class Type>
struct RPD_Max : RPD_Base<Type> {

  USING_RPD_BASE;
  
  Type operator()(const vector<Type> &x) {
    Type logFbar = x(0);
    PERREC_t<Type> r = getPerRec(logFbar);
    return -r.logYPR;
  }
};

MAKE_REFPOINT_D(Max);


//////// %dYPR(0) Reference point (e.g., F~0.1~) ////////


template<class Type>
struct RPD_xdYPR : RPD_Base<Type> {

 // dataSet<Type> dat;
 //  confSet conf;
 //  paraSet<Type> par;
 //  referencepointSet<Type> rp;
  Type dYPR0;

  RPD_xdYPR() = default;
  
  RPD_xdYPR(const dataSet<Type>& dat_,
	    const confSet& conf_,
	    const paraSet<Type>& par_,
	    const referencepointSet<Type>& rp_) : RPD_Base<Type>(dat_, conf_, par_, rp_), dYPR0(R_NegInf) {
    dYPR0 = dYPR(Type(SAM_NegInf), this->dat, this->conf, this->par, this->rp);
  };

  // using RPD_Base<Type>::getPerRec;
  
  Type operator()(const vector<Type> &x) {
    if(x.size() != this->rp.xVal.size())
      Rf_error("In reference point xdYPR, length of F does not match length of fractions.");
    Type kappa = 0.0;
    for(int i = 0; i < x.size(); ++i){
      Type logFbar = x(i);
      Type v = dYPR(logFbar, this->dat,this->conf,this->par,this->rp);
      Type tmp = v - this->rp.xVal(i) * dYPR0;
      kappa += tmp * tmp;
    }
    return kappa;
  }
};

MAKE_REFPOINT_D(xdYPR);



//////// %SPR(0) Reference point (e.g., F~35%SPR~) ////////


template<class Type>
struct RPD_xSPR : RPD_Base<Type> {

 // dataSet<Type> dat;
 //  confSet conf;
 //  paraSet<Type> par;
 //  referencepointSet<Type> rp;
  Type SPR0;

  RPD_xSPR() = default;
  
  RPD_xSPR(const dataSet<Type>& dat_,
	    const confSet& conf_,
	    const paraSet<Type>& par_,
	   const referencepointSet<Type>& rp_) : RPD_Base<Type>(dat_, conf_, par_, rp_), SPR0(spawnersPerRecruit_i(Type(SAM_NegInf), this->dat, this->conf, this->par, this->rp)) {};

  // using RPD_Base<Type>::getPerRec;
  
  Type operator()(const vector<Type> &x) {
    if(x.size() != this->rp.xVal.size())
      Rf_error("In reference point xSPR, length of F does not match length of fractions.");
    Type kappa = 0.0;
    for(int i = 0; i < x.size(); ++i){
      Type logFbar = x(i);
      Type v = spawnersPerRecruit_i(logFbar, this->dat,this->conf,this->par,this->rp);
      Type tmp = v - this->rp.xVal(i) * SPR0;
      kappa += tmp * tmp;
    }
    return kappa;
  }
};

MAKE_REFPOINT_D(xSPR);



//////// %B(0) Reference point (e.g., F~20%B0~) ////////


template<class Type>
struct RPD_xB0 : RPD_Base<Type> {

  Type B0;

  RPD_xB0() = default;
  
  RPD_xB0(const dataSet<Type>& dat_,
	    const confSet& conf_,
	    const paraSet<Type>& par_,
	   const referencepointSet<Type>& rp_) : RPD_Base<Type>(dat_, conf_, par_, rp_), B0(B0_i(this->dat, this->conf, this->par, this->rp)) {};
  
  Type operator()(const vector<Type> &x) {
    if(x.size() != this->rp.xVal.size())
      Rf_error("In reference point xB0, length of F does not match length of fractions.");
    Type kappa = 0.0;
    for(int i = 0; i < x.size(); ++i){
      Type logFbar = x(i);      
      Type v = equilibriumBiomass_i(logFbar, this->dat,this->conf,this->par,this->rp);
      Type tmp = v - this->rp.xVal(i) * B0;
      kappa += tmp * tmp;
    }
    return kappa;
  }
};

MAKE_REFPOINT_D(xB0);


//////// Maximum yield per year lost (v1) ////////


template<class Type>
struct RPD_MYPYLdiv : RPD_Base<Type> {

  Type logAgeRange;
  RPD_MYPYLdiv() = default;
  
  RPD_MYPYLdiv(const dataSet<Type>& dat_,
	       const confSet& conf_,
	       const paraSet<Type>& par_,
	       const referencepointSet<Type>& rp_) : RPD_Base<Type>(dat_, conf_, par_, rp_), logAgeRange(log((Type)conf_.maxAge - (Type)conf_.minAge + (Type)1.0)) {};

  using RPD_Base<Type>::getPerRec;
 
  Type operator()(const vector<Type> &x) {    
    Type logFbar = x(0);
    PERREC_t<Type> r = getPerRec(logFbar);
    // Yield / (1 + YearsLost/AgeRange)
    Type tmp = r.logYe - logspace_add2(Type(0.0), r.logYearsLost - logAgeRange);
    return -tmp;
  }
};

MAKE_REFPOINT_D(MYPYLdiv);

//////// Maximum yield per year lost (v2) ////////


template<class Type>
struct RPD_MYPYLprod : RPD_Base<Type> {

  Type logAgeRange;
  RPD_MYPYLprod() = default;
  
  RPD_MYPYLprod(const dataSet<Type>& dat_,
	       const confSet& conf_,
	       const paraSet<Type>& par_,
	       const referencepointSet<Type>& rp_) : RPD_Base<Type>(dat_, conf_, par_, rp_), logAgeRange(log((Type)conf_.maxAge - (Type)conf_.minAge + (Type)1.0)) {};

  using RPD_Base<Type>::getPerRec;
 
  Type operator()(const vector<Type> &x) {    
    Type logFbar = x(0);
    PERREC_t<Type> r = getPerRec(logFbar);
    // Yield * (1 - exp(-YearsLost/AgeRange))
    Type tmp = r.logYe - logspace_sub2(Type(0.0), r.logYearsLost - logAgeRange);
    return -tmp;
  }
};

MAKE_REFPOINT_D(MYPYLprod);

//////// Maximum (life year) Discounted Yield ////////


template<class Type>
struct RPD_MDY : RPD_Base<Type> {

  USING_RPD_BASE;
 
  Type operator()(const vector<Type> &x) {    
    Type logFbar = x(0);
    PERREC_t<Type> r = getPerRec(logFbar);
    return -r.logDiscYe;
  }
};

MAKE_REFPOINT_D(MDY);


//////// Crash ////////


template<class Type>
struct RPD_Crash : RPD_Base<Type> {

  Type logdSR0;
  RPD_Crash() = default;
  
  RPD_Crash(const dataSet<Type>& dat_,
	    const confSet& conf_,
	    const paraSet<Type>& par_,
	    const referencepointSet<Type>& rp_) : RPD_Base<Type>(dat_, conf_, par_, rp_), logdSR0(R_NegInf) {
    Recruitment<Type> rec = makeRecruitmentFunction(conf_,par_);
    logdSR0 = log(rec.dSR(Type(SAM_Zero)));
  };

  using RPD_Base<Type>::getPerRec;
  
  Type operator()(const vector<Type> &x) {    
    Type logFbar = x(0);
    PERREC_t<Type> r = getPerRec(logFbar);
    Type tmp = logdSR0 - (-r.logSPR);
    return tmp * tmp;
  }
};

MAKE_REFPOINT_D(Crash);



//////// Ext ////////


template<class Type>
struct RPD_Ext : RPD_Base<Type> {

  Type logdSR0;
  RPD_Ext() = default;
  
  RPD_Ext(const dataSet<Type>& dat_,
	    const confSet& conf_,
	    const paraSet<Type>& par_,
	    const referencepointSet<Type>& rp_) : RPD_Base<Type>(dat_, conf_, par_, rp_), logdSR0(R_NegInf) {
    Recruitment<Type> rec = makeRecruitmentFunction(conf_,par_);
    logdSR0 = log(rec.dSR(Type(SAM_Zero)));
  };

  using RPD_Base<Type>::getPerRec;
  
  Type operator()(const vector<Type> &x) {    
    Type logFbar = x(0);
    PERREC_t<Type> r = getPerRec(logFbar);
    Type tmp = r.logSe;
    return tmp * tmp;
  }
};

MAKE_REFPOINT_D(Ext);


//////// Lim?? ////////


/////////////////////////////// Function to call from main program ////////////////////////////////


template<class Type>
void reportDeterministicReferencePoints(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, array<Type> &logF, Recruitment<Type> &recruit, referencepointList<Type> &referencepoints, objective_function<Type> *of){
  
  if(referencepoints.size() == 0)
    return;
  
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
	logFsq(xi) = fbari(conf, logF, logF.cols()-1 - CppAD::Integer(rps.xVal(xi)), true);
      rp = Referencepoint_D<Type>("StatusQuo",i,logFsq, new RefPointD_FixedF<Type>(dat,conf,par,rps));
    }else if(rpt  == ReferencePointDeterministic::MSY){
      rp = Referencepoint_D<Type>("MSY",i,rps.logF0, new RefPointD_MSY<Type>(dat,conf,par,rps));
    }else if(rpt  == ReferencePointDeterministic::MSYRange){
    }else if(rpt  == ReferencePointDeterministic::Max){
      rp = Referencepoint_D<Type>("Max",i,rps.logF0, new RefPointD_Max<Type>(dat,conf,par,rps));
    }else if(rpt  == ReferencePointDeterministic::xdYPR){
      rp = Referencepoint_D<Type>("xdYPR",i,rps.logF0, new RefPointD_xdYPR<Type>(dat,conf,par,rps));
    }else if(rpt  == ReferencePointDeterministic::xSPR){
      rp = Referencepoint_D<Type>("xSPR",i,rps.logF0, new RefPointD_xSPR<Type>(dat,conf,par,rps));
    }else if(rpt  == ReferencePointDeterministic::xB0){
      rp = Referencepoint_D<Type>("xB0",i,rps.logF0, new RefPointD_xB0<Type>(dat,conf,par,rps));
    }else if(rpt  == ReferencePointDeterministic::MYPYLdiv){
      rp = Referencepoint_D<Type>("MYPYLdiv",i,rps.logF0, new RefPointD_MYPYLdiv<Type>(dat,conf,par,rps));
    }else if(rpt  == ReferencePointDeterministic::MYPYLprod){
      rp = Referencepoint_D<Type>("MYPYLprod",i,rps.logF0, new RefPointD_MYPYLprod<Type>(dat,conf,par,rps));
    }else if(rpt  == ReferencePointDeterministic::MDY){
      rp = Referencepoint_D<Type>("MDY",i,rps.logF0, new RefPointD_MDY<Type>(dat,conf,par,rps));
    }else if(rpt  == ReferencePointDeterministic::Crash){
      rp = Referencepoint_D<Type>("Crash",i,rps.logF0, new RefPointD_Crash<Type>(dat,conf,par,rps));
    }else if(rpt  == ReferencePointDeterministic::Ext){
      rp = Referencepoint_D<Type>("Ext",i,rps.logF0, new RefPointD_Ext<Type>(dat,conf,par,rps));
    }else if(rpt  == ReferencePointDeterministic::Lim){

    }else{
      Rf_error("Referencepoint type not implemented yet");
    }
    rp.report(of);
    rp.adreport(of);
  }

  // Report relative reference points
 
    return;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Stochastic //////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////



#endif
