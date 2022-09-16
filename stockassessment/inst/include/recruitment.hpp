#pragma once
#ifndef SAM_RECRUITMENT_HPP
#define SAM_RECRUITMENT_HPP

#include <memory>

#define SAMREC_TYPEDEFS(scalartype_)			\
public:						\
typedef scalartype_ scalartype;			\
typedef vector<scalartype> vectortype;		\
typedef matrix<scalartype> matrixtype;		\
typedef array<scalartype> arraytype

/*
Depensatory_A_%s models:
R^B_%s(S) = R_%s(S^g)
(SigmoidalBevertonHolt is of this type)

Depensatory_B_%s models:
R^A_%s(S) = R_%s(S) * S / (d + S)

Depensatory_C_%s models:
R^B_%s(S) = R_%s(S) / (1 + exp(-l * (S-d))) 


 */



enum RecruitmentModel {
		       ICESforecast = -2,
		       NoRecruit = -1,
		       LogRandomWalk = 0,
		       Ricker = 1,
		       BevertonHolt = 2,
		       ConstantMean = 3,
		       LogisticHockeyStick = 60,
		       HockeyStick = 61,
		       LogAR1 = 62,
		       BentHyperbola = 63,
		       Power_CMP = 64,
		       Power_NCMP = 65,
		       Shepherd = 66,
		       Hassel_Deriso = 67,
		       SailaLorda = 68,
		       SigmoidalBevertonHolt = 69,
		       Spline_CMP = 90,
		       Spline_Smooth = 91,
		       Spline_General = 92,
		       Spline_ConvexCompensatory = 93,
		       Depensatory_B_Ricker = 201,
		       Depensatory_B_BevertonHolt = 202,
		       Depensatory_B_LogisticHockeyStick = 260,
		       Depensatory_B_HockeyStick = 261,
		       Depensatory_B_BentHyperbola = 263,
		       Depensatory_B_Power = 264,
		       Depensatory_B_Shepherd = 266,
		       Depensatory_B_Hassel_Deriso = 267,
		       Depensatory_B_Spline_CMP = 290,
		       Depensatory_B_Spline_ConvexCompensatory = 293,
		       Depensatory_C_Ricker = 401,
		       Depensatory_C_BevertonHolt = 402,
		       Depensatory_C_Spline_CMP = 490,
		       Depensatory_C_Spline_ConvexCompensatory = 493,
		       Num_Ricker = 991,
		       Num_BevertonHolt = 992
		    
};


template<class Type>
struct RecruitmentWorker {

  int isAutoregressive;
  int isTimevarying;

  RecruitmentWorker() : isAutoregressive(0), isTimevarying(0) {}
  RecruitmentWorker(int isA, int isT) : isAutoregressive(isA), isTimevarying(isT) {}
  
  // Stock recruitment function logR = f(logssb)
  virtual Type operator()(Type logssb, Type lastLogR, Type year) = 0;

  Type R(Type logssb, Type lastLogR, Type year){
    return exp(operator()(logssb, lastLogR, year));
  }

  // Deterministic equilibrium biomass for lambda = 1/SPR
  virtual Type logSe(Type logLambda) = 0;
  // Derivative
  virtual Type dSR(Type logssb) = 0;

  virtual Type logSAtMaxR() = 0;
  virtual Type logMaxR() = 0;
  virtual Type maxGradient() = 0;

  virtual ~RecruitmentWorker() {};

};

template<class Type>
struct Recruitment {
public:
  // RecruitmentWorker<Type>* ptr;
  const char* name;
  std::shared_ptr<RecruitmentWorker<Type> > ptr;
  explicit Recruitment() : name("Uninitialized"), ptr(nullptr) {};
  explicit Recruitment(const char* nm, RecruitmentWorker<Type>* p=nullptr) : name(nm), ptr(p) {};
  virtual ~Recruitment() { ptr.reset(); };

  
  int isTimevarying(){
    return ptr->isTimevarying;
  }

  int isAutoregressive(){
    return ptr->isAutoregressive;
  }
  
  Type operator()(Type logssb, Type lastLogR, Type year){
    return ptr->operator()(logssb, lastLogR, year);
  }
  Type R(Type logssb, Type lastLogR, Type year){
    return ptr->R(logssb, lastLogR, year);
  }
  Type logSe(Type logLambda){
    return ptr->logSe(logLambda);
  }
  Type dSR(Type logssb){
    return ptr->dSR(logssb);
  }
  Type logSAtMaxR(){
    return ptr->logSAtMaxR();
  }
  Type logMaxR(){
    return ptr->logMaxR();
  }
  Type maxGradient(){
    return ptr->maxGradient();
  }

};

template<class Functor>
struct WrapSR {
  SAMREC_TYPEDEFS(typename Functor::scalartype);
  Functor f;
  WrapSR(const Functor& f_) : f(f_) {};
  WrapSR(const WrapSR<Functor>& wf_) : f(wf_.f) {};

  auto to_ad(){
    return WrapSR<decltype(f.to_ad())>(f.to_ad());
  }
  
  // template <template<class> class V, class T>
  // T operator()(const V<T> &x){
  scalartype operator()(const vectortype &x) {
    // T v = -f(x(0));
    // return v;
    return -f(x(0));
  }
};

template<class Functor>
struct WrapSRp {
  SAMREC_TYPEDEFS(typename Functor::scalartype);
  Functor f;
  WrapSRp(const Functor& f_) : f(f_) {};
  WrapSRp(const WrapSRp<Functor>& wf_) : f(wf_.f) {};

  auto to_ad(){
    return WrapSRp<decltype(f.to_ad())>(f.to_ad());
  }

 
  scalartype operator()(const vectortype &x) {
    return f(x(0));
  }
};

// template<class Functor>
template <class Functor>
struct Exp {
  SAMREC_TYPEDEFS(typename Functor::scalartype);
  Functor f;
  Exp(const Functor& f_) : f(f_) {};
  Exp(const Exp<Functor>& ef_) : f(ef_.f) {};

  auto to_ad(){
    return Exp<decltype(f.to_ad())>(f.to_ad());
  }

   
  scalartype operator()(const scalartype& x){
    return exp(f(x));
  }
  vectortype operator()(const vectortype& x){
    vectortype r(x.size()); r.setZero();
    for(int i = 0; i < r.size(); ++i)
      r(i) = this->operator()(x(i));
    return r;
  }
};

template<class Functor>
struct diffSR {
  SAMREC_TYPEDEFS(typename Functor::scalartype);
  Exp<Functor> f;
  
  diffSR(const Functor& f_) : f(f_) {}; //f(Exp<Functor>(f_)) {};
  diffSR(const Exp<Functor>& ef_) : f(ef_) {};
  diffSR(const diffSR<Functor>& wf_) : f(wf_.f) {};

  auto to_ad(){
    return diffSR<decltype(f.to_ad())>(f.to_ad());
  }

 
  
  // template<class T>
  // T operator()(const T& x){
  scalartype operator()(const scalartype& x){
    // vector<scalartype> v(1); v[0] = x;
    // auto f2 = f.to_ad();
    // vectortype g = autodiff::gradient(f, v);
    // TMBad::ADFun<> G(TMBad::StdWrap<decltype(f.to_ad()),vector<TMBad::ad_aug> >(f), v);
    // G = G.JacFun();
    // scalartype gg = newton::unsafe_cast<scalartype>(G(v)[0]);
    // scalartype gg = G(v)[0];
    // return gg / exp(x);
    // return g(0) / exp(x);
    scalartype h = 0.0001 * sqrt(x * x + SAM_Zero);
    scalartype x1 = x + 2.0 * h;
    scalartype v1 = f(x1);
    scalartype x2 = x + h;
    scalartype v2 = f(x2);
    scalartype x3 = x - h;
    scalartype v3 = f(x3);
    scalartype x4 = x - 2.0 * h;
    scalartype v4 = f(x4);
    scalartype g = (-v1 + 8.0 * v2 - 8.0 * v3 + v4) / (12.0 * h);
    return (scalartype)g / exp(x);
  }
};


template<class Functor>
struct WrapDiffSR {
  SAMREC_TYPEDEFS(typename Functor::scalartype);
  diffSR<Functor> dF;

  WrapDiffSR(const Functor& f_) : dF(f_) {};
  WrapDiffSR(const diffSR<Functor>& f_) : dF(f_) {};
  WrapDiffSR(const WrapDiffSR<Functor>& wf_) : dF(wf_.dF) {};

  auto to_ad(){
    return WrapDiffSR<decltype(dF.to_ad())>(dF.to_ad());
  }

   
  scalartype operator()(const vectortype &x) {
    return -dF(x(0));
  }
};

template<class Functor>
struct WrapEquiS {
  SAMREC_TYPEDEFS(typename Functor::scalartype);
  Functor f;
  scalartype logLambda;

  WrapEquiS(const Functor& f_, const scalartype& l) : f(f_), logLambda(l) {};

  auto to_ad(){
    return WrapEquiS<decltype(f.to_ad())>(f.to_ad(), TMBad::ad_aug(logLambda));
  }
 
  
  scalartype operator()(const vectortype &x) {
    scalartype xv = x(0);
    scalartype loga = logLambda + f(xv);
    scalartype logb = xv;
    // v = (a-b)^2
    //scalartype v = exp(2.0 * loga) + exp(2.0 * logb) - 2.0 * exp(loga + logb);
    scalartype v = (loga - logb) * (loga - logb);
    return v;
  }  
};




//template<class Type, class Functor>
// RecruitmentNumeric is not allowed to be autoregressive or timevarying
template<class Functor>
struct RecruitmentNumeric : RecruitmentWorker<typename Functor::scalartype> {
  SAMREC_TYPEDEFS(typename Functor::scalartype);
  Functor f;
  
  RecruitmentNumeric(const Functor& f_) : RecruitmentWorker<scalartype>(0,0), f(f_) {};  
   
  scalartype operator()(scalartype logssb, scalartype lastLogR, scalartype year){
    return f(logssb);
  }
  
  virtual scalartype logSe(scalartype logLambda){
    // WARNING: Only working if re-running with optimized parameter values!
    // WARNING: Do not use in nll!
    
    WrapEquiS<Functor> fx(f,logLambda);//(f.to_ad(), TMBad::ad_aug(logLambda));
    vectortype x0v(1); x0v(0) = 20.0;
    newton::newton_config cfg;
    cfg.simplify = false;
    cfg.on_failure_return_nan = false;
    auto fx_ad = fx.to_ad();
    vectortype v = newton::Newton(fx_ad,x0v,cfg);
    return TMBad::CondExpLt((scalartype)v(0),(scalartype)SAM_NegInf,(scalartype)SAM_NegInf,(scalartype)v(0));
  }

  virtual scalartype dSR(scalartype logssb){
    diffSR<Functor> dF(f);
    return dF(logssb);
  }
  
  virtual scalartype logSAtMaxR(){
    WrapSR<Functor> fx(f);//.to_ad());
    vectortype x0v(1); x0v(0) = 20.0;
    auto fx_ad = fx.to_ad();
    newton::newton_config cfg;
    cfg.simplify = false;
    cfg.on_failure_return_nan = false;
    vectortype v = newton::Newton(fx_ad,x0v);  
    return v(0);
  }
  virtual scalartype logMaxR(){   
    scalartype v = logSAtMaxR();    
    return f(v);
  }
  virtual scalartype maxGradient(){
    WrapDiffSR<Functor> fx(f);//.to_ad());
    vectortype x0v(1); x0v(0) = 0.0;
    auto fx_ad = fx.to_ad();
    newton::newton_config cfg;
    cfg.simplify = false;
    cfg.on_failure_return_nan = false;
    vectortype v = newton::Newton(fx_ad,x0v);      
    return dSR(v(0));
  }  
};




template <class Functor>
struct WrapDepensatoryA {
  SAMREC_TYPEDEFS(typename Functor::scalartype);
  Functor f;
  scalartype logd;

  WrapDepensatoryA(const Functor& f_, const scalartype& logd_) : f(f_), logd(logd_) {};

  auto to_ad(){
    return WrapDepensatoryA<decltype(f.to_ad())>(f.to_ad(), TMBad::ad_aug(logd));
  }

  
  // template<class T>
  // T operator()(const T& x){
  scalartype operator()(const scalartype& x){
    scalartype logssb = exp(logd) * x;
    scalartype v = f(logssb, R_NaReal, R_NaReal);
    return v;
  }



  
};

template <class Functor>
RecruitmentNumeric<WrapDepensatoryA<Functor> >* Rec_DepensatoryA(Functor SR, typename Functor::scalartype logd){
  WrapDepensatoryA<Functor> DepSR(SR, logd);
  return new RecruitmentNumeric<WrapDepensatoryA<Functor> >(DepSR);
}




template <class Functor>
struct WrapDepensatoryB {
  SAMREC_TYPEDEFS(typename Functor::scalartype);
  Functor f;
  scalartype logd;

  WrapDepensatoryB(const Functor& f_, const scalartype& logd_) : f(f_), logd(logd_) {};

  auto to_ad(){
    return WrapDepensatoryB<decltype(f.to_ad())>(f.to_ad(), TMBad::ad_aug(logd));
  }

  
  scalartype operator()(const scalartype& x){
    scalartype logssb = x;
    scalartype v = f(logssb, (scalartype)R_NaReal, (scalartype)R_NaReal) + logssb - logspace_add_SAM(logssb, (scalartype)logd);
    return v;
  }

};

template <class Functor>
RecruitmentNumeric<WrapDepensatoryB<Functor> >* Rec_DepensatoryB(Functor SR, typename Functor::scalartype logd){
  WrapDepensatoryB<Functor> DepSR(SR, logd);
  return new RecruitmentNumeric<WrapDepensatoryB<Functor> >(DepSR);
}




template <class Functor>
struct WrapDepensatoryC {
  SAMREC_TYPEDEFS(typename Functor::scalartype);
  Functor f;
  scalartype logd;
  scalartype logl;

  WrapDepensatoryC(const Functor& f_, const scalartype& logd_, const scalartype& logl_) : f(f_), logd(logd_), logl(logl_) {};

  auto to_ad(){
    return WrapDepensatoryC<decltype(f.to_ad())>(f.to_ad(), TMBad::ad_aug(logd), TMBad::ad_aug(logl));
  }
  
  scalartype operator()(const scalartype& x){
    scalartype logssb = x;
    scalartype v = f(logssb, (scalartype)R_NaReal, (scalartype)R_NaReal) - logspace_add_SAM((scalartype)0.0, -exp(logl) * (exp(logssb) - exp(logd)));
    return v;
  }

};

template <class Functor>
RecruitmentNumeric<WrapDepensatoryC<Functor> >* Rec_DepensatoryC(Functor SR, typename Functor::scalartype logd, typename Functor::scalartype logl){
  WrapDepensatoryC<Functor> DepSR(SR, logd, logl);
  return new RecruitmentNumeric<WrapDepensatoryC<Functor> >(DepSR);
}


// Recruitment function -2
// No recruitment
template<class scalartype_>
struct Rec_ICESforecast : RecruitmentWorker<scalartype_> {
  SAMREC_TYPEDEFS(scalartype_);
  scalartype logMean;
  scalartype logSd;

  Rec_ICESforecast(scalartype lm, scalartype lsd) : RecruitmentWorker<scalartype>(0,0), logMean(lm), logSd(lsd) {}
  
  scalartype operator()(scalartype logssb, scalartype lastLogR, scalartype year){
    return logMean;
  }
  scalartype logSe(scalartype logLambda){
    return logLambda + logMean;
  }
  scalartype dSR(scalartype logssb){
    return 0.0;
  }
  scalartype logSAtMaxR(){
    return R_NaReal;
  }
  scalartype logMaxR(){
    return logMean;
  }
  scalartype maxGradient(){
    return 0.0;
  }

  Rec_ICESforecast<TMBad::ad_aug> to_ad(){
    return Rec_ICESforecast<TMBad::ad_aug>(TMBad::ad_aug(logMean),TMBad::ad_aug(logSd));
  }

};


// Recruitment function -1
// No recruitment
template<class scalartype_>
struct Rec_None : RecruitmentWorker<scalartype_> {
  SAMREC_TYPEDEFS(scalartype_);

  Rec_None() : RecruitmentWorker<scalartype>(0,0) {};
  
  scalartype operator()(scalartype logssb, scalartype lastLogR, scalartype year){
    return R_NegInf;
  }
  scalartype logSe(scalartype logLambda){
    return R_NegInf;
  }
  scalartype dSR(scalartype logssb){
    return R_NaReal;
  }
  scalartype logSAtMaxR(){
    return R_NaReal;
  }
   scalartype logMaxR(){
    return R_NegInf;
  }
  scalartype maxGradient(){
    return R_NaReal;
  }
  Rec_None<TMBad::ad_aug> to_ad(){
    return Rec_None<TMBad::ad_aug>();
  }

};

// Recruitment function 0
// Random walk
template<class scalartype_>
struct Rec_LogRW : RecruitmentWorker<scalartype_> {
  SAMREC_TYPEDEFS(scalartype_);

  Rec_LogRW() : RecruitmentWorker<scalartype>(1,0) {};
  
  scalartype operator()(scalartype logssb, scalartype lastLogR, scalartype year){
    return lastLogR;
  }
  scalartype logSe(scalartype logLambda){
    return R_NaReal;
  }
  scalartype dSR(scalartype logssb){
    return R_NaReal;
  }
  scalartype logSAtMaxR(){
    return R_NaReal;
  }
   scalartype logMaxR(){
    return R_NaReal;
  }
  scalartype maxGradient(){
    return R_NaReal;
  }
  Rec_LogRW<TMBad::ad_aug> to_ad(){
    return Rec_LogRW<TMBad::ad_aug>();
  }

};

// Recruitment function 1
// Ricker
template<class scalartype_>
struct Rec_Ricker : RecruitmentWorker<scalartype_> {

  SAMREC_TYPEDEFS(scalartype_);
  
  scalartype loga;
  scalartype logb;

  Rec_Ricker(scalartype la, scalartype lb) : RecruitmentWorker<scalartype>(0,0), loga(la), logb(lb) {};
  template<class T>
  Rec_Ricker(const Rec_Ricker<T>& other) : RecruitmentWorker<scalartype>(0,0), loga(other.loga), logb(other.logb) {}
  
  scalartype operator()(scalartype logssb, scalartype lastLogR, scalartype year){
    return loga + logssb - exp(logb + logssb);
  }
  
  scalartype logSe(scalartype logLambda){    
    scalartype v = TMBad::CondExpLe(loga+logLambda, scalartype(0.0), exp(loga+logLambda), (loga + logLambda) * exp(-logb));
    return log(v + SAM_Zero);
  }
  scalartype dSR(scalartype logssb){
    scalartype a = exp(loga);
    scalartype ee = -exp(logb + logssb);
    return -a * exp(ee) * (-ee - 1.0);
  }
  scalartype logSAtMaxR(){
    return -logb;
  }
  scalartype logMaxR(){
    return loga - 1.0 - logb;
  }
  scalartype maxGradient(){
    return exp(loga);
  }

  Rec_Ricker<TMBad::ad_aug> to_ad(){
    return Rec_Ricker<TMBad::ad_aug>(TMBad::ad_aug(loga),TMBad::ad_aug(logb));
  }

  
};


// Recruitment function 2
// Beverton Holt
template<class scalartype_>
struct Rec_BevertonHolt : RecruitmentWorker<scalartype_> {
  SAMREC_TYPEDEFS(scalartype_);

  scalartype loga;
  scalartype logb;

  Rec_BevertonHolt(scalartype la, scalartype lb) : RecruitmentWorker<scalartype>(0,0), loga(la), logb(lb) {};
  template<class T>
  Rec_BevertonHolt(const Rec_BevertonHolt<T>& other) : RecruitmentWorker<scalartype>(0,0), loga(other.loga), logb(other.logb) {}

  scalartype operator()(scalartype logssb, scalartype lastR, scalartype year){
    return loga + logssb - logspace_add_SAM(scalartype(0.0),logb + logssb);
  }
  
  scalartype logSe(scalartype logLambda){
    scalartype v = TMBad::CondExpLe(loga+logLambda, scalartype(0.0), scalartype(0.0), (exp(loga + logLambda) - 1.0) * exp(-logb));
    return log(v + 1.0e-16);
  }
  scalartype dSR(scalartype logssb){
    scalartype a = exp(loga);
    scalartype ee = exp(logb + logssb) + 1.0;
    return a / (ee * ee);
  }
  scalartype logSAtMaxR(){
    return R_PosInf;
  }
  scalartype logMaxR(){
    return loga - logb;
  }
  scalartype maxGradient(){
    return exp(loga);
  }

  Rec_BevertonHolt<TMBad::ad_aug> to_ad(){
    return Rec_BevertonHolt<TMBad::ad_aug>(TMBad::ad_aug(loga),TMBad::ad_aug(logb));
  }

};

// Recruitment function 3
// Constant mean
template<class scalartype_>
struct Rec_ConstantMean : RecruitmentWorker<scalartype_> {
 SAMREC_TYPEDEFS(scalartype_);

  vectortype logRvalue;
  vectortype constRecBreaks;

  Rec_ConstantMean(const vectortype& logRv, const vectortype& crb) : RecruitmentWorker<scalartype>(0,1), logRvalue(logRv), constRecBreaks(crb) {};

  template<class T>
  Rec_ConstantMean(const vector<T>& logRv, const vector<T>& crb) : RecruitmentWorker<scalartype>(0,1), logRvalue(logRv), constRecBreaks(crb) {}
  
  scalartype operator()(scalartype logssb, scalartype lastLogR, scalartype year){
    int usepar=0;
    for(int ii=0; ii<constRecBreaks.size(); ++ii){
      if(year>(scalartype)constRecBreaks(ii)){usepar++;}
    }
    return logRvalue(usepar);
  }
  
  scalartype logSe(scalartype logLambda){
    return  logLambda + logRvalue(logRvalue.size()-1);
  }
  scalartype dSR(scalartype logssb){
    return 0.0;
  }
  scalartype logSAtMaxR(){
    return R_NaReal;
  }
  scalartype logMaxR(){
    return R_NaReal;
  }
  scalartype maxGradient(){
    return 0.0;
  }

  Rec_ConstantMean<TMBad::ad_aug> to_ad(){
    return Rec_ConstantMean<TMBad::ad_aug>(logRvalue, constRecBreaks);
  }

};



// Recruitment function 60
// Logistic Hockey Stick

//     // 0: alpha
//     // 1: mu
//     // 2: theta
//     predN = rec_pars(0) + rec_pars(1) + rec_pars(2) + log(1.0 + exp(-exp(-rec_pars(2)))) + log(exp(log(thisSSB)-rec_pars(1) - rec_pars(2)) - log(1.0 + exp((thisSSB-exp(rec_pars(1)))/exp(rec_pars(1) + rec_pars(2)))) + log(1.0 + exp(-exp(-rec_pars(2)))));

// mdsr = exp(par.rec_pars(0));


template<class scalartype_>
struct RF_LogisticHockeyStick_t {
 SAMREC_TYPEDEFS(scalartype_);
  scalartype loga;			// alpha
  scalartype logm;			// mu
  scalartype logt;			// theta

  RF_LogisticHockeyStick_t(scalartype la, scalartype lm, scalartype lt) :
    loga(la), logm(lm), logt(lt) {}

  template<class T>
  RF_LogisticHockeyStick_t(const RF_LogisticHockeyStick_t<T>& other) :
    loga(other.loga), logm(other.logm), logt(other.logt) {}

  RF_LogisticHockeyStick_t<TMBad::ad_aug> to_ad(){
    return RF_LogisticHockeyStick_t<TMBad::ad_aug>(TMBad::ad_aug(loga),TMBad::ad_aug(logm),TMBad::ad_aug(logt));
  }
  
  // template <template<class> class V, class T>
  // T operator()(const V<T> &logssb0){
  // template<class T>
  // T operator()(const T& logssb){
  scalartype operator()(const scalartype& logssb){
    scalartype thisSSB = exp(logssb);
    scalartype v= loga + logm + logt + log(1.0 + exp(-exp(-logt))) + log(exp(logssb - logm - logt) - log(1.0 + exp((thisSSB-exp(logm))/exp(logm + logt))) + log(1.0 + exp(-exp(-logt))));
    return v;
  }
  
};
  
template<class scalartype_>
struct Rec_LogisticHockeyStick : RecruitmentNumeric<RF_LogisticHockeyStick_t<scalartype_> >  {
 SAMREC_TYPEDEFS(scalartype_);

  scalartype loga;
  scalartype logm;
  scalartype logt;
  
  Rec_LogisticHockeyStick(scalartype la, scalartype lm, scalartype lt) :
    RecruitmentNumeric<RF_LogisticHockeyStick_t<scalartype> >(RF_LogisticHockeyStick_t<scalartype>(la, lm, lt)), loga(la), logm(lm), logt(lt) {};

  template<class T>
  Rec_LogisticHockeyStick(const Rec_LogisticHockeyStick<T>& other) :
    RecruitmentNumeric<RF_LogisticHockeyStick_t<scalartype> >(other.f), loga(other.loga), logm(other.logm), logt(other.logt) {}

  
  scalartype logSAtMaxR(){
    return R_PosInf;
  }
   scalartype logMaxR(){
     return loga + logm + logt + log(1.0 + exp(-exp(-logt))) + log(exp(-logt) - log(1.0 + exp(-exp(-logt))));
  }
  scalartype maxGradient(){
    return exp(loga);
  }

  Rec_LogisticHockeyStick<TMBad::ad_aug> to_ad(){
    return Rec_LogisticHockeyStick<TMBad::ad_aug>(TMBad::ad_aug(loga),TMBad::ad_aug(logm),TMBad::ad_aug(logt));
  }

  
};


// Recruitment function 61
// Hockey Stick

template<class scalartype_>
struct Rec_HockeyStick : RecruitmentWorker<scalartype_> {
 SAMREC_TYPEDEFS(scalartype_);

  scalartype loglevel;
  scalartype logblim;

  Rec_HockeyStick(scalartype ll, scalartype lbl) :
    RecruitmentWorker<scalartype>(0,0), loglevel(ll), logblim(lbl) {};

  template<class T>
  Rec_HockeyStick(const Rec_HockeyStick<T>& other) :
    RecruitmentWorker<scalartype>(0,0), loglevel(other.loglevel), logblim(other.logblim) {}

  
  scalartype operator()(scalartype logssb, scalartype lastLogR, scalartype year){
    scalartype thisSSB = exp(logssb);
    return loglevel - logblim +
      log(thisSSB - (0.5 * ((thisSSB - exp(logblim))+scalartype(0.0)+TMBad::fabs((thisSSB - exp(logblim))-scalartype(0.0)))));
  }
  
  scalartype logSe(scalartype logLambda){
    return exp(logLambda + loglevel);
  }
  scalartype dSR(scalartype logssb){
    return TMBad::CondExpLt(logssb, logblim, exp(loglevel-logblim), scalartype(0.0));
  }
  scalartype logSAtMaxR(){
    return R_NaReal;
  }
  scalartype logMaxR(){
    return loglevel;
  }
  scalartype maxGradient(){
    return exp(loglevel-logblim);
  }

  Rec_HockeyStick<TMBad::ad_aug> to_ad(){
    return Rec_HockeyStick<TMBad::ad_aug>(TMBad::ad_aug(loglevel),TMBad::ad_aug(logblim));
  }
  
};


// Recruitment function 62
// AR(1) on log-recruitment

template<class scalartype_>
struct Rec_LogAR1 : RecruitmentWorker<scalartype_> {
 SAMREC_TYPEDEFS(scalartype_);

  scalartype loglevel;
  scalartype logitPhi;

  Rec_LogAR1(scalartype ll, scalartype lp) : RecruitmentWorker<scalartype>(1,0), loglevel(ll), logitPhi(lp) {};
  
  scalartype operator()(scalartype logssb, scalartype lastLogR, scalartype year){
    return loglevel + (2.0 / (1.0 + exp(-logitPhi)) - 1.0) * (lastLogR - loglevel);
  }
  
  scalartype logSe(scalartype logLambda){
    return exp(logLambda + loglevel);
  }
  scalartype dSR(scalartype logssb){
    return 0.0;
  }
  scalartype logSAtMaxR(){
    return R_NegInf;
  }
  scalartype logMaxR(){
    return loglevel;
  }
  scalartype maxGradient(){
    return 0.0;
  }
  Rec_LogAR1<TMBad::ad_aug> to_ad(){
    return Rec_LogAR1<TMBad::ad_aug>(TMBad::ad_aug(loglevel),TMBad::ad_aug(logitPhi));
  }

};


// Recruitment function 63
// Bent hyperbola

//     /*
//       Source: e.g., DOI:10.1093/icesjms/fsq055
//       rec_pars(0): log-Blim
//       rec_pars(1): log of half the slope from 0 to Blim
//       rec_pars(2): log-Smoothness parameter
//      */
//     predN = rec_pars(1) +
//       log(thisSSB + sqrt(exp(2.0 * rec_pars(0)) + (exp(2.0 * rec_pars(2)) / 4.0)) -
//     	  sqrt(pow(thisSSB-exp(rec_pars(0)),2) + (exp(2.0 * rec_pars(2)) / 4.0)));

//   Se = (2.0 * sqrt(exp(2.0 * newPar.rec_pars(0)) + exp(2.0 * newPar.rec_pars(2)) / 4.0) / (lambda * exp(newPar.rec_pars(1) + logRecCorrection)) - 2.0 * exp(newPar.rec_pars(0)) - 2.0 * sqrt(exp(2.0 * newPar.rec_pars(0)) + exp(2.0 * newPar.rec_pars(2)) / 4.0)) / ( 1.0 / ((lambda * lambda * exp(2.0 * (newPar.rec_pars(1) + logRecCorrection)))) - 2.0 / (lambda * exp(newPar.rec_pars(1) + logRecCorrection))  );  

//    mdsr = 2.0 * exp(par.rec_pars(1));


#define THIS_RF_TO_AD(NAME)			\
  RF_##NAME##_t<TMBad::ad_aug> to_ad(){		\
    return RF_##NAME##_t<TMBad::ad_aug>(*this);	\
  }

#define THIS_REC_TO_AD(NAME)			\
  Rec_##NAME<TMBad::ad_aug> to_ad(){		\
    return Rec_##NAME<TMBad::ad_aug>(*this);	\
  }


template<class scalartype_>
struct RF_BentHyperbola_t {
  SAMREC_TYPEDEFS(scalartype_);
 scalartype logBlim;
  scalartype logHalfSlope;
  scalartype logSmooth;

  RF_BentHyperbola_t(scalartype loga_, scalartype logb_, scalartype logg_) :
    logBlim(loga_), logHalfSlope(logb_), logSmooth(logg_) {}

  template<class T>
  RF_BentHyperbola_t(const RF_BentHyperbola_t<T>& other) :
    logBlim(other.logBlim), logHalfSlope(other.logHalfSlope), logSmooth(other.logSmooth) {}

  THIS_RF_TO_AD(BentHyperbola);
  
  
  // template <template<class> class V, class T>
  // T operator()(const V<T> &logssb0){
    // T logssb = logssb0(0);
  // template<class T>
  // T operator()(const T& logssb){
  scalartype operator()(const scalartype& logssb){
    scalartype thisSSB = exp(logssb);
    scalartype v = logHalfSlope + log(thisSSB + sqrt(exp(2.0 * logBlim) + (exp(2.0 * logSmooth) / 4.0)) -
			     sqrt(pow(thisSSB-exp(logBlim),2) + (exp(2.0 * logSmooth) / 4.0)));
    return v;
  }
 
};
  
template<class scalartype_>
struct Rec_BentHyperbola : RecruitmentNumeric<RF_BentHyperbola_t<scalartype_> >  {
  SAMREC_TYPEDEFS(scalartype_);
 // Implement with known values when time permits!
  Rec_BentHyperbola(scalartype logBlim, scalartype logHalfSlope, scalartype logSmooth) :
    RecruitmentNumeric<RF_BentHyperbola_t<scalartype> >(RF_BentHyperbola_t<scalartype>(logBlim, logHalfSlope, logSmooth)) {};

  template<class T>
  Rec_BentHyperbola(const Rec_BentHyperbola<T>& other) :
    RecruitmentNumeric<RF_BentHyperbola_t<scalartype> >(other.f) {}

  THIS_REC_TO_AD(BentHyperbola);
  
};




// Recruitment function 64
// Power function with compensatory mortality property


//     predN = rec_pars(0) + invlogit(rec_pars(1)) * log(thisSSB);
//     break;

 //   Se = exp(1.0 / (1.0 - invlogit(newPar.rec_pars(1))) * (newPar.rec_pars(0) + logRecCorrection + log(lambda)));

//    mdsr = R_PosInf;


template<class scalartype_>
struct RF_PowerCMP_t {
  SAMREC_TYPEDEFS(scalartype_);
 scalartype loga;
  scalartype logb;

  RF_PowerCMP_t(scalartype loga_, scalartype logb_) :
    loga(loga_), logb(logb_) {}
  
  template<class T>
  RF_PowerCMP_t(const RF_PowerCMP_t<T>& other) :
    loga(other.loga), logb(other.logb) {}

  // template <template<class> class V, class T>
  // T operator()(const V<T> &logssb0){
  //   T logssb = logssb0(0);
  // template<class T>
  // T operator()(const T& logssb){
  scalartype operator()(const scalartype& logssb){
    return loga + invlogit(logb) * logssb;
    // return v;
  }

   THIS_RF_TO_AD(PowerCMP);

 
};
  
template<class scalartype_>
struct Rec_PowerCMP : RecruitmentNumeric<RF_PowerCMP_t<scalartype_> >  {
  SAMREC_TYPEDEFS(scalartype_);
 // Implement with known values when time permits!
  Rec_PowerCMP(scalartype loga, scalartype logb) :
    RecruitmentNumeric<RF_PowerCMP_t<scalartype> >(RF_PowerCMP_t<scalartype>(loga, logb)) {};

  template<class T>
  Rec_PowerCMP(const Rec_PowerCMP<T>& other) :
    RecruitmentNumeric<RF_PowerCMP_t<scalartype> >(other.f) {}

  THIS_REC_TO_AD(PowerCMP);
  
};



// Recruitment function 65
// Power function without compensatory mortality property

//     predN = rec_pars(0) + (exp(rec_pars(1))+1.0001) * log(thisSSB);
//     break;

 //   Se = exp(1.0 / (1.0 - (exp(newPar.rec_pars(1)) + 1.0001)) * (newPar.rec_pars(0) + logRecCorrection + log(lambda)));

//   mdsr = R_PosInf;


template<class scalartype_>
struct RF_PowerNCMP_t {
 SAMREC_TYPEDEFS(scalartype_);
  scalartype loga;
  scalartype logb;

  RF_PowerNCMP_t(scalartype loga_, scalartype logb_) :
    loga(loga_), logb(logb_) {}
  
  template<class T>
  RF_PowerNCMP_t(const RF_PowerNCMP_t<T>& other) :
    loga(other.loga), logb(other.logb) {}

  THIS_RF_TO_AD(PowerNCMP);
  
  // template <template<class> class V, class T>
  // T operator()(const V<T> &logssb0){
  //   T logssb = logssb0(0);
  // template<class T>
  // T operator()(const T& logssb){
  scalartype operator()(const scalartype& logssb){
    return loga + (exp(logb)+1.0001) * logssb;
    // return v;
  }
 
};
  
template<class scalartype_>
struct Rec_PowerNCMP : RecruitmentNumeric<RF_PowerNCMP_t<scalartype_> >  {
 SAMREC_TYPEDEFS(scalartype_);
  // Implement with known values when time permits!
  Rec_PowerNCMP(scalartype loga, scalartype logb) :
    RecruitmentNumeric<RF_PowerNCMP_t<scalartype> >(RF_PowerNCMP_t<scalartype>(loga, logb)) {};

  template<class T>
  Rec_PowerNCMP(const Rec_PowerNCMP<T>& other) :
    RecruitmentNumeric<RF_PowerNCMP_t<scalartype> >(other.f) {}

  THIS_REC_TO_AD(PowerNCMP);

};



// Recruitment function 66
// Shepherd

//     predN = rec_pars(0) + log(thisSSB) - log(1.0 + exp(exp(rec_pars(2)) * (log(thisSSB) - rec_pars(1))));

//   // Se = CppAD::CondExpGt((newPar.rec_pars(0) + logSPR), T(SAM_Zero),
  //   // 			  exp( newPar.rec_pars(1) + 1.0 / exp(newPar.rec_pars(2)) * log(fabs(exp(newPar.rec_pars(0)) * lambda - 1.0))),
  //   // 			  T(SAM_Zero));
  //   Se = exp( newPar.rec_pars(1) + 1.0 / exp(newPar.rec_pars(2)) * log(softmax(exp(newPar.rec_pars(0) + logRecCorrection) * lambda - 1.0,(T)SAM_Zero, (T)100.0)) );


  // xx = CppAD::CondExpGt(par.rec_pars(2),Type(0.0),exp(log((1.0+exp(par.rec_pars(2))) / (exp(par.rec_pars(2)) - 1.0)) / exp(par.rec_pars(2))), Type(1e-10));
  //   mdsr = -(-1.0 + (exp(par.rec_pars(2)) - 1.0)*pow(xx/exp(par.rec_pars(1)),exp(par.rec_pars(2))))*exp(par.rec_pars(0))/pow(pow(xx/exp(par.rec_pars(1)),exp(par.rec_pars(2))) + 1.0,2.0);


template<class scalartype_>
struct RF_Shepherd_t {
 SAMREC_TYPEDEFS(scalartype_);
  scalartype loga;
  scalartype logb;
  scalartype logg;

  RF_Shepherd_t(scalartype loga_, scalartype logb_, scalartype logg_) :
    loga(loga_), logb(logb_), logg(logg_) {}

  template<class T>
  RF_Shepherd_t(const RF_Shepherd_t<T>& other) :
    loga(other.loga), logb(other.logb), logg(other.logg) {}

    THIS_RF_TO_AD(Shepherd);
  
  // template <template<class> class V, class T>
  // T operator()(const V<T> &logssb0){
  //   T logssb = logssb0(0);
  // template<class T>
  // T operator()(const T& logssb){
  scalartype operator()(const scalartype& logssb){
    return loga + logssb - logspace_add_SAM(scalartype(0.0), exp(logg) * (logssb - logb));
    // return v;
  }
 
};
  
template<class scalartype_>
struct Rec_Shepherd : RecruitmentNumeric<RF_Shepherd_t<scalartype_> >  {
 SAMREC_TYPEDEFS(scalartype_);

  Rec_Shepherd(scalartype loga, scalartype logb, scalartype logg) :
    RecruitmentNumeric<RF_Shepherd_t<scalartype> >(RF_Shepherd_t<scalartype>(loga, logb, logg)) {};

  template<class T>
  Rec_Shepherd(const Rec_Shepherd<T>& other) :
    RecruitmentNumeric<RF_Shepherd_t<scalartype> >(other.f) {}

  THIS_REC_TO_AD(Shepherd);
  
};




// Recruitment function 67
// Hassel/Deriso

// //     predN = rec_pars(0)+log(thisSSB)-exp(rec_pars(2)) * log(1.0+exp(rec_pars(1) + rec_pars(2))*thisSSB); 

 //   // Handle negative values below!
  //   // Se = CppAD::CondExpGt((newPar.rec_pars(0) + logSPR), T(SAM_Zero),
  //   // 			  fabs(exp(exp(-newPar.rec_pars(2)) * (newPar.rec_pars(0) + logSPR)) - 1.0) * exp(-newPar.rec_pars(1)),
  //   // 			   T(SAM_Zero));
  //   Se = (exp(exp(-newPar.rec_pars(2)) * (newPar.rec_pars(0) + logRecCorrection + logSPR)) - 1.0) * exp(-newPar.rec_pars(1) - newPar.rec_pars(2));

// mdsr: Numeric



template<class scalartype_>
struct RF_HasselDeriso_t {
  SAMREC_TYPEDEFS(scalartype_);
 scalartype loga;
  scalartype logb;
  scalartype logg;

  RF_HasselDeriso_t(scalartype loga_, scalartype logb_, scalartype logg_) :
    loga(loga_), logb(logb_), logg(logg_) {}

  template<class T>
  RF_HasselDeriso_t(const RF_HasselDeriso_t<T>& other) :
    loga(other.loga), logb(other.logb), logg(other.logg) {}

  THIS_RF_TO_AD(HasselDeriso);
  
  // template <template<class> class V, class T>
  // T operator()(const V<T> &logssb0){
  //   T logssb = logssb0(0);
  // template<class T>
  // T operator()(const T& logssb){
  scalartype operator()(const scalartype& logssb){
    scalartype thisSSB = exp(logssb);
    return loga+logssb-exp(logg) * log(1.0+exp(logb + logg)*thisSSB);
    // return v;
  }
 
};
  
template<class scalartype_>
struct Rec_HasselDeriso : RecruitmentNumeric<RF_HasselDeriso_t<scalartype_> >  {
 SAMREC_TYPEDEFS(scalartype_);

  Rec_HasselDeriso(scalartype loga, scalartype logb, scalartype logg) :
    RecruitmentNumeric<RF_HasselDeriso_t<scalartype> >(RF_HasselDeriso_t<scalartype>(loga, logb, logg)) {};

  template<class T>
  Rec_HasselDeriso(const Rec_HasselDeriso<T>& other) :
    RecruitmentNumeric<RF_HasselDeriso_t<scalartype> >(other.f) {}

  THIS_REC_TO_AD(HasselDeriso);
  
};




// Recruitment function 68
// Saila-Lorda

//     predN = rec_pars(0)+exp(rec_pars(2)) * log(thisSSB) - exp(rec_pars(1))*thisSSB;

 //   Se = Se_sl(lambda, exp(newPar.rec_pars(0) + logRecCorrection), exp(newPar.rec_pars(1)), exp(newPar.rec_pars(2)));



template<class scalartype_>
struct RF_SailaLorda_t {
  SAMREC_TYPEDEFS(scalartype_);
 scalartype loga;
  scalartype logb;
  scalartype logg;

  RF_SailaLorda_t(scalartype loga_, scalartype logb_, scalartype logg_) :
    loga(loga_), logb(logb_), logg(logg_) {}

  template<class T>
  RF_SailaLorda_t(const RF_SailaLorda_t<T>& other) :
    loga(other.loga), logb(other.logb), logg(other.logg) {}

    THIS_RF_TO_AD(SailaLorda);
  
  // template <template<class> class V, class T>
  // T operator()(const V<T> &logssb0){
  //   T logssb = logssb0(0);
  // template<class T>
  // T operator()(const T& logssb){
  scalartype operator()(const scalartype& logssb){
    scalartype thisSSB = exp(logssb);
    return loga+exp(logg) * logssb - exp(logb)*thisSSB;
    // return v;
  }
 
};
  
template<class scalartype_>
struct Rec_SailaLorda : RecruitmentNumeric<RF_SailaLorda_t<scalartype_> >  {
 SAMREC_TYPEDEFS(scalartype_);

  Rec_SailaLorda(scalartype loga, scalartype logb, scalartype logg) :
    RecruitmentNumeric<RF_SailaLorda_t<scalartype> >(RF_SailaLorda_t<scalartype>(loga, logb, logg)) {};

  template<class T>
  Rec_SailaLorda(const Rec_SailaLorda<T>& other) :
    RecruitmentNumeric<RF_SailaLorda_t<scalartype> >(other.f) {}

  THIS_REC_TO_AD(SailaLorda);
  
};



// Recruitment function 69
// Sigmoidal Beverton-Holt

//     predN = rec_pars(0)+exp(rec_pars(2)) * log(thisSSB)-log(1.0+exp(rec_pars(1))*exp(exp(rec_pars(2)) * log(thisSSB)));

 //   Se = Se_sbh(lambda, exp(newPar.rec_pars(0) + logRecCorrection), exp(newPar.rec_pars(1)), exp(newPar.rec_pars(2)));
    // if(g < 1.0){
    //   return (1.0 - g) / b * lambertW_raw( b / (1.0 - g) * pow(a * l, 1 / (1.0 - g) ));
    // }
 


template<class scalartype_>
struct RF_SigmoidalBevHolt_t {
 SAMREC_TYPEDEFS(scalartype_);
  scalartype loga;
  scalartype logb;
  scalartype logg;

  RF_SigmoidalBevHolt_t(scalartype loga_, scalartype logb_, scalartype logg_) :
    loga(loga_), logb(logb_), logg(logg_) {}

  template<class T>
  RF_SigmoidalBevHolt_t(const RF_SigmoidalBevHolt_t<T>& other) :
    loga(other.loga), logb(other.logb), logg(other.logg) {}

  THIS_RF_TO_AD(SigmoidalBevHolt);
  
  // template <template<class> class V, class T>
  // T operator()(const V<T> &logssb0){
  //   T logssb = logssb0(0);
  // template<class T>
  // T operator()(const T& logssb){
  scalartype operator()(const scalartype& logssb){
    return loga+exp(logg) * logssb-log(1.0+exp(logb)*exp(exp(logg) * logssb));
    // return v;
  }
 
};
  
template<class scalartype_>
struct Rec_SigmoidalBevHolt : RecruitmentNumeric<RF_SigmoidalBevHolt_t<scalartype_> >  {
 SAMREC_TYPEDEFS(scalartype_);

  Rec_SigmoidalBevHolt(scalartype loga, scalartype logb, scalartype logg) :
    RecruitmentNumeric<RF_SigmoidalBevHolt_t<scalartype> >(RF_SigmoidalBevHolt_t<scalartype>(loga, logb, logg)) {};

  template<class T>
  Rec_SigmoidalBevHolt(const Rec_SigmoidalBevHolt<T>& other) :
    RecruitmentNumeric<RF_SigmoidalBevHolt_t<scalartype> >(other.f) {}

  THIS_REC_TO_AD(SigmoidalBevHolt);

  
};




// Recruitment function 90
// Spline with compensatory mortality property (non-increasing on log(R/SSB))

 //   predN(0) = log(thisSSB) + ibcdspline(log(thisSSB),
  // 					 (vector<Type>)(conf.constRecBreaks.template cast<Type>()),
  // 					 par.rec_pars);


template<class scalartype_>
struct RF_SplineCMP_t {
  SAMREC_TYPEDEFS(scalartype_);
 vectortype pars;
  vectortype knots;

  RF_SplineCMP_t(vectortype pars_, vectortype knots_) :
    pars(pars_), knots(knots_) {}

  template<class T>
  RF_SplineCMP_t(const RF_SplineCMP_t<T>& other) :
    pars(other.pars), knots(other.knots) {}

  THIS_RF_TO_AD(SplineCMP);
  
  // template <template<class> class V, class T>
  // T operator()(const V<T> &logssb0){
  //   T logssb = logssb0(0);
  //  template<class T>
  // T operator()(const T& logssb){
  //  vector<T> k2(knots);
  //   vector<T> p2(pars);
  scalartype operator()(const scalartype& logssb){
    return logssb + ibcdspline(logssb,knots,pars);
    // return v;
  }
 
};
  
template<class scalartype_>
struct Rec_SplineCMP : RecruitmentNumeric<RF_SplineCMP_t<scalartype_> >  {
 SAMREC_TYPEDEFS(scalartype_);

  Rec_SplineCMP(vectortype pars, vectortype knots) :
    RecruitmentNumeric<RF_SplineCMP_t<scalartype> >(RF_SplineCMP_t<scalartype>(pars,knots)) {}

  template<class T>
  Rec_SplineCMP(const Rec_SplineCMP<T>& other) :
    RecruitmentNumeric<RF_SplineCMP_t<scalartype> >(other.f) {}

  THIS_REC_TO_AD(SplineCMP);

  
};


template<class scalartype_>
struct RF_SplineConvexCompensatory_t {
 SAMREC_TYPEDEFS(scalartype_);
  vectortype pars;
  vectortype knots;		// These knots are for SSB, input knots are for logssb to be consistent with other splines.

  RF_SplineConvexCompensatory_t(vectortype pars_, vectortype knots_) :
    pars(pars_), knots(exp(knots_)) {}

  template<class T>
  RF_SplineConvexCompensatory_t(const RF_SplineConvexCompensatory_t<T>& other) :
    pars(other.pars), knots(other.knots) {}

  THIS_RF_TO_AD(SplineConvexCompensatory);
  
  // template <template<class> class V, class T>
  // T operator()(const V<T> &logssb0){
  //   T logssb = logssb0(0);
  // template<class T>
  // T operator()(const T& logssb){
  //   vector<T> k2(knots);
  //   vector<T> p2(pars);
  scalartype operator()(const scalartype& logssb){
    return logssb + iibcispline(exp(logssb),knots,pars);
    // return v;
  }
 
};
  
template<class scalartype_>
struct Rec_SplineConvexCompensatory : RecruitmentNumeric<RF_SplineConvexCompensatory_t<scalartype_> >  {
 SAMREC_TYPEDEFS(scalartype_);

  Rec_SplineConvexCompensatory(vectortype pars, vectortype knots) :
    RecruitmentNumeric<RF_SplineConvexCompensatory_t<scalartype> >(RF_SplineConvexCompensatory_t<scalartype>(pars,knots)) {}

  template<class T>
  Rec_SplineConvexCompensatory(const Rec_SplineConvexCompensatory<T>& other) :
    RecruitmentNumeric<RF_SplineConvexCompensatory_t<scalartype> >(other.f) {}

  THIS_REC_TO_AD(SplineConvexCompensatory);
  
};




// Recruitment function 91
// Smooth spline (integrated spline on log(R/SSB))

  //   predN(0) = log(thisSSB) + ibcspline(log(thisSSB),
  // 					(vector<Type>)(conf.constRecBreaks.template cast<Type>()),
  // 					par.rec_pars);


template<class scalartype_>
struct RF_SplineSmooth_t {
  SAMREC_TYPEDEFS(scalartype_);
 vectortype pars;
  vectortype knots;

  RF_SplineSmooth_t(vectortype pars_, vectortype knots_) :
    pars(pars_), knots(knots_) {}

  template<class T>
  RF_SplineSmooth_t(const RF_SplineSmooth_t<T>& other) :
    pars(other.pars), knots(other.knots) {}

  THIS_RF_TO_AD(SplineSmooth);
  
  // template <template<class> class V, class T>
  // T operator()(const V<T> &logssb0){
  //   T logssb = logssb0(0);
  // template<class T>
  // T operator()(const T& logssb){
  //   vector<T> k2(knots);
  //   vector<T> p2(pars.segment(0, pars.size()-1));
  scalartype operator()(const scalartype& logssb){
    scalartype mu(pars(pars.size()-1));
    return logssb + ibcspline(logssb,knots,(vectortype)pars.segment(0, pars.size()-1)) + mu;
    // return v;
  }
 
};
  
template<class scalartype_>
struct Rec_SplineSmooth : RecruitmentNumeric<RF_SplineSmooth_t<scalartype_> >  {
 SAMREC_TYPEDEFS(scalartype_);

  Rec_SplineSmooth(vectortype pars, vectortype knots) :
    RecruitmentNumeric<RF_SplineSmooth_t<scalartype> >(RF_SplineSmooth_t<scalartype>(pars,knots)) {}

  template<class T>
  Rec_SplineSmooth(const Rec_SplineSmooth<T>& other) :
    RecruitmentNumeric<RF_SplineSmooth_t<scalartype> >(other.f) {}

  THIS_REC_TO_AD(SplineSmooth);
  
};


// Recruitment function 92
// General spline (on log(R/SSB))

 //   predN(0) = log(thisSSB) + bcspline(log(thisSSB),
  // 				       (vector<Type>)(conf.constRecBreaks.template cast<Type>()),
  // 				       par.rec_pars);



template<class scalartype_>
struct RF_SplineGeneral_t {
  SAMREC_TYPEDEFS(scalartype_);
 vectortype pars;
  vectortype knots;

  RF_SplineGeneral_t(vectortype pars_, vectortype knots_) :
    pars(pars_), knots(knots_) {}

  template<class T>
  RF_SplineGeneral_t(const RF_SplineGeneral_t<T>& other) :
    pars(other.pars), knots(other.knots) {}

THIS_RF_TO_AD(SplineGeneral);
  
  // template <template<class> class V, class T>
  // T operator()(const V<T> &logssb0){
  //   T logssb = logssb0(0);
  // template<class T>
  // T operator()(const T& logssb){
  scalartype operator()(const scalartype& logssb){
    // vector k2(knots);
    //   vector<T> p2(pars.segment(0, pars.size()-1));
    scalartype mu(pars(pars.size()-1));
    return logssb + bcspline(logssb,knots,(vectortype)pars.segment(0, pars.size()-1)) + mu;
    // return v;
  }
 
};
  
template<class scalartype_>
struct Rec_SplineGeneral : RecruitmentNumeric<RF_SplineGeneral_t<scalartype_> >  {
 SAMREC_TYPEDEFS(scalartype_);

  Rec_SplineGeneral(vectortype pars, vectortype knots) :
    RecruitmentNumeric<RF_SplineGeneral_t<scalartype> >(RF_SplineGeneral_t<scalartype>(pars,knots)) {}

  template<class T>
  Rec_SplineGeneral(const Rec_SplineGeneral<T>& other) :
    RecruitmentNumeric<RF_SplineGeneral_t<scalartype> >(other.f) {}

  THIS_REC_TO_AD(SplineGeneral);
};



//////////////////////////////////////////////////////////////////////////////////////////

// Numerical Ricker for testing


template<class scalartype_>
struct RF_NumRicker_t {
   SAMREC_TYPEDEFS(scalartype_);
  scalartype loga;
  scalartype logb;

  RF_NumRicker_t(scalartype loga_, scalartype logb_) :
    loga(loga_), logb(logb_) {}

  template<class T>
  RF_NumRicker_t(const RF_NumRicker_t<T>& other) :
    loga(other.loga), logb(other.logb) {}

  THIS_RF_TO_AD(NumRicker);
  
  // template <template<class> class V, class T>
  // T operator()(const V<T> &logssb0){
  //   T logssb = logssb0(0);
  // template<class T>
  // T operator()(const T& logssb){
  scalartype operator()(const scalartype& logssb){
    return loga + logssb - exp(logb + logssb);
  }
 
};
  
template<class scalartype_>
struct Rec_NumRicker : RecruitmentNumeric<RF_NumRicker_t<scalartype_> >  {
 SAMREC_TYPEDEFS(scalartype_);

  Rec_NumRicker(scalartype loga, scalartype logb) :
    RecruitmentNumeric<RF_NumRicker_t<scalartype> >(RF_NumRicker_t<scalartype>(loga, logb)) {};

  template<class T>
  Rec_NumRicker(const Rec_NumRicker<T>& other) :
    RecruitmentNumeric<RF_NumRicker_t<scalartype> >(other.f) {}

  THIS_REC_TO_AD(NumRicker);
  
};

// Numerical Beverton Holt


template<class scalartype_>
struct RF_NumBevHolt_t {
 SAMREC_TYPEDEFS(scalartype_);
  scalartype loga;
  scalartype logb;

  RF_NumBevHolt_t(scalartype loga_, scalartype logb_) :
    loga(loga_), logb(logb_) {}

  template<class T>
  RF_NumBevHolt_t(const RF_NumBevHolt_t<T>& other) :
    loga(other.loga), logb(other.logb) {}

  THIS_RF_TO_AD(NumBevHolt);
  
  // template <template<class> class V, class T>
  // T operator()(const V<T> &logssb0){
    // T logssb = logssb0(0);
  // template<class T>
  // T operator()(const T& logssb){

  scalartype operator()(const scalartype& logssb){
    return loga + logssb - logspace_add_SAM(scalartype(0.0), logb + logssb);
  }
 
};
  
template<class scalartype_>
struct Rec_NumBevHolt : RecruitmentNumeric<RF_NumBevHolt_t<scalartype_> >  {
 SAMREC_TYPEDEFS(scalartype_);

  Rec_NumBevHolt(scalartype loga, scalartype logb) :
    RecruitmentNumeric<RF_NumBevHolt_t<scalartype> >(RF_NumBevHolt_t<scalartype>(loga, logb)) {};

  template<class T>
  Rec_NumBevHolt(const Rec_NumBevHolt<T>& other) :
    RecruitmentNumeric<RF_NumBevHolt_t<scalartype> >(other.f) {}

  THIS_REC_TO_AD(NumBevHolt);
  
};



///////////////////////////////////////////////////////////////////////////////////////////////////

template<class Type>
Recruitment<Type> makeRecruitmentFunction(const confSet& conf, const paraSet<Type>& par){
  RecruitmentModel rm = static_cast<RecruitmentModel>(conf.stockRecruitmentModelCode);
  Recruitment<Type> r;

  if(rm == RecruitmentModel::NoRecruit){
    r = Recruitment<Type>("zero",new Rec_None<Type>());
    
////////////////////////////////////////// The Beginning //////////////////////////////////////////
    
  }else if(rm == RecruitmentModel::LogRandomWalk){
    if(par.rec_pars.size() != 0)
      Rf_error("The random walk recruitment should not have any parameters.");
    r = Recruitment<Type>("log-random walk",new Rec_LogRW<Type>());
  }else if(rm == RecruitmentModel::Ricker){
    if(par.rec_pars.size() != 2)
      Rf_error("The Ricker recruitment must have two parameters.");
    r = Recruitment<Type>("Ricker",new Rec_Ricker<Type>(par.rec_pars(0), par.rec_pars(1)));
  }else if(rm == RecruitmentModel::BevertonHolt){
    if(par.rec_pars.size() != 2)
      Rf_error("The Beverton Holt recruitment must have two parameters.");
    r = Recruitment<Type>("Beverton-Holt",new Rec_BevertonHolt<Type>(par.rec_pars(0), par.rec_pars(1)));
  }else if(rm == RecruitmentModel::ConstantMean){
    if(par.rec_pars.size() != conf.constRecBreaks.size() + 1)
      Rf_error("The constant mean recruitment should have one more parameter than constRecBreaks.");
    r = Recruitment<Type>("constant mean", new Rec_ConstantMean<Type>(par.rec_pars, conf.constRecBreaks));
  }else if(rm == RecruitmentModel::LogisticHockeyStick){
    if(par.rec_pars.size() != 3)
      Rf_error("The logistic hockey stick recruitment should have three parameters.");
    r = Recruitment<Type>("logistic hockey stick",new Rec_LogisticHockeyStick<Type>(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)));
  }else if(rm == RecruitmentModel::HockeyStick){
   if(par.rec_pars.size() != 2)
      Rf_error("The hockey stick recruitment should have two parameters.");
   r = Recruitment<Type>("hockey stick",new Rec_HockeyStick<Type>(par.rec_pars(0), par.rec_pars(1)));
  }else if(rm == RecruitmentModel::LogAR1){
   if(par.rec_pars.size() != 2)
      Rf_error("The log-AR(1) recruitment should have two parameters.");
   r = Recruitment<Type>("log-AR(1)",new Rec_LogAR1<Type>(par.rec_pars(0), par.rec_pars(1)));
  }else if(rm == RecruitmentModel::BentHyperbola){
   if(par.rec_pars.size() != 3)
      Rf_error("The bent hyperbola recruitment should have three parameters.");
   r = Recruitment<Type>("bent hyperbola",new Rec_BentHyperbola<Type>(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)));
  }else if(rm == RecruitmentModel::Power_CMP){
   if(par.rec_pars.size() != 2)
      Rf_error("The power law recruitment should have two parameters.");
   r = Recruitment<Type>("CMP power law",new Rec_PowerCMP<Type>(par.rec_pars(0), par.rec_pars(1)));
  }else if(rm == RecruitmentModel::Power_NCMP){
   if(par.rec_pars.size() != 2)
      Rf_error("The power law recruitment should have two parameters.");
   r = Recruitment<Type>("non-CMP power law",new Rec_PowerNCMP<Type>(par.rec_pars(0), par.rec_pars(1)));

//////////////////////////////////////// 3 parameter models ///////////////////////////////////////
    
  }else if(rm == RecruitmentModel::Shepherd){
   if(par.rec_pars.size() != 3)
      Rf_error("The Shepherd recruitment should have three parameters.");
   r = Recruitment<Type>("Shepherd",new Rec_Shepherd<Type>(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)));
  }else if(rm == RecruitmentModel::Hassel_Deriso){
   if(par.rec_pars.size() != 3)
      Rf_error("The Hassel/Deriso recruitment should have three parameters.");
   r = Recruitment<Type>("Hassel/Deriso",new Rec_HasselDeriso<Type>(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)));
  }else if(rm == RecruitmentModel::SailaLorda){
   if(par.rec_pars.size() != 3)
      Rf_error("The Saila-Lorda recruitment should have three parameters.");
   r = Recruitment<Type>("Saila-Lorda",new Rec_SailaLorda<Type>(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)));
  }else if(rm == RecruitmentModel::SigmoidalBevertonHolt){
   if(par.rec_pars.size() != 3)
     Rf_error("The sigmoidal Beverton-Holt recruitment should have three parameters.");
   r = Recruitment<Type>("sigmoidal Beverton-Holt",new Rec_SigmoidalBevHolt<Type>(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)));

//////////////////////////////////////// Spline recruitment ///////////////////////////////////////
    
  }else if(rm == RecruitmentModel::Spline_CMP){
    if(par.rec_pars.size() != conf.constRecBreaks.size() + 1)
      Rf_error("The spline recruitment should have one more parameter than constRecBreaks.");
    r = Recruitment<Type>("CMP spline",new Rec_SplineCMP<Type>(par.rec_pars, conf.constRecBreaks.template cast<Type>()));
  }else if(rm == RecruitmentModel::Spline_Smooth){
    if(par.rec_pars.size() != conf.constRecBreaks.size() + 1)
      Rf_error("The spline recruitment should have one more parameter than constRecBreaks.");
    r = Recruitment<Type>("smooth spline",new Rec_SplineSmooth<Type>(par.rec_pars, conf.constRecBreaks.template cast<Type>()));
  }else if(rm == RecruitmentModel::Spline_General){
    if(par.rec_pars.size() != conf.constRecBreaks.size() + 1)
      Rf_error("The spline recruitment should have one more parameter than constRecBreaks.");
    r = Recruitment<Type>("unrestricted spline",new Rec_SplineGeneral<Type>(par.rec_pars, conf.constRecBreaks.template cast<Type>()));

  }else if(rm == RecruitmentModel::Spline_ConvexCompensatory){
    if(par.rec_pars.size() != conf.constRecBreaks.size() + 1)
      Rf_error("The spline recruitment should have one more parameter than constRecBreaks.");
    r = Recruitment<Type>("convex-compensatory spline",new Rec_SplineConvexCompensatory<Type>(par.rec_pars, conf.constRecBreaks.template cast<Type>()));


//////////////////////// S/(d+S) type Depensatory recruitment models //////////////////////////////

  }else if(rm == RecruitmentModel::Depensatory_B_Ricker){
   if(par.rec_pars.size() != 3)
     Rf_error("The depensatory B Ricker recruitment should have three parameters.");
   r = Recruitment<Type>("type B depensatory Ricker",Rec_DepensatoryB(Rec_Ricker<Type>(par.rec_pars(0), par.rec_pars(1)),par.rec_pars(2)));
  }else if(rm == RecruitmentModel::Depensatory_B_BevertonHolt){
    if(par.rec_pars.size() != 3)
     Rf_error("The depensatory B Beverton-Holt recruitment should have three parameters.");
    r = Recruitment<Type>("type B depensatory Beverton-Holt",Rec_DepensatoryB(Rec_BevertonHolt<Type>(par.rec_pars(0), par.rec_pars(1)),par.rec_pars(2)));
  }else if(rm == RecruitmentModel::Depensatory_B_LogisticHockeyStick){
     if(par.rec_pars.size() != 4)
     Rf_error("The depensatory B logistic hockey stick recruitment should have four parameters.");
     r = Recruitment<Type>("type B depensatory logistic hockey stick",Rec_DepensatoryB(Rec_LogisticHockeyStick<Type>(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)),par.rec_pars(3)));
  }else if(rm == RecruitmentModel::Depensatory_B_HockeyStick){
    if(par.rec_pars.size() != 3)
      Rf_error("The depensatory B hockey stick recruitment should have three parameters.");
    r = Recruitment<Type>("type B depensatory hockey stick",Rec_DepensatoryB(Rec_HockeyStick<Type>(par.rec_pars(0), par.rec_pars(1)),par.rec_pars(2)));
  }else if(rm == RecruitmentModel::Depensatory_B_BentHyperbola){
    if(par.rec_pars.size() != 4)
     Rf_error("The depensatory B bent hyperbola recruitment should have four parameters.");
    r = Recruitment<Type>("type B depensatory bent hyperbola",Rec_DepensatoryB(Rec_BentHyperbola<Type>(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)),par.rec_pars(3)));
  }else if(rm == RecruitmentModel::Depensatory_B_Power){
    if(par.rec_pars.size() != 3)
     Rf_error("The depensatory B power law recruitment should have three parameters.");
    r = Recruitment<Type>("type B depensatory CMP power law",Rec_DepensatoryB(Rec_PowerCMP<Type>(par.rec_pars(0), par.rec_pars(1)),par.rec_pars(2)));
  }else if(rm == RecruitmentModel::Depensatory_B_Shepherd){
    if(par.rec_pars.size() != 4)
      Rf_error("The depensatory B Shepherd recruitment should have four parameters.");
    r = Recruitment<Type>("type B depensatory Shepherd",Rec_DepensatoryB(Rec_Shepherd<Type>(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)),par.rec_pars(3)));
  }else if(rm == RecruitmentModel::Depensatory_B_Hassel_Deriso){
       if(par.rec_pars.size() != 4)
     Rf_error("The depensatory B Hassel/Deriso recruitment should have four parameters.");
       r = Recruitment<Type>("type B depensatory Hassel/Deriso",Rec_DepensatoryB(Rec_HasselDeriso<Type>(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)),par.rec_pars(3)));
  }else if(rm == RecruitmentModel::Depensatory_B_Spline_CMP){
    if(par.rec_pars.size() != conf.constRecBreaks.size() + 2 && conf.constRecBreaks.size() >= 3)
      Rf_error("The depensatory B spline recruitment should have two parameters more than the number of knots which should be at least 3.");
    r = Recruitment<Type>("type B depensatory CMP spline",Rec_DepensatoryB(Rec_SplineCMP<Type>(par.rec_pars.segment(0,par.rec_pars.size()-1), conf.constRecBreaks),par.rec_pars(par.rec_pars.size()-1)));

 }else if(rm == RecruitmentModel::Depensatory_B_Spline_ConvexCompensatory){
    if(par.rec_pars.size() != conf.constRecBreaks.size() + 2 && conf.constRecBreaks.size() >= 3)
      Rf_error("The depensatory B spline recruitment should have two parameters more than the number of knots which should be at least 3.");
    r = Recruitment<Type>("type B depensatory convex compensatory spline",Rec_DepensatoryB(Rec_SplineConvexCompensatory<Type>(par.rec_pars.segment(0,par.rec_pars.size()-1), conf.constRecBreaks),par.rec_pars(par.rec_pars.size()-1)));

    
//////////////////////// 1/(1+exp(-e * (S-d))) type Depensatory recruitment models //////////////////////////////

  }else if(rm == RecruitmentModel::Depensatory_C_Ricker){
    if(par.rec_pars.size() != 4)
      Rf_error("The depensatory C Ricker recruitment should have four parameters.");
    r = Recruitment<Type>("type C depensatory Ricker",Rec_DepensatoryC(Rec_Ricker<Type>(par.rec_pars(0), par.rec_pars(1)),par.rec_pars(2), par.rec_pars(3)));
    
  }else if(rm == RecruitmentModel::Depensatory_C_BevertonHolt){
    if(par.rec_pars.size() != 4)
      Rf_error("The depensatory C Beverton-Holt recruitment should have four parameters.");
    r = Recruitment<Type>("type C depensatory Beverton-Holt",Rec_DepensatoryC(Rec_BevertonHolt<Type>(par.rec_pars(0), par.rec_pars(1)),par.rec_pars(2), par.rec_pars(3)));

   }else if(rm == RecruitmentModel::Depensatory_C_Spline_CMP){
    if(par.rec_pars.size() != conf.constRecBreaks.size() + 3 && conf.constRecBreaks.size() >= 3)
      Rf_error("The depensatory C spline recruitment should have three parameters more than the number of knots which should be at least 3.");
    r = Recruitment<Type>("type C depensatory CMP spline",Rec_DepensatoryC(Rec_SplineCMP<Type>(par.rec_pars.segment(0,par.rec_pars.size()-2), conf.constRecBreaks),par.rec_pars(par.rec_pars.size()-2),par.rec_pars(par.rec_pars.size()-1)));

  }else if(rm == RecruitmentModel::Depensatory_C_Spline_ConvexCompensatory){
    if(par.rec_pars.size() != conf.constRecBreaks.size() + 3 && conf.constRecBreaks.size() >= 3)
      Rf_error("The depensatory C spline recruitment should have three parameters more than the number of knots which should be at least 3.");
    r = Recruitment<Type>("type C depensatory convex compensatory spline",Rec_DepensatoryC(Rec_SplineConvexCompensatory<Type>(par.rec_pars.segment(0,par.rec_pars.size()-2), conf.constRecBreaks),par.rec_pars(par.rec_pars.size()-2),par.rec_pars(par.rec_pars.size()-1)));

    
///////////////////////////////////////////// For testing /////////////////////////////////////////////
  }else if(rm == RecruitmentModel::Num_Ricker){
    if(par.rec_pars.size() != 2)
      Rf_error("The Numeric Ricker recruitment must have two parameters.");
    r = Recruitment<Type>("numeric Ricker",new Rec_NumRicker<Type>(par.rec_pars(0), par.rec_pars(1)));

  }else if(rm == RecruitmentModel::Num_BevertonHolt){
    if(par.rec_pars.size() != 2)
      Rf_error("The Numeric Beverton Holt recruitment must have two parameters.");
    r = Recruitment<Type>("numeric Beverton-Holt",new Rec_NumBevHolt<Type>(par.rec_pars(0), par.rec_pars(1)));

///////////////////////////////////////////// The End /////////////////////////////////////////////
    
  }else{
    Rf_error("Stock-recruitment model code not implemented.");
  }

  return r;
  
}

  
  
#endif
