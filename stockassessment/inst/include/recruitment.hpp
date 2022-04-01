#pragma once
#ifndef SAM_RECRUITMENT_HPP
#define SAM_RECRUITMENT_HPP

#include <memory>

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
		       ICESforecast = -2, // Known, Implemented
		       NoRecruit = -1, // Known, Implemented
		       LogRandomWalk = 0, // Known, Implemented
		       Ricker = 1,	  // Known, Implemented
		       BevertonHolt = 2,  // Known, Implemented
		       ConstantMean = 3,  // Known, Implemented
		       LogisticHockeyStick = 60,     // Mix, Implemented
		       HockeyStick = 61,	     // Known, Implemented
		       LogAR1 = 62,		     // Known, Implemented
		       BentHyperbola = 63,	     // Mix
		       Power_CMP = 64,		     // Known
		       Power_NCMP = 65,		     // Known
		       Shepherd = 66,		     // Mix		       
		       Hassel_Deriso = 67,	     // Mix
		       SailaLorda = 68,		     // Mix
		       SigmoidalBevertonHolt = 69,   // Mix
		       Spline_CMP = 90,		     // Numeric
		       Spline_Smooth = 91,	     // Numeric
		       Spline_General = 92,	     // Numeric
		       Spline_ConvexCompensatory = 93,
		       Depensatory_B_Ricker = 201,	     // Numeric, Implemented
		       Depensatory_B_BevertonHolt = 202,	     // Numeric, Implemented
		       Depensatory_B_LogisticHockeyStick = 260, // Numeric, Implemented
		       Depensatory_B_HockeyStick = 261,	     // Numeric, Implemented
		       Depensatory_B_BentHyperbola = 263,	     // Numeric, Implemented
		       Depensatory_B_Power = 264,	     // Numeric, Implemented
		       Depensatory_B_Shepherd = 266,	     // Numeric, Implemented
		       Depensatory_B_Hassel_Deriso = 267,	     // Numeric, Implemented
		       Depensatory_B_Spline_CMP = 290,	     // Numeric
		       Depensatory_C_Ricker = 401,	     // Numeric, Implemented
		       Depensatory_C_BevertonHolt = 402,	     // Numeric, Implemented
		    
};

// Convert newton::vector to utils::vector
template<class Type>
vector<Type> n2u(const newton::vector<Type>& x){
  vector<Type> r(x);
  return r;
}

template<class Type>
vector<Type> n2u(vector<Type>& x){
  return x;
}
template<class Type>
vector<Type> n2u(const vector<Type>& x){
  return x;
}

template<class Type>
struct RecruitmentWorker {
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
class Recruitment {
public:
  // RecruitmentWorker<Type>* ptr;
  std::shared_ptr<RecruitmentWorker<Type> > ptr;
  explicit Recruitment(RecruitmentWorker<Type>* p=NULL) : ptr(p) {};
  virtual ~Recruitment() { ptr.reset(); };

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


/*
The functor f should return logR for a given logSSB
 */

template<class Functor>
struct WrapSR {
  Functor f;
  WrapSR(Functor f_) : f(f_) {};
  
  // template<class T>
  // T operator()(vector<T> x){
  //    return -f(x);
  // }
  template <template<class> class V, class T>
  T operator()(const V<T> &x){
  //  template<class T>
  // T operator()(const vector<T> &x) {
    return -f(x);
  }
};

template<class Functor>
struct Exp {
  Functor f;

  Exp(Functor f_) : f(f_) {};

  // template<class T>
  // T operator()(vector<T> x){
  //   return exp(f(x));
  // }
  template <template<class> class V, class T>
  T operator()(const V<T> &x){
  // template<class T>
  // T operator()(const vector<T> &x) {
    return exp(f(x));
  }
};

template<class Functor>
struct diffSR {
  Exp<Functor> f;

  diffSR(Functor f_) : f(f_) {};
  
  // template<class T>
  // T operator()(vector<T> x){
  //    vector<T> g = autodiff::gradient(f, x);
  //   return (T)g(0) / exp(x(0));
  // }
  template <template<class> class V, class T>
  T operator()(const V<T> &x){
  // template<class T>
  // T operator()(const vector<T> &x) {
    vector<T> g = autodiff::gradient(f, n2u(x));
    return (T)g(0) / exp(x(0));
  }
};


template<class Functor>
struct WrapDiffSR {
  diffSR<Functor> dF;

  WrapDiffSR(Functor f_) : dF(f_) {};
  
  // template<class T>
  // T operator()(vector<T> x){
  //    return -dF(x);
  // }
  template <template<class> class V, class T>
  T operator()(const V<T> &x){
  // template<class T>
  // T operator()(const vector<T> &x) {
     return -dF(x);
  }
};

template<class Type, class Functor>
struct WrapEquiS {
  Functor f;
  Type logLambda;

  WrapEquiS(Functor f_, const Type& l) : f(f_), logLambda(l) {};

  // template<class T>
  // T operator()(vector<T> x){
  //   T v = (T)logLambda + f(x) - x(0);
  //   return v * v;
  // }
  template <template<class> class V, class T>
  T operator()(const V<T> &x){
  // template<class T>
  // T operator()(const vector<T> &x) {
    T v = (T)logLambda + f(x) - x(0);
    return v * v;
  }  
};

template<class Type, class Functor>
struct RecruitmentNumeric : RecruitmentWorker<Type> {

  Functor f;
  Type x0;
  
  RecruitmentNumeric(Functor f_) : f(f_), x0(10.0) {};
  RecruitmentNumeric(Functor f_, Type x0_) : f(f_), x0(x0_) {};  
  
  Type operator()(Type logssb, Type lastR, Type year){
    vector<Type> ls(1); ls(0) = logssb;
    return f(ls);
  }
  
  virtual Type logSe(Type logLambda){
    WrapEquiS<Type, Functor> fx(f, logLambda);
    vector<Type> x0v(1); x0v(0) = x0;
    vector<Type> v = newton::Newton(fx,x0v);  
    return (Type)v(0);
  }
  virtual Type dSR(Type logssb){
    vector<Type> ls(1); ls(0) = logssb;
    diffSR<Functor> dF(f);
    Type r = dF(ls);
    return r;
  }
  virtual Type logSAtMaxR(){
    WrapSR<Functor> fx(f);
    vector<Type> x0v(1); x0v(0) = x0;
    vector<Type> v = newton::Newton(fx,x0v);  
    return v(0);
  }
  virtual Type logMaxR(){   
    vector<Type> v(1); v(0) = logSAtMaxR();    
    return f(v);
  }
  virtual Type maxGradient(){
    WrapDiffSR<Functor> fx(f);
    vector<Type> x0v(1); x0v(0) = x0;
    vector<Type> v = newton::Newton(fx,x0v);      
    return dSR(v(0));
  }  
};




template <template<class> class Functor, class Type>
struct WrapDepensatoryA {
  Functor<Type> f;
  Type logd;

  WrapDepensatoryA(Functor<Type> f_, const Type& logd_) : f(f_), logd(logd_) {};

  // template<class T>
  // T operator()(const vector<T> &x) {
  template <template<class> class V, class T>
  T operator()(const V<T> &x){
    Functor<T>f2(f);
    T logssb = exp((T)logd) * x(0);
    T v = f2(logssb, (T)R_NaReal, (T)R_NaReal);
    return v;
  }  
};

template <template<class> class Functor, class Type>
RecruitmentNumeric<Type, WrapDepensatoryA<Functor, Type> >* Rec_DepensatoryA(Functor<Type> SR, Type logd){
  WrapDepensatoryA<Functor, Type> DepSR(SR, logd);
  return new RecruitmentNumeric<Type, WrapDepensatoryA<Functor, Type> >(DepSR, logd + 2.0);
}




template <template<class> class Functor, class Type>
struct WrapDepensatoryB {
  Functor<Type> f;
  Type logd;

  WrapDepensatoryB(Functor<Type> f_, const Type& logd_) : f(f_), logd(logd_) {};

  // template<class T>
  // T operator()(const vector<T> &x) {
  template <template<class> class V, class T>
  T operator()(const V<T> &x){
    Functor<T>f2(f);
    T logssb = x(0);
    T v = f2(logssb, (T)R_NaReal, (T)R_NaReal) + logssb - logspace_add2(logssb, (T)logd);
    return v;
  }  
};

template <template<class> class Functor, class Type>
RecruitmentNumeric<Type, WrapDepensatoryB<Functor, Type> >* Rec_DepensatoryB(Functor<Type> SR, Type logd){
  WrapDepensatoryB<Functor, Type> DepSR(SR, logd);
  return new RecruitmentNumeric<Type, WrapDepensatoryB<Functor, Type> >(DepSR, logd + 2.0);
}




template <template<class> class Functor, class Type>
struct WrapDepensatoryC {
  Functor<Type> f;
  Type logd;
  Type logl;

  WrapDepensatoryC(Functor<Type> f_, const Type& logd_, const Type& logl_) : f(f_), logd(logd_), logl(logl_) {};

  // template<class T>
  // T operator()(const vector<T> &x) {
  template <template<class> class V, class T>
  T operator()(const V<T> &x){
    Functor<T>f2(f);
    T logssb = x(0);
    T v = f2(logssb, (T)R_NaReal, (T)R_NaReal) - logspace_add2((T)0.0, -exp((T)logl) * (exp(logssb) - exp((T)logd)));
    return v;
  }  
};

template <template<class> class Functor, class Type>
RecruitmentNumeric<Type, WrapDepensatoryC<Functor, Type> >* Rec_DepensatoryC(Functor<Type> SR, Type logd, Type logl){
  WrapDepensatoryC<Functor, Type> DepSR(SR, logd, logl);
  return new RecruitmentNumeric<Type, WrapDepensatoryC<Functor, Type> >(DepSR, logd + 2.0);
}



// Recruitment function -2
// No recruitment
template<class Type>
struct Rec_ICESforecast : RecruitmentWorker<Type> {
  Type logMean;
  Type logSd;
  Type operator()(Type logssb, Type lastLogR, Type year){
    return logMean;
  }
  Type logSe(Type logLambda){
    return logLambda + logMean;
  }
  Type dSR(Type logssb){
    return 0.0;
  }
  Type logSAtMaxR(){
    return R_NaReal;
  }
   Type logMaxR(){
    return logMean;
  }
  Type maxGradient(){
    return 0.0;
  }
};


// Recruitment function -1
// No recruitment
template<class Type>
struct Rec_None : RecruitmentWorker<Type> {
  Type operator()(Type logssb, Type lastLogR, Type year){
    return R_NegInf;
  }
  Type logSe(Type logLambda){
    return R_NegInf;
  }
  Type dSR(Type logssb){
    return R_NaReal;
  }
  Type logSAtMaxR(){
    return R_NaReal;
  }
   Type logMaxR(){
    return R_NegInf;
  }
  Type maxGradient(){
    return R_NaReal;
  }
};

// Recruitment function 0
// Random walk
template<class Type>
struct Rec_LogRW : RecruitmentWorker<Type> {
  Type operator()(Type logssb, Type lastLogR, Type year){
    return lastLogR;
  }
  Type logSe(Type logLambda){
    return R_NaReal;
  }
  Type dSR(Type logssb){
    return R_NaReal;
  }
  Type logSAtMaxR(){
    return R_NaReal;
  }
   Type logMaxR(){
    return R_NaReal;
  }
  Type maxGradient(){
    return R_NaReal;
  }
};

// Recruitment function 1
// Ricker
template<class Type>
struct Rec_Ricker : RecruitmentWorker<Type> {

  Type loga;
  Type logb;

  Rec_Ricker(Type la, Type lb) : loga(la), logb(lb) {};
  template<class T>
  Rec_Ricker(const Rec_Ricker<T>& other) : loga(other.loga), logb(other.logb) {}
  
  Type operator()(Type logssb, Type lastR, Type year){
    return loga + logssb - exp(logb + logssb);
  }
  
  Type logSe(Type logLambda){
    Type v = TMBad::CondExpLe(loga+logLambda, Type(0.0), Type(0.0), (loga + logLambda) * exp(-logb));
    return log(v);
  }
  Type dSR(Type logssb){
    Type a = exp(loga);
    Type ee = -exp(logb + logssb);
    return -a * exp(ee) * (-ee - 1.0);
  }
  Type logSAtMaxR(){
    return -logb;
  }
  Type logMaxR(){
    return loga - 1.0 - logb;
  }
  Type maxGradient(){
    return exp(loga);
  }
};


// Recruitment function 2
// Beverton Holt
template<class Type>
struct Rec_BevertonHolt : RecruitmentWorker<Type> {

  Type loga;
  Type logb;

  Rec_BevertonHolt(Type la, Type lb) : loga(la), logb(lb) {};
  template<class T>
  Rec_BevertonHolt(const Rec_BevertonHolt<T>& other) : loga(other.loga), logb(other.logb) {}

  Type operator()(Type logssb, Type lastR, Type year){
    return loga + logssb - logspace_add2(Type(0.0),logb + logssb);
  }
  
  Type logSe(Type logLambda){
    Type v = TMBad::CondExpLe(loga+logLambda, Type(0.0), Type(0.0), (exp(loga + logLambda) - 1.0) * exp(-logb));
    return log(v);
  }
  Type dSR(Type logssb){
    Type a = exp(loga);
    Type ee = exp(logb + logssb) + 1.0;
    return a / (ee * ee);
  }
  Type logSAtMaxR(){
    return R_PosInf;
  }
  Type logMaxR(){
    return loga - logb;
  }
  Type maxGradient(){
    return exp(loga);
  }
};

// Recruitment function 3
// Constant mean
template<class Type>
struct Rec_ConstantMean : RecruitmentWorker<Type> {

  vector<Type> logRvalue;

  vector<Type> constRecBreaks;

  Rec_ConstantMean(vector<Type> logRv, vector<Type> crb) :
    logRvalue(logRv), constRecBreaks(crb) {};
  
  Type operator()(Type logssb, Type lastR, Type year){
    int usepar=0;
    for(int ii=0; ii<constRecBreaks.size(); ++ii){
      if(year>(Type)constRecBreaks(ii)){usepar++;}
    }
    return logRvalue(usepar);
  }
  
  Type logSe(Type logLambda){
    return  logLambda + logRvalue(logRvalue.size()-1);
  }
  Type dSR(Type logssb){
    return 0.0;
  }
  Type logSAtMaxR(){
    return R_NaReal;
  }
  Type logMaxR(){
    return R_NaReal;
  }
  Type maxGradient(){
    return R_NaReal;
  }
};



// Recruitment function 60
// Logistic Hockey Stick

//     // 0: alpha
//     // 1: mu
//     // 2: theta
//     predN = rec_pars(0) + rec_pars(1) + rec_pars(2) + log(1.0 + exp(-exp(-rec_pars(2)))) + log(exp(log(thisSSB)-rec_pars(1) - rec_pars(2)) - log(1.0 + exp((thisSSB-exp(rec_pars(1)))/exp(rec_pars(1) + rec_pars(2)))) + log(1.0 + exp(-exp(-rec_pars(2)))));

// mdsr = exp(par.rec_pars(0));


template<class Type>
struct RF_LogisticHockeyStick_t {
  Type loga;			// alpha
  Type logm;			// mu
  Type logt;			// theta

  RF_LogisticHockeyStick_t(Type la, Type lm, Type lt) :
    loga(la), logm(lm), logt(lt) {}

  template<class T>
  RF_LogisticHockeyStick_t(const RF_LogisticHockeyStick_t<T>& other) :
    loga(other.loga), logm(other.logm), logt(other.logt) {}
  
  template <template<class> class V, class T>
  T operator()(const V<T> &logssb0){
    T logssb = logssb0(0);
    T thisSSB = exp(logssb);
    T v= loga + logm + logt + log(1.0 + exp(-exp(-logt))) + log(exp(logssb - logm - logt) - log(1.0 + exp((thisSSB-exp(logm))/exp(logm + logt))) + log(1.0 + exp(-exp(-logt))));
    return v;
  } 
};
  
template<class Type>
struct Rec_LogisticHockeyStick : RecruitmentNumeric<Type, RF_LogisticHockeyStick_t<Type> >  {

  Type loga;
  Type logm;
  Type logt;
  
  Rec_LogisticHockeyStick(Type la, Type lm, Type lt) :
    RecruitmentNumeric<Type, RF_LogisticHockeyStick_t<Type> >(RF_LogisticHockeyStick_t<Type>(la, lm, lt), la), loga(la), logm(lm), logt(lt) {};

  template<class T>
  Rec_LogisticHockeyStick(const Rec_LogisticHockeyStick<T>& other) :
    RecruitmentNumeric<Type, RF_LogisticHockeyStick_t<Type> >(other.f, other.x0), loga(other.loga), logm(other.logm), logt(other.logt) {}

  
  Type logSAtMaxR(){
    return R_PosInf;
  }
   Type logMaxR(){
     return loga + logm + logt + log(1.0 + exp(-exp(-logt))) + log(exp(-logt) - log(1.0 + exp(-exp(-logt))));
  }
  Type maxGradient(){
    return exp(loga);
  }
};


// Recruitment function 61
// Hockey Stick

template<class Type>
struct Rec_HockeyStick : RecruitmentWorker<Type> {

  Type loglevel;
  Type logblim;

  Rec_HockeyStick(Type ll, Type lbl) :
    loglevel(ll), logblim(lbl) {};

  template<class T>
  Rec_HockeyStick(const Rec_HockeyStick<T>& other) :
    loglevel(other.loglevel), logblim(other.logblim) {}

  
  Type operator()(Type logssb, Type lastR, Type year){
    Type thisSSB = exp(logssb);
    return loglevel - logblim +
      log(thisSSB - (0.5 * ((thisSSB - exp(logblim))+Type(0.0)+TMBad::fabs((thisSSB - exp(logblim))-Type(0.0)))));
  }
  
  Type logSe(Type logLambda){
    return exp(logLambda + loglevel);
  }
  Type dSR(Type logssb){
    return TMBad::CondExpLt(logssb, logblim, exp(loglevel-logblim), Type(0.0));
  }
  Type logSAtMaxR(){
    return R_NaReal;
  }
  Type logMaxR(){
    return loglevel;
  }
  Type maxGradient(){
    return exp(loglevel-logblim);
  }
};


// Recruitment function 62
// AR(1) on log-recruitment

template<class Type>
struct Rec_LogAR1 : RecruitmentWorker<Type> {

  Type loglevel;
  Type logitPhi;

  Rec_LogAR1(Type ll, Type lp) : loglevel(ll), logitPhi(lp) {};
  
  Type operator()(Type logssb, Type lastR, Type year){
    return loglevel + (2.0 / (1.0 + exp(-logitPhi)) - 1.0) * (lastR - loglevel);
  }
  
  Type logSe(Type logLambda){
    return exp(logLambda + loglevel);
  }
  Type dSR(Type logssb){
    return 0.0;
  }
  Type logSAtMaxR(){
    return R_NegInf;
  }
  Type logMaxR(){
    return loglevel;
  }
  Type maxGradient(){
    return R_PosInf;
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


template<class Type>
struct RF_BentHyperbola_t {
  Type logBlim;
  Type logHalfSlope;
  Type logSmooth;

  RF_BentHyperbola_t(Type loga_, Type logb_, Type logg_) :
    logBlim(loga_), logHalfSlope(logb_), logSmooth(logg_) {}

  template<class T>
  RF_BentHyperbola_t(const RF_BentHyperbola_t<T>& other) :
    logBlim(other.logBlim), logHalfSlope(other.logHalfSlope), logSmooth(other.logSmooth) {}

  
  
  template <template<class> class V, class T>
  T operator()(const V<T> &logssb0){
    T logssb = logssb0(0);
    T thisSSB = exp(logssb);
    T v = logHalfSlope + log(thisSSB + sqrt(exp(2.0 * logBlim) + (exp(2.0 * logSmooth) / 4.0)) -
			     sqrt(pow(thisSSB-exp(logBlim),2) + (exp(2.0 * logSmooth) / 4.0)));
    return v;
  }
 
};
  
template<class Type>
struct Rec_BentHyperbola : RecruitmentNumeric<Type, RF_BentHyperbola_t<Type> >  {
  // Implement with known values when time permits!
  Rec_BentHyperbola(Type logBlim, Type logHalfSlope, Type logSmooth) :
    RecruitmentNumeric<Type, RF_BentHyperbola_t<Type> >(RF_BentHyperbola_t<Type>(logBlim, logHalfSlope, logSmooth), logBlim) {};

  template<class T>
  Rec_BentHyperbola(const Rec_BentHyperbola<T>& other) :
    RecruitmentNumeric<Type, RF_BentHyperbola_t<Type> >(other.f, other.x0) {}

};




// Recruitment function 64
// Power function with compensatory mortality property


//     predN = rec_pars(0) + invlogit(rec_pars(1)) * log(thisSSB);
//     break;

 //   Se = exp(1.0 / (1.0 - invlogit(newPar.rec_pars(1))) * (newPar.rec_pars(0) + logRecCorrection + log(lambda)));

//    mdsr = R_PosInf;


template<class Type>
struct RF_PowerCMP_t {
  Type loga;
  Type logb;

  RF_PowerCMP_t(Type loga_, Type logb_) :
    loga(loga_), logb(logb_) {}
  
  template<class T>
  RF_PowerCMP_t(const RF_PowerCMP_t<T>& other) :
    loga(other.loga), logb(other.logb) {}

  template <template<class> class V, class T>
  T operator()(const V<T> &logssb0){
    T logssb = logssb0(0);
    T v = loga + invlogit(logb) * logssb;
    return v;
  }
 
};
  
template<class Type>
struct Rec_PowerCMP : RecruitmentNumeric<Type, RF_PowerCMP_t<Type> >  {
  // Implement with known values when time permits!
  Rec_PowerCMP(Type loga, Type logb) :
    RecruitmentNumeric<Type, RF_PowerCMP_t<Type> >(RF_PowerCMP_t<Type>(loga, logb), loga) {};

  template<class T>
  Rec_PowerCMP(const Rec_PowerCMP<T>& other) :
    RecruitmentNumeric<Type, RF_PowerCMP_t<Type> >(other.f, other.x0) {}
};



// Recruitment function 65
// Power function without compensatory mortality property

//     predN = rec_pars(0) + (exp(rec_pars(1))+1.0001) * log(thisSSB);
//     break;

 //   Se = exp(1.0 / (1.0 - (exp(newPar.rec_pars(1)) + 1.0001)) * (newPar.rec_pars(0) + logRecCorrection + log(lambda)));

//   mdsr = R_PosInf;


template<class Type>
struct RF_PowerNCMP_t {
  Type loga;
  Type logb;

  RF_PowerNCMP_t(Type loga_, Type logb_) :
    loga(loga_), logb(logb_) {}
  
  template<class T>
  RF_PowerNCMP_t(const RF_PowerNCMP_t<T>& other) :
    loga(other.loga), logb(other.logb) {}

  template <template<class> class V, class T>
  T operator()(const V<T> &logssb0){
    T logssb = logssb0(0);
    T v = loga + (exp(logb)+1.0001) * logssb;
    return v;
  }
 
};
  
template<class Type>
struct Rec_PowerNCMP : RecruitmentNumeric<Type, RF_PowerNCMP_t<Type> >  {
  // Implement with known values when time permits!
  Rec_PowerNCMP(Type loga, Type logb) :
    RecruitmentNumeric<Type, RF_PowerNCMP_t<Type> >(RF_PowerNCMP_t<Type>(loga, logb), loga) {};

  template<class T>
  Rec_PowerNCMP(const Rec_PowerNCMP<T>& other) :
    RecruitmentNumeric<Type, RF_PowerNCMP_t<Type> >(other.f, other.x0) {}
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


template<class Type>
struct RF_Shepherd_t {
  Type loga;
  Type logb;
  Type logg;

  RF_Shepherd_t(Type loga_, Type logb_, Type logg_) :
    loga(loga_), logb(logb_), logg(logg_) {}

  template<class T>
  RF_Shepherd_t(const RF_Shepherd_t<T>& other) :
    loga(other.loga), logb(other.logb), logg(other.logg) {}

  
  template <template<class> class V, class T>
  T operator()(const V<T> &logssb0){
    T logssb = logssb0(0);
    T v = loga + logssb - logspace_add2(T(0.0), exp(logg) * (logssb - logb));
    return v;
  }
 
};
  
template<class Type>
struct Rec_Shepherd : RecruitmentNumeric<Type, RF_Shepherd_t<Type> >  {

  Rec_Shepherd(Type loga, Type logb, Type logg) :
    RecruitmentNumeric<Type, RF_Shepherd_t<Type> >(RF_Shepherd_t<Type>(loga, logb, logg), loga) {};

  template<class T>
  Rec_Shepherd(const Rec_Shepherd<T>& other) :
    RecruitmentNumeric<Type, RF_Shepherd_t<Type> >(other.f, other.x0) {}

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



template<class Type>
struct RF_HasselDeriso_t {
  Type loga;
  Type logb;
  Type logg;

  RF_HasselDeriso_t(Type loga_, Type logb_, Type logg_) :
    loga(loga_), logb(logb_), logg(logg_) {}

  template<class T>
  RF_HasselDeriso_t(const RF_HasselDeriso_t<T>& other) :
    loga(other.loga), logb(other.logb), logg(other.logg) {}

  template <template<class> class V, class T>
  T operator()(const V<T> &logssb0){
    T logssb = logssb0(0);
    T thisSSB = exp(logssb);
    T v = loga+logssb-exp(logg) * log(1.0+exp(logb + logg)*thisSSB);
    return v;
  }
 
};
  
template<class Type>
struct Rec_HasselDeriso : RecruitmentNumeric<Type, RF_HasselDeriso_t<Type> >  {

  Rec_HasselDeriso(Type loga, Type logb, Type logg) :
    RecruitmentNumeric<Type, RF_HasselDeriso_t<Type> >(RF_HasselDeriso_t<Type>(loga, logb, logg), loga) {};

  template<class T>
  Rec_HasselDeriso(const Rec_HasselDeriso<T>& other) :
    RecruitmentNumeric<Type, RF_HasselDeriso_t<Type> >(other.f, other.x0) {}

};




// Recruitment function 68
// Saila-Lorda

//     predN = rec_pars(0)+exp(rec_pars(2)) * log(thisSSB) - exp(rec_pars(1))*thisSSB;

 //   Se = Se_sl(lambda, exp(newPar.rec_pars(0) + logRecCorrection), exp(newPar.rec_pars(1)), exp(newPar.rec_pars(2)));



template<class Type>
struct RF_SailaLorda_t {
  Type loga;
  Type logb;
  Type logg;

  RF_SailaLorda_t(Type loga_, Type logb_, Type logg_) :
    loga(loga_), logb(logb_), logg(logg_) {}

  template<class T>
  RF_SailaLorda_t(const RF_SailaLorda_t<T>& other) :
    loga(other.loga), logb(other.logb), logg(other.logg) {}

  
  template <template<class> class V, class T>
  T operator()(const V<T> &logssb0){
    T logssb = logssb0(0);
    T thisSSB = exp(logssb);
    T v = loga+exp(logg) * logssb - exp(logb)*thisSSB;
    return v;
  }
 
};
  
template<class Type>
struct Rec_SailaLorda : RecruitmentNumeric<Type, RF_SailaLorda_t<Type> >  {

  Rec_SailaLorda(Type loga, Type logb, Type logg) :
    RecruitmentNumeric<Type, RF_SailaLorda_t<Type> >(RF_SailaLorda_t<Type>(loga, logb, logg), loga) {};

  template<class T>
  Rec_SailaLorda(const Rec_SailaLorda<T>& other) :
    RecruitmentNumeric<Type, RF_SailaLorda_t<Type> >(other.f, other.x0) {}

};



// Recruitment function 69
// Sigmoidal Beverton-Holt

//     predN = rec_pars(0)+exp(rec_pars(2)) * log(thisSSB)-log(1.0+exp(rec_pars(1))*exp(exp(rec_pars(2)) * log(thisSSB)));

 //   Se = Se_sbh(lambda, exp(newPar.rec_pars(0) + logRecCorrection), exp(newPar.rec_pars(1)), exp(newPar.rec_pars(2)));
    // if(g < 1.0){
    //   return (1.0 - g) / b * lambertW_raw( b / (1.0 - g) * pow(a * l, 1 / (1.0 - g) ));
    // }
 


template<class Type>
struct RF_SigmoidalBevHolt_t {
  Type loga;
  Type logb;
  Type logg;

  RF_SigmoidalBevHolt_t(Type loga_, Type logb_, Type logg_) :
    loga(loga_), logb(logb_), logg(logg_) {}

  template<class T>
  RF_SigmoidalBevHolt_t(const RF_SigmoidalBevHolt_t<T>& other) :
    loga(other.loga), logb(other.logb), logg(other.logg) {}

  
  template <template<class> class V, class T>
  T operator()(const V<T> &logssb0){
    T logssb = logssb0(0);
    T v = loga+exp(logg) * logssb-log(1.0+exp(logb)*exp(exp(logg) * logssb));
    return v;
  }
 
};
  
template<class Type>
struct Rec_SigmoidalBevHolt : RecruitmentNumeric<Type, RF_SigmoidalBevHolt_t<Type> >  {

  Rec_SigmoidalBevHolt(Type loga, Type logb, Type logg) :
    RecruitmentNumeric<Type, RF_SigmoidalBevHolt_t<Type> >(RF_SigmoidalBevHolt_t<Type>(loga, logb, logg), loga) {};

  template<class T>
  Rec_SigmoidalBevHolt(const Rec_SigmoidalBevHolt<T>& other) :
    RecruitmentNumeric<Type, RF_SigmoidalBevHolt_t<Type> >(other.f, other.x0) {}

};




// Recruitment function 90
// Spline with compensatory mortality property (non-increasing on log(R/SSB))

 //   predN(0) = log(thisSSB) + ibcdspline(log(thisSSB),
  // 					 (vector<Type>)(conf.constRecBreaks.template cast<Type>()),
  // 					 par.rec_pars);


template<class Type>
struct RF_SplineCMP_t {
  vector<Type> pars;
  vector<Type> knots;

  RF_SplineCMP_t(vector<Type> pars_, vector<Type> knots_) :
    pars(pars_), knots(knots_) {}

  template<class T>
  RF_SplineCMP_t(const RF_SplineCMP_t<T>& other) :
    pars(other.pars), knots(other.knots) {}
  
  template <template<class> class V, class T>
  T operator()(const V<T> &logssb0){
    T logssb = logssb0(0);
    vector<T> k2(knots);
    vector<T> p2(pars);
    T v = logssb + ibcdspline(logssb,k2,p2);
    return v;
  }
 
};
  
template<class Type>
struct Rec_SplineCMP : RecruitmentNumeric<Type, RF_SplineCMP_t<Type> >  {

  Rec_SplineCMP(vector<Type> pars, vector<Type> knots) :
    RecruitmentNumeric<Type, RF_SplineCMP_t<Type> >(RF_SplineCMP_t<Type>(pars,knots), pars(pars.size()-1)) {}

  template<class T>
  Rec_SplineCMP(const Rec_SplineCMP<T>& other) :
    RecruitmentNumeric<Type, RF_SplineCMP_t<Type> >(other.f, other.x0) {}
};


template<class Type>
struct RF_SplineConvexCompensatory_t {
  vector<Type> pars;
  vector<Type> knots;

  RF_SplineConvexCompensatory_t(vector<Type> pars_, vector<Type> knots_) :
    pars(pars_), knots(knots_) {}

  template<class T>
  RF_SplineConvexCompensatory_t(const RF_SplineCMP_t<T>& other) :
    pars(other.pars), knots(other.knots) {}
  
  template <template<class> class V, class T>
  T operator()(const V<T> &logssb0){
    T logssb = logssb0(0);
    vector<T> k2(knots);
    vector<T> p2(pars);
    T v = logssb + iibcispline(logssb,k2,p2);
    return v;
  }
 
};
  
template<class Type>
struct Rec_SplineConvexCompensatory : RecruitmentNumeric<Type, RF_SplineConvexCompensatory_t<Type> >  {

  Rec_SplineConvexCompensatory(vector<Type> pars, vector<Type> knots) :
    RecruitmentNumeric<Type, RF_SplineConvexCompensatory_t<Type> >(RF_SplineConvexCompensatory_t<Type>(pars,knots), pars(pars.size()-1)) {}

  template<class T>
  Rec_SplineConvexCompensatory(const Rec_SplineConvexCompensatory<T>& other) :
    RecruitmentNumeric<Type, RF_SplineConvexCompensatory_t<Type> >(other.f, other.x0) {}
};




// Recruitment function 91
// Smooth spline (integrated spline on log(R/SSB))

  //   predN(0) = log(thisSSB) + ibcspline(log(thisSSB),
  // 					(vector<Type>)(conf.constRecBreaks.template cast<Type>()),
  // 					par.rec_pars);


template<class Type>
struct RF_SplineSmooth_t {
  vector<Type> pars;
  vector<Type> knots;

  RF_SplineSmooth_t(vector<Type> pars_, vector<Type> knots_) :
    pars(pars_), knots(knots_) {}

  template<class T>
  RF_SplineSmooth_t(const RF_SplineSmooth_t<T>& other) :
    pars(other.pars), knots(other.knots) {}
  
  template <template<class> class V, class T>
  T operator()(const V<T> &logssb0){
    T logssb = logssb0(0);
    vector<T> k2(knots);
    vector<T> p2(pars.segment(0, pars.size()-1));
    T mu(pars(pars.size()-1));
    T v = logssb + ibcspline(logssb,k2,p2) + mu;
    return v;
  }
 
};
  
template<class Type>
struct Rec_SplineSmooth : RecruitmentNumeric<Type, RF_SplineSmooth_t<Type> >  {

  Rec_SplineSmooth(vector<Type> pars, vector<Type> knots) :
    RecruitmentNumeric<Type, RF_SplineSmooth_t<Type> >(RF_SplineSmooth_t<Type>(pars,knots), pars(pars.size()-1)) {}

  template<class T>
  Rec_SplineSmooth(const Rec_SplineSmooth<T>& other) :
    RecruitmentNumeric<Type, RF_SplineSmooth_t<Type> >(other.f, other.x0) {}
};


// Recruitment function 92
// General spline (on log(R/SSB))

 //   predN(0) = log(thisSSB) + bcspline(log(thisSSB),
  // 				       (vector<Type>)(conf.constRecBreaks.template cast<Type>()),
  // 				       par.rec_pars);



template<class Type>
struct RF_SplineGeneral_t {
  vector<Type> pars;
  vector<Type> knots;

  RF_SplineGeneral_t(vector<Type> pars_, vector<Type> knots_) :
    pars(pars_), knots(knots_) {}

  template<class T>
  RF_SplineGeneral_t(const RF_SplineGeneral_t<T>& other) :
    pars(other.pars), knots(other.knots) {}
  
  template <template<class> class V, class T>
  T operator()(const V<T> &logssb0){
    T logssb = logssb0(0);
    vector<T> k2(knots);
    vector<T> p2(pars.segment(0, pars.size()-1));
    T mu(pars(pars.size()-1));
    T v = logssb + bcspline(logssb,k2,p2) + mu;
    return v;
  }
 
};
  
template<class Type>
struct Rec_SplineGeneral : RecruitmentNumeric<Type, RF_SplineGeneral_t<Type> >  {

  Rec_SplineGeneral(vector<Type> pars, vector<Type> knots) :
    RecruitmentNumeric<Type, RF_SplineGeneral_t<Type> >(RF_SplineGeneral_t<Type>(pars,knots), pars(pars.size()-1)) {}

  template<class T>
  Rec_SplineGeneral(const Rec_SplineGeneral<T>& other) :
    RecruitmentNumeric<Type, RF_SplineGeneral_t<Type> >(other.f, other.x0) {}
};


///////////////////////////////////////////////////////////////////////////////////////////////////

template<class Type>
Recruitment<Type> makeRecruitmentFunction(const confSet& conf, const paraSet<Type>& par){
  RecruitmentModel rm = static_cast<RecruitmentModel>(conf.stockRecruitmentModelCode);
  Recruitment<Type> r;

  if(rm == RecruitmentModel::NoRecruit){
    r = Recruitment<Type>(new Rec_None<Type>());
    
////////////////////////////////////////// The Beginning //////////////////////////////////////////
    
  }else if(rm == RecruitmentModel::LogRandomWalk){
    if(par.rec_pars.size() != 0)
      Rf_error("The random walk recruitment should not have any parameters.");
    r = Recruitment<Type>(new Rec_LogRW<Type>());
  }else if(rm == RecruitmentModel::Ricker){
    if(par.rec_pars.size() != 2)
      Rf_error("The Ricker recruitment must have two parameters.");
    r = Recruitment<Type>(new Rec_Ricker<Type>(par.rec_pars(0), par.rec_pars(1)));
  }else if(rm == RecruitmentModel::BevertonHolt){
    if(par.rec_pars.size() != 2)
      Rf_error("The Beverton Holt recruitment must have two parameters.");
    r = Recruitment<Type>(new Rec_BevertonHolt<Type>(par.rec_pars(0), par.rec_pars(1)));
  }else if(rm == RecruitmentModel::ConstantMean){
    if(par.rec_pars.size() != conf.constRecBreaks.size() + 1)
      Rf_error("The constant mean recruitment should have one more parameter than constRecBreaks.");
    r = Recruitment<Type>(new Rec_ConstantMean<Type>(par.rec_pars, conf.constRecBreaks));
  }else if(rm == RecruitmentModel::LogisticHockeyStick){
    if(par.rec_pars.size() != 3)
      Rf_error("The logistic hockey stick recruitment should have three parameters.");
    r = Recruitment<Type>(new Rec_LogisticHockeyStick<Type>(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)));
  }else if(rm == RecruitmentModel::HockeyStick){
   if(par.rec_pars.size() != 2)
      Rf_error("The hockey stick recruitment should have two parameters.");
     r = Recruitment<Type>(new Rec_HockeyStick<Type>(par.rec_pars(0), par.rec_pars(1)));
  }else if(rm == RecruitmentModel::LogAR1){
   if(par.rec_pars.size() != 2)
      Rf_error("The log-AR(1) recruitment should have two parameters.");
    r = Recruitment<Type>(new Rec_LogAR1<Type>(par.rec_pars(0), par.rec_pars(1)));
  }else if(rm == RecruitmentModel::BentHyperbola){
   if(par.rec_pars.size() != 3)
      Rf_error("The bent hyperbola recruitment should have three parameters.");
    r = Recruitment<Type>(new Rec_BentHyperbola<Type>(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)));
  }else if(rm == RecruitmentModel::Power_CMP){
   if(par.rec_pars.size() != 2)
      Rf_error("The power law recruitment should have two parameters.");
    r = Recruitment<Type>(new Rec_PowerCMP<Type>(par.rec_pars(0), par.rec_pars(1)));
  }else if(rm == RecruitmentModel::Power_NCMP){
   if(par.rec_pars.size() != 2)
      Rf_error("The power law recruitment should have two parameters.");
    r = Recruitment<Type>(new Rec_PowerNCMP<Type>(par.rec_pars(0), par.rec_pars(1)));

//////////////////////////////////////// 3 parameter models ///////////////////////////////////////
    
  }else if(rm == RecruitmentModel::Shepherd){
   if(par.rec_pars.size() != 3)
      Rf_error("The Shepherd recruitment should have three parameters.");
    r = Recruitment<Type>(new Rec_Shepherd<Type>(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)));
  }else if(rm == RecruitmentModel::Hassel_Deriso){
   if(par.rec_pars.size() != 3)
      Rf_error("The Hassel/Deriso recruitment should have three parameters.");
    r = Recruitment<Type>(new Rec_HasselDeriso<Type>(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)));
  }else if(rm == RecruitmentModel::SailaLorda){
   if(par.rec_pars.size() != 3)
      Rf_error("The Saila-Lorda recruitment should have three parameters.");
    r = Recruitment<Type>(new Rec_SailaLorda<Type>(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)));
  }else if(rm == RecruitmentModel::SigmoidalBevertonHolt){
   if(par.rec_pars.size() != 3)
     Rf_error("The sigmoidal Beverton-Holt recruitment should have three parameters.");
    r = Recruitment<Type>(new Rec_SigmoidalBevHolt<Type>(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)));

//////////////////////////////////////// Spline recruitment ///////////////////////////////////////
    
  }else if(rm == RecruitmentModel::Spline_CMP){
    if(par.rec_pars.size() != conf.constRecBreaks.size() + 1)
      Rf_error("The spline recruitment should have one more parameter than constRecBreaks.");
    r = Recruitment<Type>(new Rec_SplineCMP<Type>(par.rec_pars, conf.constRecBreaks.template cast<Type>()));
  }else if(rm == RecruitmentModel::Spline_Smooth){
    if(par.rec_pars.size() != conf.constRecBreaks.size() + 1)
      Rf_error("The spline recruitment should have one more parameter than constRecBreaks.");
    r = Recruitment<Type>(new Rec_SplineSmooth<Type>(par.rec_pars, conf.constRecBreaks.template cast<Type>()));
  }else if(rm == RecruitmentModel::Spline_General){
    if(par.rec_pars.size() != conf.constRecBreaks.size() + 1)
      Rf_error("The spline recruitment should have one more parameter than constRecBreaks.");
    r = Recruitment<Type>(new Rec_SplineGeneral<Type>(par.rec_pars, conf.constRecBreaks.template cast<Type>()));

  }else if(rm == RecruitmentModel::Spline_ConvexCompensatory){
    if(par.rec_pars.size() != conf.constRecBreaks.size() + 1)
      Rf_error("The spline recruitment should have one more parameter than constRecBreaks.");
    r = Recruitment<Type>(new Rec_SplineConvexCompensatory<Type>(par.rec_pars, conf.constRecBreaks.template cast<Type>()));


//////////////////////// S/(d+S) type Depensatory recruitment models //////////////////////////////

  }else if(rm == RecruitmentModel::Depensatory_B_Ricker){
   if(par.rec_pars.size() != 3)
     Rf_error("The depensatory B Ricker recruitment should have three parameters.");
    r = Recruitment<Type>(Rec_DepensatoryB(Rec_Ricker<Type>(par.rec_pars(0), par.rec_pars(1)),par.rec_pars(2)));
  }else if(rm == RecruitmentModel::Depensatory_B_BevertonHolt){
    if(par.rec_pars.size() != 3)
     Rf_error("The depensatory B Beverton-Holt recruitment should have three parameters.");
    r = Recruitment<Type>(Rec_DepensatoryB(Rec_BevertonHolt<Type>(par.rec_pars(0), par.rec_pars(1)),par.rec_pars(2)));
  }else if(rm == RecruitmentModel::Depensatory_B_LogisticHockeyStick){
     if(par.rec_pars.size() != 4)
     Rf_error("The depensatory B logistic hockey stick recruitment should have four parameters.");
  r = Recruitment<Type>(Rec_DepensatoryB(Rec_LogisticHockeyStick<Type>(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)),par.rec_pars(3)));
  }else if(rm == RecruitmentModel::Depensatory_B_HockeyStick){
    if(par.rec_pars.size() != 3)
     Rf_error("The depensatory B hockey stick recruitment should have three parameters.");
r = Recruitment<Type>(Rec_DepensatoryB(Rec_HockeyStick<Type>(par.rec_pars(0), par.rec_pars(1)),par.rec_pars(2)));
  }else if(rm == RecruitmentModel::Depensatory_B_BentHyperbola){
    if(par.rec_pars.size() != 4)
     Rf_error("The depensatory B bent hyperbola recruitment should have four parameters.");
    r = Recruitment<Type>(Rec_DepensatoryB(Rec_BentHyperbola<Type>(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)),par.rec_pars(3)));
  }else if(rm == RecruitmentModel::Depensatory_B_Power){
    if(par.rec_pars.size() != 3)
     Rf_error("The depensatory B power law recruitment should have three parameters.");
   r = Recruitment<Type>(Rec_DepensatoryB(Rec_PowerCMP<Type>(par.rec_pars(0), par.rec_pars(1)),par.rec_pars(2)));
  }else if(rm == RecruitmentModel::Depensatory_B_Shepherd){
    if(par.rec_pars.size() != 4)
      Rf_error("The depensatory B Shepherd recruitment should have four parameters.");
    r = Recruitment<Type>(Rec_DepensatoryB(Rec_Shepherd<Type>(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)),par.rec_pars(3)));
  }else if(rm == RecruitmentModel::Depensatory_B_Hassel_Deriso){
       if(par.rec_pars.size() != 4)
     Rf_error("The depensatory B Hassel/Deriso recruitment should have four parameters.");
r = Recruitment<Type>(Rec_DepensatoryB(Rec_HasselDeriso<Type>(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)),par.rec_pars(3)));
  }else if(rm == RecruitmentModel::Depensatory_B_Spline_CMP){
    if(par.rec_pars.size() != conf.constRecBreaks.size() + 2)
      Rf_error("The depensatory B spline recruitment should have two parameters more than constRecBreaks.");
    r = Recruitment<Type>(Rec_DepensatoryB(Rec_SplineCMP<Type>(par.rec_pars.segment(0,par.rec_pars.size()-1), conf.constRecBreaks),par.rec_pars(par.rec_pars.size()-1)));

//////////////////////// 1/(1+exp(-e * (S-d))) type Depensatory recruitment models //////////////////////////////

  }else if(rm == RecruitmentModel::Depensatory_C_Ricker){
    if(par.rec_pars.size() != 4)
      Rf_error("The depensatory C Ricker recruitment should have four parameters.");
    r = Recruitment<Type>(Rec_DepensatoryC(Rec_Ricker<Type>(par.rec_pars(0), par.rec_pars(1)),par.rec_pars(2), par.rec_pars(3)));
  }else if(rm == RecruitmentModel::Depensatory_C_BevertonHolt){
    if(par.rec_pars.size() != 4)
      Rf_error("The depensatory C Beverton-Holt recruitment should have four parameters.");
    r = Recruitment<Type>(Rec_DepensatoryC(Rec_BevertonHolt<Type>(par.rec_pars(0), par.rec_pars(1)),par.rec_pars(2), par.rec_pars(3)));


///////////////////////////////////////////// The End /////////////////////////////////////////////
    
  }else{
    Rf_error("Stock-recruitment model code not implemented.");
  }

  return r;
  
}

  
  
#endif
