#ifndef SAM_RECRUITMENT_HPP
#define SAM_RECRUITMENT_HPP

#include <memory>

enum RecruitmentModel {
		       ICESforecast = -2, // Known, Implemented
		       NoRecruit = -1, // Known, Implemented
		       LogRandomWalk = 0, // Known, Implemented
		       Ricker = 1,	  // Known, Implemented
		       BevertonHolt = 2,  // Known, Implemented
		       ConstantMean = 3,  // Known, Implemented
		       T2D_LogisticHockeyStick = 50, // Numeric, Implemented
		       T2D_Ricker = 51,		     // Numeric, Implemented
		       T2D_BevertonHolt = 52,	     // Numeric, Implemented
		       T2D_HockeyStick = 53,	     // Numeric, Implemented
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
		       T2D_Spline_CMP = 93	     // Numeric
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
  // template <template<class> class V, class T>
  // T operator()(V<T> &x){
  template<class T>
  T operator()(const vector<T> &x) {
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
  // template <template<class> class V, class T>
  // T operator()(const V<T> &x){
  template<class T>
  T operator()(const vector<T> &x) {
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
  // template <template<class> class V, class T>
  // T operator()(V<T> &x){
  template<class T>
  T operator()(const vector<T> &x) {
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
  // template <template<class> class V, class T>
  // T operator()(V<T> &x){
  template<class T>
  T operator()(const vector<T> &x) {
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
  // template <template<class> class V, class T>
  // T operator()(V<T> &x){
  template<class T>
  T operator()(const vector<T> &x) {
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

  
// Recruitment function 50
// Type 2 depensatory logistic hockey stick
template<class Type>
struct RF_T2D_LHS_t {
  Type loga;
  Type logb;
  Type logc;
  Type logd;

  RF_T2D_LHS_t(Type la, Type lb, Type lc, Type ld) :
    loga(la), logb(lb), logc(lc), logd(ld) {};
  
  template <template<class> class V, class T>
  T operator()(V<T> &logssb0){
    T logssb = logssb0(0);
    T v = (T)loga + (T)logb + (T)logc + log(1.0 + exp(-exp(-(T)logc))) + log(exp(logssb-(T)logb - (T)logc) - log(1.0 + exp((exp(logssb)-exp((T)logb))/exp((T)logb + (T)logc))) + log(1.0 + exp(-exp(-logb)))) +
      logssb - logspace_add2(logssb,(T)logd);
    return v;
  }
 
};
  
template<class Type>
struct Rec_T2D_logisticHockeyStick : RecruitmentNumeric<Type, RF_T2D_LHS_t<Type> >  {

  Rec_T2D_logisticHockeyStick(Type loga, Type logb, Type logc, Type logd) :
    RecruitmentNumeric<Type, RF_T2D_LHS_t<Type> >(RF_T2D_LHS_t<Type>(loga, logb, logc, logd), logd + 2.0) {};
};

// Recruitment function 51
// Type 2 depensatory Ricker

template<class Type>
struct RF_T2D_Ricker_t {
  Type loga;
  Type logb;
  Type logd;

  RF_T2D_Ricker_t(Type la, Type lb, Type ld) :
    loga(la), logb(lb), logd(ld) {};
  
  template <template<class> class V, class T>
  T operator()(V<T> &logssb0){
    T logssb = logssb0(0);
    T v = (T)loga + logssb - exp((T)logb + logssb) +
      logssb - logspace_add2(logssb,(T)logd);
    return v;
  }
 
};
  
template<class Type>
struct Rec_T2D_Ricker : RecruitmentNumeric<Type, RF_T2D_LHS_t<Type> >  {

  Rec_T2D_Ricker(Type loga, Type logb, Type logc, Type logd) :
    RecruitmentNumeric<Type, RF_T2D_Ricker_t<Type> >(RF_T2D_Ricker_t<Type>(loga, logb, logc, logd), logd + 2.0) {};
};


// Recruitment function 52
// Type 2 depensatory Beverton-Holt


template<class Type>
struct RF_T2D_BevHolt_t {
  Type loga;
  Type logb;
  Type logd;

  RF_T2D_BevHolt_t(Type la, Type lb, Type ld) :
    loga(la), logb(lb), logd(ld) {};
  
  template <template<class> class V, class T>
  T operator()(V<T> &logssb0){
    T logssb = logssb0(0);
    T v = (T)loga + logssb - logspace_add2(T(0.0),(T)logb + logssb) +
      logssb - logspace_add2(logssb,(T)logd);
    return v;
  }
 
};
  
template<class Type>
struct Rec_T2D_BevHolt : RecruitmentNumeric<Type, RF_T2D_BevHolt_t<Type> >  {

  Rec_T2D_BevHolt(Type loga, Type logb, Type logc, Type logd) :
    RecruitmentNumeric<Type, RF_T2D_BevHolt_t<Type> >(RF_T2D_BevHolt_t<Type>(loga, logb, logc, logd), logd + 2.0) {};
};


// Recruitment function 53
// Type 2 depensatory Hockey Stick

template<class Type>
struct RF_T2D_HockeyStick_t {
  Type loglevel;
  Type logblim;
  Type logd;

  RF_T2D_HockeyStick_t(Type ll, Type lbl, Type ld) :
    loglevel(ll), logblim(lbl), logd(ld) {};
  
  template <template<class> class V, class T>
  T operator()(V<T> &logssb0){
    T logssb = logssb0(0);
    T thisSSB = exp(logssb);
    T v= loglevel - logblim +
      log(thisSSB - (0.5 * ((thisSSB - exp(logblim))+Type(0.0)+CppAD::abs((thisSSB - exp(logblim))-Type(0.0))))) +
      logssb - logspace_add2(logssb,(T)logd);
    return v;
  }
 
};
  
template<class Type>
struct Rec_T2D_HockeyStick : RecruitmentNumeric<Type, RF_T2D_HockeyStick_t<Type> >  {

  Rec_T2D_HockeyStick(Type loglevel, Type logblim, Type logd) :
    RecruitmentNumeric<Type, RF_T2D_HockeyStick_t<Type> >(RF_T2D_HockeyStick_t<Type>(loglevel, logblim, logd), logd + 2.0) {};
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
    loga(la), logm(lm), logt(lt) {};
  
  template <template<class> class V, class T>
  T operator()(V<T> &logssb0){
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
    RecruitmentNumeric<Type, RF_T2D_HockeyStick_t<Type> >(RF_LogisticHockeyStick_t<Type>(la, lm, lt), la), loga(la), logm(lm), logt(lt) {};

  
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
struct Rec_logAR1 : RecruitmentWorker<Type> {

  Type loglevel;
  Type logitPhi;

  Rec_logAR1(Type ll, Type lp) : loglevel(ll), logitPhi(lp) {};
  
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

// Recruitment function 64
// Power function with compensatory mortality property


//     predN = rec_pars(0) + invlogit(rec_pars(1)) * log(thisSSB);
//     break;

 //   Se = exp(1.0 / (1.0 - invlogit(newPar.rec_pars(1))) * (newPar.rec_pars(0) + logRecCorrection + log(lambda)));

//    mdsr = R_PosInf;

// Recruitment function 65
// Power function without compensatory mortality property

//     predN = rec_pars(0) + (exp(rec_pars(1))+1.0001) * log(thisSSB);
//     break;

 //   Se = exp(1.0 / (1.0 - (exp(newPar.rec_pars(1)) + 1.0001)) * (newPar.rec_pars(0) + logRecCorrection + log(lambda)));


//   mdsr = R_PosInf;

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
    loga(loga_), logb(logb_), logg(logg_) {};
  
  template <template<class> class V, class T>
  T operator()(V<T> &logssb0){
    T logssb = logssb0(0);
    T thisSSB = exp(logssb);
    T v = loga + logssb - logspace_add2(Type(0.0), exp(logg) * (logssb - logb));
    return v;
  }
 
};
  
template<class Type>
struct Rec_Shepherd : RecruitmentNumeric<Type, RF_Shepherd_t<Type> >  {

  Rec_Shepherd(Type loga, Type logb, Type logg) :
    RecruitmentNumeric<Type, RF_Shepherd_t<Type> >(RF_Shepherd_t<Type>(loga, logb, logg), loga) {};
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

// Recruitment function 68
// Saila-Lorda

//     predN = rec_pars(0)+exp(rec_pars(2)) * log(thisSSB) - exp(rec_pars(1))*thisSSB;

 //   Se = Se_sl(lambda, exp(newPar.rec_pars(0) + logRecCorrection), exp(newPar.rec_pars(1)), exp(newPar.rec_pars(2)));


// Recruitment function 69
// Sigmoidal Beverton-Holt

//     predN = rec_pars(0)+exp(rec_pars(2)) * log(thisSSB)-log(1.0+exp(rec_pars(1))*exp(exp(rec_pars(2)) * log(thisSSB)));

 //   Se = Se_sbh(lambda, exp(newPar.rec_pars(0) + logRecCorrection), exp(newPar.rec_pars(1)), exp(newPar.rec_pars(2)));
    // if(g < 1.0){
    //   return (1.0 - g) / b * lambertW_raw( b / (1.0 - g) * pow(a * l, 1 / (1.0 - g) ));
    // }
 


// Recruitment function 90
// Spline with compensatory mortality property (non-increasing on log(R/SSB))

 //   predN(0) = log(thisSSB) + ibcdspline(log(thisSSB),
  // 					 (vector<Type>)(conf.constRecBreaks.template cast<Type>()),
  // 					 par.rec_pars);
  

// Recruitment function 91
// Smooth spline (integrated spline on log(R/SSB))

  //   predN(0) = log(thisSSB) + ibcspline(log(thisSSB),
  // 					(vector<Type>)(conf.constRecBreaks.template cast<Type>()),
  // 					par.rec_pars);


// Recruitment function 92
// General spline (on log(R/SSB))

 //   predN(0) = log(thisSSB) + bcspline(log(thisSSB),
  // 				       (vector<Type>)(conf.constRecBreaks.template cast<Type>()),
  // 				       par.rec_pars);


// Recruitment function 
// Type 2 depensatory CMP spline (recruitment function 90)

 //    predN(0) = log(thisSSB) + ibcdspline(log(thisSSB),
  // 					 (vector<Type>)(conf.constRecBreaks.template cast<Type>()),
  // 					  (vector<Type>)par.rec_pars.segment(0,par.rec_pars.size()-1)) +
  //      log(thisSSB) - logspace_add2(log(thisSSB),par.rec_pars(par.rec_pars.size()-1));



///////////////////////////////////////////////////////////////////////////////////////////////////

template<class Type>
Recruitment<Type> makeRecruitmentFunction(dataSet<Type>& dat, confSet& conf, paraSet<Type>& par){
  RecruitmentModel rm = static_cast<RecruitmentModel>(conf.stockRecruitmentModelCode);
  Recruitment<Type> r;

  if(rm == RecruitmentModel::NoRecruit){
    r = Recruitment<Type>(new Rec_None<Type>());
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
      Rf_error("The Ricker recruitment must have two parameters.");
    r = Recruitment<Type>(new Rec_BevertonHolt<Type>(par.rec_pars(0), par.rec_pars(1)));
  }else if(rm == RecruitmentModel::ConstantMean){

  }else if(rm == RecruitmentModel::T2D_LogisticHockeyStick){

  }else if(rm == RecruitmentModel::T2D_Ricker){

  }else if(rm == RecruitmentModel::T2D_BevertonHolt){

  }else if(rm == RecruitmentModel::T2D_HockeyStick){

  }else if(rm == RecruitmentModel::LogisticHockeyStick){

  }else if(rm == RecruitmentModel::HockeyStick){

  }else if(rm == RecruitmentModel::LogAR1){

  }else if(rm == RecruitmentModel::BentHyperbola){

  }else if(rm == RecruitmentModel::Power_CMP){

  }else if(rm == RecruitmentModel::Power_NCMP){

  }else if(rm == RecruitmentModel::Shepherd){

  }else if(rm == RecruitmentModel::Hassel_Deriso){

  }else if(rm == RecruitmentModel::SailaLorda){

  }else if(rm == RecruitmentModel::SigmoidalBevertonHolt){

  }else if(rm == RecruitmentModel::Spline_CMP){

  }else if(rm == RecruitmentModel::Spline_Smooth){

  }else if(rm == RecruitmentModel::Spline_General){

  }else if(rm == RecruitmentModel::T2D_Spline_CMP){

  }else{
	Rf_error("Stock-recruitment model code not implemented.");
  }

  return r;
  
}

  
  
#endif
