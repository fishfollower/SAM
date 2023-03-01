SAM_DEPENDS(spline)
SAM_DEPENDS(define)
SAM_DEPENDS(logspace)
SAM_DEPENDS(newton)


#include <memory>


HEADER(
       template<class Type>
       struct RecruitmentWorker {

	 int isAutoregressive;
	 int isTimevarying;

	 RecruitmentWorker(): isAutoregressive(0), isTimevarying(0) {};
	 RecruitmentWorker(int isA, int isT): isAutoregressive(isA), isTimevarying(isT) {};
  
	 // Stock recruitment function logR = f(logssb)
	 virtual Type operator()(Type logssb, Type lastLogR, Type year) = 0;

	 Type R(Type logssb, Type lastLogR, Type year){
	   return exp(operator()(logssb, lastLogR, year));
	 };

	 // Deterministic equilibrium biomass for lambda = 1/SPR
	 virtual Type logSe(Type logLambda) = 0;
	 // Derivative
	 virtual Type dSR(Type logssb) = 0;

	 virtual Type logSAtMaxR() = 0;
	 virtual Type logMaxR() = 0;
	 virtual Type maxGradient() = 0;

	 virtual ~RecruitmentWorker() {};

       });


SAM_SPECIALIZATION(struct RecruitmentWorker<double>);
SAM_SPECIALIZATION(struct RecruitmentWorker<TMBad::ad_aug>);


HEADER(
       template<class Type>
       struct Recruitment {
       public:
	 const char* name;
	 std::shared_ptr<RecruitmentWorker<Type> > ptr;
	 explicit Recruitment() : name("Uninitialized"), ptr(nullptr) {};
	 explicit Recruitment(const char* nm, std::shared_ptr<RecruitmentWorker<Type> > p=nullptr) : name(nm), ptr(p) {};
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

       });

SAM_SPECIALIZATION(struct Recruitment<double>);
SAM_SPECIALIZATION(struct Recruitment<TMBad::ad_aug>);



#ifndef WITH_SAM_LIB

namespace RecruitmentConvenience {

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

  typedef TMBad::ad_aug ad;
  
  struct RecruitmentFunctor { // : NewtonFunctor {
    virtual ad operator()(const ad& x){
      return R_NaReal;
    }

    ad safe_eval(const ad& x){
      return operator()(x);
    }
    double safe_eval(const double& x){
      TMBad::global dummy;
      dummy.ad_start();
      double ans = asDouble(operator()(x));
      dummy.ad_stop();
      return ans;
    }
  };

#define USING_RECFUN				\
  using RecruitmentFunctor::safe_eval

  struct WrapSR : NewtonFunctor {
    std::shared_ptr<RecruitmentFunctor> ptr;
    WrapSR(const std::shared_ptr<RecruitmentFunctor>& p_) : ptr(p_) {};
    
    ad operator()(const vector<ad> &x) {
      return -ptr->operator()(x(0));
    }
  };

  struct WrapSRp : NewtonFunctor {
    std::shared_ptr<RecruitmentFunctor> ptr;
    WrapSRp(const std::shared_ptr<RecruitmentFunctor>& p_) : ptr(p_) {};
    
    ad operator()(const vector<ad> &x) {
      return ptr->operator()(x(0));
    }
  };
  
  struct Exp : RecruitmentFunctor {
    std::shared_ptr<RecruitmentFunctor> ptr;
    Exp(const std::shared_ptr<RecruitmentFunctor>& p_) : ptr(p_) {};  
    ad operator()(const ad &logssb) {
      return exp(ptr->operator()(logssb));
    }
    USING_RECFUN;
  };
  
  struct diffSR : RecruitmentFunctor {
    std::shared_ptr<RecruitmentFunctor> ptr;
    diffSR(const std::shared_ptr<RecruitmentFunctor>& p_) : ptr(p_) {};  
    ad operator()(const ad &logssb) {
      ad h = 0.0001 * sqrt(logssb * logssb + (ad)SAM_Zero);
      ad x1 = logssb + 2.0 * h;
      ad v1 = ptr->operator()(x1);
      ad x2 = logssb + h;
      ad v2 = ptr->operator()(x2);
      ad x3 = logssb - h;
      ad v3 = ptr->operator()(x3);
      ad x4 = logssb - 2.0 * h;
      ad v4 = ptr->operator()(x4);
      ad g = (-v1 + 8.0 * v2 - 8.0 * v3 + v4) / (12.0 * h);
      return (ad)g / exp(logssb);
    };

    USING_RECFUN;
  };


  struct WrapDiffSR : NewtonFunctor {
    std::shared_ptr<diffSR> ptr;

    WrapDiffSR(const std::shared_ptr<RecruitmentFunctor>& p_) : ptr(new diffSR(p_)) {};  
   
    ad operator()(const vector<ad> &x) {
      return -ptr->operator()(x(0));
    }
  };
  
  struct WrapEquiS : NewtonFunctor {
    std::shared_ptr<RecruitmentFunctor> ptr;
    ad logLambda;

    WrapEquiS(const std::shared_ptr<RecruitmentFunctor>& p_,
	      const ad& l) : ptr(p_), logLambda(l) {};
 
    ad operator()(const vector<ad> &x) {
      ad xv = x(0);
      ad loga = logLambda + ptr->operator()(xv);
      ad logb = xv;
      ad v = (loga - logb) * (loga - logb);
      return v;
    }  
  };




  // RecruitmentNumeric is not allowed to be autoregressive or timevarying
  template<class Type>
  struct RecruitmentNumeric : RecruitmentWorker<Type> {

    std::shared_ptr<RecruitmentFunctor> ptr;
  
    RecruitmentNumeric(const std::shared_ptr<RecruitmentFunctor>& p_) : RecruitmentWorker<Type>(0,0), ptr(p_) {};  
   
    Type operator()(Type logssb, Type lastLogR, Type year){
      return ptr->safe_eval(logssb);
    }
  
    virtual Type logSe(Type logLambda){
      std::shared_ptr<WrapEquiS> pWrap = std::make_shared<WrapEquiS>(ptr, logLambda);
      std::shared_ptr<NewtonFunctor> pUse = std::static_pointer_cast<NewtonFunctor>(pWrap);
      vector<Type> x0v(1); x0v(0) = 20.0;
      newton::newton_config cfg;
      cfg.on_failure_return_nan = false;
      vector<Type> v = SAM_Newton(pUse,x0v,cfg);
      return TMBad::CondExpLt((Type)v(0),(Type)SAM_NegInf,(Type)SAM_NegInf,(Type)v(0));
    }

    virtual Type dSR(Type logssb){
      diffSR dF(ptr);
      return dF.safe_eval(logssb);
    }
  
    virtual Type logSAtMaxR(){
      std::shared_ptr<WrapSR> pWrap = std::make_shared<WrapSR>(ptr);
      std::shared_ptr<NewtonFunctor> pUse = std::static_pointer_cast<NewtonFunctor>(pWrap);
      vector<Type> x0v(1); x0v(0) = 20.0;
      newton::newton_config cfg;
      cfg.on_failure_return_nan = false;
      vector<Type> v = SAM_Newton(pUse,x0v,cfg);  
      return v(0);
    }
  
    virtual Type logMaxR(){   
      Type v = logSAtMaxR();    
      return ptr->safe_eval(v);
    }
    virtual Type maxGradient(){
      std::shared_ptr<WrapDiffSR> pWrap = std::make_shared<WrapDiffSR>(ptr);
      std::shared_ptr<NewtonFunctor> pUse = std::static_pointer_cast<NewtonFunctor>(pWrap);
      vector<Type> x0v(1); x0v(0) = 0.0;
      newton::newton_config cfg;
      // cfg.simplify = false;
      cfg.on_failure_return_nan = false;
      vector<Type> v = SAM_Newton(pUse,x0v, cfg);
      return pWrap->ptr->safe_eval(v(0));
    }  
  };



  //////////////////////////////////////////////////////////////////////////////////////////
  // Depensatory wrappers
  //////////////////////////////////////////////////////////////////////////////////////////

  /*
    Form A: R_D(S; d) = R(S^d)
  */
  struct WrapDepensatoryA : RecruitmentFunctor {
    std::shared_ptr<RecruitmentWorker<ad> > ptr;
    ad logd;

    WrapDepensatoryA(std::shared_ptr<RecruitmentWorker<ad> >& p_,
		     ad& logd_) : ptr(p_), logd(logd_) {};

    ad operator()(const ad& x){
      ad logssb = exp(logd) * x;
      ad v = ptr->operator()(logssb, (ad)R_NaReal, (ad)R_NaReal);
      return v;
    }

    USING_RECFUN;

  };


  /*
    Form B: R_D(S; d) = R(S) * S / (S + d)
  */
  template <class Type>
  std::shared_ptr<RecruitmentNumeric<Type> > Rec_DepensatoryA(std::shared_ptr<RecruitmentWorker<ad> > SR,
							      ad logd){
    return std::make_shared<RecruitmentNumeric<Type> >(std::make_shared<WrapDepensatoryA>(SR, logd));
  }


  struct WrapDepensatoryB : RecruitmentFunctor {
    std::shared_ptr<RecruitmentWorker<ad> > ptr;
    ad logd;

    WrapDepensatoryB(std::shared_ptr<RecruitmentWorker<ad> >& p_,
		     ad& logd_) : ptr(p_), logd(logd_) {};

    ad operator()(const ad& x){
      ad logssb = x;
      ad v = ptr->operator()(logssb, (ad)R_NaReal, (ad)R_NaReal) + logssb - logspace_add_SAM(logssb, logd);
      return v;
    }

    USING_RECFUN;

  };

  template <class Type>
  std::shared_ptr<RecruitmentNumeric<Type> > Rec_DepensatoryB(std::shared_ptr<RecruitmentWorker<ad> > SR,
							      ad logd){
    return std::make_shared<RecruitmentNumeric<Type> >(std::make_shared<WrapDepensatoryB>(SR, logd));
  }




  /*
    Form C: R_D(S; d) = R(S) * 1 / (1 + exp(-l * (log(S) - log(d))))
  */
  struct WrapDepensatoryC : RecruitmentFunctor {
    std::shared_ptr<RecruitmentWorker<ad> > ptr;
    ad logd;
    ad logl;

    WrapDepensatoryC(std::shared_ptr<RecruitmentWorker<ad> >& p_,
		     ad& logd_,
		     ad& logl_) : ptr(p_), logd(logd_), logl(logl_) {};

  
    ad operator()(const ad& x){
      ad logssb = x;
      ad v = ptr->operator()(logssb, (ad)R_NaReal, (ad)R_NaReal) - logspace_add_SAM((ad)0.0, -exp(logl) * (exp(logssb) - exp(logd)));
      return v;
    }

    USING_RECFUN;

  };

  template <class Type>
  std::shared_ptr<RecruitmentNumeric<Type> > Rec_DepensatoryC(std::shared_ptr<RecruitmentWorker<ad> > SR,
							      ad logd,
							      ad logl){
    return std::make_shared<RecruitmentNumeric<Type> >(std::make_shared<WrapDepensatoryC>(SR, logd, logl));
  }


  /*
    Form D: R_D(S; d) = R(S * (1 - exp(log(0.5) * S/d)))
  */
  struct WrapDepensatoryD : RecruitmentFunctor {
    std::shared_ptr<RecruitmentWorker<ad> > ptr;
    ad logd;
    ad logl;

    WrapDepensatoryD(std::shared_ptr<RecruitmentWorker<ad> >& p_,
		     ad& logd_,
		     ad& logl_) : ptr(p_), logd(logd_), logl(logl_) {};

  
    ad operator()(const ad& x){
      ad logssb = x + logspace_sub_SAM((ad)0.0, log(0.5 * exp(x - logd)));
      ad v = ptr->operator()(logssb, (ad)R_NaReal, (ad)R_NaReal);
      return v;
    }

    USING_RECFUN;

  };

  template <class Type>
  std::shared_ptr<RecruitmentNumeric<Type> > Rec_DepensatoryD(std::shared_ptr<RecruitmentWorker<ad> > SR,
							      ad logd,
							      ad logl){
    return std::make_shared<RecruitmentNumeric<Type> >(std::make_shared<WrapDepensatoryD>(SR, logd, logl));
  }


  //////////////////////////////////////////////////////////////////////////////////////////
  // Implementation of recruitment functions
  //////////////////////////////////////////////////////////////////////////////////////////
  

  // Recruitment function -2
  // No recruitment
  template<class Type>
  struct Rec_ICESforecast : RecruitmentWorker<Type> {
    Type logMean;
    Type logSd;

    Rec_ICESforecast(Type lm, Type lsd) : RecruitmentWorker<Type>(0,0), logMean(lm), logSd(lsd) {}
  
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

    Rec_None() : RecruitmentWorker<Type>(0,0) {};
  
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

    Rec_LogRW() : RecruitmentWorker<Type>(1,0) {};
  
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

    Rec_Ricker(Type la, Type lb) : RecruitmentWorker<Type>(0,0), loga(la), logb(lb) {};
  
    Type operator()(Type logssb, Type lastLogR, Type year){
      return loga + logssb - exp(logb + logssb);
    }
  
    Type logSe(Type logLambda){    
      Type v = TMBad::CondExpLe(loga+logLambda, Type(0.0), exp(loga+logLambda), (loga + logLambda) * exp(-logb));
      return log(v + SAM_Zero);
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

    Rec_BevertonHolt(Type la, Type lb) : RecruitmentWorker<Type>(0,0), loga(la), logb(lb) {};

    Type operator()(Type logssb, Type lastR, Type year){
      return loga + logssb - logspace_add_SAM(Type(0.0),logb + logssb);
    }
  
    Type logSe(Type logLambda){
      Type v = TMBad::CondExpLe(loga+logLambda, Type(0.0), Type(0.0), (exp(loga + logLambda) - 1.0) * exp(-logb));
      return log(v + 1.0e-16);
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

    Rec_ConstantMean(const vector<Type>& logRv, const vector<Type>& crb) : RecruitmentWorker<Type>(0,1), logRvalue(logRv), constRecBreaks(crb) {};

    Type operator()(Type logssb, Type lastLogR, Type year){
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
      return 0.0;
    }

  };


  // Recruitment function 60
  // Logistic hockey stick
  struct RF_LogisticHockeyStick_t : RecruitmentFunctor {
    ad loga;			// alpha
    ad logm;			// mu
    ad logt;			// theta

    RF_LogisticHockeyStick_t(ad la, ad lm, ad lt) :
      loga(la), logm(lm), logt(lt) {}

    ad operator()(const ad& logssb){
      ad thisSSB = exp(logssb);
      ad v= loga + logm + logt + log(1.0 + exp(-exp(-logt))) + log(exp(logssb - logm - logt) - log(1.0 + exp((thisSSB-exp(logm))/exp(logm + logt))) + log(1.0 + exp(-exp(-logt))));
      return v;
    }

    USING_RECFUN;
  
  };
  
  template<class Type>
  struct Rec_LogisticHockeyStick : RecruitmentNumeric<Type>  {

    Type loga;
    Type logm;
    Type logt;
  
    Rec_LogisticHockeyStick(Type la, Type lm, Type lt) :
      RecruitmentNumeric<Type>(std::make_shared<RF_LogisticHockeyStick_t>(la, lm, lt)), loga(la), logm(lm), logt(lt) {};

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
      RecruitmentWorker<Type>(0,0), loglevel(ll), logblim(lbl) {};
  
    Type operator()(Type logssb, Type lastLogR, Type year){
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

    Rec_LogAR1(Type ll, Type lp) : RecruitmentWorker<Type>(1,0), loglevel(ll), logitPhi(lp) {};
  
    Type operator()(Type logssb, Type lastLogR, Type year){
      return loglevel + (2.0 / (1.0 + exp(-logitPhi)) - 1.0) * (lastLogR - loglevel);
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
      return 0.0;
    }

  };


#define TO_REC_2PAR(NAME)						\
  template<class Type>							\
  struct Rec_##NAME : RecruitmentNumeric<Type>  {			\
    Rec_##NAME(Type p1, Type p2) :					\
    RecruitmentNumeric<Type>(std::make_shared<RF_##NAME##_t>(p1, p2)) {}; \
  };

#define TO_REC_3PAR(NAME)						\
  template<class Type>							\
  struct Rec_##NAME : RecruitmentNumeric<Type>  {			\
    Rec_##NAME(Type p1, Type p2, Type p3) :				\
    RecruitmentNumeric<Type>(std::make_shared<RF_##NAME##_t>(p1, p2, p3)) {}; \
  };

#define TO_REC_SPLINE(NAME)						\
  template<class Type>							\
  struct Rec_##NAME : RecruitmentNumeric<Type>  {			\
    Rec_##NAME(vector<Type> knots, vector<Type> pars) :			\
    RecruitmentNumeric<Type>(std::make_shared<RF_##NAME##_t>(knots, pars)) {}; \
  };

  

  // Recruitment function 63
  // Bent hyperbola
  /*
    Source: e.g., DOI:10.1093/icesjms/fsq055
    rec_pars(0): log-Blim
    rec_pars(1): log of half the slope from 0 to Blim
    rec_pars(2): log-Smoothness parameter
  */
  //   Se = (2.0 * sqrt(exp(2.0 * newPar.rec_pars(0)) + exp(2.0 * newPar.rec_pars(2)) / 4.0) / (lambda * exp(newPar.rec_pars(1) + logRecCorrection)) - 2.0 * exp(newPar.rec_pars(0)) - 2.0 * sqrt(exp(2.0 * newPar.rec_pars(0)) + exp(2.0 * newPar.rec_pars(2)) / 4.0)) / ( 1.0 / ((lambda * lambda * exp(2.0 * (newPar.rec_pars(1) + logRecCorrection)))) - 2.0 / (lambda * exp(newPar.rec_pars(1) + logRecCorrection))  );  

  struct RF_BentHyperbola_t : RecruitmentFunctor {
    ad logBlim;
    ad logHalfSlope;
    ad logSmooth;

    RF_BentHyperbola_t(ad loga_, ad logb_, ad logg_) :
      logBlim(loga_), logHalfSlope(logb_), logSmooth(logg_) {}

    ad operator()(const ad& logssb){
      ad thisSSB = exp(logssb);
      ad v = logHalfSlope + log(thisSSB + sqrt(exp(2.0 * logBlim) + (exp(2.0 * logSmooth) / 4.0)) -
				sqrt(pow(thisSSB-exp(logBlim),2) + (exp(2.0 * logSmooth) / 4.0)));
      return v;
    }


    USING_RECFUN;
    
  };

  TO_REC_3PAR(BentHyperbola);


  // Recruitment function 64
  // Power function with compensatory mortality property
  //   Se = exp(1.0 / (1.0 - invlogit(newPar.rec_pars(1))) * (newPar.rec_pars(0) + logRecCorrection + log(lambda)));
  //    mdsr = R_PosInf;


  struct RF_PowerCMP_t : RecruitmentFunctor{
    ad loga;
    ad logb;

    RF_PowerCMP_t(ad loga_, ad logb_) :
      loga(loga_), logb(logb_) {}
  
    ad operator()(const ad& logssb){
      return loga + invlogit(logb) * logssb;
    }

    USING_RECFUN;
 
  };

  TO_REC_2PAR(PowerCMP);

  // Recruitment function 65
  // Power function without compensatory mortality property
  //   Se = exp(1.0 / (1.0 - (exp(newPar.rec_pars(1)) + 1.0001)) * (newPar.rec_pars(0) + logRecCorrection + log(lambda)));
  //   mdsr = R_PosInf;


  struct RF_PowerNCMP_t : RecruitmentFunctor {
    ad loga;
    ad logb;

    RF_PowerNCMP_t(ad loga_, ad logb_) :
      loga(loga_), logb(logb_) {}
  
    ad operator()(const ad& logssb){
      return loga + (exp(logb)+1.0001) * logssb;
    }

    USING_RECFUN;
 
  };

  TO_REC_2PAR(PowerNCMP);
  


  // Recruitment function 66
  // Shepherd
  // Se = CppAD::CondExpGt((newPar.rec_pars(0) + logSPR), T(SAM_Zero),
  // 			  exp( newPar.rec_pars(1) + 1.0 / exp(newPar.rec_pars(2)) * log(fabs(exp(newPar.rec_pars(0)) * lambda - 1.0))),
  // 			  T(SAM_Zero));
  //   Se = exp( newPar.rec_pars(1) + 1.0 / exp(newPar.rec_pars(2)) * log(softmax(exp(newPar.rec_pars(0) + logRecCorrection) * lambda - 1.0,(T)SAM_Zero, (T)100.0)) );


  // xx = CppAD::CondExpGt(par.rec_pars(2),Type(0.0),exp(log((1.0+exp(par.rec_pars(2))) / (exp(par.rec_pars(2)) - 1.0)) / exp(par.rec_pars(2))), Type(1e-10));
  //   mdsr = -(-1.0 + (exp(par.rec_pars(2)) - 1.0)*pow(xx/exp(par.rec_pars(1)),exp(par.rec_pars(2))))*exp(par.rec_pars(0))/pow(pow(xx/exp(par.rec_pars(1)),exp(par.rec_pars(2))) + 1.0,2.0);


  struct RF_Shepherd_t : RecruitmentFunctor {
    ad loga;
    ad logb;
    ad logg;

    RF_Shepherd_t(ad loga_, ad logb_, ad logg_) :
      loga(loga_), logb(logb_), logg(logg_) {}

    ad operator()(const ad& logssb){
      return loga + logssb - logspace_add_SAM(ad(0.0), exp(logg) * (logssb - logb));
    }

    USING_RECFUN;
 
  };

  TO_REC_3PAR(Shepherd);

  // Recruitment function 67
  // Hassel/Deriso

  // Se = CppAD::CondExpGt((newPar.rec_pars(0) + logSPR), T(SAM_Zero),
  // 			  fabs(exp(exp(-newPar.rec_pars(2)) * (newPar.rec_pars(0) + logSPR)) - 1.0) * exp(-newPar.rec_pars(1)),
  // 			   T(SAM_Zero));
  //   Se = (exp(exp(-newPar.rec_pars(2)) * (newPar.rec_pars(0) + logRecCorrection + logSPR)) - 1.0) * exp(-newPar.rec_pars(1) - newPar.rec_pars(2));


  struct RF_HasselDeriso_t : RecruitmentFunctor {
    ad loga;
    ad logb;
    ad logg;

    RF_HasselDeriso_t(ad loga_, ad logb_, ad logg_) :
      loga(loga_), logb(logb_), logg(logg_) {}

    ad operator()(const ad& logssb){
      ad thisSSB = exp(logssb);
      return loga+logssb-exp(logg) * log(1.0+exp(logb + logg)*thisSSB);
    }

    USING_RECFUN;
 
  };

  TO_REC_3PAR(HasselDeriso);



  // Recruitment function 68
  // Saila-Lorda

  //   Se = Se_sl(lambda, exp(newPar.rec_pars(0) + logRecCorrection), exp(newPar.rec_pars(1)), exp(newPar.rec_pars(2)));



  struct RF_SailaLorda_t : RecruitmentFunctor {
    ad loga;
    ad logb;
    ad logg;

    RF_SailaLorda_t(ad loga_, ad logb_, ad logg_) :
      loga(loga_), logb(logb_), logg(logg_) {}

    ad operator()(const ad& logssb){
      ad thisSSB = exp(logssb);
      return loga+exp(logg) * logssb - exp(logb)*thisSSB;
    }

    USING_RECFUN;
 
  };

  TO_REC_3PAR(SailaLorda);


  // Recruitment function 69
  // Sigmoidal Beverton-Holt

  struct RF_SigmoidalBevHolt_t : RecruitmentFunctor{
    ad loga;
    ad logb;
    ad logg;

    RF_SigmoidalBevHolt_t(ad loga_, ad logb_, ad logg_) :
      loga(loga_), logb(logb_), logg(logg_) {}

    ad operator()(const ad& logssb){
      return loga+exp(logg) * logssb-log(1.0+exp(logb)*exp(exp(logg) * logssb));
    }

    USING_RECFUN;
  
  };

  TO_REC_3PAR(SigmoidalBevHolt);
  

  // Recruitment function 90
  // Spline with compensatory mortality property (non-increasing on log(R/SSB))

  struct RF_SplineCMP_t : RecruitmentFunctor {
    vector<ad> pars;
    vector<ad> knots;

    RF_SplineCMP_t(vector<ad> pars_, vector<ad> knots_) :
      pars(pars_), knots(knots_) {}

    ad operator()(const ad& logssb){
      return logssb + ibcdspline(logssb,knots,pars);
    }

    USING_RECFUN;
 
  };

  TO_REC_SPLINE(SplineCMP);

  
  struct RF_SplineConvexCompensatory_t : RecruitmentFunctor {
    vector<ad> pars;
    vector<ad> knots;		// These knots are for SSB, input knots are for logssb to be consistent with other splines.

    RF_SplineConvexCompensatory_t(vector<ad> pars_, vector<ad> knots_) :
      pars(pars_), knots(exp(knots_)) {}

    ad operator()(const ad& logssb){
      return logssb + iibcispline(exp(logssb),knots,pars);
    }

    USING_RECFUN;
  };

  TO_REC_SPLINE(SplineConvexCompensatory);


  // Recruitment function 91
  // Smooth spline (integrated spline on log(R/SSB))

  struct RF_SplineSmooth_t : RecruitmentFunctor {
    vector<ad> pars;
    vector<ad> knots;

    RF_SplineSmooth_t(vector<ad> pars_, vector<ad> knots_) :
      pars(pars_), knots(knots_) {}

    ad operator()(const ad& logssb){
      ad mu(pars(pars.size()-1));
      return logssb + ibcspline(logssb,knots,(vector<ad>)pars.segment(0, pars.size()-1)) + mu;
    }

    USING_RECFUN;
 
  };

  TO_REC_SPLINE(SplineSmooth)


  // Recruitment function 92
  // General spline (on log(R/SSB))

  struct RF_SplineGeneral_t : RecruitmentFunctor {
    vector<ad> pars;
    vector<ad> knots;

    RF_SplineGeneral_t(vector<ad> pars_, vector<ad> knots_) :
      pars(pars_), knots(knots_) {}

    ad operator()(const ad& logssb){
      ad mu(pars(pars.size()-1));
      return logssb + bcspline(logssb,knots,(vector<ad>)pars.segment(0, pars.size()-1)) + mu;
    }

    USING_RECFUN;
 
  };

  TO_REC_SPLINE(SplineGeneral);



  //////////////////////////////////////////////////////////////////////////////////////////

  // Numerical Ricker for testing
  struct RF_NumRicker_t : RecruitmentFunctor{
    ad loga;
    ad logb;

    RF_NumRicker_t(ad loga_, ad logb_) :
      loga(loga_), logb(logb_) {}

    ad operator()(const ad& logssb){
      return loga + logssb - exp(logb + logssb);
    }

    USING_RECFUN;
 
  };

  TO_REC_2PAR(NumRicker)
  
  
  // Numerical Beverton Holt
  struct RF_NumBevHolt_t : RecruitmentFunctor{
    ad loga;
    ad logb;

    RF_NumBevHolt_t(ad loga_, ad logb_) :
      loga(loga_), logb(logb_) {}

    ad operator()(const ad& logssb){
      return loga + logssb - logspace_add_SAM(ad(0.0), logb + logssb);
    }

    USING_RECFUN;
 
  };

  TO_REC_2PAR(NumBevHolt)

}
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////

template<class Type>
Recruitment<Type> makeRecruitmentFunction(const confSet& conf, const paraSet<Type>& par)SOURCE({
    RecruitmentConvenience::RecruitmentModel rm = static_cast<RecruitmentConvenience::RecruitmentModel>(conf.stockRecruitmentModelCode);
    typedef TMBad::ad_aug ad;
    Recruitment<Type> r;

    if(rm == RecruitmentConvenience::RecruitmentModel::NoRecruit){
      r = Recruitment<Type>("zero",std::make_shared<RecruitmentConvenience::Rec_None<Type> >());
   
      ////////////////////////////////////////// The Beginning //////////////////////////////////////////
    
    }else if(rm == RecruitmentConvenience::RecruitmentModel::LogRandomWalk){
      if(par.rec_pars.size() != 0)
	Rf_error("The random walk recruitment should not have any parameters.");
      r = Recruitment<Type>("log-random walk",std::make_shared<RecruitmentConvenience::Rec_LogRW<Type> >());
    }else if(rm == RecruitmentConvenience::RecruitmentModel::Ricker){
      if(par.rec_pars.size() != 2)
	Rf_error("The Ricker recruitment must have two parameters.");
      r = Recruitment<Type>("Ricker",std::make_shared<RecruitmentConvenience::Rec_Ricker<Type> >(par.rec_pars(0), par.rec_pars(1)));
    }else if(rm == RecruitmentConvenience::RecruitmentModel::BevertonHolt){
      if(par.rec_pars.size() != 2)
	Rf_error("The Beverton Holt recruitment must have two parameters.");
      r = Recruitment<Type>("Beverton-Holt",std::make_shared<RecruitmentConvenience::Rec_BevertonHolt<Type> >(par.rec_pars(0), par.rec_pars(1)));
    }else if(rm == RecruitmentConvenience::RecruitmentModel::ConstantMean){
      if(par.rec_pars.size() != conf.constRecBreaks.size() + 1)
	Rf_error("The constant mean recruitment should have one more parameter than constRecBreaks.");
      r = Recruitment<Type>("constant mean", std::make_shared<RecruitmentConvenience::Rec_ConstantMean<Type> >(par.rec_pars, conf.constRecBreaks));
    }else if(rm == RecruitmentConvenience::RecruitmentModel::LogisticHockeyStick){
      if(par.rec_pars.size() != 3)
	Rf_error("The logistic hockey stick recruitment should have three parameters.");
      r = Recruitment<Type>("logistic hockey stick",std::make_shared<RecruitmentConvenience::Rec_LogisticHockeyStick<Type> >(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)));
    }else if(rm == RecruitmentConvenience::RecruitmentModel::HockeyStick){
      if(par.rec_pars.size() != 2)
	Rf_error("The hockey stick recruitment should have two parameters.");
      r = Recruitment<Type>("hockey stick",std::make_shared<RecruitmentConvenience::Rec_HockeyStick<Type> >(par.rec_pars(0), par.rec_pars(1)));
    }else if(rm == RecruitmentConvenience::RecruitmentModel::LogAR1){
      if(par.rec_pars.size() != 2)
	Rf_error("The log-AR(1) recruitment should have two parameters.");
      r = Recruitment<Type>("log-AR(1)",std::make_shared<RecruitmentConvenience::Rec_LogAR1<Type> >(par.rec_pars(0), par.rec_pars(1)));
    }else if(rm == RecruitmentConvenience::RecruitmentModel::BentHyperbola){
      if(par.rec_pars.size() != 3)
	Rf_error("The bent hyperbola recruitment should have three parameters.");
      r = Recruitment<Type>("bent hyperbola",std::make_shared<RecruitmentConvenience::Rec_BentHyperbola<Type> >(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)));
    }else if(rm == RecruitmentConvenience::RecruitmentModel::Power_CMP){
      if(par.rec_pars.size() != 2)
	Rf_error("The power law recruitment should have two parameters.");
      r = Recruitment<Type>("CMP power law",std::make_shared<RecruitmentConvenience::Rec_PowerCMP<Type> >(par.rec_pars(0), par.rec_pars(1)));
    }else if(rm == RecruitmentConvenience::RecruitmentModel::Power_NCMP){
      if(par.rec_pars.size() != 2)
	Rf_error("The power law recruitment should have two parameters.");
      r = Recruitment<Type>("non-CMP power law",std::make_shared<RecruitmentConvenience::Rec_PowerNCMP<Type> >(par.rec_pars(0), par.rec_pars(1)));

      //////////////////////////////////////// 3 parameter models ///////////////////////////////////////
    
    }else if(rm == RecruitmentConvenience::RecruitmentModel::Shepherd){
      if(par.rec_pars.size() != 3)
	Rf_error("The Shepherd recruitment should have three parameters.");
      r = Recruitment<Type>("Shepherd",std::make_shared<RecruitmentConvenience::Rec_Shepherd<Type> >(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)));
    }else if(rm == RecruitmentConvenience::RecruitmentModel::Hassel_Deriso){
      if(par.rec_pars.size() != 3)
	Rf_error("The Hassel/Deriso recruitment should have three parameters.");
      r = Recruitment<Type>("Hassel/Deriso",std::make_shared<RecruitmentConvenience::Rec_HasselDeriso<Type> >(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)));
    }else if(rm == RecruitmentConvenience::RecruitmentModel::SailaLorda){
      if(par.rec_pars.size() != 3)
	Rf_error("The Saila-Lorda recruitment should have three parameters.");
      r = Recruitment<Type>("Saila-Lorda",std::make_shared<RecruitmentConvenience::Rec_SailaLorda<Type> >(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)));
    }else if(rm == RecruitmentConvenience::RecruitmentModel::SigmoidalBevertonHolt){
      if(par.rec_pars.size() != 3)
	Rf_error("The sigmoidal Beverton-Holt recruitment should have three parameters.");
      r = Recruitment<Type>("sigmoidal Beverton-Holt",std::make_shared<RecruitmentConvenience::Rec_SigmoidalBevHolt<Type> >(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)));

      //////////////////////////////////////// Spline recruitment ///////////////////////////////////////
    
    }else if(rm == RecruitmentConvenience::RecruitmentModel::Spline_CMP){
      if(par.rec_pars.size() != conf.constRecBreaks.size() + 1)
	Rf_error("The spline recruitment should have one more parameter than constRecBreaks.");
      r = Recruitment<Type>("CMP spline",std::make_shared<RecruitmentConvenience::Rec_SplineCMP<Type> >(par.rec_pars, conf.constRecBreaks.template cast<Type>()));
    }else if(rm == RecruitmentConvenience::RecruitmentModel::Spline_Smooth){
      if(par.rec_pars.size() != conf.constRecBreaks.size() + 1)
	Rf_error("The spline recruitment should have one more parameter than constRecBreaks.");
      r = Recruitment<Type>("smooth spline",std::make_shared<RecruitmentConvenience::Rec_SplineSmooth<Type> >(par.rec_pars, conf.constRecBreaks.template cast<Type>()));
    }else if(rm == RecruitmentConvenience::RecruitmentModel::Spline_General){
      if(par.rec_pars.size() != conf.constRecBreaks.size() + 1)
	Rf_error("The spline recruitment should have one more parameter than constRecBreaks.");
      r = Recruitment<Type>("unrestricted spline",std::make_shared<RecruitmentConvenience::Rec_SplineGeneral<Type> >(par.rec_pars, conf.constRecBreaks.template cast<Type>()));

    }else if(rm == RecruitmentConvenience::RecruitmentModel::Spline_ConvexCompensatory){
      if(par.rec_pars.size() != conf.constRecBreaks.size() + 1)
	Rf_error("The spline recruitment should have one more parameter than constRecBreaks.");
      r = Recruitment<Type>("convex-compensatory spline",std::make_shared<RecruitmentConvenience::Rec_SplineConvexCompensatory<Type> >(par.rec_pars, conf.constRecBreaks.template cast<Type>()));


      //////////////////////// S/(d+S) type Depensatory recruitment models //////////////////////////////

    }else if(rm == RecruitmentConvenience::RecruitmentModel::Depensatory_B_Ricker){
      if(par.rec_pars.size() != 3)
	Rf_error("The depensatory B Ricker recruitment should have three parameters.");
      r = Recruitment<Type>("type B depensatory Ricker",RecruitmentConvenience::Rec_DepensatoryB<Type>(std::make_shared<RecruitmentConvenience::Rec_Ricker<ad> >(par.rec_pars(0), par.rec_pars(1)),par.rec_pars(2)));
    }else if(rm == RecruitmentConvenience::RecruitmentModel::Depensatory_B_BevertonHolt){
      if(par.rec_pars.size() != 3)
	Rf_error("The depensatory B Beverton-Holt recruitment should have three parameters.");
      r = Recruitment<Type>("type B depensatory Beverton-Holt",RecruitmentConvenience::Rec_DepensatoryB<Type>(std::make_shared<RecruitmentConvenience::Rec_BevertonHolt<ad> >(par.rec_pars(0), par.rec_pars(1)),par.rec_pars(2)));
    }else if(rm == RecruitmentConvenience::RecruitmentModel::Depensatory_B_LogisticHockeyStick){
      if(par.rec_pars.size() != 4)
	Rf_error("The depensatory B logistic hockey stick recruitment should have four parameters.");
      r = Recruitment<Type>("type B depensatory logistic hockey stick",RecruitmentConvenience::Rec_DepensatoryB<Type>(std::make_shared<RecruitmentConvenience::Rec_LogisticHockeyStick<ad> >(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)),par.rec_pars(3)));
    }else if(rm == RecruitmentConvenience::RecruitmentModel::Depensatory_B_HockeyStick){
      if(par.rec_pars.size() != 3)
	Rf_error("The depensatory B hockey stick recruitment should have three parameters.");
      r = Recruitment<Type>("type B depensatory hockey stick",RecruitmentConvenience::Rec_DepensatoryB<Type>(std::make_shared<RecruitmentConvenience::Rec_HockeyStick<ad> >(par.rec_pars(0), par.rec_pars(1)),par.rec_pars(2)));
    }else if(rm == RecruitmentConvenience::RecruitmentModel::Depensatory_B_BentHyperbola){
      if(par.rec_pars.size() != 4)
	Rf_error("The depensatory B bent hyperbola recruitment should have four parameters.");
      r = Recruitment<Type>("type B depensatory bent hyperbola",RecruitmentConvenience::Rec_DepensatoryB<Type>(std::make_shared<RecruitmentConvenience::Rec_BentHyperbola<ad> >(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)),par.rec_pars(3)));
    }else if(rm == RecruitmentConvenience::RecruitmentModel::Depensatory_B_Power){
      if(par.rec_pars.size() != 3)
	Rf_error("The depensatory B power law recruitment should have three parameters.");
      r = Recruitment<Type>("type B depensatory CMP power law",RecruitmentConvenience::Rec_DepensatoryB<Type>(std::make_shared<RecruitmentConvenience::Rec_PowerCMP<ad> >(par.rec_pars(0), par.rec_pars(1)),par.rec_pars(2)));
    }else if(rm == RecruitmentConvenience::RecruitmentModel::Depensatory_B_Shepherd){
      if(par.rec_pars.size() != 4)
	Rf_error("The depensatory B Shepherd recruitment should have four parameters.");
      r = Recruitment<Type>("type B depensatory Shepherd",RecruitmentConvenience::Rec_DepensatoryB<Type>(std::make_shared<RecruitmentConvenience::Rec_Shepherd<ad> >(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)),par.rec_pars(3)));

    }else if(rm == RecruitmentConvenience::RecruitmentModel::Depensatory_B_Hassel_Deriso){
      if(par.rec_pars.size() != 4)
	Rf_error("The depensatory B Hassel/Deriso recruitment should have four parameters.");
      r = Recruitment<Type>("type B depensatory Hassel/Deriso",RecruitmentConvenience::Rec_DepensatoryB<Type>(std::make_shared<RecruitmentConvenience::Rec_HasselDeriso<ad> >(par.rec_pars(0), par.rec_pars(1), par.rec_pars(2)),par.rec_pars(3)));

    }else if(rm == RecruitmentConvenience::RecruitmentModel::Depensatory_B_Spline_CMP){
      if(par.rec_pars.size() != conf.constRecBreaks.size() + 2 && conf.constRecBreaks.size() >= 3)
	Rf_error("The depensatory B spline recruitment should have two parameters more than the number of knots which should be at least 3.");
      r = Recruitment<Type>("type B depensatory CMP spline",RecruitmentConvenience::Rec_DepensatoryB<Type>(std::make_shared<RecruitmentConvenience::Rec_SplineCMP<ad> >(par.rec_pars.segment(0,par.rec_pars.size()-1), conf.constRecBreaks),par.rec_pars(par.rec_pars.size()-1)));

    }else if(rm == RecruitmentConvenience::RecruitmentModel::Depensatory_B_Spline_ConvexCompensatory){
      if(par.rec_pars.size() != conf.constRecBreaks.size() + 2 && conf.constRecBreaks.size() >= 3)
	Rf_error("The depensatory B spline recruitment should have two parameters more than the number of knots which should be at least 3.");
      r = Recruitment<Type>("type B depensatory convex compensatory spline",RecruitmentConvenience::Rec_DepensatoryB<Type>(std::make_shared<RecruitmentConvenience::Rec_SplineConvexCompensatory<ad> >(par.rec_pars.segment(0,par.rec_pars.size()-1), conf.constRecBreaks),par.rec_pars(par.rec_pars.size()-1)));

    
      //////////////////////// 1/(1+exp(-e * (S-d))) type Depensatory recruitment models //////////////////////////////

    }else if(rm == RecruitmentConvenience::RecruitmentModel::Depensatory_C_Ricker){
      if(par.rec_pars.size() != 4)
	Rf_error("The depensatory C Ricker recruitment should have four parameters.");
      r = Recruitment<Type>("type C depensatory Ricker",RecruitmentConvenience::Rec_DepensatoryC<Type>(std::make_shared<RecruitmentConvenience::Rec_Ricker<ad> >(par.rec_pars(0), par.rec_pars(1)),par.rec_pars(2), par.rec_pars(3)));
    
    }else if(rm == RecruitmentConvenience::RecruitmentModel::Depensatory_C_BevertonHolt){
      if(par.rec_pars.size() != 4)
	Rf_error("The depensatory C Beverton-Holt recruitment should have four parameters.");
      r = Recruitment<Type>("type C depensatory Beverton-Holt",RecruitmentConvenience::Rec_DepensatoryC<Type>(std::make_shared<RecruitmentConvenience::Rec_BevertonHolt<ad> >(par.rec_pars(0), par.rec_pars(1)),par.rec_pars(2), par.rec_pars(3)));

    }else if(rm == RecruitmentConvenience::RecruitmentModel::Depensatory_C_Spline_CMP){
      if(par.rec_pars.size() != conf.constRecBreaks.size() + 3 && conf.constRecBreaks.size() >= 3)
	Rf_error("The depensatory C spline recruitment should have three parameters more than the number of knots which should be at least 3.");
      r = Recruitment<Type>("type C depensatory CMP spline",RecruitmentConvenience::Rec_DepensatoryC<Type>(std::make_shared<RecruitmentConvenience::Rec_SplineCMP<ad> >(par.rec_pars.segment(0,par.rec_pars.size()-2), conf.constRecBreaks),par.rec_pars(par.rec_pars.size()-2),par.rec_pars(par.rec_pars.size()-1)));

    }else if(rm == RecruitmentConvenience::RecruitmentModel::Depensatory_C_Spline_ConvexCompensatory){
      if(par.rec_pars.size() != conf.constRecBreaks.size() + 3 && conf.constRecBreaks.size() >= 3)
	Rf_error("The depensatory C spline recruitment should have three parameters more than the number of knots which should be at least 3.");
      r = Recruitment<Type>("type C depensatory convex compensatory spline",
			    RecruitmentConvenience::Rec_DepensatoryC<Type>(std::make_shared<RecruitmentConvenience::Rec_SplineConvexCompensatory<ad> >(par.rec_pars.segment(0,par.rec_pars.size()-2), conf.constRecBreaks),par.rec_pars(par.rec_pars.size()-2),par.rec_pars(par.rec_pars.size()-1)));

    
      ///////////////////////////////////////////// For testing /////////////////////////////////////////////
    }else if(rm == RecruitmentConvenience::RecruitmentModel::Num_Ricker){
      if(par.rec_pars.size() != 2)
	Rf_error("The Numeric Ricker recruitment must have two parameters.");
      r = Recruitment<Type>("numeric Ricker",std::make_shared<RecruitmentConvenience::Rec_NumRicker<Type> >(par.rec_pars(0), par.rec_pars(1)));

    }else if(rm == RecruitmentConvenience::RecruitmentModel::Num_BevertonHolt){
      if(par.rec_pars.size() != 2)
	Rf_error("The Numeric Beverton Holt recruitment must have two parameters.");
      r = Recruitment<Type>("numeric Beverton-Holt",std::make_shared<RecruitmentConvenience::Rec_NumBevHolt<Type> >(par.rec_pars(0), par.rec_pars(1)));

      ///////////////////////////////////////////// The End /////////////////////////////////////////////
    
    }else{
      Rf_error("Stock-recruitment model code not implemented.");
    }

    return r;
  
  });

SAM_SPECIALIZATION(Recruitment<double> makeRecruitmentFunction<double>(const confSet&,const paraSet<double>&));
SAM_SPECIALIZATION(Recruitment<TMBad::ad_aug> makeRecruitmentFunction<TMBad::ad_aug>(const confSet&,const paraSet<TMBad::ad_aug>&));

  
 
template<class Type>
Recruitment<Type> makeICESrecruitment(Type lm, Type sd)SOURCE({
  return Recruitment<Type>("ICES",std::make_shared<RecruitmentConvenience::Rec_ICESforecast<Type> >(lm,sd));
  })

  SAM_SPECIALIZATION(Recruitment<double> makeICESrecruitment<double>(double, double));
SAM_SPECIALIZATION(Recruitment<TMBad::ad_aug> makeICESrecruitment<TMBad::ad_aug>(TMBad::ad_aug, TMBad::ad_aug));
