#ifndef WITH_SAM_LIB
namespace hcr_fun {


  template <class Type>
  Type hcr_min(Type a, Type b){
    return 0.5 * (a + b - sqrt((a-b) * (a-b)));
  }

  SAM_SPECIALIZATION(double hcr_min(double, double));
  SAM_SPECIALIZATION(TMBad::ad_aug hcr_min(TMBad::ad_aug, TMBad::ad_aug));
  
  template <class Type>
  Type hcr_max(Type a, Type b){
    return 0.5 * (a + b + sqrt((a-b) * (a-b)));
  }

  SAM_SPECIALIZATION(double hcr_max(double, double));
  SAM_SPECIALIZATION(TMBad::ad_aug hcr_max(TMBad::ad_aug, TMBad::ad_aug));

} // end of namespace hcr_fun
#endif



template <class Type>
Type hcr(Type ssb, vector<Type> hcrConf)SOURCE({
    Type Ftarget = hcrConf(0);
    Type Forigin = hcrConf(1);
    Type Fcap = hcrConf(2);
    Type Borigin = hcrConf(3);
    Type Bcap = hcrConf(4);
    Type Btrigger = hcrConf(5);

    Type newF = Ftarget;
    if(fabs(Btrigger - Borigin) < 1e-12){
      newF = TMBad::CondExpLt(ssb,
			      Bcap,
			      Fcap,
			      Ftarget);			    
    }else{
      newF = TMBad::CondExpLt(ssb,
			      Bcap,
			      Fcap,
			      hcr_fun::hcr_min(Ftarget, hcr_fun::hcr_max(Forigin, Forigin + (ssb - Borigin) * (Ftarget - Forigin) / (Btrigger - Borigin))));
    }
    return log(newF);	  
  });

SAM_SPECIALIZATION(double hcr(double, vector<double>));
SAM_SPECIALIZATION(TMBad::ad_aug hcr(TMBad::ad_aug, vector<TMBad::ad_aug>));

