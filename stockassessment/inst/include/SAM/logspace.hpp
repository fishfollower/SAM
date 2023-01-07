
#ifndef WITH_SAM_LIB

namespace rec_atomic {

  
  template<class T>
  T logspace_add2_raw (const T &logx, const T &logy) {
    // Was:
    //  fmax2 (logx, logy) + log1p (exp (-fabs (logx - logy)));
    //if(logx == logy)
    if(fabs(logx - logy) < 1e-16)
      return 2.0 + logx;
    if(logx == R_NegInf)
      return(logy);
    if(logy == R_NegInf)
      return(logx);
    return ( logx < logy ?
             logy + log1p (exp (logx - logy)) :
             logx + log1p (exp (logy - logx)) );
  }

  
  TMB_BIND_ATOMIC(logspace_add2,
  		  11,
  		  logspace_add2_raw(x[0], x[1]) )

  template<class T>
  T logspace_sub2_raw (const T &logx, const T &logy) {
    if(fabs(logx - logy) <= 1e-8)
      return R_NegInf;
    if(logy == R_NegInf)
      return(logx);
    if(logx < logy+1e-8){
      Rf_warning("logx < logy in logspace_sub2 (logx=%.9f, logy=%.9f, fabs=%.9f)",logx,logy,(logx-logy));
      return(R_NaReal);
    }
    return logx + atomic::robust_utils::R_Log1_Exp(logy - logx);
  }

  
  TMB_BIND_ATOMIC(logspace_sub2,
		  11,
		  logspace_sub2_raw(x[0], x[1]) )

 
}
#endif

template<class Type>
Type logspace_add2(Type logx, Type logy) SOURCE({  
  if ( !CppAD::Variable(logx) && logx == Type(R_NegInf) )
    return logy;
  if ( !CppAD::Variable(logy) && logy == Type(R_NegInf) )
    return logx;
  CppAD::vector<Type> tx(3);
  tx[0] = logx;
  tx[1] = logy;
  tx[2] = 0; // order
  return rec_atomic::logspace_add2(tx)[0];
  })

SAM_SPECIALIZATION(double logspace_add2(double, double));
SAM_SPECIALIZATION(TMBad::ad_aug logspace_add2(TMBad::ad_aug, TMBad::ad_aug));

template<class Type>
Type logspace_sub2(Type logx, Type logy) SOURCE({
  if ( !CppAD::Variable(logy) && logy == Type(R_NegInf) )
    return logx;
  CppAD::vector<Type> tx(3);
  tx[0] = logx;
  tx[1] = logy;
  tx[2] = 0; // order
  return rec_atomic::logspace_sub2(tx)[0];
  })

SAM_SPECIALIZATION(double logspace_sub2(double, double));
SAM_SPECIALIZATION(TMBad::ad_aug logspace_sub2(TMBad::ad_aug, TMBad::ad_aug));


// tiny_ad versions of logspace_add gives segmentation fault for Fcrash reference point
// template<class Type>
// Type logspace_add3(Type logx, Type logy){
//   return TMBad::CondExpLt(logx,logy,logy + log1p(exp(logx - logy)),logx + log1p(exp(logy - logx)));
// }

// template<class Type>
// Type logspace_sub3(Type logx, Type logy){
//   Type xx = logy-logx;
//   // xx = TMBad::CondExpGt(xx,(Type)0.0,(Type)0.0,xx);
//   // ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x))) ;
//   Type log1_exp = TMBad::CondExpGt(xx, Type(-0.6931472),log(-expm1(xx)), log1p(-exp(xx)));
//   return logx + log1_exp;
// }


template<class Type>
Type logspace_sum(vector<Type> logx)SOURCE({
  Type r = R_NegInf;
  for(int i = 0; i < logx.size(); ++i)
    r = logspace_add_SAM(r, logx(i));
  return r;
  })

SAM_SPECIALIZATION(double logspace_sum(vector<double>));
SAM_SPECIALIZATION(TMBad::ad_aug logspace_sum(vector<TMBad::ad_aug>));
