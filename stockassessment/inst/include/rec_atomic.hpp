#pragma once
#ifndef SAM_REC_ATOMIC_HPP
#define SAM_REC_ATOMIC_HPP


namespace rec_atomic {

  
  template<class T>
  T logspace_add2_raw (const T &logx, const T &logy) {
    // Was:
    //  fmax2 (logx, logy) + log1p (exp (-fabs (logx - logy)));
    if(logx == logy)
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
    if(logx == logy)
      return R_NegInf;
    if(logy == R_NegInf)
      return(logx);
    if(logx < logy)
      Rf_error("logx < logy in logspace_sub2");
    return logx + atomic::robust_utils::R_Log1_Exp(logy - logx);
  }

  
  TMB_BIND_ATOMIC(logspace_sub2,
		  11,
		  logspace_sub2_raw(x[0], x[1]) )


  // struct LOGSPACE_CUMSUM_t {
  //   template<class Type>
  //   vector<Type> operator()(const vector<Type>& xx){
  //     vector<Type> r(xx.size());
  //     r.setZero();
  //     r(0) = xx(0);
  //     for(int i = 1; i < xx.size(); ++i){
  // 	Type logx = r(i-1);
  // 	Type logy = xx(i);
  // 	r(i) = TMBad::CondExpLt(logx,logy,logy + log1p(exp(logx - logy)),logx + log1p(exp(logy - logx)));
  //     }
  //     return r;
  //   }
  // };

    
  // struct LOGSPACE_SUB_t {
  //   template<class Type>
  //   vector<Type> operator()(const vector<Type>& xx){
  //     Type logx = xx(0);
  //     Type logy = xx(1);
  //     Type yy = logy - logx;
  //      yy = TMBad::CondExpGt(yy,(Type)0.0,(Type)0.0,yy);
  //     // ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x))) ;
  //     Type log1_exp = TMBad::CondExpGt(yy, Type(-0.6931472),log(-expm1(yy)), log1p(-exp(yy)));
  //     vector<Type> r(1);
  //     r(0) = logx + log1_exp;
  //     return r;
  //   }
  // };
  
}

// template<class Type>
// Type logspace_add_SAM(Type logx, Type logy){
//   atomic::AtomicLocal<rec_atomic::LOGSPACE_CUMSUM_t> F((rec_atomic::LOGSPACE_CUMSUM_t()));
//   vector<Type> xx(2);
//   xx(0) = logx;
//   xx(1) = logy;
//   return F(xx)(1);
// }



// template<class Type>
// Type logspace_sum_SAM(vector<Type> logx){
//   atomic::AtomicLocal<rec_atomic::LOGSPACE_CUMSUM_t> F((rec_atomic::LOGSPACE_CUMSUM_t()));
//   return F(logx)(logx.size()-1);
// }

// template<class Type>
// vector<Type> logspace_cumsum_SAM(vector<Type> logx){
//   atomic::AtomicLocal<rec_atomic::LOGSPACE_CUMSUM_t> F((rec_atomic::LOGSPACE_CUMSUM_t()));
//   return F(logx);
// }


// template<class Type>
// Type logspace_sub_SAM(Type logx, Type logy){
//   atomic::AtomicLocal<rec_atomic::LOGSPACE_SUB_t> F((rec_atomic::LOGSPACE_SUB_t()));
//   vector<Type> xx(2);
//   xx(0) = logx;
//   xx(1) = logy;
//   return F(xx)(0);
// }


template<class Type>
Type logspace_add2(Type logx, Type logy) {  
  if ( !CppAD::Variable(logx) && logx == Type(R_NegInf) )
    return logy;
  if ( !CppAD::Variable(logy) && logy == Type(R_NegInf) )
    return logx;
  CppAD::vector<Type> tx(3);
  tx[0] = logx;
  tx[1] = logy;
  tx[2] = 0; // order
  return rec_atomic::logspace_add2(tx)[0];
}

template<class Type>
Type logspace_sub2(Type logx, Type logy) {
  if ( !CppAD::Variable(logy) && logy == Type(R_NegInf) )
    return logx;
  CppAD::vector<Type> tx(3);
  tx[0] = logx;
  tx[1] = logy;
  tx[2] = 0; // order
  return rec_atomic::logspace_sub2(tx)[0];
}


// tiny_ad versions of logspace_add gives segmentation fault for Fcrash reference point
template<class Type>
Type logspace_add3(Type logx, Type logy){
  return TMBad::CondExpLt(logx,logy,logy + log1p(exp(logx - logy)),logx + log1p(exp(logy - logx)));
}

#ifndef logspace_add_SAM 
#define logspace_add_SAM logspace_add
#endif

template<class Type>
Type logspace_sub3(Type logx, Type logy){
  Type xx = logy-logx;
  // xx = TMBad::CondExpGt(xx,(Type)0.0,(Type)0.0,xx);
  // ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x))) ;
  Type log1_exp = TMBad::CondExpGt(xx, Type(-0.6931472),log(-expm1(xx)), log1p(-exp(xx)));
  return logx + log1_exp;
}

#ifndef logspace_sub_SAM
#define logspace_sub_SAM logspace_sub
#endif




template<class Type>
Type logspace_sum(vector<Type> logx){
  Type r = R_NegInf;
  for(int i = 0; i < logx.size(); ++i)
    r = logspace_add_SAM(r, logx(i));
  return r;
}


#endif
