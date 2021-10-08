#ifndef SAM_REC_ATOMIC_HPP
#define SAM_REC_ATOMIC_HPP


namespace rec_atomic {

  
  template<class T>
  T logspace_add2_raw (const T &logx, const T &logy) {
    // Was:
    //  fmax2 (logx, logy) + log1p (exp (-fabs (logx - logy)));
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
    if(logy == R_NegInf)
      return(logx);
    return logx + atomic::robust_utils::R_Log1_Exp(logy - logx);
  }

  
  TMB_BIND_ATOMIC(logspace_sub2,
		  11,
		  logspace_sub2_raw(x[0], x[1]) )

  
}



template<class Type>
Type logspace_add2(Type logx, Type logy) {
  if ( !CppAD::Variable(logx) && logx == Type(-INFINITY) )
    return logy;
  if ( !CppAD::Variable(logy) && logy == Type(-INFINITY) )
    return logx;
  CppAD::vector<Type> tx(3);
  tx[0] = logx;
  tx[1] = logy;
  tx[2] = 0; // order
  return rec_atomic::logspace_add2(tx)[0];
}

template<class Type>
Type logspace_sub2(Type logx, Type logy) {
  if ( !CppAD::Variable(logy) && logy == Type(-INFINITY) )
    return logx;
  CppAD::vector<Type> tx(3);
  tx[0] = logx;
  tx[1] = logy;
  tx[2] = 0; // order
  return rec_atomic::logspace_sub2(tx)[0];
}


#endif
