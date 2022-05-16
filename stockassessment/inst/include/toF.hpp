#pragma once
#ifndef SAM_TOF_HPP
#define SAM_TOF_HPP

// Faster versions

template<class Type>
struct CATCH2F_QUICK {
  vector<Type> Flast;
  vector<Type> M;
  vector<Type> N;
  vector<Type> w;
  vector<Type> frac;
  Type catchval;
    
  Type operator()(Type logFScale){
    vector<Type> Fa = exp(logFScale) * Flast;
    vector<Type> Z = Fa + M;
    vector<Type> C = Fa * (Type(1.0) - exp(-Z)) * N / Z * w * frac;
    return log(catchval) - log(sum(C) + 1e-5) + logFScale; //(sum(C) / catchval) * FScale;
  }    
};



template<class Type>
Type catch2F_quick(Type catchval, vector<Type> lastF, vector<Type> M, vector<Type> N, vector<Type> w) {
  Type sv = 0.0;
  int maxAge = lastF.size();
  vector<Type> frac(maxAge);
  frac.setZero();
  frac += 1.0;
  CATCH2F_QUICK<Type> f = {lastF, M, N, w, frac, catchval};
  for(int i = 0; i < 100; ++i){
    Type tmp = f(sv);
    sv = tmp;
  }
  return exp(sv);
}


template<class Type>
Type landing2F_quick(Type landingval, vector<Type> lastF, vector<Type> M, vector<Type> N, vector<Type> w, vector<Type> frac) {
  Type sv = 0.0;
  CATCH2F_QUICK<Type> f = {lastF, M, N, w, frac, landingval};
  for(int i = 0; i < 100; ++i){
    Type tmp = f(sv);
    sv = tmp;
  }
  return exp(sv);
}

template<class Type>
struct SSB2F_QUICK {
  vector<Type> logFlast;
    
  array<Type> logN;
  array<Type> logF;

  MortalitySet<Type> mort;
  
  confSet cf;
  dataSet<Type> ds;
  paraSet<Type> ps;
  int i;

  Type ssbval;
  Type rec_mean;

  Type f(Type logFScale){
    Recruitment<Type> recruit = Recruitment<Type>(new Rec_None<Type>());
    array<Type> lN = logN;
    array<Type> lF = logF;
    lF += logFScale;
    vector<Type> nextN = predNFun(ds, cf, ps, lN, lF, recruit, mort, i);
    if(!isNA(rec_mean))
      nextN(0) = rec_mean;
    int jj = i;
    if(i + 1 < lN.cols())
      jj = i+1;
    Type newSSB = 0.0;
    for(int q = 0; q < nextN.size(); ++q)
      newSSB += exp(nextN(q)) * ds.propMat(jj,q) * ds.stockMeanWeight(jj,q);    
    // if(i + 1 < logN.cols()){      
    //   lN.col(i+1) = nextN;
    //   lF.col(i+1) = logFScale + logFlast;
    //   newSSB = ssbi(ds, cf, lN, lF, i+1);
    // }else{
    //   lN.col(i) = nextN;
    //   lF.col(i) = logFScale + logFlast;
    //   newSSB = ssbi(ds, cf, lN, lF, i);
    // }
    return log(ssbval) - log(newSSB + 1e-5);
  }

  Type softmax(Type x, Type y, Type k = 1.0){
    return logspace_add2(k * x, k * y) / k;
  }
  Type sign0(Type x){
    return x / (fabs(x) + 1e-8);
  }
  Type numnewt(Type logs){
    Type h = 0.0001 * softmax(fabs(logs),Type(0.01), Type(100.0));
    Type a = f(logs);
    Type g = (-f(logs + 2.0 * h) + 8.0 * f(logs + h) - 8.0 * f(logs-h) + f(logs-2.0*h)) / (12 * h);
    // Type g = (f(logs + h) - f(logs - h)) / (2.0 * h);
    Type s = sign0(a) * sign0(g);
    Type y = log(fabs(a)) - log(softmax(fabs(g), (Type)0.001, (Type)1000.0)); // Damp the gradient
    return logs - 0.9 * s * exp(y); //softmax(exp(y), 0.5 * fabs(logs), Type(1000.0));
  }
  
};



template<class Type>
Type ssb2F_quick(Type ssbval, vector<Type> logFlast, dataSet<Type> dat, confSet conf, paraSet<Type> par, array<Type> logF, array<Type> logN, MortalitySet<Type> mort, int i, Type rec_mean) {
  Type sv = 0;
  SSB2F_QUICK<Type> f = {logFlast, logN, logF, mort, conf, dat, par, i, ssbval, rec_mean};
  for(int j = 0; j < 30; ++j){
    sv = f.numnewt(sv);
  }
  return exp(sv);
}


#endif
