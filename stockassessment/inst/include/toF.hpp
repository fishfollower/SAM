
namespace toF_atomic {

  // Catch to F
  
  template<class Float>
  struct CATCH2F : public UNIROOT<Float> {
    vector<Float> Flast;
    vector<Float> M;
    vector<Float> N;
    vector<Float> w;
    vector<Float> frac;
    Float catchval;
    
    CATCH2F(vector<Float> Flast_,
	    vector<Float> M_,
	    vector<Float> N_,
	    vector<Float> w_,
	    vector<Float> frac_,
	    Float catchval_) :
      UNIROOT<Float>(),
      Flast(Flast_),
      M(M_),
      N(N_),
      w(w_),
      frac(frac_),
      catchval(catchval_) {};

    Float operator()(Float FScale){
      vector<Float> Fa = FScale * Flast;
      vector<Float> Z = Fa;
      Z += M;
      vector<Float> C = Fa;
      C *= (Float(1.0) - exp(-Z));
      C *= N;
      C /= Z;
      C *= w;
      C *= frac;
      return catchval - sum(C);
    }
  };
 
  template<class Float>
  Float getCatch2F1(CppAD::vector<Float> val) {
    vector<Float> valIn(val);
    //Float order = valIn(valIn.size() - 1);
    int maxAge = (valIn.size() - 2) / 5;

    if(valIn.size() != 5 * maxAge + 2)
      Rf_error("Wrong lengths in catch2F");
    
    Float catchval = valIn(valIn.size() - 2);
    vector<Float> N = valIn.segment(0 * maxAge,maxAge);
    vector<Float> lastF = valIn.segment(1 * maxAge,maxAge);
    vector<Float> M = valIn.segment(2 * maxAge,maxAge);
    vector<Float> w = valIn.segment(3 * maxAge,maxAge);
    vector<Float> frac = valIn.segment(4 * maxAge,maxAge);
    CATCH2F<Float> f(lastF, M, N, w, frac, catchval);
    return f.R_zeroin2(0.0,100.0);
  }
  
  TMB_BIND_ATOMIC_FLEX_PART(getCatch2F, getCatch2F1(x), (tx.size()-2)/5*2-1) // Only derivative of N and lastF



  // Next SSB to F
}

template<class Type>
Type catch2F(Type catchval, vector<Type> lastF, vector<Type> M, vector<Type> N, vector<Type> w) {
  int maxAge = lastF.size();
  vector<Type> args(maxAge * 5 + 2);
  args.setZero();
  args.segment(0 * maxAge,maxAge) = N;
  args.segment(1 * maxAge,maxAge) = lastF;
  args.segment(2 * maxAge,maxAge) = M;
  args.segment(3 * maxAge,maxAge) = w;
  args.segment(4 * maxAge,maxAge) = 1.0;
  args(args.size() - 2) = catchval;
  args(args.size() - 1) = 0; // Last index reserved for derivative order
  return toF_atomic::getCatch2F(CppAD::vector<Type>(args))[0];
}


template<class Type>
Type landing2F(Type landingval, vector<Type> lastF, vector<Type> M, vector<Type> N, vector<Type> w, vector<Type> frac) {
  int maxAge = lastF.size();
  vector<Type> args(maxAge * 5 + 2);
  args.setZero();
  args.segment(0 * maxAge,maxAge) = N;
  args.segment(1 * maxAge,maxAge) = lastF;
  args.segment(2 * maxAge,maxAge) = M;
  args.segment(3 * maxAge,maxAge) = w;
  args.segment(4 * maxAge,maxAge) = frac;
  args(args.size() - 2) = landingval;
  args(args.size() - 1) = 0; // Last index reserved for derivative order
  return toF_atomic::getCatch2F(CppAD::vector<Type>(args))[0];
}

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

