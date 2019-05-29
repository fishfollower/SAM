
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

