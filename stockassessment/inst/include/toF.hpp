
namespace toF_atomic {

  // Catch to F
  
  template<class Float>
  struct CATCH2F : public UNIROOT<Float> {
    vector<Float> sel;
    vector<Float> M;
    vector<Float> N;
    vector<Float> w;
    Float catchval;

    CATCH2F(vector<Float> sel_,
	    vector<Float> M_,
	    vector<Float> N_,
	    vector<Float> w_,
	    Float catchval_) :
      UNIROOT<Float>(),
      sel(sel_),
      M(M_),
      N(N_),
      w(w_),
      catchval(catchval_) {};

    Float operator()(Float Fbar){
      vector<Float> Fa = Fbar * sel;
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
    int maxAge = (valIn.size() - 2) / 4;

    if(valIn.size() != 4 * maxAge + 2)
      Rf_error("Wrong lengths in catch2F");
    
    Float catchval = valIn(valIn.size() - 2);
    vector<Float> N = valIn.segment(0 * maxAge,maxAge);
    vector<Float> M = valIn.segment(1 * maxAge,maxAge);
    vector<Float> sel = valIn.segment(2 * maxAge,maxAge);
    vector<Float> w = valIn.segment(3 * maxAge,maxAge);
    CATCH2F<Float> f(sel, M, N, w, catchval);
    return f.R_zeroin2(0.0,100.0);
  }
  
  // TMB_BIND_ATOMIC_FLEX(getCatch2F, getCatch2F1(x))
  TMB_BIND_ATOMIC_FLEX_PART(getCatch2F, getCatch2F1(x), (tx.size()-2)/4-1)

  // Landings to F


  // Next SSB to F
}

template<class Type>
Type catch2F(Type catchval, vector<Type> sel, vector<Type> M, vector<Type> N, vector<Type> w) {
  int maxAge = sel.size();
  vector<Type> args(maxAge * 4 + 2);
  args.setZero();
  args.segment(0 * maxAge,maxAge) = N;
  args.segment(1 * maxAge,maxAge) = M;
  args.segment(2 * maxAge,maxAge) = sel;
  args.segment(3 * maxAge,maxAge) = w;
  args(args.size() - 2) = catchval;
  args(args.size() - 1) = 0; // Last index reserved for derivative order
  return toF_atomic::getCatch2F(CppAD::vector<Type>(args))[0];
}

