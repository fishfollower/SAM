/* WITH_SAM_LIB: Get headers */
/* SAM_COMPILEUNITS:  Use for separate compilation units*/

#ifndef SAM_DEPENDS
#define SAM_DEPENDS(x) ;
#endif

#ifdef WITH_SAM_LIB		/* Only get headers */
#undef SOURCE
#define SOURCE(...) ;
#else
#undef SOURCE
#define SOURCE(...) __VA_ARGS__;
#endif

#ifdef SAM_COMPILEUNITS
#ifdef WITH_SAM_LIB
#undef DEFARG
#define DEFARG(...) __VA_ARGS__
#undef SAM_SPECIALIZATION
#define SAM_SPECIALIZATION(...) extern template __VA_ARGS__
#undef HEADER
#define HEADER(...) __VA_ARGS__;
#else 
#undef DEFARG
#define DEFARG(...)
#undef SAM_SPECIALIZATION
#define SAM_SPECIALIZATION(...) template __VA_ARGS__
#undef HEADER
#define HEADER(...);
#endif
#else  /* Use as header only library */
#undef DEFARG
#define DEFARG(...) __VA_ARGS__
#undef SAM_SPECIALIZATION
#define SAM_SPECIALIZATION(...) ;
#undef HEADER
#define HEADER(...) __VA_ARGS__;
#endif


#ifndef SAM_NegInf
#define SAM_NegInf -20.0
#endif

#ifndef SAM_NIZero
#define SAM_NIZero -20.0
#endif

#ifndef SAM_Zero
#define SAM_Zero exp(SAM_NIZero)
#endif

#ifndef SAM_ASSERT
#define SAM_ASSERT(x,msg) if(!(x)) Rf_error(msg);
#endif

#ifndef REPORT_F
#define REPORT_F(name,F)					\
  if(isDouble<Type>::value && F->current_parallel_region<0) {		\
    Rf_defineVar(Rf_install(#name),					\
	      PROTECT(asSEXP(name)),F->report);			\
    UNPROTECT(1);						\
  }
#endif

#ifndef ADREPORT_F
#define ADREPORT_F(name,F) F->reportvector.push(name,#name);
#endif

#ifndef SIMULATE_F
#define SIMULATE_F(F)				\
  if(isDouble<Type>::value && F->do_simulate)
#endif

#ifndef NOT_SIMULATE_F
#define NOT_SIMULATE_F(F)				\
  if(!(isDouble<Type>::value && F->do_simulate))
#endif

#ifndef logspace_add_SAM
#define logspace_add_SAM logspace_add2
#endif


#ifndef logspace_sub_SAM
#define logspace_sub_SAM logspace_sub2
#endif
