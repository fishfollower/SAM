#define REPORT_F(name,F)					\
if(isDouble<Type>::value && F->current_parallel_region<0) {     \
  defineVar(install(#name),                                     \
            asSEXP_protect(name),F->report);                    \
  UNPROTECT(1);                                                 \
}

#define ADREPORT_F(name,F) F->reportvector.push(name,#name);

#define SIMULATE_F(F)						\
if(isDouble<Type>::value && F->do_simulate)

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

bool isNAINT(int x){
  return NA_INTEGER==x;
}
