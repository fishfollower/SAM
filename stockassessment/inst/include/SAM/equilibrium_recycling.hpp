

HEADER(
template<class Type>
struct PERREC_t {
  Type logFbar;
  Type logYPR;
  Type logSPR;
  Type logSe;
  Type logRe;
  Type logYe;
  Type dSR0;
  Type logLifeExpectancy;
  Type logYearsLost;
  Type logDiscYPR;
  Type logDiscYe;  
});

SAM_SPECIALIZATION(struct PERREC_t<double>);
SAM_SPECIALIZATION(struct PERREC_t<TMBad::ad_aug>);

template<class Type>
SEXP asSEXP(const PERREC_t<Type> &x) SOURCE({
    const char *resNms[] = {"logF", "logYPR", "logSPR", "logSe", "logRe", "logYe", "dSR0", "logLifeExpectancy", "logYearsLost","logDiscYe","logDiscYPR", ""}; // Must end with ""
    SEXP res;
    PROTECT(res = Rf_mkNamed(VECSXP, resNms));
    SET_VECTOR_ELT(res, 0, asSEXP(x.logFbar));
    SET_VECTOR_ELT(res, 1, asSEXP(x.logYPR));
    SET_VECTOR_ELT(res, 2, asSEXP(x.logSPR));
    SET_VECTOR_ELT(res, 3, asSEXP(x.logSe));
    SET_VECTOR_ELT(res, 4, asSEXP(x.logRe));
    SET_VECTOR_ELT(res, 5, asSEXP(x.logYe));
    SET_VECTOR_ELT(res, 6, asSEXP(x.dSR0));
    SET_VECTOR_ELT(res, 7, asSEXP(x.logLifeExpectancy));
    SET_VECTOR_ELT(res, 8, asSEXP(x.logYearsLost));
    SET_VECTOR_ELT(res, 9, asSEXP(x.logDiscYe));
    SET_VECTOR_ELT(res, 10, asSEXP(x.logDiscYPR));
    
    UNPROTECT(1);    
    return res;
  })


SAM_SPECIALIZATION(SEXP asSEXP(const PERREC_t<double>&));
SAM_SPECIALIZATION(SEXP asSEXP(const PERREC_t<TMBad::ad_aug>&));



HEADER(
template<class Type>
struct STOCHASTIC_PERREC_t {
  Type E_logFbar;
  
  Type E_logYPR;
  Type E_logSPR;
  Type E_logSe;
  Type E_logRe;
  Type E_logYe;
  Type E_logLifeExpectancy;
  Type E_logYearsLost;

  Type V_logFbar;
  Type V_logYPR;
  Type V_logSPR;
  Type V_logSe;
  Type V_logRe;
  Type V_logYe;
  Type V_logLifeExpectancy;
  Type V_logYearsLost;

  vector<Type> E_logN;
  matrix<Type> V_logN;

  vector<Type> lastLogNDiff;

  Type dSR0;
  
});

SAM_SPECIALIZATION(struct STOCHASTIC_PERREC_t<double>);
SAM_SPECIALIZATION(struct STOCHASTIC_PERREC_t<TMBad::ad_aug>);



template<class Type>
SEXP asSEXP(const STOCHASTIC_PERREC_t<Type> &x)
  SOURCE({
      const char *names[] = {
	"E_logFbar",
	"E_logYPR",
	"E_logSPR",
	"E_logSe",
	"E_logRe",
	"E_logYe",
	"E_logLifeExpectancy",
	"E_logYearsLost",
	"V_logFbar",
	"V_logYPR",
	"V_logSPR",
	"V_logSe",
	"V_logRe",
	"V_logYe",
	"V_logLifeExpectancy",
	"V_logYearsLost",
	"E_logN",
	"V_logN",
	"lastLogNDiff",
	"dSR0",
	""
      };  
      SEXP res = PROTECT(Rf_mkNamed(VECSXP, names));
      SET_VECTOR_ELT(res, 0, asSEXP(x.E_logFbar));
      SET_VECTOR_ELT(res, 1, asSEXP(x.E_logYPR));
      SET_VECTOR_ELT(res, 2, asSEXP(x.E_logSPR));
      SET_VECTOR_ELT(res, 3, asSEXP(x.E_logSe));
      SET_VECTOR_ELT(res, 4, asSEXP(x.E_logRe));
      SET_VECTOR_ELT(res, 5, asSEXP(x.E_logYe));
      SET_VECTOR_ELT(res, 6, asSEXP(x.E_logLifeExpectancy));
      SET_VECTOR_ELT(res, 7, asSEXP(x.E_logYearsLost));
      SET_VECTOR_ELT(res, 8, asSEXP(x.V_logFbar));
      SET_VECTOR_ELT(res, 9, asSEXP(x.V_logYPR));
      SET_VECTOR_ELT(res, 10, asSEXP(x.V_logSPR));
      SET_VECTOR_ELT(res, 11, asSEXP(x.V_logSe));
      SET_VECTOR_ELT(res, 12, asSEXP(x.V_logRe));
      SET_VECTOR_ELT(res, 13, asSEXP(x.V_logYe));
      SET_VECTOR_ELT(res, 14, asSEXP(x.V_logLifeExpectancy));
      SET_VECTOR_ELT(res, 15, asSEXP(x.V_logYearsLost));
      SET_VECTOR_ELT(res, 16, asSEXP(x.E_logN));
      SET_VECTOR_ELT(res, 17, asSEXP(x.V_logN));
      SET_VECTOR_ELT(res, 18, asSEXP(x.lastLogNDiff));
      SET_VECTOR_ELT(res, 19, asSEXP(x.dSR0));

   UNPROTECT(1);
   return res;
    })

SAM_SPECIALIZATION(SEXP asSEXP(const STOCHASTIC_PERREC_t<double>&));
SAM_SPECIALIZATION(SEXP asSEXP(const STOCHASTIC_PERREC_t<TMBad::ad_aug>&));



HEADER(
template<class Type>
struct EquilibriumRecycler {
  virtual PERREC_t<Type> operator()(Type logFbar) = 0;
}
       )
