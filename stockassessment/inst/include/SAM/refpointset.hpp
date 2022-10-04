SAM_DEPENDS(define)
SAM_DEPENDS(logspace)

HEADER(
template<class Type>
struct referencepointSet {

  enum CatchType {
		  totalCatch,
		  landings,
		  discard
  };

  int nYears;
  int rpType;
  vector<int> aveYears;
  vector<int> selYears;
  vector<Type> logCustomSel;
  vector<Type> xVal;
  CatchType catchType;
  vector<Type> logF0;		// Starting value for optimization
  vector<Type> logSel;

  inline referencepointSet() = default;

  referencepointSet(int nYears_, int CT, int i, array<Type> logF, confSet conf);
  
  referencepointSet(SEXP x);

  template<class T>
  inline referencepointSet(const referencepointSet<T>& other) :
    nYears(other.nYears),
    rpType(other.rpType),
    aveYears(other.aveYears),
    selYears(other.selYears),
    logCustomSel(other.logCustomSel),
    xVal(other.xVal),
    catchType(static_cast<typename referencepointSet<Type>::CatchType>((int)other.catchType)),
    logF0(other.logF0),
    logSel(other.logSel)
  {}    
 
  vector<Type> getLogSelectivity();
  
  void setLogSelectivity(array<Type>& logF, confSet& conf);

  Type logFbar(array<Type>& logF, confSet& conf);
});


SOURCE(
	 template<class Type>
	 referencepointSet<Type>::referencepointSet(int nYears_,
						    int CT,
						    int i,
						    array<Type> logF,
						    confSet conf) :
	 nYears(nYears_), rpType(-99), aveYears(1), selYears(1), logCustomSel(0), xVal(0), catchType(static_cast<typename referencepointSet<Type>::CatchType>(CT)), logF0(0), logSel(0) {
	   aveYears(0) = i;
	   selYears(0) = i;
	   setLogSelectivity(logF,conf);
	 }
	 )

SOURCE(
	 template<class Type>
	 referencepointSet<Type>::referencepointSet(SEXP x) {
	   if(!Rf_isNull(getListElement(x,"nYears"))){
	     nYears = Rf_asInteger(getListElement(x,"nYears"));
	   }else{
	     nYears = 0;
	   }
	   if(!Rf_isNull(getListElement(x,"rpType"))){
	     rpType = Rf_asInteger(getListElement(x,"rpType"));
	   }else{
	     rpType = -99;
	   }
	   if(!Rf_isNull(getListElement(x,"aveYears"))){
	     aveYears = asVector<int>(getListElement(x,"aveYears"));
	   }else{
	     aveYears = vector<int>(0);
	   }
	   if(!Rf_isNull(getListElement(x,"selYears"))){
	     selYears = asVector<int>(getListElement(x,"selYears"));
	   }else{
	     selYears = vector<int>(0);
	   }
	   if(!Rf_isNull(getListElement(x,"logCustomSel"))){
	     logCustomSel = asVector<Type>(getListElement(x,"logCustomSel"));
	   }else{
	     logCustomSel = vector<Type>(0);
	   }
	   if(!Rf_isNull(getListElement(x,"xVal"))){
	     xVal = asVector<Type>(getListElement(x,"xVal"));
	   }else{
	     xVal = vector<Type>(0);
	   }
	   if(!Rf_isNull(getListElement(x,"catchType"))){
	     catchType = static_cast<typename referencepointSet<Type>::CatchType>(Rf_asInteger(getListElement(x,"catchType")));
	   }else{
	     catchType = static_cast<typename referencepointSet<Type>::CatchType>(0);
	   }
	   if(!Rf_isNull(getListElement(x,"logF0"))){
	     logF0 = asVector<Type>(getListElement(x,"logF0"));
	   }else{
	     logF0 = vector<Type>(0);
	   }
	   logSel = vector<Type>(0);

	 }
	 )

SOURCE(
	 template<class Type>
	 vector<Type> referencepointSet<Type>::getLogSelectivity(){
	   if(logSel.size() == 0)
	     Rf_error("logSelectivity not set. Call with logF and conf.");
	   return logSel;
	 }
	 )

SOURCE(
	 template<class Type>
	 void referencepointSet<Type>::setLogSelectivity(array<Type>& logF, confSet& conf){
	   if(logCustomSel.size() == logF.rows()){
	     logSel = logCustomSel;
	     return;
	   }else if(logCustomSel.size() > 0){
	     Rf_error("Size of logCustomSel does not match logF");
	   }
	   vector<Type> logfbartmp(selYears.size());
	   logfbartmp.setConstant(R_NegInf);
	   Type logfsum = R_NegInf;
	   logSel = vector<Type>(logF.rows());
	   logSel.setConstant(R_NegInf);
	   for(int y = 0; y < selYears.size(); ++y){
	     for(int i = 0; i < logSel.size(); ++i){
	       logSel(i) = logspace_add_SAM(logSel(i), logF(i,selYears(y)));
	     }
	     for(int a = conf.fbarRange(0); a <= conf.fbarRange(1); a++){
	       for(int f = 0; f < conf.keyLogFsta.dim(0); ++f)
		 if(conf.keyLogFsta(f,a-conf.minAge) > (-1))
		   logfbartmp(y) = logspace_add_SAM(logfbartmp(y), logF(conf.keyLogFsta(f,a-conf.minAge),selYears(y)));
	     }
	     logfbartmp(y) -= log(Type(conf.fbarRange(1)-conf.fbarRange(0)+1));
	     logfsum = logspace_add_SAM(logfsum, logfbartmp(y));
	   }
	   logSel -= logfsum;
	   return;
	 }
	 )

SOURCE(
template<class Type>
Type referencepointSet<Type>::logFbar(array<Type>& logF, confSet& conf){
    vector<Type> logfbartmp(selYears.size());
    logfbartmp.setConstant(R_NegInf);
    Type logfsum = R_NegInf;
     for(int y = 0; y < selYears.size(); ++y){
       for(int a = conf.fbarRange(0); a <= conf.fbarRange(1); a++){  
	 for(int f = 0; f < conf.keyLogFsta.dim(0); ++f)
	   if(conf.keyLogFsta(f,a-conf.minAge) > (-1))
	     logfbartmp(y) = logspace_add_SAM(logfbartmp(y), logF(conf.keyLogFsta(f,a-conf.minAge),selYears(y)));
       }
       logfbartmp(y) -= log(Type(conf.fbarRange(1)-conf.fbarRange(0)+1));
       logfsum = logspace_add_SAM(logfsum, logfbartmp(y));
     }
     return logfsum;
  }
	 )

SAM_SPECIALIZATION(struct referencepointSet<double>);
SAM_SPECIALIZATION(struct referencepointSet<TMBad::ad_aug>);


HEADER(
template<class Type>
struct referencepointList : vector<referencepointSet<Type> > {

  newton::newton_config cfg;

  referencepointList();
  referencepointList(int i);
  referencepointList(SEXP x);

  template<class T>
  inline referencepointList(const listMatrixFromR<T>& other) : vector<referencepointSet<Type> >(other.size()), cfg(other.cfg) {
    for(int i = 0; i < other.size(); ++i)
      (*this)(i) = referencepointSet<T>(other(i));
  }
  
});


SOURCE(
	 template<class Type>
	 referencepointList<Type>::referencepointList() : vector<referencepointSet<Type> >(0), cfg(){}
	 )

SOURCE(
	 template<class Type>
	 referencepointList<Type>::referencepointList(int i) :
	 vector<referencepointSet<Type> >(i), cfg(){}
	 )

SOURCE(
	 template<class Type>
	 referencepointList<Type>::referencepointList(SEXP x) :
	 vector<referencepointSet<Type> >(Rf_length(x)), cfg(Rf_getAttrib(x,Rf_install("newton_config"))) {
	   for(int i = 0; i < Rf_length(x); ++i){
	     (*this)(i) = referencepointSet<Type>(VECTOR_ELT(x,i));
	   }
	 }
	 )


SAM_SPECIALIZATION(struct tmbutils::vector<referencepointSet<double> >);
SAM_SPECIALIZATION(struct tmbutils::vector<referencepointSet<TMBad::ad_aug> >);

SAM_SPECIALIZATION(struct referencepointList<double>);
SAM_SPECIALIZATION(struct referencepointList<TMBad::ad_aug>);

