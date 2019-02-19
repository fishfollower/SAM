#include <TMB.hpp>

template <class Type>
Type itrans(Type x){ // scaled ilogit
  return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(year)
  DATA_VECTOR(age)
  DATA_ARRAY(Fobs)
  DATA_INTEGER(cormode)

  PARAMETER_VECTOR(logsdU);
  PARAMETER(logsdV);
  PARAMETER(logsd);
  PARAMETER(betaU);
  PARAMETER(betaV);
  PARAMETER_VECTOR(alphaU);
  PARAMETER(alphaV);
  PARAMETER_ARRAY(logU);
  PARAMETER_VECTOR(logV);
  vector<Type> sdU = exp(logsdU);
  Type sdV = exp(logsdV);
  Type sd = exp(logsd);  

  int nrow=Fobs.rows();
  int ncol=Fobs.cols();
  
  Type jnll=0;
 
  Type rhoU = itrans(betaU);
  Type rhoV = itrans(betaV);
 
  matrix<Type> SigmaU(ncol-1,ncol-1);
  SigmaU.setZero();

  switch(cormode){
    case 0:
      SigmaU.diagonal() = sdU*sdU;
    break;
      
    default:
      std::cout<<"Error: This cormode not implemented yet."<<std::endl;
      exit(EXIT_FAILURE);
    break;
  }
  
  density::MVNORM_t<Type> nldens(SigmaU);
  for(int y=1; y<nrow; ++y){
    jnll += nldens(logU.col(y)-rhoU* logU.col(y-1) -alphaU );
  }
  for(int y=1; y<nrow; ++y){
    jnll += -dnorm(logV(y),rhoV* logV(y-1) - alphaV ,sdV,true);
  }
  
  array<Type> logF = Fobs*0;
  for(int y=0; y<nrow; ++y){
    for(int a=0; a<(ncol-1); ++a){
      logF(y,a) = logU(a,y) + logV(y);
      jnll += -dnorm(log(Fobs(y,a)), logF(y,a), sd, true);
    }
  }
  
  int a = ncol-1;
  for(int y=0; y<nrow; ++y){ //Sum to zero constraint
      logF(y,a) = 1- sum(logU.col(y)) + logV(y);
      jnll += -dnorm(log(Fobs(y,a)), logF(y,a), sd, true);
  }
  return jnll;
}

