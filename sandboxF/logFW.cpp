#include <TMB.hpp>

template <class Type>
Type itrans(Type x){ // scaled ilogit
  return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);
}

template <class Type>
Type jacobiUVtrans(matrix<Type> logF){
  int nr=logF.rows();
  int nc=logF.cols();
  matrix<Type> A(nc,nc);
  for(int i=0; i<nc; ++i){
    for(int j=0; j<nc; ++j){
      A(i,j) = -1;
    }
  }
  for(int i=0; i<nc; ++i){
    A(0,i)=1;
  }
  for(int i=1; i<nc; ++i){
    A(i,i-1)=nc-1;
  }
  A/=nc;
  return nr*log(A.determinant());
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
  //PARAMETER_ARRAY(logU);
  //PARAMETER_VECTOR(logV);
  PARAMETER_MATRIX(logF);
  PARAMETER(logsdW);
  PARAMETER_MATRIX(logW);  
  vector<Type> sdU = exp(logsdU);
  Type sdV = exp(logsdV);
  Type sd = exp(logsd);  

  int nrow=Fobs.rows();
  int ncol=Fobs.cols();

  matrix<Type> logU(nrow,ncol-1);
  logU.setZero();
  vector<Type> logV(nrow);
  logV.setZero();
  for(int i=0; i<nrow; ++i){
    logV(i)=(logF-logW).row(i).mean();
  }
  for(int i=0; i<nrow; ++i){
    for(int j=0; j<ncol-1; ++j){
      logU(i,j)=logF(i,j)-logW(i,j)-logV(i);
    }
  }
  
  Type jnll=0;
 
  Type rhoU = itrans(betaU);
  Type rhoV = itrans(betaV);
  
  for(int y=0; y<nrow; ++y){
    for(int a=0; a<ncol; ++a){
      jnll+= -dnorm(logW(y,a),Type(0),exp(logsdW),true);
    }
  }
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
    vector<Type> diff=vector<Type>(logU.row(y))-rhoU*vector<Type>(logU.row(y-1))-alphaU;
    jnll += nldens(diff);
  }
  for(int y=1; y<nrow; ++y){
    jnll += -dnorm(logV(y),rhoV* logV(y-1) - alphaV ,sdV,true);
  }
  jnll += -jacobiUVtrans(logF);
  
  for(int y=0; y<nrow; ++y){
    for(int a=0; a<ncol; ++a){
      jnll += -dnorm(log(Fobs(y,a)), logF(y,a), sd, true);
    }
  }
  REPORT(logU)
  REPORT(logV)    
  return jnll;
}

