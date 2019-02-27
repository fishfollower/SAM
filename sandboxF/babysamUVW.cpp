#include <TMB.hpp>

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

bool isNAINT(int x){
  return NA_INTEGER==x;
}

template <class Type>
Type itrans(Type x){ // scaled ilogit
  return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);
}

template<class Type>
vector<Type> ssbFUN(matrix<Type> logN, matrix<Type> logFF, matrix<Type> M, matrix<Type> SW, matrix<Type> MO, matrix<Type> PF, matrix<Type> PM){
  int nrow = logN.rows();
  int ncol = logN.cols();
  vector<Type> ret(nrow);
  ret.setZero();
  for(int y=0; y<nrow; ++y){
    for(int a=0; a<ncol; ++a){
      ret(y) += SW(y,a)*MO(y,a)*exp(logN(y,a))*exp(-PF(y,a)*exp(logFF(y,a))-PM(y,a)*M(y,a));
    }
  }
  return ret;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(obs);
  DATA_IMATRIX(aux);
  DATA_IMATRIX(idx1);
  DATA_IMATRIX(idx2);
  DATA_INTEGER(minYear);
  DATA_INTEGER(minAge);    
  DATA_IVECTOR(fleetTypes);
  DATA_VECTOR(sampleTimes);     
  DATA_VECTOR(year);
  DATA_VECTOR(age);
  DATA_MATRIX(M);
  DATA_MATRIX(SW);
  DATA_MATRIX(MO);
  DATA_MATRIX(PF);
  DATA_MATRIX(PM);
  
  DATA_INTEGER(srmode);
  DATA_INTEGER(fcormode);
  DATA_IVECTOR(keyF);       
  DATA_IMATRIX(keyQ);   
  DATA_IMATRIX(keySd);
  DATA_IVECTOR(covType);
  DATA_IMATRIX(keyIGAR);
  DATA_IVECTOR(noParUS);
  DATA_INTEGER(nrowF);
  DATA_INTEGER(ncolF);
  
  int nobs=obs.size();    
  int nrow=M.rows();
  int ncol=M.cols();
  
  PARAMETER(logsdR);
  PARAMETER(logsdS);
  PARAMETER_VECTOR(rickerpar);
  PARAMETER_VECTOR(transRhoF); 
  PARAMETER_VECTOR(bhpar);    
  PARAMETER_VECTOR(logQ);
  PARAMETER_VECTOR(logsd);
  PARAMETER_VECTOR(logIGARdist);
  PARAMETER_VECTOR(parUS);
  PARAMETER_MATRIX(logN);
  PARAMETER_VECTOR(missing);
  
  PARAMETER_VECTOR(logsdU);
  PARAMETER(logsdV);
  PARAMETER(logsdW);
  PARAMETER(betaU);
  PARAMETER(betaV);
  PARAMETER_VECTOR(alphaU);
  PARAMETER(alphaV);
  PARAMETER_MATRIX(logU);
  PARAMETER_VECTOR(logV);
  PARAMETER_MATRIX(logW);
  
  Type sdR = exp(logsdR);
  Type sdS = exp(logsdS);
  vector<Type> sd=exp(logsd);

  vector<Type> sdU = exp(logsdU);
  Type sdV = exp(logsdV);
  Type sdW = exp(logsdW);  
  Type rhoU = itrans(betaU);
  Type rhoV = itrans(betaV);
  
  //patch missing 
  int idxmis=0; 
  for(int i=0; i<nobs; i++){
    if(isNA(obs(i))){
      obs(i)=exp(missing(idxmis++));
    }    
  }
  
  //Set logF equal logU + logV + logW
  matrix<Type> logF = logN*0;
  for(int y=0; y<nrowF; ++y){
    for(int a=0; a<(ncolF-1); ++a){
      logF(y,a) = logU(y,a) + logV(y) + logW(y,a);
    }
  }
  int aMax = ncolF-1;
  for(int y=0; y<nrowF; ++y){ //Sum to zero constraint
    logF(y,aMax) = 1- logU.row(y).sum() + logV(y)+ logW(y,aMax);
  }
  
  matrix<Type> logFF(nrow,ncol);
  for(int i=0; i<nrow; ++i){
    for(int j=0; j<ncol; ++j){
      logFF(i,j)=logF(i,keyF(j));
    }
  }
  
  Type jnll=0;
  
  vector<Type> ssb = ssbFUN(logN,logFF,M,SW,MO,PF,PM);
  
  // ############################################# N part
  Type predN;
  for(int y=1; y<nrow; ++y){
    switch(srmode){
    case 0:
      predN = logN(y-1,0);
      break;

    default:
      std::cout<<"Stock recruitment code not implemented yet."<<std::endl;
    exit(EXIT_FAILURE);
    break;
    }
    jnll += -dnorm(logN(y,0),predN,sdR,true);
  }
  
  for(int y=1; y<nrow; ++y){
    for(int a=1; a<ncol; ++a){
      predN=logN(y-1,a-1)-exp(logFF(y-1,a-1))-M(y-1,a-1);
      if(a==(ncol-1)){
        predN=log(exp(predN)+exp(logN(y-1,a)-exp(logFF(y-1,a))-M(y-1,a)));
      }
      jnll += -dnorm(logN(y,a),predN,sdS,true);
    }
  }
  
  // ############################################# F part
  matrix<Type> SigmaU(ncolF-1,ncolF-1);
  SigmaU.setZero();
  
  switch(fcormode){
  case 0:
    SigmaU.diagonal() = sdU*sdU;
    break;
    
  default:
    std::cout<<"Error: This fcormode not implemented yet."<<std::endl;
  exit(EXIT_FAILURE);
  break;
  }  
  
  density::MVNORM_t<Type> nldens(SigmaU);
  for(int y=1; y<nrowF; ++y){
    jnll += nldens(vector<Type>(logU.row(y))-rhoU* vector<Type>(logU.row(y-1) )-alphaU );
  }

  for(int y=1; y<nrowF; ++y){
    jnll += -dnorm(logV(y),rhoV* logV(y-1) - alphaV ,sdV,true);
  }
  
  for(int y=0; y<nrowF; ++y){
    for(int a=0; a<ncolF; ++a){
      jnll += -dnorm(logW(y,a),Type(0), sdW, true); // Add noize which does not come from U or V.
    }
  }

  /// obs part
  vector<Type> logPred(nobs);
  Type Z;
  int y, a, f;
  for(int i=0; i<nobs; ++i){
    y=aux(i,0)-minYear;
    f=aux(i,1)-1;
    a=aux(i,2)-minAge;
    Z=exp(logFF(y,a))+M(y,a);
    switch(fleetTypes(f)){
    case 0:
      logPred(i) = logN(y,a)-log(Z)+log(1-exp(-Z))+logFF(y,a);
      break;
      
    case 2:
      logPred(i) = logQ(keyQ(f,a))+logN(y,a)-Z*sampleTimes(f);
      break;
      
    default:
      std::cout<<"Error: This fleet type not implemented yet."<<std::endl;
    exit(EXIT_FAILURE);
    break;
    }
  }
  
  vector< density::MVNORM_t<Type> >  nllVec(fleetTypes.size());
  vector< matrix<Type> >  Svec(fleetTypes.size());
  
  for(int f=0; f<idx1.rows(); ++f){
    int d=0;
    for(int a=0; a<keySd.cols(); ++a){
      if(!isNAINT(keySd(f,a)))d++;
    };
    int thisdim=d;
    matrix<Type> S(thisdim,thisdim);
    S.setZero();
    if(covType(f)==0){
      d=0;
      for(int a=0; a<keySd.cols(); ++a){
        if(!isNAINT(keySd(f,a))){
          S(d,d)=sd(keySd(f,a))*sd(keySd(f,a));
          d++;
        }
      };
    }
    if(covType(f)==1){
      vector<Type> dist(thisdim);
      dist.setZero();
      d=1;
      for(int a=0; a<keyIGAR.cols(); ++a){
        if(!isNAINT(keyIGAR(f,a))){
          dist(d)=dist(d-1)+exp(logIGARdist(keyIGAR(f,a)));
          d++;
        }
      }
      vector<Type> sdvec(thisdim);
      d=0;
      for(int a=0; a<keySd.cols(); ++a){
        if(!isNAINT(keySd(f,a))){
          sdvec(d++)=sd(keySd(f,a));
        }
      }
      for(int i=0; i<S.rows(); ++i){
        for(int j=0; j<i; ++j){
          S(i,j) = pow(0.5,dist(i)-dist(j))*sdvec(i)*sdvec(j);
          S(j,i) = S(i,j);
        }
        S(i,i)=sdvec(i)*sdvec(i);
      }
    }
    if(covType(f)==2){
      vector<Type> sdvec(thisdim);
      d=0;
      for(int a=0; a<keySd.cols(); ++a){
        if(!isNAINT(keySd(f,a))){
          sdvec(d++)=sd(keySd(f,a));
        }
      }
      int from=0;
      for(int i=0; i<f; ++i) from+=noParUS(i);
      vector<Type> thispar=parUS.segment(from,noParUS(f));
      density::UNSTRUCTURED_CORR_t<Type> USDENS(thispar);
      matrix<Type> DD(S.rows(),S.cols());
      DD.setZero();
      for(int i=0; i<DD.rows(); ++i)DD(i,i)=sdvec(i);
      S = DD*(USDENS.cov())*DD;
    }
    
    nllVec(f).setSigma(S);
    Svec(f)=S;
  }
  
  for(int f=0; f<idx1.rows(); ++f){
    for(int y=0; y<idx1.cols(); ++y){
      if(!isNAINT(idx1(f,y))){
        vector<Type> o=obs.segment(idx1(f,y),idx2(f,y)-idx1(f,y)+1);
        vector<Type> p=logPred.segment(idx1(f,y),idx2(f,y)-idx1(f,y)+1);
        jnll += nllVec(f)(log(o)-p);
      }
    }
  }
  
  ADREPORT(ssb);
  return jnll;
}

