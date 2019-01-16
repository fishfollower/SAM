
template <class Type>
vector<Type> predNFun(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, array<Type> &logF, int i){
  int stateDimN=logN.dim[0];

  vector<Type> predN(stateDimN); 
  Type thisSSB=Type(0);
<<<<<<< HEAD
  if((i-conf.minAge)>=0){
    thisSSB=ssbi(dat,conf,logN,logF,i-conf.minAge);
  }else{
    thisSSB=ssbi(dat,conf,logN,logF,0); // use first in beginning       
  } 

  int ii, usepar;
  switch(conf.stockRecruitmentModelCode){
    case 0: // straight RW 
      predN(0)=logN(0,i-1);
    break;
    case 1: //ricker
      predN(0)=par.rec_loga(0)+log(thisSSB)-exp(par.rec_logb(0))*thisSSB;
    break;
    case 2:  //BH
      predN(0)=par.rec_loga(0)+log(thisSSB)-log(1.0+exp(par.rec_logb(0))*thisSSB); 
    break;
    case 3: //Constant mean
      usepar=0;
      for(ii=0; ii<conf.constRecBreaks.size(); ++ii){
        if(dat.years(i)>conf.constRecBreaks(ii)){usepar++;}
      }
      predN(0)=par.rec_loga(usepar);
    break;
  case 61: // Hockey stick
    predN(0) = par.rec_loga(0) - par.rec_logb(0) + log(thisSSB- (0.5 * (a+b+CppAD::abs((thisSSB - exp(log_blim))-Type(0.0)))));
    break;
  case 62: // AR1 (on log-scale)
    predN(0) = par.rec_loga(0) + (2.0 / (1.0 + exp(-par.rec_logb(0))) - 1.0) * (logN(0,i-1) - par.rec_loga(0));
    break;
  case 63: //Bent hypoerbola
    predN(0) = par.rec_logb(0) +
      log(thisSSB + sqrt(exp(2.0 * par.rec_loga(0)) + (20.0 * 20.0 / 4.0)) -
	  sqrt(pow(thisSSB-exp(par.rec_loga(0)),2) + (20.0 * 20.0 / 4.0)));
    break;
    default:
      error("SR model code not recognized");
    break;
  }
  
  for(int j=1; j<stateDimN; ++j){
    if(conf.keyLogFsta(0,j-1)>(-1)){
      predN(j)=logN(j-1,i-1)-exp(logF(conf.keyLogFsta(0,j-1),i-1))-dat.natMor(i-1,j-1); 
    }else{
      predN(j)=logN(j-1,i-1)-dat.natMor(i-1,j-1); 
    }
  }
  if(conf.maxAgePlusGroup(0)==1){// plusgroup adjustment if catches need them 
    predN(stateDimN-1)=log(exp(logN(stateDimN-2,i-1)-exp(logF(conf.keyLogFsta(0,stateDimN-2),i-1))-dat.natMor(i-1,stateDimN-2))+
                           exp(logN(stateDimN-1,i-1)-exp(logF(conf.keyLogFsta(0,stateDimN-1),i-1))-dat.natMor(i-1,stateDimN-1)));
  }
  return predN;  
}

