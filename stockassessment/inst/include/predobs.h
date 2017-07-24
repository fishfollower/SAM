template <class Type>
vector<Type> predObsFun(array<Type> &logF,
                        array<Type> &logN,
                        vector<Type> &logFpar,
                        vector<Type> &logScale,
                        vector<Type> &logQpow,
                        int nobs,
                        int minAge,
                        int maxAge,
                        int noScaledYears,
                        vector<int> &fleetTypes,
                        vector<int> &keyScaledYears,
                        array<int> &keyQpow,
                        vector<int> &keyBiomassTreat,
                        array<int> &aux,
                        array<int> &keyLogFsta,
                        array<int> &keyLogFpar,
                        array<int> &keyParScaledYA,
                        array<Type> &natMor,
                        vector<Type> &sampleTimes,
                        vector<Type> &logssb,
                        vector<Type> &logfsb,
                        vector<Type> &logCatch,
                        vector<Type> &releaseSurvivalVec){
  vector<Type> pred(nobs);
  pred.setZero();
  // Calculate predicted observations
  int f, ft, a, y, yy, scaleIdx;  // a is no longer just ages, but an attribute (e.g. age or length) 
  int minYear=aux(0,0);
  Type zz=Type(0);
  for(int i=0;i<nobs;i++){
    y=aux(i,0)-minYear;
    f=aux(i,1);
    ft=fleetTypes(f-1);
    a=aux(i,2)-minAge;
    if(ft==3){a=0;}
    if(ft<3){ 
      zz = natMor(y,a);
      if(keyLogFsta(0,a)>(-1)){
        zz+=exp(logF(keyLogFsta(0,a),y));
      }
    }    

    switch(ft){
      case 0:
        pred(i)=logN(a,y)-log(zz)+log(1-exp(-zz));
        if(keyLogFsta(f-1,a)>(-1)){
          pred(i)+=logF(keyLogFsta(0,a),y);
        }
        scaleIdx=-1;
        yy=aux(i,0);
        for(int j=0; j<noScaledYears; ++j){
          if(yy==keyScaledYears(j)){
            scaleIdx=keyParScaledYA(j,a);
            if(scaleIdx>=0){
              pred(i)-=logScale(scaleIdx);
            }
            break;
          }
        }
      break;
  
      case 1:
  	error("Unknown fleet code");
        return(0);
      break;
  
      case 2:
        pred(i)=logN(a,y)-zz*sampleTimes(f-1);
        if(keyQpow(f-1,a)>(-1)){
          pred(i)*=exp(logQpow(keyQpow(f-1,a))); 
        }
        if(keyLogFpar(f-1,a)>(-1)){
          pred(i)+=logFpar(keyLogFpar(f-1,a));
        }
        
      break;
  
      case 3:// biomass survey
        if(keyBiomassTreat(f-1)==0){
          pred(i) = logssb(y)+logFpar(keyLogFpar(f-1,a));
        }
        if(keyBiomassTreat(f-1)==1){
          pred(i) = logCatch(y)+logFpar(keyLogFpar(f-1,a));
        }
        if(keyBiomassTreat(f-1)==2){
          pred(i) = logfsb(y)+logFpar(keyLogFpar(f-1,a));
        }
      break;
  
      case 4:
  	error("Unknown fleet code");
        return 0;
      break;
  
      case 5:// tags  
        if((a+minAge)>maxAge){a=maxAge-minAge;} 
	pred(i)=exp(log(aux(i,6))+log(aux(i,5))-logN(a,y)-log(1000))*releaseSurvivalVec(i);
      break;
  
      case 6:
  	error("Unknown fleet code");
        return 0;
      break;
  
      case 7:
  	error("Unknown fleet code");
        return 0;
      break;
  
      default:
  	error("Unknown fleet code");
        return 0 ;
      break;
    }    
  }
  return pred;
}
