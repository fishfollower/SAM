template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x,eps,Type(0.01)*pow(x-eps,2),Type(0));
  return CppAD::CondExpGe(x,eps,x,eps/(Type(2)-x/eps));
}

template <class Type>
vector<Type> predObsFun(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, array<Type> &logF, vector<Type> &logssb, vector<Type> &logfsb, vector<Type> &logCatch, vector<Type> &logLand, objective_function<Type> *of){
  vector<Type> pred(dat.nobs);
  pred.setZero();

  
  vector<Type> releaseSurvival(par.logitReleaseSurvival.size());
  vector<Type> releaseSurvivalVec(dat.nobs);
  if(par.logitReleaseSurvival.size()>0){
    releaseSurvival=invlogit(par.logitReleaseSurvival);
    for(int j=0; j<dat.nobs; ++j){
      if(!isNAINT(dat.aux(j,7))){
        releaseSurvivalVec(j)=releaseSurvival(dat.aux(j,7)-1);
      }
    }
  }

  array<Type> logNdev(logN.dim[0],logN.dim[1]); 
  logNdev.setZero();
  if(conf.assignProcessNoiseToM==1){
    for(int y=0; y<(logNdev.dim[1]-1); ++y){
      for(int a=0; a<(logNdev.dim[0]-1); ++a){
        if(a==(logNdev.dim[0]-2)){ // plus group
	  Type FAm1y=0;	  
          if(conf.keyLogFsta(0,a)>(-1)){
            FAm1y=exp(logF(conf.keyLogFsta(0,a),y));
          }
          Type FAy=0;
          if(conf.keyLogFsta(0,a+1)>(-1)){
            FAy=exp(logF(conf.keyLogFsta(0,a+1),y));
          }	  
          Type pen1=0;
          //Type tmpDev=logN(a,y)-log(posfun(exp(logN(a+1,y+1))-exp(logN(a+1,y))*exp(-dat.natMor(y,a+1)-FAy),Type(1.0e-6),pen1))-dat.natMor(y,a);
	  Type tmpDev=logN(a+1,y+1)-log(exp(logN(a,y)-FAm1y-dat.natMor(y,a))+exp(logN(a+1,y)-FAy-dat.natMor(y,a+1)));
          logNdev(a,y)=tmpDev;
          logNdev(a+1,y)=tmpDev;
        }else{
          logNdev(a,y)=logN(a+1,y+1)-logN(a,y)+dat.natMor(y,a);
          if(conf.keyLogFsta(0,a)>(-1)){
            logNdev(a,y)+=exp(logF(conf.keyLogFsta(0,a),y));
          }
        }
      }
    }
  }
  REPORT_F(logNdev,of)


  // Calculate predicted observations
  int f, ft, a, y, yy, scaleIdx;  // a is no longer just ages, but an attribute (e.g. age or length) 
  int minYear=dat.aux(0,0);
  Type zz=Type(0);
  Type pen=0;
  for(int i=0;i<dat.nobs;i++){
    y=dat.aux(i,0)-minYear;
    f=dat.aux(i,1);
    ft=dat.fleetTypes(f-1);
    a=dat.aux(i,2)-conf.minAge;
    if(ft==3){a=0;}
    if(ft<3){ 
      zz = posfun(dat.natMor(y,a)-logNdev(a,y),Type(1.0e-6),pen);
      if(conf.keyLogFsta(0,a)>(-1)){
        zz+=exp(logF(conf.keyLogFsta(0,a),y));
      }
    }    

    switch(ft){
      case 0:
        pred(i)=logN(a,y)-log(zz)+log(1-exp(-zz));
        if(conf.keyLogFsta(f-1,a)>(-1)){
          pred(i)+=logF(conf.keyLogFsta(0,a),y);
        }
        scaleIdx=-1;
        yy=dat.aux(i,0);
        for(int j=0; j<conf.noScaledYears; ++j){
          if(yy==conf.keyScaledYears(j)){
            scaleIdx=conf.keyParScaledYA(j,a);
            if(scaleIdx>=0){
              pred(i)-=par.logScale(scaleIdx);
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
        pred(i)=logN(a,y)-zz*dat.sampleTimes(f-1);
        if(conf.keyQpow(f-1,a)>(-1)){
          pred(i)*=exp(par.logQpow(conf.keyQpow(f-1,a))); 
        }
        if(conf.keyLogFpar(f-1,a)>(-1)){
          pred(i)+=par.logFpar(conf.keyLogFpar(f-1,a));
        }
        
      break;
  
      case 3:// biomass or catch survey
        if(conf.keyBiomassTreat(f-1)==0){
          pred(i) = logssb(y)+par.logFpar(conf.keyLogFpar(f-1,a));
        }
        if(conf.keyBiomassTreat(f-1)==1){
          pred(i) = logCatch(y)+par.logFpar(conf.keyLogFpar(f-1,a));
        }
        if(conf.keyBiomassTreat(f-1)==2){
          pred(i) = logfsb(y)+par.logFpar(conf.keyLogFpar(f-1,a));
        }
        if(conf.keyBiomassTreat(f-1)==3){
          pred(i) = logCatch(y);
        }
        if(conf.keyBiomassTreat(f-1)==4){
          pred(i) = logLand(y);
        }
	break;
  
      case 4:
  	error("Unknown fleet code");
        return 0;
      break;
  
      case 5:// tags  
        if((a+conf.minAge)>conf.maxAge){a=conf.maxAge-conf.minAge;} 
	pred(i)=exp(log(dat.aux(i,6))+log(dat.aux(i,5))-logN(a,y)-log(1000))*releaseSurvivalVec(i);
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
