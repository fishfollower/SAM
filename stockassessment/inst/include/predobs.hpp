template <class Type>
vector<Type> predObsFun(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, array<Type> &logF, array<Type> &logP, vector<Type> &logssb, vector<Type> &logfsb, vector<Type> &logCatch, objective_function<Type> *of){
  vector<Type> pred(dat.nobs);
  pred.setZero();
  array<Type> totF=totFFun(conf, logF);
  //LAI parts
  array<Type> logPs = scalePFun(conf,dat,logP);
  vector<Type> varAlphaSCB = scaleWeekFun(par,dat,logP,of);
  int noYearsLAI = yearsPFun(conf,dat);
  //END LAI PARTS
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

  // Calculate predicted observations
  int f, ft, a, y, yy, scaleIdx, LAIf, lyr, alpha;  // a is no longer just ages, but an attribute (e.g. age or length) 
  int minYear=dat.aux(0,0);
  Type zz;
  Type sumF=Type(0); 
  for(int i=0;i<dat.nobs;i++){
    y=dat.aux(i,0)-minYear;
    f=dat.aux(i,1);
    
    //Get the LAI component
    LAIf = -1;
    for(int lf=0;lf<dat.noFleets;lf++){
  	  if(dat.fleetTypes(lf)==6){
  	  	++LAIf;
      }
      if(lf == f){
      	break;	
	  }
    }
    
    ft=dat.fleetTypes(f-1);
    a=dat.aux(i,2)-conf.minAge;
    if(ft==3){a=0;}
    zz=0;
    if(ft<3 || ft==7){ 
      zz = dat.natMor(y,a)+totF(a,y);
    }    

    switch(ft){
      case 0:// residual fleets (without effort information)
        pred(i)=logN(a,y)-log(zz)+log(1-exp(-zz));
        if(conf.keyLogFsta(f-1,a)>(-1)){
          pred(i)+=logF(conf.keyLogFsta(f-1,a),y);
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
  
      case 2:// standard survey 
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
		lyr = y - (dat.noYears - noYearsLAI);
		alpha = dat.minWeek(LAIf) + a;	
	    pred(i)=logssb(y) + par.logFpar(conf.keyLogFpar(f-1,0)) + logPs(LAIf,lyr) + varAlphaSCB(alpha);
      break;
  
      case 7:// sum residual fleets 
	pred(i)=logN(a,y)-log(zz)+log(1-exp(-zz));
        sumF=0;
        for(int ff=1; ff<=dat.noFleets; ++ff){
          if(dat.sumKey(f-1,ff-1)==1){
            if(conf.keyLogFsta(ff-1,a)>(-1)){
              sumF+=exp(logF(conf.keyLogFsta(ff-1,a),y));
            }
          }
        }
        pred(i)+=log(sumF);
      break;
  
      default:
  	error("Unknown fleet code");
        return 0 ;
      break;
    }    
  }
  
  
  REPORT_F(logPs,of);
  REPORT_F(varAlphaSCB,of);
  return pred;
}
