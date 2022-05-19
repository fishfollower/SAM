template <class Type>
vector<Type> predObsFun(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, array<Type> &logF, vector<Type> &logssb, vector<Type> &logtsb, vector<Type> &logfsb, vector<Type> &logCatch, vector<Type> &logLand){
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

  // Calculate predicted observations
  int f, ft, a, y, yy, scaleIdx, ma, pg;  // a is no longer just ages, but an attribute (e.g. age or length) 
  int minYear=dat.aux(0,0);
  Type zz=Type(0);
  for(int i=0;i<dat.nobs;i++){
    y=dat.aux(i,0)-minYear;
    f=dat.aux(i,1);
    ft=dat.fleetTypes(f-1);
    a=dat.aux(i,2)-conf.minAge;
    if(dat.aux(i,2)==dat.maxAgePerFleet(f-1)){ma=1;}else{ma=0;}
    pg=conf.maxAgePlusGroup(f-1);
    if(ft==3){a=0;}
    if(ft<3){ 
      zz = dat.natMor(y,a);
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
  	Rf_error("Unknown fleet code");
        return(0);
      break;
      
      case 2:
 	if((pg!=conf.maxAgePlusGroup(0))&&(a==(conf.maxAge-conf.minAge))){
          Rf_error("When maximum age for the fleet is the same as maximum age in the assessment it must be treated the same way as catches w.r.t. plusgroup configuration");
  	}

	if((ma==1) && (pg==1)){
	  pred(i)=0;
	  for(int aa=a; aa<=(conf.maxAge-conf.minAge); aa++){
	    zz = dat.natMor(y,aa);
            if(conf.keyLogFsta(0,aa)>(-1)){
              zz+=exp(logF(conf.keyLogFsta(0,aa),y));
	    }
	    pred(i)+=exp(logN(aa,y)-zz*dat.sampleTimes(f-1));
	  }
	  pred(i)=log(pred(i));
	}else{
          pred(i)=logN(a,y)-zz*dat.sampleTimes(f-1);
	}
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
	if(conf.keyBiomassTreat(f-1)==5){
	  Type tsb = 0;
	  for(int aa=a; aa<=(conf.maxAge-conf.minAge); aa++){
	    zz = dat.natMor(y,aa);
	    if(conf.keyLogFsta(0,aa)>(-1)){
	      zz+=exp(logF(conf.keyLogFsta(0,aa),y));
	    }
	    tsb +=  exp(logN(aa,y)-zz*dat.sampleTimes(f-1))*dat.stockMeanWeight(y,aa);
	  }
	  pred(i) = log(tsb) +par.logFpar(conf.keyLogFpar(f-1,a));
	}
        if(conf.keyBiomassTreat(f-1)==6){
          Type N = 0;
          for(int aa=a; aa<=(conf.maxAge-conf.minAge); aa++){
            zz = dat.natMor(y,aa);
            if(conf.keyLogFsta(0,aa)>(-1)){
              zz+=exp(logF(conf.keyLogFsta(0,aa),y));
            }
            N +=  exp(logN(aa,y)-zz*dat.sampleTimes(f-1));
          }
          pred(i) = log(N) +par.logFpar(conf.keyLogFpar(f-1,a));
        }
	break;
  
      case 4:
  	Rf_error("Unknown fleet code");
        return 0;
      break;
  
      case 5:// tags  
        if((a+conf.minAge)>conf.maxAge){a=conf.maxAge-conf.minAge;} 
	pred(i)=exp(log(dat.aux(i,6))+log(dat.aux(i,5))-logN(a,y)-log(1000))*releaseSurvivalVec(i);
      break;
  
      case 6:
  	Rf_error("Unknown fleet code");
        return 0;
      break;
  
      case 7:
  	Rf_error("Unknown fleet code");
        return 0;
      break;
  
      default:
  	Rf_error("Unknown fleet code");
        return 0 ;
      break;
    }    
  }
  return pred;
}
