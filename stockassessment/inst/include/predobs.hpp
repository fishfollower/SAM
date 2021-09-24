#ifndef SAM_PREDOBS_HPP
#define SAM_PREDOBS_HPP

template<class Type>
Type predOneObs(int fleet,	// obs.aux(i,1)
		int fleetType,	// obs.fleetTypes(f-1)
		int age,	// obs.aux(i,2)-confA(s).minAge
		int year, // obs.aux(i,0)
		int minYear, // obs.aux(0,0)
	        dataSet<Type> dat,
		confSet conf,
		paraSet<Type> par,
		array<Type> logF,
		array<Type> logN,
		Type logssb,
		Type logtsb,
		Type logfsb,
		Type logCatch,
		Type logLand,
		Type tagv1,	     // dat.aux(i,5)
		Type tagv2,	     // dat.aux(i,6)
		Type releaseSurvival // releaseSurvivalVec(i)
		){
  int f, ft, a, y, yy, scaleIdx, ma, pg;  // a is no longer just ages, but an attribute (e.g. age or length) 
  y=year - minYear;
  f=fleet;
  ft=fleetType;
  a=age-conf.minAge;

  Type pred = 0.0;
  Type logzz = R_NegInf;

   if(age==dat.maxAgePerFleet(f-1)){ma=1;}else{ma=0;}
    pg=conf.maxAgePlusGroup(f-1);
    if(ft==3){a=0;}
    if(ft<3){ 
      logzz = log(dat.natMor(y,a));
      if(conf.keyLogFsta(0,a)>(-1)){
        logzz = logspace_add2(logzz, logF(conf.keyLogFsta(0,a),y));
      }
    }    

    switch(ft){
      case 0:
        //pred(i)=logN(a,y)-logzz+log(1-exp(-exp(logzz)));
	pred=logN(a,y)-logzz+logspace_sub2(Type(0.0),-exp(logzz));
        if(conf.keyLogFsta(f-1,a)>(-1)){
          pred+=logF(conf.keyLogFsta(0,a),y);
        }
        scaleIdx=-1;
        yy=year;
        for(int j=0; j<conf.noScaledYears; ++j){
          if(yy==conf.keyScaledYears(j)){
            scaleIdx=conf.keyParScaledYA(j,a);
            if(scaleIdx>=0){
              pred-=par.logScale(scaleIdx);
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
	  pred=0;
	  for(int aa=a; aa<=(conf.maxAge-conf.minAge); aa++){
	    logzz = log(dat.natMor(y,aa));
            if(conf.keyLogFsta(0,aa)>(-1)){
              logzz = logspace_add2(logzz, logF(conf.keyLogFsta(0,aa),y));
	    }
	    pred+=exp(logN(aa,y)-exp(logzz)*dat.sampleTimes(f-1));
	  }
	  pred=log(pred);
	}else{
          pred=logN(a,y)-exp(logzz)*dat.sampleTimes(f-1);
	}
        if(conf.keyQpow(f-1,a)>(-1)){
          pred*=exp(par.logQpow(conf.keyQpow(f-1,a))); 
        }
        if(conf.keyLogFpar(f-1,a)>(-1)){
          pred+=par.logFpar(conf.keyLogFpar(f-1,a));
        }
        
      break;
  
      case 3:// biomass or catch survey
        if(conf.keyBiomassTreat(f-1)==0){
          pred = logssb+par.logFpar(conf.keyLogFpar(f-1,a));
        }
        if(conf.keyBiomassTreat(f-1)==1){
          pred = logCatch+par.logFpar(conf.keyLogFpar(f-1,a));
        }
        if(conf.keyBiomassTreat(f-1)==2){
          pred = logfsb+par.logFpar(conf.keyLogFpar(f-1,a));
        }
        if(conf.keyBiomassTreat(f-1)==3){
          pred = logCatch;
        }
        if(conf.keyBiomassTreat(f-1)==4){
          pred = logLand;
        }
        if(conf.keyBiomassTreat(f-1)==5){
          pred = logtsb+par.logFpar(conf.keyLogFpar(f-1,a));
        }
	break;
  
      case 4:
  	Rf_error("Unknown fleet code");
        return 0;
      break;
  
      case 5:// tags  
        if((a+conf.minAge)>conf.maxAge){a=conf.maxAge-conf.minAge;} 
	pred=exp(log(tagv2)+log(tagv1)-logN(a,y)-log(1000))*releaseSurvival;
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
    return pred;
}


template <class Type>
vector<Type> predObsFun(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, array<Type> &logF, vector<Type> &logssb, vector<Type> &logtsb, vector<Type> &logfsb, vector<Type> &logCatch, vector<Type> &logLand){
  vector<Type> pred(dat.nobs);
  pred.setConstant(R_NegInf);

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
  // int f, ft, a, y, yy, scaleIdx, ma, pg;  // a is no longer just ages, but an attribute (e.g. age or length) 
  // int minYear=dat.aux(0,0);
  // Type logzz=Type(R_NegInf);
  for(int i=0;i<dat.nobs;i++){
    Type tagv1 = 0.0;
    Type tagv2 = 0.0;
    if(dat.aux.cols() >= 7){
      tagv1 = dat.aux(i,5);
      tagv2 = dat.aux(i,6);
    }
    pred(i) = predOneObs(dat.aux(i,1), // Fleet
			 dat.fleetTypes(dat.aux(i,1)-1), // FleetType
			 dat.aux(i,2),	      // Age
			 dat.aux(i,0),	      // Year
			 dat.aux(0,0),	      // minYear
			 dat,
			 conf,
			 par,
			 logF,
			 logN,
			 logssb(dat.aux(i,0)-dat.aux(0,0)),
			 logtsb(dat.aux(i,0)-dat.aux(0,0)),
			 logfsb(dat.aux(i,0)-dat.aux(0,0)),
			 logCatch(dat.aux(i,0)-dat.aux(0,0)),
			 logLand(dat.aux(i,0)-dat.aux(0,0)),
			 tagv1,
			 tagv2,	     
			 releaseSurvivalVec(i) // releaseSurvival
			 );
    
    // y=dat.aux(i,0)-minYear;
    // f=dat.aux(i,1);
    // ft=dat.fleetTypes(f-1);
    // a=dat.aux(i,2)-conf.minAge;
    
    // if(dat.aux(i,2)==dat.maxAgePerFleet(f-1)){ma=1;}else{ma=0;}
    // pg=conf.maxAgePlusGroup(f-1);
    // if(ft==3){a=0;}
    // if(ft<3){ 
    //   logzz = log(dat.natMor(y,a));
    //   if(conf.keyLogFsta(0,a)>(-1)){
    //     logzz = logspace_add2(logzz, logF(conf.keyLogFsta(0,a),y));
    //   }
    // }    

    // switch(ft){
    //   case 0:
    //     //pred(i)=logN(a,y)-logzz+log(1-exp(-exp(logzz)));
    // 	pred(i)=logN(a,y)-logzz+logspace_sub2(Type(0.0),-exp(logzz));
    //     if(conf.keyLogFsta(f-1,a)>(-1)){
    //       pred(i)+=logF(conf.keyLogFsta(0,a),y);
    //     }
    //     scaleIdx=-1;
    //     yy=dat.aux(i,0);
    //     for(int j=0; j<conf.noScaledYears; ++j){
    //       if(yy==conf.keyScaledYears(j)){
    //         scaleIdx=conf.keyParScaledYA(j,a);
    //         if(scaleIdx>=0){
    //           pred(i)-=par.logScale(scaleIdx);
    //         }
    //         break;
    //       }
    //     }
    //   break;
  
    //   case 1:
    // 	Rf_error("Unknown fleet code");
    //     return(0);
    //   break;
      
    //   case 2:
    // 	if((pg!=conf.maxAgePlusGroup(0))&&(a==(conf.maxAge-conf.minAge))){
    //       Rf_error("When maximum age for the fleet is the same as maximum age in the assessment it must be treated the same way as catches w.r.t. plusgroup configuration");
    // 	}

    // 	if((ma==1) && (pg==1)){
    // 	  pred(i)=0;
    // 	  for(int aa=a; aa<=(conf.maxAge-conf.minAge); aa++){
    // 	    logzz = log(dat.natMor(y,aa));
    //         if(conf.keyLogFsta(0,aa)>(-1)){
    //           logzz = logspace_add2(logzz, logF(conf.keyLogFsta(0,aa),y));
    // 	    }
    // 	    pred(i)+=exp(logN(aa,y)-exp(logzz)*dat.sampleTimes(f-1));
    // 	  }
    // 	  pred(i)=log(pred(i));
    // 	}else{
    //       pred(i)=logN(a,y)-exp(logzz)*dat.sampleTimes(f-1);
    // 	}
    //     if(conf.keyQpow(f-1,a)>(-1)){
    //       pred(i)*=exp(par.logQpow(conf.keyQpow(f-1,a))); 
    //     }
    //     if(conf.keyLogFpar(f-1,a)>(-1)){
    //       pred(i)+=par.logFpar(conf.keyLogFpar(f-1,a));
    //     }
        
    //   break;
  
    //   case 3:// biomass or catch survey
    //     if(conf.keyBiomassTreat(f-1)==0){
    //       pred(i) = logssb(y)+par.logFpar(conf.keyLogFpar(f-1,a));
    //     }
    //     if(conf.keyBiomassTreat(f-1)==1){
    //       pred(i) = logCatch(y)+par.logFpar(conf.keyLogFpar(f-1,a));
    //     }
    //     if(conf.keyBiomassTreat(f-1)==2){
    //       pred(i) = logfsb(y)+par.logFpar(conf.keyLogFpar(f-1,a));
    //     }
    //     if(conf.keyBiomassTreat(f-1)==3){
    //       pred(i) = logCatch(y);
    //     }
    //     if(conf.keyBiomassTreat(f-1)==4){
    //       pred(i) = logLand(y);
    //     }
    //     if(conf.keyBiomassTreat(f-1)==5){
    //       pred(i) = logtsb(y)+par.logFpar(conf.keyLogFpar(f-1,a));
    //     }
    // 	break;
  
    //   case 4:
    // 	Rf_error("Unknown fleet code");
    //     return 0;
    //   break;
  
    //   case 5:// tags  
    //     if((a+conf.minAge)>conf.maxAge){a=conf.maxAge-conf.minAge;} 
    // 	pred(i)=exp(log(dat.aux(i,6))+log(dat.aux(i,5))-logN(a,y)-log(1000))*releaseSurvivalVec(i);
    //   break;
  
    //   case 6:
    // 	Rf_error("Unknown fleet code");
    //     return 0;
    //   break;
  
    //   case 7:
    // 	Rf_error("Unknown fleet code");
    //     return 0;
    //   break;
  
    //   default:
    // 	Rf_error("Unknown fleet code");
    //     return 0 ;
    //   break;
    // }    
  }
  return pred;
}

#endif
