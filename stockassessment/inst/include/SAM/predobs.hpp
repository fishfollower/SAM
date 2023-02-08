SAM_DEPENDS(convenience)
SAM_DEPENDS(define)
SAM_DEPENDS(incidence)


#ifndef COMPILING_FROM_MSAM
void msam_reserved()SOURCE({
  return;
  })
#else
void msam_reserved()SOURCE({
  Rf_warning("Reserved for multiStockassessment");
  return;
  })
#endif


template<class Type>
Type predOneObs(int fleet,	// obs.aux(i,1)
		int fleetType,	// obs.fleetTypes(f-1)
		int age,	// obs.aux(i,2)-confA(s).minAge
		int year, // obs.aux(i,0)
		int minYear, // obs.aux(0,0)
		int noYearsLAI,
	        dataSet<Type>& dat,
		confSet& conf,
		paraSet<Type>& par,
		array<Type>& logF,
		array<Type>& logN,
		array<Type>& logPs,
		array<Type>& logitFseason,
		vector<Type>& varAlphaSCB,
		MortalitySet<Type>& mort,
		Type logssb,
		Type logtsb,
		Type logfsb,
		Type logCatch,
		Type logLand,
		Type logfbar,
		// Type tagv1,	     // dat.aux(i,5)
		// Type tagv2,	     // dat.aux(i,6)
		Type releaseSurvival, // releaseSurvivalVec(i)
		vector<Type> auxData // dat.aux.block(i,5,1,?)
		)SOURCE({
		    int f, ft, a, y, yy, scaleIdx, ma, pg;  // a is no longer just ages, but an attribute (e.g. age or length) 
		    y=year - minYear;
		    f=fleet;
		    ft=fleetType;
		    a=age-conf.minAge;

		    Type pred = 0.0;
		    // Type logzz = R_NegInf;
		    // Type zz = 0.0;
		    Type sumF = 0.0;
		    // Type vv = 0.0;
		    int LAIf = -1;
		    int lyr = -1;
		    int diffyears = 0;
		    if(ft == 6){
		      for(int lf = 0; lf < dat.noFleets; ++lf){
			if(dat.fleetTypes(lf) == 6)
			  ++LAIf;
			if(lf == (f-1))
			  break;
		      }
		      lyr = y - (dat.noYears - noYearsLAI);
		      if(lyr < 0){
			diffyears = - 1 * lyr;
			lyr = lyr + diffyears;
		      }
		    }
  
		    if(age==dat.maxAgePerFleet(f-1)){ma=1;}else{ma=0;}
		    pg=conf.maxAgePlusGroup(f-1);
		    if(ft==3){a=0;}
		    // if(ft<3){ 
		    //   logzz = log(dat.natMor(y,a));
		    //   // for(int fx = 0; fx < conf.keyLogFsta.dim[0]; ++fx)
		    //   // 	if(conf.keyLogFsta(fx,a)>(-1)){
		    //   // 	  logzz = logspace_add2(logzz, logF(conf.keyLogFsta(fx,a),y));
		    //   // 	}
		    //   //logzz = log(mort.totalZ(a,y));
		    // }    

		    switch(ft){
		    case 0:
		      //pred(i)=logN(a,y)-logzz+log(1-exp(-exp(logzz)));
		      // vv=-logzz+logspace_sub2(Type(0.0),-exp(logzz));
		      // if(conf.keyLogFsta(f-1,a)>(-1)){
		      //   vv+=logF(conf.keyLogFsta(f-1,a),y);
		      // }
		      // N * Survival until start * cumulative incidence
		      pred = logN(a,y) + mort.logFleetSurvival_before(a,y,f-1) + log(mort.fleetCumulativeIncidence(a,y,f-1));
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
			  // logzz = log(dat.natMor(y,aa));
			  // for(int fx = 0; fx < conf.keyLogFsta.dim[0]; ++fx)
			  //   if(conf.keyLogFsta(fx,aa)>(-1)){
			  //   logzz = logspace_add2(logzz, logF(conf.keyLogFsta(fx,aa),y));
			  // }
			  // pred+=exp(logN(aa,y)-exp(logzz)*dat.sampleTimes(f-1));
			  pred += exp(logN(aa,y) + mort.logFleetSurvival_before(aa,y,f-1));
			}
			pred=log(pred);
		      }else{
			//pred=logN(a,y)-exp(logzz)*dat.sampleTimes(f-1);
			pred = logN(a,y) + mort.logFleetSurvival_before(a,y,f-1);
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
			pred = logssb;
  		        if(conf.keyQpow(f-1,a)>(-1)){
			  pred*=exp(par.logQpow(conf.keyQpow(f-1,a))); 
  		        }
		        if(conf.keyLogFpar(f-1,a)>(-1)){
			  pred+=par.logFpar(conf.keyLogFpar(f-1,a));
		        }
		      }
		      if(conf.keyBiomassTreat(f-1)==1){
			pred = logCatch+par.logFpar(conf.keyLogFpar(f-1,a));
		      }
		      if(conf.keyBiomassTreat(f-1)==2){
			pred = logfsb+par.logFpar(conf.keyLogFpar(f-1,a));
			if(conf.keyQpow(f-1,0)>(-1)){
			  pred = logfsb*exp(par.logQpow(conf.keyQpow(f-1,0)))+par.logFpar(conf.keyLogFpar(f-1,a));
			}
		      }
		      if(conf.keyBiomassTreat(f-1)==3){
			pred = logCatch;
		      }
		      if(conf.keyBiomassTreat(f-1)==4){
			pred = logLand;
		      }
		      if(conf.keyBiomassTreat(f-1)==5){
		        Type tsbPred = 0;
		        for(int aa=a; aa<=(conf.maxAge-conf.minAge); aa++){
		          tsbPred += exp(logN(aa,y))*dat.stockMeanWeight(y,aa) * mort.fleetSurvival_before(aa,y,f-1);
		        }
			pred = log(tsbPred)+par.logFpar(conf.keyLogFpar(f-1,a));
		      }
		      if(conf.keyBiomassTreat(f-1)==6){
			Type N = 0;
			for(int aa=a; aa<=(conf.maxAge-conf.minAge); aa++){
			  // zz = dat.natMor(y,aa);
			  // for(int fx = 0; fx < conf.keyLogFsta.dim[0]; ++fx)
			  //   if(conf.keyLogFsta(fx,aa)>(-1)){
			  // 	zz+=exp(logF(conf.keyLogFsta(fx,aa),y));
			  //   }
			  // N +=  exp(logN(aa,y)-zz*dat.sampleTimes(f-1));
			  N += exp(logN(aa,y) + mort.logFleetSurvival_before(aa,y,f-1));
			}
			pred = log(N) + par.logFpar(conf.keyLogFpar(f-1,a));
		      }
		      if(conf.keyBiomassTreat(f-1)==10){
			pred = logfbar+par.logFpar(conf.keyLogFpar(f-1,a));
		      }		      
		      break;
  
		    case 4:
		      Rf_error("Unknown fleet code");
		      return 0;
		      break;
  
		    case 5:// tags
		      // (0) RecaptureY, (1) Yearclass, (2) Nscan, (3) R, (4) Type
		      SAM_ASSERT(auxData.size() >= 5,"aux is not large enough for fleet type 5");
		      if((a+conf.minAge)>conf.maxAge){a=conf.maxAge-conf.minAge;} 
		      pred=exp(log(auxData(3))+log(auxData(2))-logN(a,y)-log(1000.0))*releaseSurvival;
		      break;
  
		    case 6:
		      pred=logssb + par.logFpar(conf.keyLogFpar(f-1,0)) + logPs(LAIf,lyr) + varAlphaSCB(dat.minWeek(LAIf) + a);
		      if(conf.keyQpow(f-1,0)>(-1)){
			pred*=exp(par.logQpow(conf.keyQpow(f-1,0)));
		      }
		      break;
  
		    case 7:// sum residual fleets 
		      //pred=logN(a,y)-log(zz)+log(1-exp(-zz));
		      pred = logN(a,y);
		      sumF=0;
		      for(int ff=1; ff<=dat.noFleets; ++ff){
			if(dat.sumKey(f-1,ff-1)==1){
			  sumF += mort.fleetCumulativeIncidence(a,y,ff-1);
			  // if(conf.keyLogFsta(ff-1,a)>(-1)){
			  //   sumF+=exp(logF(conf.keyLogFsta(ff-1,a),y));
			  // }
			}
		      }
		      pred+=log(sumF);
		      break;

		    case 80:	// Catch/Landing proportion in part of year for total stock
		      {
		      SAM_ASSERT(auxData.size() >= 6,"aux is not large enough for fleet type 80");
		      Type Ctotal = 0.0;
		      Type Cseason = 0.0;
		      int flt0 = CppAD::Integer(auxData(0))-1; // TODO: allow auxData(0)==0 to sum over all fleets
		      int flt1 = flt0;
		      if(flt0 < 0){
			flt0 = 0;
			flt1 = dat.fleetTypes.size()-1;
		      }
		      for(int flt = flt0; flt <= flt1; ++flt){
			if(dat.fleetTypes(flt) == 0){
			int aMin = dat.minAgePerFleet(flt);
			int aMax = dat.maxAgePerFleet(flt);
			if(age > -1){
			  aMin = age;
			  aMax = age;
			}
			for(int aa = aMin - conf.minAge; aa < aMax - conf.minAge; ++aa){
			  //Type Cttmp = exp(logN(aa,y)) * mort.CIF(flt,aa,y,dat.sampleTimesStart(flt),dat.sampleTimesEnd(flt));
			  Type Cttmp = exp(logN(aa,y) + mort.logFleetSurvival_before(aa,y,flt) + log(mort.fleetCumulativeIncidence(aa,y,flt)));
			  Type Cstmp = exp(logN(aa,y)) * mort.partialCIF(flt,aa,y, auxData(1), auxData(2));
			  // 0: Catch numbers
			  if(CppAD::Integer(auxData(3)) == 1){ // 1: Catch weight
			    Cttmp *= dat.catchMeanWeight(y,aa, flt);
			    Cstmp *= dat.catchMeanWeight(y,aa, flt);
			  }else if(CppAD::Integer(auxData(3)) == 2){ // 2: Landing numbers
			    Cttmp *= dat.landFrac(y,aa, flt);
			    Cstmp *= dat.landFrac(y,aa, flt);
			  }else if(CppAD::Integer(auxData(3)) == 3){ // 3: Landing weight
			    Cttmp *= dat.landFrac(y,aa, flt) * dat.landMeanWeight(y,aa, flt);
			    Cstmp *= dat.landFrac(y,aa, flt) * dat.landMeanWeight(y,aa, flt);
			  }
			  Ctotal += Cttmp;
			  Cseason += Cstmp;
			}
			}
		      }
		      pred = log(Cseason) - log(Ctotal);		      
		      break;
		      }		 
		    case 90:	// Total Stock composition in catch/landing
		      msam_reserved();
		      pred = 0.0;
		      break;
		    case 92:	// Total area/Stock composition in catch/landing
		      msam_reserved();
		      pred = 0.0;
		      break;
		      
		    default:
		      Rf_error("Unknown fleet code");
		      return 0 ;
		      break;
		    }
		    return pred;
		  });
		  
SAM_SPECIALIZATION(double predOneObs(int, int,int,int,int, int,dataSet<double>&,confSet&,paraSet<double>&,array<double>&, array<double>&, array<double>&,array<double>&, vector<double>& ,MortalitySet<double>&,double,double,double,double,double,double,double,vector<double>));
SAM_SPECIALIZATION(TMBad::ad_aug predOneObs(int, int,int,int,int, int,dataSet<TMBad::ad_aug>&,confSet&,paraSet<TMBad::ad_aug>&,array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&,array<TMBad::ad_aug>&, vector<TMBad::ad_aug>&,MortalitySet<TMBad::ad_aug>&,TMBad::ad_aug,TMBad::ad_aug,TMBad::ad_aug,TMBad::ad_aug,TMBad::ad_aug,TMBad::ad_aug,TMBad::ad_aug,vector<TMBad::ad_aug>));


template <class Type>
vector<Type> predObsFun(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, array<Type> &logF, array<Type>& logPs, array<Type>& logitFseason, vector<Type>& varAlphaSCB, MortalitySet<Type>& mort, vector<Type> &logssb, vector<Type> &logtsb, vector<Type> &logfsb, vector<Type> &logCatch, vector<Type> &logLand, vector<Type> &logfbar, int noYearsLAI)SOURCE({
    vector<Type> pred(dat.nobs);
    pred.setConstant(R_NegInf);

    vector<Type> releaseSurvival(par.logitReleaseSurvival.size());
    vector<Type> releaseSurvivalVec(dat.nobs);
    if(par.logitReleaseSurvival.size()>0){
      releaseSurvival=invlogit(par.logitReleaseSurvival);
      for(int j=0; j<dat.nobs; ++j){
	int indx = CppAD::Integer(dat.auxData(j,4));
	if(!isNAINT(indx) && dat.fleetTypes(dat.aux(j,1)-1) == 5){
	  releaseSurvivalVec(j)=releaseSurvival(indx-1);
	}
      }
    }

    // Calculate predicted observations
    // int f, ft, a, y, yy, scaleIdx, ma, pg;  // a is no longer just ages, but an attribute (e.g. age or length) 
    // int minYear=dat.aux(0,0);
    // Type logzz=Type(R_NegInf);
    for(int i=0;i<dat.nobs;i++){

      vector<Type> auxData(dat.auxData.cols());
      for(int q = 0; q < auxData.size(); ++q)
	auxData(q) = dat.auxData(i,q);
      
      int y = dat.aux(i,0)-dat.aux(0,0);
      Type lfsb = R_NaReal;
      if(y < logfsb.size())
	lfsb = logfsb(y);
      Type lctch = R_NaReal;
      if(y < logCatch.size())
	lctch = logCatch(y);
      Type lland = R_NaReal;
      if(y < logLand.size())
	lland = logLand(y);
      pred(i) = predOneObs(dat.aux(i,1), // Fleet
			   dat.fleetTypes(dat.aux(i,1)-1), // FleetType
			   dat.aux(i,2),	      // Age
			   dat.aux(i,0),	      // Year
			   dat.aux(0,0),	      // minYear
			   noYearsLAI,
			   dat,
			   conf,
			   par,
			   logF,
			   logN,
			   logPs,
			   logitFseason,
			   varAlphaSCB,
			   mort,
			   logssb(y),
			   logtsb(y),
			   lfsb,
			   lctch,
			   lland,
			   logfbar(y),
			   // tagv1,
			   // tagv2,	     
			   releaseSurvivalVec(i), // releaseSurvival
			   auxData
			   );

    }

    return pred;
  });

SAM_SPECIALIZATION(vector<double> predObsFun(dataSet<double>&, confSet&, paraSet<double>&, array<double>&, array<double>&, array<double>&, array<double>&, vector<double>&, MortalitySet<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, int));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> predObsFun(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, array<TMBad::ad_aug>&, vector<TMBad::ad_aug>&, MortalitySet<TMBad::ad_aug>&, vector<TMBad::ad_aug>&, vector<TMBad::ad_aug>&, vector<TMBad::ad_aug>&, vector<TMBad::ad_aug>&, vector<TMBad::ad_aug>&, vector<TMBad::ad_aug>&, int));
