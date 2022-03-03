#ifndef SAM_PREDN_HPP
#define SAM_PREDN_HPP



template<class Type>
Type functionalStockRecruitment(Type thisSSB, vector<Type> rec_pars, int stockRecruitmentModelCode){
  Type predN = R_NegInf;
  switch(stockRecruitmentModelCode){
    // Old recruitment models
    case 0: // straight RW 
      Rf_error("Not a functional recruitment");
    break;
    case 1: //ricker
      predN=rec_pars(0)+log(thisSSB)-exp(rec_pars(1))*thisSSB;
    break;
    case 2:  //BH
      predN=rec_pars(0)+log(thisSSB)-log(1.0+exp(rec_pars(1))*thisSSB); 
    break;
    case 3: //Constant mean
      Rf_error("Not a functional recruitment");
    break;
    // Type 2 depensatory
    // Parameterization 1: S/(S+d) => log(thisSSB) - logspace_add(log(thisSSB),d)
    // Parameterization 2: (S/d)/((S/d)+1) => log(thisSSB) - log(d) - logspace_add(log(thisSSB)-log(d),0.0)
    // Parameterization 3:
  case 50: // Type 2 depensatory logistic hockey stick (R(S) * S/(S+d))
    predN = rec_pars(0) + rec_pars(1) + rec_pars(2) + log(1.0 + exp(-exp(-rec_pars(2)))) + log(exp(log(thisSSB)-rec_pars(1) - rec_pars(2)) - log(1.0 + exp((thisSSB-exp(rec_pars(1)))/exp(rec_pars(1) + rec_pars(2)))) + log(1.0 + exp(-exp(-rec_pars(2))))) +
      log(thisSSB) - logspace_add2(log(thisSSB),rec_pars(3));
      //log(thisSSB) - log(rec_pars(3)) - logspace_add2(log(thisSSB) - rec_pars(3), (Type)0.0);
    break;
  case 51:			// Type 2 depensatory Ricker (R(S) * S/(S+d))
    predN = rec_pars(0)+log(thisSSB)-exp(rec_pars(1))*thisSSB +
      log(thisSSB) - logspace_add2(log(thisSSB),rec_pars(2));
      //log(thisSSB) - log(rec_pars(2)) - logspace_add2(log(thisSSB) - rec_pars(2), (Type)0.0);
    break;
  case 52:			// Type 2 depensatory Beverton-Holt (BH * S/(S+d))
    predN = rec_pars(0)+log(thisSSB)-log(1.0+exp(rec_pars(1))*thisSSB) + 
      log(thisSSB) - logspace_add2(log(thisSSB),rec_pars(2));
      // log(thisSSB) - log(rec_pars(2)) - logspace_add2(log(thisSSB) - rec_pars(2), (Type)0.0);
    break;
  case 53:			// Type 2 depensatory hockey stick
     predN = rec_pars(0) - rec_pars(1) +
      log(thisSSB - (0.5 * ((thisSSB - exp(rec_pars(1)))+Type(0.0)+CppAD::abs((thisSSB - exp(rec_pars(1)))-Type(0.0))))) +
       log(thisSSB) - logspace_add2(log(thisSSB),rec_pars(2));
     break;
  case 60:	// logistic Hockey stick
    // 0: alpha
    // 1: mu
    // 2: theta
    predN = rec_pars(0) + rec_pars(1) + rec_pars(2) + log(1.0 + exp(-exp(-rec_pars(2)))) + log(exp(log(thisSSB)-rec_pars(1) - rec_pars(2)) - log(1.0 + exp((thisSSB-exp(rec_pars(1)))/exp(rec_pars(1) + rec_pars(2)))) + log(1.0 + exp(-exp(-rec_pars(2)))));
    break;
  case 61: // Hockey stick
    // Type log_level = rec_pars(0);
    // Type log_blim = rec_pars(1);
    // Type a = thisSSB - exp(log_blim);
    // Type b = 0.0;
    // Type cut = 0.5 * (a+b+CppAD::abs(a-b)); // max(a,b)
    predN = rec_pars(0) - rec_pars(1) +
      log(thisSSB - (0.5 * ((thisSSB - exp(rec_pars(1)))+Type(0.0)+CppAD::abs((thisSSB - exp(rec_pars(1)))-Type(0.0)))));
    break;
  case 62: // AR1 (on log-scale)
    Rf_error("Not a functional recruitment");
    break;
  case 63: //Bent hyperbola / Hockey-stick-like
    /*
      Source: e.g., DOI:10.1093/icesjms/fsq055
      rec_pars(0): log-Blim
      rec_pars(1): log of half the slope from 0 to Blim
      rec_pars(2): log-Smoothness parameter
     */
    predN = rec_pars(1) +
      log(thisSSB + sqrt(exp(2.0 * rec_pars(0)) + (exp(2.0 * rec_pars(2)) / 4.0)) -
    	  sqrt(pow(thisSSB-exp(rec_pars(0)),2) + (exp(2.0 * rec_pars(2)) / 4.0)));
    break;
  case 64: // Power CMP
    predN = rec_pars(0) + invlogit(rec_pars(1)) * log(thisSSB);
    break;
  case 65: // Power Non-CMP
    predN = rec_pars(0) + (exp(rec_pars(1))+1.0001) * log(thisSSB);
    break;
  case 66: // Shepherd
    predN = rec_pars(0) + log(thisSSB) - log(1.0 + exp(exp(rec_pars(2)) * (log(thisSSB) - rec_pars(1))));
    break;
  case 67: // Deriso
    predN = rec_pars(0)+log(thisSSB)-exp(rec_pars(2)) * log(1.0+exp(rec_pars(1) + rec_pars(2))*thisSSB); 
    break;
  case 68: // Saila-Lorda (Iles 1994)
    predN = rec_pars(0)+exp(rec_pars(2)) * log(thisSSB) - exp(rec_pars(1))*thisSSB;
    break;
  case 69: // Sigmoidal Beverton-Holt
    predN = rec_pars(0)+exp(rec_pars(2)) * log(thisSSB)-log(1.0+exp(rec_pars(1))*exp(exp(rec_pars(2)) * log(thisSSB)));
    break;
  case 90: // Non-increasing spline on log(R/S)
    Rf_error("Needs more info");
    break;
  case 91: // Integrated spline on log(R/S)
    Rf_error("Needs more info");
    break;
  case 92: // Spline on log(R/S)
    Rf_error("Needs more info");
    break;
  case 93: // Type 2 depensatory Non-increasing spline on log(R/S)
    Rf_error("Needs more info");
    break;
  default:    
      Rf_error("SR model code not recognized");
    break;
  }
  return predN;
}


template <class Type>
vector<Type> predNFun(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, array<Type> &logF, int i){
  int stateDimN=logN.dim[0];

  vector<Type> predN(stateDimN);
  predN.setZero();
  Type thisSSB=Type(0);

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
  case 3: //Constant mean
    usepar=0;
    for(ii=0; ii<conf.constRecBreaks.size(); ++ii){
      if(dat.years(i)>conf.constRecBreaks(ii)){usepar++;}
    }
    predN(0)=par.rec_pars(usepar);
    break;
  case 62: // AR1 (on log-scale)
    predN(0) = par.rec_pars(0) + (2.0 / (1.0 + exp(-par.rec_pars(1))) - 1.0) * (logN(0,i-1) - par.rec_pars(0));
    break;
  case 90: // Non-increasing spline on log R/S
    predN(0) = log(thisSSB) + ibcdspline(log(thisSSB),
					 (vector<Type>)(conf.constRecBreaks.template cast<Type>()),
					 par.rec_pars);
    break;
  case 91: // integrated spline on log R/S
    predN(0) = log(thisSSB) + ibcspline(log(thisSSB),
					(vector<Type>)(conf.constRecBreaks.template cast<Type>()),
					par.rec_pars);
    break;
  case 92: // spline on log R/S
    predN(0) = log(thisSSB) + bcspline(log(thisSSB),
				       (vector<Type>)(conf.constRecBreaks.template cast<Type>()),
				       par.rec_pars);
    break;
  case 93: 			// Type 2 depensatory Non-increasing spline on log R/S
     predN(0) = log(thisSSB) + ibcdspline(log(thisSSB),
					 (vector<Type>)(conf.constRecBreaks.template cast<Type>()),
					  (vector<Type>)par.rec_pars.segment(0,par.rec_pars.size()-1)) +
       log(thisSSB) - logspace_add2(log(thisSSB),par.rec_pars(par.rec_pars.size()-1));
    break;
  default:
    predN(0)=functionalStockRecruitment(thisSSB, par.rec_pars, conf.stockRecruitmentModelCode);
    break;
  }

  switch(conf.logNMeanCorrection(0)){
  case 1:			// Mean on natural scale
    predN(0) -= 0.5 * exp(2.0 * par.logSdLogN(conf.keyVarLogN(0)));
    break;
  case 2:			// Mode on natural scale
    predN(0) += exp(2.0 * par.logSdLogN(conf.keyVarLogN(0)));
    break;
  default:			// Median on natural scale
    predN(0) += 0.0;
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

  for(int j=1; j<stateDimN; ++j){
    switch(conf.logNMeanCorrection(1)){
    case 1:			// Mean on natural scale
      predN(j) -= 0.5 * exp(2.0 * par.logSdLogN(conf.keyVarLogN(j)));
      break;
    case 2:			// Mode on natural scale
      predN(j) += exp(2.0 * par.logSdLogN(conf.keyVarLogN(j)));
      break;
    default:			// Median on natural scale
      predN(j) += 0.0;
      break;    
    }
  }
  
  return predN;  
}

#endif

