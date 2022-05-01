template <class Type>
Type nllBioProcess(array<Type> P, vector<Type> meanVec, vector<int> keyMeanVec, vector<Type> logPhi, Type logSdP){
    int nrow=P.dim[0];
    int ncol=P.dim[1];
    int n=nrow*ncol;
    vector<int> r(n); r.setZero();   
    vector<int> c(n); c.setZero();
    int idx=0;
    for(int j=0; j<ncol; ++j){
      for(int i=0; i<nrow; ++i){
        r(idx)=i;
	c(idx)=j;
	++idx;
      }
    }
    matrix<Type> Wc(n,n); Wc.setZero();
    matrix<Type> Wd(n,n); Wd.setZero();
    for(int i=0; i<n; ++i){
      for(int j=0; j<n; ++j){
	if((c(i)==c(j))&&(abs(r(i)-r(j))==1)){
	  Wc(i,j)=1;
	  Wc(i,i)-=1;
        }   
 	if( (((r(i)-r(j))==1)&&((c(i)-c(j))==1))||(((r(i)-r(j)==(-1)))&&((c(i)-c(j))==(-1))) ){
       	  Wd(i,j)=1;
	  Wd(i,i)-=1;
	}
      }
    }
    
    array<Type> mP(nrow,ncol);
    for(int i=0; i<nrow; ++i){
      for(int j=0; j<ncol; ++j){
	mP(i,j)=meanVec(keyMeanVec(j));
      }
    }
    
    vector<Type> phi=exp(logPhi);
    matrix<Type> I(Wc.rows(),Wc.cols());
    I.setIdentity();
    matrix<Type> Q=I-phi(0)*Wc-phi(1)*Wd;
    
    using namespace density;
    return SCALE(GMRF(asSparseMatrix(Q)),exp(logSdP))((P-mP).vec());
}



template <class Type>
Type nllSW(array<Type> &logSW, dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, objective_function<Type> *of){
  if(conf.stockWeightModel==1){
    Type nll=0;
    array<Type> sw=dat.stockMeanWeight;
    nll += nllBioProcess(logSW, par.meanLogSW, conf.keyStockWeightMean, par.logPhiSW, par.logSdProcLogSW(0));
    for(int i=0; i<sw.dim[0]; ++i){
      for(int j=0; j<sw.dim[1]; ++j){
        if(!isNA(sw(i,j))){
          nll += -dnorm(log(sw(i,j)),logSW(i,j),exp(par.logSdLogSW(conf.keyStockWeightObsVar(j))),true);
        }
	dat.stockMeanWeight(i,j)=exp(logSW(i,j));
      }
    }
    ADREPORT_F(logSW,of);
    return nll;
  }
  return Type(0);
}


template <class Type>
Type nllCW(array<Type> &logCW, dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, objective_function<Type> *of){
  if(conf.catchWeightModel==1){
    Type nll=0;
    array<Type> cw=dat.catchMeanWeight;
    nll += nllBioProcess(logCW, par.meanLogCW, conf.keyCatchWeightMean, par.logPhiCW, par.logSdProcLogCW(0));
    for(int i=0; i<cw.dim[0]; ++i){
      for(int j=0; j<cw.dim[1]; ++j){
        if(!isNA(cw(i,j))){
          nll += -dnorm(log(cw(i,j)),logCW(i,j),exp(par.logSdLogCW(conf.keyCatchWeightObsVar(j))),true);
        }
	dat.catchMeanWeight(i,j)=exp(logCW(i,j));
      }
    }
    ADREPORT_F(logCW,of);
    return nll;
  }
  return Type(0);
}


template<class Type>
Type squash(Type u){
  Type eps = 1.0e-6;
  u = (1.0 - eps) * (u - .5) + .5;
  return u;
}

template <class Type>
Type nllMO(array<Type> &logitMO, dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, objective_function<Type> *of){
  if(conf.matureModel==1){
    Type nll=0;
    array<Type> mo=dat.propMat;
    nll += nllBioProcess(logitMO, par.meanLogitMO, conf.keyMatureMean, par.logPhiMO, par.logSdProcLogitMO(0));
    Type m,a,b,v, prec;
    for(int i=0; i<mo.dim[0]; ++i){
      for(int j=0; j<mo.dim[1]; ++j){
	m = invlogit(logitMO(i,j));
        if(!isNA(mo(i,j))){

 	  //Francisco Cribari-Neto, Achim Zeileis (2010). "Beta Regression in R."
	  prec=invlogit(par.logSdMO(0))*1.0e3;
	  a = m*prec;                 //
          b = (Type(1)-m)*prec;       //v=mu*(1-mu)/(1+precision)
	  //v=m*(1-m)/(1+exp(logPrecision));
	
	  if(mo(i,j)>0.5){
            nll += -dbeta(Type(squash(1-mo(i,j))),b,a,true);
	  }else{
	    nll += -dbeta(Type(squash(mo(i,j))),a,b,true);
	  }
        }
	dat.propMat(i,j)=m;
      }
    }
    ADREPORT_F(logitMO,of);
    return nll;
  }
  return Type(0);
}



template <class Type>
Type nllNM(array<Type> &logNM, dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, objective_function<Type> *of){
  if(conf.mortalityModel==1){
    Type nll=0;
    array<Type> nm=dat.natMor;
    nll += nllBioProcess(logNM, par.meanLogNM, conf.keyMortalityMean, par.logPhiNM, par.logSdProcLogNM(0));
    for(int i=0; i<nm.dim[0]; ++i){
      for(int j=0; j<nm.dim[1]; ++j){
        if(!isNA(nm(i,j))){
          nll += -dnorm(log(nm(i,j)),logNM(i,j),exp(par.logSdLogNM(conf.keyMortalityObsVar(j))),true);
        }
	dat.natMor(i,j)=exp(logNM(i,j));
      }
    }
    ADREPORT_F(logNM,of);
    return nll;
  }
  return Type(0);
}

template <class Type>
Type nllASD(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, objective_function<Type> *of){
  if(par.logAlphaASD.size()>0){
    vector<Type> alpha=exp(par.logAlphaASD);
    vector<Type> beta=exp(par.logBetaASD);
    vector<Type> logAge=par.logAgeASD;
    Type nll=0, m, sd;
    matrix<int> obs=dat.agesampledata.matrix();
    int mg=0, pg=conf.maxAge;//mg=obs(0,0), pg=obs(0,0);
    //for(int i=1; i<obs.rows(); ++i){
    //  if(obs(i,0)<mg){mg=obs(i,0);}
    //  if(obs(i,0)>pg){pg=obs(i,0);}
    //}
    int aa,ss,dd;
    for(int i=0; i<obs.rows(); ++i){
      aa=obs(i,0);
      ss=obs(i,1);
      dd=obs(i,2);
      m=logAge(ss);
      sd=sqrt(log(alpha(0)*exp(m*(beta(0)-2))+1)); // corresponds to v=alpha*mu^beta on natural scale ()
      if(pg==aa){
	nll += -pnorm5(log(Type(aa)),m,sd,Type(0),Type(1));
      }else{
        if(mg==aa){
	  nll += -pnorm5(log(Type(aa)+Type(1.0)),m,sd,Type(1),Type(1));
        }else{
	  nll += -logspace_sub(pnorm5(log(Type(aa)+Type(1.0)),m,sd,Type(1),Type(1)),pnorm5(log(Type(aa)),m,sd,Type(1),Type(1)));
        }
      }
    }
  
    for(int f=0; f<dat.noFleets; ++f){
      if(!isNAINT(conf.keyAgeSampleData(f))){
        matrix<Type> confusion=dat.ageConfusion(f); confusion.setZero();
	if(conf.offsetAgeSampleData(f)>0){
          for(int i=0; i<confusion.rows(); ++i){
            for(int j=0; j<confusion.cols(); ++j){ 
      	      m=log(Type(i)+dat.minAgePerFleet(f)+conf.offsetAgeSampleData(f));
              sd=sqrt(log(alpha(0)*exp(m*(beta(0)-2))+1));
  	      if(pg==(j+dat.minAgePerFleet(f))){
                confusion(i,j)=pnorm5(log(Type(j+dat.minAgePerFleet(f))),m,sd,Type(0),Type(0)); 
  	      }else{
                if(mg==(j+dat.minAgePerFleet(f))){
  	          confusion(i,j)=pnorm5(log(Type(j+dat.minAgePerFleet(f))+Type(1.0)),m,sd,Type(1),Type(0));
  	        }else{
                  confusion(i,j)=pnorm5(log(Type(j+dat.minAgePerFleet(f))+Type(1.0)),m,sd,Type(1),Type(0))-pnorm5(log(Type(j+dat.minAgePerFleet(f))),m,sd,Type(1),Type(0));
  	        }
  	      }
            }
  	    Type rs=(confusion.row(i)).sum();
  	    for(int j=0; j<confusion.cols(); ++j){
              confusion(i,j)=confusion(i,j)/rs; 
  	    }
          }
	}else{
          int nsteps=12;
	  Type dt=Type(1)/Type(nsteps);
	  for(int k=0; k<nsteps; ++k){
            for(int i=0; i<confusion.rows(); ++i){
              for(int j=0; j<confusion.cols(); ++j){ 
      	        m=log(Type(i)+dat.minAgePerFleet(f)+dt*k+0.5*dt);
                sd=sqrt(log(alpha(0)*exp(m*(beta(0)-2))+1));
  	        if(pg==(j+dat.minAgePerFleet(f))){
                  confusion(i,j)+=dt*(pnorm5(log(Type(j+dat.minAgePerFleet(f))),m,sd,Type(0),Type(0))); 
  	        }else{
                  if(mg==(j+dat.minAgePerFleet(f))){
  	            confusion(i,j)+=dt*(pnorm5(log(Type(j+dat.minAgePerFleet(f))+Type(1.0)),m,sd,Type(1),Type(0)));
  	          }else{
                    confusion(i,j)+=dt*(pnorm5(log(Type(j+dat.minAgePerFleet(f))+Type(1.0)),m,sd,Type(1),Type(0))-
					pnorm5(log(Type(j+dat.minAgePerFleet(f))),m,sd,Type(1),Type(0)));
  	          }
  	        }
              }
	    }
	  }
	}
	for(int i=0; i<confusion.rows(); ++i){
          for(int j=0; j<confusion.cols(); ++j){
	    Type rs=(confusion.row(i)).sum();
  	    for(int j=0; j<confusion.cols(); ++j){
              confusion(i,j)=confusion(i,j)/rs; 
  	    }
	  }
	}
	std::cout<<confusion<<std::endl;
        dat.ageConfusion(f)=confusion;
      }
    }
    //ADREPORT_F(logNM,of);
    return nll;
  }else{
    return Type(0);
  }
}
