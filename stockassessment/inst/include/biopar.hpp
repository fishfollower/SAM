template <class Type>
Type nllSW(array<Type> &logSW, dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, objective_function<Type> *of){
  if(conf.stockWeightModel==1){
    Type nll=0;
    array<Type> sw=dat.stockMeanWeight;
    int nrow=sw.dim[0];
    int ncol=sw.dim[1];
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
    
    array<Type> mLogSW(nrow,ncol);
    for(int i=0; i<nrow; ++i){
      for(int j=0; j<ncol; ++j){
	mLogSW(i,j)=par.meanLogSW(conf.keyStockWeightMean(j));
      }
    }
    
    vector<Type> phi=exp(par.logPhiSW);
    matrix<Type> I(Wc.rows(),Wc.cols());
    I.setIdentity();
    matrix<Type> Q=I-phi(0)*Wc-phi(1)*Wd;
    
    using namespace density;
    nll += SCALE(GMRF(asSparseMatrix(Q)),exp(par.logSdProcLogSW(0)))((logSW-mLogSW).vec());
        
    for(int i=0; i<nrow; ++i){
      for(int j=0; j<ncol; ++j){
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
    int nrow=mo.dim[0];
    int ncol=mo.dim[1];
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
    
    array<Type> mLogitMO(nrow,ncol);
    for(int i=0; i<nrow; ++i){
      for(int j=0; j<ncol; ++j){
	mLogitMO(i,j)=par.meanLogitMO(conf.keyMatureMean(j));
      }
    }
    
    vector<Type> phi=exp(par.logPhiMO);
    matrix<Type> I(Wc.rows(),Wc.cols());
    I.setIdentity();
    matrix<Type> Q=I-phi(0)*Wc-phi(1)*Wd;
    
    using namespace density;
    nll += SCALE(GMRF(asSparseMatrix(Q)),exp(par.logSdProcLogitMO(0)))((logitMO-mLogitMO).vec());

    Type m,a,b,v, prec;
    for(int i=0; i<nrow; ++i){
      for(int j=0; j<ncol; ++j){
        if(!isNA(mo(i,j))){
	  m = invlogit(logitMO(i,j));
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
        	  
          //nll += -dnorm(log(sw(i,j)),logSW(i,j),exp(par.logSdLogSW(conf.keyStockWeightObsVar(j))),true);
        }
	dat.propMat(i,j)=m;
      }
    }
    ADREPORT_F(logitMO,of);
    return nll;
  }
  return Type(0);
}
