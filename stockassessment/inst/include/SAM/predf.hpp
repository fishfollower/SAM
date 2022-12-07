SAM_DEPENDS(convenience)
SAM_DEPENDS(define)


template<class Type>
vector<Type> get_fmu(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logF)SOURCE({    
  int stateDimF=logF.dim[0];

  vector<Type> mu(logF.dim(0));
  mu.setZero();

  if(par.muF.size() == 0)
    return mu;

  
  int stateDimN=conf.keyLogFsta.dim[1];
  int noFleets=conf.keyLogFsta.dim[0];
 vector<Type> statesFleets(stateDimF);
   //Fill statesFleets:
  for(int f=0; f<noFleets;f++){
    for(int i=0;i<stateDimF;i++){
      for(int j=0;j<stateDimN;j++){
	if(conf.keyLogFsta(f,j)==i){
	  statesFleets(i)=f;
	  break;
	}  
      }
    }
  }
  

  vector<bool> done(mu.size());
  done.setConstant(false);
  
  for(int f=0;f<noFleets;f++){
    for(int a = 0; a < conf.keyLogFsta.dim(1); ++a){
      int Findx = conf.keyLogFsta(f,a);
      int muIndx = conf.keyLogFmu(f,a);
      if(Findx > (-1) && !done(Findx)){
	if(muIndx > (-1)){
	  mu(Findx) = par.muF(muIndx);
	}else{
	  mu(Findx) = 0.0;
	}
	done(Findx) = true;
      }
    }
  }
  return mu;
  })


SAM_SPECIALIZATION(vector<double> get_fmu(dataSet<double>&, confSet&, paraSet<double>&, array<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> get_fmu(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&));


template<class Type>
vector<Type> get_frho(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logF)SOURCE({
  int stateDimF=logF.dim[0];

  vector<Type> rho(logF.dim(0));
  rho.setConstant(1.0);

  if(par.trans_rho_F.size() == 0)
    return rho;

  
  int noFleets=conf.keyLogFsta.dim[0];
  int stateDimN=conf.keyLogFsta.dim[1];
 vector<Type> statesFleets(stateDimF);
   //Fill statesFleets:
  for(int f=0; f<noFleets;f++){
    for(int i=0;i<stateDimF;i++){
      for(int j=0;j<stateDimN;j++){
	if(conf.keyLogFsta(f,j)==i){
	  statesFleets(i)=f;
	  break;
	}  
      }
    }
  }
  
  vector<bool> done(rho.size());
  done.setConstant(false);
  
  for(int f=0;f<noFleets;f++){
    for(int a = 0; a < conf.keyLogFsta.dim(1); ++a){
      int Findx = conf.keyLogFsta(f,a);
      int muIndx = conf.keyLogFrho(f,a);
      if(Findx > (-1) && !done(Findx)){
	if(muIndx > (-1)){
	  rho(Findx) = toInterval((Type)par.trans_rho_F(muIndx),Type(0.0),Type(1.0),Type(1.0));
	}else{
	  rho(Findx) = 1.0;
	}
	done(Findx) = true;
      }
    }
  }
  return rho;
  })



SAM_SPECIALIZATION(vector<double> get_frho(dataSet<double>&, confSet&, paraSet<double>&, array<double>&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> get_frho(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&));

  


template<class Type>
matrix<Type> get_fvar(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logF)SOURCE({
  // using CppAD::abs;
  int stateDimF=logF.dim[0];
  int timeSteps=logF.dim[1];
  int noFleets=conf.keyLogFsta.dim[0];
  int stateDimN=conf.keyLogFsta.dim[1];
  vector<Type> sdLogFsta=exp(par.logSdLogFsta);
  array<Type> resF(stateDimF,timeSteps-1);
  matrix<Type> fvar(stateDimF,stateDimF);
  matrix<Type> fcor(stateDimF,stateDimF);
  vector<Type> fsd(stateDimF);  
  vector<Type> statesFleets(stateDimF);
  
  //Fill statesFleets:
  for(int f=0; f<noFleets;f++){
    for(int i=0;i<stateDimF;i++){
      for(int j=0;j<stateDimN;j++){
	if(conf.keyLogFsta(f,j)==i){
	  statesFleets(i)=f;
	  break;
	}  
      }
    }
  }
  
  fcor.setZero();
  // for(int i=0; i<stateDimF; ++i){
  //   fcor(i,i)=1.0;
  // }
  
  int count=0; //if corFlag varies between 0-2, itrans_rho is shorter than comm fleet length
  for(int f=0;f<noFleets;f++){
    bool nxtPar = false;
    for(int i=0; i<stateDimF; ++i){
      fcor(i,i)=1.0;
      for(int j=0; j<i; ++j){
        if(statesFleets(i)==f && statesFleets(j)==f){
	  switch(conf.corFlag(f)){
	  case 0:		// Independent
	    fcor(j,i) = 0.0;
	    break;
	  case 1:		// Compound symmetry
	    fcor(j,i)=toInterval((Type)par.itrans_rho(count),Type(-1.0),Type(1.0),Type(2.0));
	    nxtPar = true;
	    break;
	  case 2:		// AR(1) structure
	    fcor(j,i)=pow(toInterval((Type)par.itrans_rho(count),Type(-1.0),Type(1.0),Type(2.0)),abs(i-j));
	    nxtPar = true;
	    break;
	    // case 3: separable structure
	  case 4:		// (almost) Perfect correlation
            fcor(j,i)= 0.99;
	    break;
	  default:
	    Rf_error("F correlation not implemented");
	    break;
	  }
	  fcor(i,j) = fcor(j,i);
        }
      }
    }
    if(nxtPar)
      count++;
  }

  for(int i=0; i<stateDimF; ++i){
    bool stop = false;
    int ff = 0;
    int j = 0;
    for(ff=0; ff<noFleets; ff++){
      for(j=0; j<stateDimN; j++){
        if(conf.keyLogFsta(ff,j)==i){
          stop=true;
          break;
        } 
      }
      if(stop)break;
    }
    fsd(i)=sdLogFsta(conf.keyVarF(ff,j));
  }
 
  for(int i=0; i<stateDimF; ++i){
    for(int j=0; j<stateDimF; ++j){
      fvar(i,j)=fsd(i)*fsd(j)*fcor(i,j);
    }
  }
  return fvar;
  })

SAM_SPECIALIZATION(matrix<double> get_fvar(dataSet<double>&, confSet&, paraSet<double>&, array<double>&));
SAM_SPECIALIZATION(matrix<TMBad::ad_aug> get_fvar(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, array<TMBad::ad_aug>&));

  
