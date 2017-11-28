
template <class Type>
vector<Type> predNFun(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, array<Type> &logF, int i){
  int stateDimN=logN.dim[0];
  int timeSteps=logN.dim[1];
  array<Type> totF=totFFun(conf, logF);
  vector<Type> predN(stateDimN); 
  Type thisSSB=Type(0);
    if(conf.stockRecruitmentModelCode==0){ // straight RW 
      predN(0)=logN(0,i-1);
    }else{
      if((i-conf.minAge)>=0){
        thisSSB=ssbi(dat,conf,logN,totF,i-conf.minAge,logF);
      }else{
        thisSSB=ssbi(dat,conf,logN,totF,0,logF); // use first in beginning       
      } 
      if(conf.stockRecruitmentModelCode==1){//ricker
        predN(0)=par.rec_loga(0)+log(thisSSB)-exp(par.rec_logb(0))*thisSSB;
      }else{
        if(conf.stockRecruitmentModelCode==2){//BH
          predN(0)=par.rec_loga(0)+log(thisSSB)-log(1.0+exp(par.rec_logb(0))*thisSSB); 
        }else{
          error("SR model code not recognized");
        }
      }
    }
  
    for(int j=1; j<stateDimN; ++j){
      predN(j)=logN(j-1,i-1)-totF(j-1,i-1)-dat.natMor(i-1,j-1);
    }  
    if(conf.maxAgePlusGroup==1){
      predN(stateDimN-1)=log(exp(logN(stateDimN-2,i-1)-totF(stateDimN-2,i-1)-dat.natMor(i-1,stateDimN-2))+
                             exp(logN(stateDimN-1,i-1)-totF(stateDimN-1,i-1)-dat.natMor(i-1,stateDimN-1))); 
    }
  return predN;  
}
