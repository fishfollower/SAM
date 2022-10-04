SAM_DEPENDS(define)

template<class Type>
Type nllSplinePenalty(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, objective_function<Type> *of)SOURCE({
  Type ans = 0.0;
  if(CppAD::Variable(par.splinePenalty)){
    if(conf.stockRecruitmentModelCode == 90 ||
       conf.stockRecruitmentModelCode == 91 ||
       conf.stockRecruitmentModelCode == 92){
      ans -= dnorm(par.rec_pars(0), Type(0.0), exp(par.splinePenalty), true);
      for(int i = 1; i < par.rec_pars.size() - 1; ++i)
	ans -= dnorm(par.rec_pars(i), par.rec_pars(i-1), exp(par.splinePenalty), true);
      ans -= dnorm(par.rec_pars(par.rec_pars.size()-1), Type(0.0), Type(10.0), true);
    }
  }
  vector<Type> rec_pars = par.rec_pars;
  ADREPORT_F(rec_pars, of);
  return ans;
  });

SAM_SPECIALIZATION(double nllSplinePenalty(dataSet<double>&, confSet&, paraSet<double>&, objective_function<double>*));
SAM_SPECIALIZATION(TMBad::ad_aug nllSplinePenalty(dataSet<TMBad::ad_aug>&, confSet&, paraSet<TMBad::ad_aug>&, objective_function<TMBad::ad_aug>*));

