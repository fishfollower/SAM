SAM_DEPENDS(logspace)


/* 
Implementation of

Trijoulet et al. (2023) Model validation for compositional data in stock assessment models: Calculating residuals with correct properties. Fish Res. DOI: 10.1016/j.fishres.2022.106487

*/

// The Dirichlet distribution is parameterized to have proportions p (summing to one) 

template<class Type>
vector<Type> rdirichlet(vector<Type> log_p, Type log_s)SOURCE({
    vector<Type> r = rgamma((vector<Type>)exp(log_p + log_s), Type(1.0));
    return r / r.sum();
  })

SAM_SPECIALIZATION(vector<double> rdirichlet(vector<double>, double));



template <class Type>
Type ddirichlet_vtri(vector<Type> x, vector<Type> alpha, int give_log)
  SOURCE({
  Type logB = lgamma(alpha).sum() - lgamma(alpha.sum());
  Type logres=((alpha-Type(1))*log(x)).sum() - logB;
  if(give_log){
    return logres;
  }else{
    return exp(logres);
  }
    })

template<class Type>
Type ddirichlet(vector<Type> log_x, vector<Type> log_p, Type log_s, data_indicator<vector<Type>, Type> keep, int give_log DEFARG(=0))SOURCE({
  SAM_ASSERT(log_x.size() == log_p.size(), "ddirichlet: x and p must have the same length");
  vector<Type> log_alpha = log_p + log_s;
  Type log_as = logspace_sum(log_alpha);
  Type log_1mcx = 0.0;
  Type logLik = 0.0;
  for(int i = 0; i < log_x.size()-1; ++i){ // Last is fixed given previous
    log_as = logspace_sub_SAM(log_as,log_alpha(i));
    logLik += keep(i) * (dbeta(exp(log_x(i) - log_1mcx), exp(log_alpha(i)), exp(log_as), true) - log_1mcx);
    if(i < log_x.size()-2)
      log_1mcx = logspace_sub_SAM(log_1mcx,log_x(i));
  }
  if(give_log)
    return logLik;
  return exp(logLik);
  })

SAM_SPECIALIZATION(double ddirichlet(vector<double>, vector<double>, double, data_indicator<vector<double>, double>, int));
SAM_SPECIALIZATION(TMBad::ad_aug ddirichlet(vector<TMBad::ad_aug>, vector<TMBad::ad_aug>, TMBad::ad_aug, data_indicator<vector<TMBad::ad_aug>, TMBad::ad_aug>, int));


template<class Type>
Type ddirichlet(vector<Type> log_x, vector<Type> log_p, Type log_s, vector<Type> keep, int give_log DEFARG(=0))SOURCE({
  SAM_ASSERT(log_x.size() == log_p.size(), "ddirichlet: x and p must have the same length");
  vector<Type> log_alpha = log_p + log_s;
  Type log_as = logspace_sum(log_alpha);
  Type log_1mcx = 0.0;
  Type logLik = 0.0;
  for(int i = 0; i < log_x.size()-1; ++i){ // Last is fixed given previous
    log_as = logspace_sub_SAM(log_as,log_alpha(i));
    logLik += keep(i) * (dbeta(exp(log_x(i) - log_1mcx), exp(log_alpha(i)), exp(log_as), true) - log_1mcx);
    if(i < log_x.size()-2)
      log_1mcx = logspace_sub_SAM(log_1mcx,log_x(i));
  }
  if(give_log)
    return logLik;
  return exp(logLik);
  })

SAM_SPECIALIZATION(double ddirichlet(vector<double>, vector<double>, double, vector<double>, int));
SAM_SPECIALIZATION(TMBad::ad_aug ddirichlet(vector<TMBad::ad_aug>, vector<TMBad::ad_aug>, TMBad::ad_aug, vector<TMBad::ad_aug>, int));
