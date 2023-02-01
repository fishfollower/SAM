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



// From https://github.com/vtrijoulet/OSA_multivariate_dists/blob/main/distr.hpp
template<class Type>
vector<int> order_keep(vector<Type> k)SOURCE({
  int n=k.size();
  vector<int> o(n);
  o.setZero();
  int at=-1;
  for(int i=0; i<n;++i){
    if(k(i)>0.5){o(++at) = i;}  
  }
  at=n;  
  for(int i=n-1; i>=0;--i){
    if(k(i)<0.5){o(--at) = i;}  
  }
  return o;
  })

SAM_SPECIALIZATION(vector<int> order_keep(vector<double>));
SAM_SPECIALIZATION(vector<int> order_keep(vector<TMBad::ad_aug>));
  

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

SAM_SPECIALIZATION(double ddirichlet_vtri(vector<double>, vector<double>,  int));
SAM_SPECIALIZATION(TMBad::ad_aug ddirichlet_vtri(vector<TMBad::ad_aug>, vector<TMBad::ad_aug>, int));

  
template<class Type>
Type ddirichlet(vector<Type> log_x, vector<Type> log_p, Type log_s, data_indicator<vector<Type>, Type> keep, int give_log DEFARG(=0))SOURCE({
  SAM_ASSERT(log_x.size() == log_p.size(), "ddirichlet: x and p must have the same length");
  vector<Type> l=keep.cdf_lower;
  vector<Type> u=keep.cdf_upper;
  vector<int> o = order_keep(keep);
  vector<Type> log_alpha = log_p + log_s;
  log_x = log_x(o); log_alpha = log_alpha(o); keep = keep(o); l = l(o); u = u(o);
  Type log_as = logspace_sum(log_alpha);
  Type log_1mcx = 0.0;
  Type logLik = 0.0;
  if((keep.sum() < 1e-8) && (l.sum() < 1e-8) && (u.sum()  < 1e-8)){
    if(give_log)
      return 0.0;
    return 1.0;
  }
  for(int i = 0; i < log_x.size()-1; ++i){ // Last is fixed given previous
    log_as = logspace_sub_SAM(log_as,log_alpha(i));
    logLik += keep(i) * (dbeta(exp(log_x(i) - log_1mcx), exp(log_alpha(i)), exp(log_as), true) - log_1mcx);
    Type cdf = squeeze(pbeta(exp(log_x(i) - log_1mcx), exp(log_alpha(i)), exp(log_as)));
    logLik += l(i) * log(cdf);
    logLik += u(i) * log(1.0 - cdf);
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
Type ddirichlet(vector<Type> log_x, vector<Type> log_p, Type log_s, vector<Type> keep, vector<Type> low_cdf, vector<Type> high_cdf, vector<int> order, int give_log DEFARG(=0))SOURCE({
  SAM_ASSERT(log_x.size() == log_p.size(), "ddirichlet: x and p must have the same length");
  vector<Type> log_alpha = log_p + log_s;
  log_x = log_x(order); log_alpha = log_alpha(order); keep = keep(order); low_cdf = low_cdf(order); high_cdf = high_cdf(order);
  Type log_as = logspace_sum(log_alpha);
  Type log_1mcx = 0.0;
  Type logLik = 0.0;
  for(int i = 0; i < log_x.size()-1; ++i){ // Last is fixed given previous
    log_as = logspace_sub_SAM(log_as,log_alpha(i));
    logLik += keep(i) * (dbeta(exp(log_x(i) - log_1mcx), exp(log_alpha(i)), exp(log_as), true) - log_1mcx);
    Type cdf = squeeze(pbeta(exp(log_x(i) - log_1mcx), exp(log_alpha(i)), exp(log_as)));
    logLik += low_cdf(i) * log(cdf);
    logLik += high_cdf(i) * log(1.0 - cdf);
    if(i < log_x.size()-2)
      log_1mcx = logspace_sub_SAM(log_1mcx,log_x(i));
  }
  if(give_log)
    return logLik;
  return exp(logLik);
  })

SAM_SPECIALIZATION(double ddirichlet(vector<double>, vector<double>, double, vector<double>,vector<double>,vector<double>,vector<int>, int));
SAM_SPECIALIZATION(TMBad::ad_aug ddirichlet(vector<TMBad::ad_aug>, vector<TMBad::ad_aug>, TMBad::ad_aug, vector<TMBad::ad_aug>,vector<TMBad::ad_aug>,vector<TMBad::ad_aug>,vector<int>, int));


template <class Type>
Type ddirichlet_osa(vector<Type>& x, vector<Type>& alpha, data_indicator<vector<Type>, Type>& keep, int give_log DEFARG(=0))
  SOURCE({
  vector<Type> k=keep;
  vector<Type> l=keep.cdf_lower;
  vector<Type> h=keep.cdf_upper;
  vector<int> o=order_keep(k);
  x=x(o); alpha=alpha(o); k=k(o); l=l(o); h=h(o);
  
  int n = alpha.size();
  Type cdf;
  Type sx = 1; // was: x.sum();
  Type sa = alpha.sum();
  sa -= alpha(0);
  Type logres=k(0)*dbeta(x(0),alpha(0),sa,true);
  cdf = pbeta(x(0),alpha(0),sa);
  cdf = squeeze(cdf);
  logres += l(0) * log( cdf );       
  logres += h(0) * log( 1.0 - cdf ); 
  
  for(int i=1; i<(n-1); ++i){
    sx -= x(i-1);
    sa -= alpha(i);
    logres += k(i)*(dbeta(x(i)/sx,alpha(i),sa,true)-log(sx));
    cdf = pbeta(x(i)/sx,alpha(i),sa);
    cdf = squeeze(cdf);
    logres += l(i) * log( cdf );       
    logres += h(i) * log( 1.0 - cdf ); 
  }
  logres += k(n-1)*Type(0);
  cdf=Type(1);
  cdf = squeeze(cdf);
  logres += l(n-1) * log( cdf );       
  logres += h(n-1) * log( 1.0 - cdf ); 
  
  if(give_log){
    return logres;
  }else{
    return exp(logres);
  }
}
)


SAM_SPECIALIZATION(double ddirichlet_osa(vector<double>&, vector<double>&, data_indicator<vector<double>, double>&, int));
SAM_SPECIALIZATION(TMBad::ad_aug ddirichlet_osa(vector<TMBad::ad_aug>&, vector<TMBad::ad_aug>&, data_indicator<vector<TMBad::ad_aug>, TMBad::ad_aug>&, int));
