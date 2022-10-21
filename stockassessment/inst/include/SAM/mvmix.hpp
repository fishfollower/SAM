SAM_DEPENDS(convenience)


template <class Type>
matrix<Type> diagonalMatrix(Type v, int n) SOURCE({
  matrix<Type> r(n,n);
  r.setZero();
  for(int i = 0; i < n; ++i)
    r(i,i) = v;
  return r;
  })

SAM_SPECIALIZATION(matrix<double> diagonalMatrix(double, int));
SAM_SPECIALIZATION(matrix<TMBad::ad_aug> diagonalMatrix(TMBad::ad_aug, int));


HEADER(
template <class Type>
class MVMIX_t{
  Type halfLogDetS;         
  vector<Type> p1;                  /*fraction t3*/
  matrix<Type> Sigma;       
  vector<Type> sd;
  matrix<Type> L_Sigma;
  matrix<Type> inv_L_Sigma;
public:
  MVMIX_t();
  MVMIX_t(matrix<Type> Sigma_, Type p1_);
  MVMIX_t(matrix<Type> Sigma_, vector<Type> p1_);
  MVMIX_t(matrix<Type> Sigma_, Type p1_, bool useAtomic);
  matrix<Type> cov();
  void setSigma(matrix<Type> Sigma_, bool useAtomic = true);
  void setSigma(matrix<Type> Sigma_, Type p1_);
  /** \brief Evaluate the negative log density */
  Type operator()(vector<Type> x);
  Type operator()(vector<Type> x, vector<Type> keep);

  vector<Type> simulate();
});

SOURCE(
	 template<class Type>
	 MVMIX_t<Type>::MVMIX_t(){}
	 );

SOURCE(
	 template<class Type>
	 MVMIX_t<Type>::MVMIX_t(matrix<Type> Sigma_, Type p1_){
	   setSigma(Sigma_);
	   p1=vector<Type>(Sigma_.rows());
	   p1.setConstant(p1_);
	 }
	 );

SOURCE(
	 template<class Type>
	 MVMIX_t<Type>::MVMIX_t(matrix<Type> Sigma_, vector<Type> p1_){
	   setSigma(Sigma_);
	   p1=p1_;
	 }
	 );

SOURCE(
	 template<class Type>
	 MVMIX_t<Type>::MVMIX_t(matrix<Type> Sigma_, Type p1_, bool useAtomic){
	   setSigma(Sigma_, useAtomic);
	   p1=p1_;
	 }
	 );

SOURCE(
	 template<class Type>
	 matrix<Type> MVMIX_t<Type>::cov(){return Sigma;}
	 );

SOURCE(
	 template<class Type>
	 void MVMIX_t<Type>::setSigma(matrix<Type> Sigma_, bool useAtomic){
	   Sigma = Sigma_;
	   p1 = vector<Type>(Sigma_.rows());
	   p1.setZero();
	   sd = sqrt(vector<Type>(Sigma.diagonal()));
	   Eigen::LLT<Eigen::Matrix<Type,Eigen::Dynamic,Eigen::Dynamic> > llt(Sigma);
	   L_Sigma = llt.matrixL();
	   vector<Type> D=L_Sigma.diagonal();
	   halfLogDetS = sum(log(D));
	   // if(useAtomic){
	   //inv_L_Sigma = atomic::matinv(L_Sigma);
	   // }else{
	   inv_L_Sigma = L_Sigma.inverse();
	   // }
	 }
	 );

SOURCE(
	 template<class Type>
	 void MVMIX_t<Type>::setSigma(matrix<Type> Sigma_, Type p1_){
	   setSigma(Sigma_);
	   p1.setConstant(p1_);
	 }
	 )

SOURCE(
	 template<class Type>
	 Type MVMIX_t<Type>::operator()(vector<Type> x){
	   vector<Type> z = inv_L_Sigma*x;
	   return -sum(logdrobust(z,p1))+halfLogDetS;
	 }
	 )

SOURCE(
	 template<class Type>
	 Type MVMIX_t<Type>::operator()(vector<Type> x, vector<Type> keep){
	   matrix<Type> S = Sigma;
	   vector<Type> p = p1; 
	   vector<Type> not_keep = Type(1.0) - keep;
	   for(int i = 0; i < S.rows(); i++){
	     for(int j = 0; j < S.cols(); j++){
	       S(i,j) = S(i,j) * keep(i) * keep(j);
	     }
	     S(i,i) += not_keep(i) * pow((Type(1)-p(i))*sqrt(Type(0.5)/M_PI)+p(i)*(Type(2)/(M_PI*sqrt(Type(3)))),2);
	   }
	   return MVMIX_t<Type>(S,p)(x * keep);
	 }
	 )

SOURCE(
	 template<class Type>
	 vector<Type> MVMIX_t<Type>::simulate(){
	   int siz = Sigma.rows();
	   vector<Type> x(siz);
	   for(int i=0; i<siz; ++i){
	     Type u = runif(0.0,1.0);
	     if(u<p1(i)){
	       x(i) = rt(3.0);
	     }else{
	       x(i) = rnorm(0.0,1.0);
	     }
	   }
	   x = L_Sigma*x;
	   return x;
	 }
	 )

SAM_SPECIALIZATION(struct MVMIX_t<double>);
SAM_SPECIALIZATION(struct MVMIX_t<TMBad::ad_aug>);

SAM_SPECIALIZATION(struct tmbutils::vector<MVMIX_t<double> >);
SAM_SPECIALIZATION(struct tmbutils::vector<MVMIX_t<TMBad::ad_aug> >);

template <class Type>
MVMIX_t<Type> MVMIX(matrix<Type> Sigma, Type p1)SOURCE({
  return MVMIX_t<Type>(Sigma,p1);
  });

SAM_SPECIALIZATION(MVMIX_t<double> MVMIX(matrix<double>,double));
SAM_SPECIALIZATION(MVMIX_t<TMBad::ad_aug> MVMIX(matrix<TMBad::ad_aug>,TMBad::ad_aug));

