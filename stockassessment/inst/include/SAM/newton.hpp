#include <memory>

HEADER(
struct NewtonFunctor {
  virtual TMBad::ad_aug operator()(const vector<TMBad::ad_aug>& x);
};
       )

SOURCE(
       TMBad::ad_aug NewtonFunctor::operator()(const vector<TMBad::ad_aug>& x){
	 return R_NaReal;
       };
       )

HEADER(
struct NewtonWrapper {
  std::shared_ptr<NewtonFunctor> ptr;
  NewtonWrapper();
  NewtonWrapper(std::shared_ptr<NewtonFunctor> p);
  ~NewtonWrapper();
  TMBad::ad_aug operator()(const vector<TMBad::ad_aug>& x);
};
       )

SOURCE(
NewtonWrapper::NewtonWrapper() : ptr(nullptr) {};
       )
SOURCE(
       NewtonWrapper::NewtonWrapper(std::shared_ptr<NewtonFunctor> p) : ptr(p) {};
       )
SOURCE(
       NewtonWrapper::~NewtonWrapper() {
	 if(ptr != nullptr)
	   ptr.reset();
       };
       )
SOURCE(
TMBad::ad_aug NewtonWrapper::operator()(const vector<TMBad::ad_aug>& x){
  if(ptr != nullptr)
    return ptr->operator()(x);
  return R_NaReal;
};
       )
  
SAM_SPECIALIZATION(newton::vector<double> newton::Newton(NewtonWrapper&, Eigen::Array<double, Eigen::Dynamic, 1>, newton::newton_config));
SAM_SPECIALIZATION(newton::vector<TMBad::ad_aug> newton::Newton(NewtonWrapper&, Eigen::Array<TMBad::ad_aug, Eigen::Dynamic, 1>, newton::newton_config));

template<class Type>
vector<Type> SAM_Newton(std::shared_ptr<NewtonFunctor>& fn, vector<Type>& start, newton::newton_config& cfg DEFARG(= newton::newton_config())) SOURCE({
  NewtonWrapper nw(fn);
  return newton::Newton(nw, start, cfg);  
  });

SAM_SPECIALIZATION(vector<double> SAM_Newton(std::shared_ptr<NewtonFunctor>&, vector<double>&, newton::newton_config&));
SAM_SPECIALIZATION(vector<TMBad::ad_aug> SAM_Newton(std::shared_ptr<NewtonFunctor>&, vector<TMBad::ad_aug>&, newton::newton_config&));
