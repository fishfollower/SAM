// Copyright (C) 2016 Kasper Kristensen
// License: GPL-2

/*
 * Modified to allow unknown number of parameters at compile time
 * 2019 Christoffer Moesgaard Albertsen
 */


#ifndef TINY_ADFLEX_H
#define TINY_ADFLEX_H

/* Standalone ? */
#ifndef R_RCONFIG_H
#include <cmath>
#include <iostream>
#define CSKIP(x) x
#endif

template<class Type>
vector<Type> resize(vector<Type> x, int n){
  if(n == x.size())
    return x;
  vector<Type> r(n);
  r.setZero();
  if(x.size() == 0)
    return r;
  int nuse = std::min((int)x.size(), n);
  r.segment(0,nuse) = x.segment(0,nuse);
  return r;
}


template<class Type>
struct adflex {
  typedef vector<Type> Vector;
  Type value;
  Vector deriv;
  adflex(){}
  adflex(Type v, Vector d){value = v; deriv = d;}
  adflex(double v)        {value = v; deriv.setZero();}

  
  adflex operator+ (const adflex &other) const{
    int newn = std::max((int)deriv.size(),(int)other.deriv.size());
    Vector d1 = resize(deriv,newn);
    Vector d2 = resize(other.deriv,newn);
    d1 += d2;
     return adflex(value + other.value,
		  d1);
  }
  adflex operator+ () const{
    return *this;
  }
  adflex operator- (const adflex &other) const{
    int newn = std::max((int)deriv.size(),(int)other.deriv.size());
    Vector d1 = resize(deriv,newn);
    Vector d2 = resize(other.deriv,newn);
    d1 -= d2;
    return adflex(value - other.value,
		  d1);
  }
  adflex operator- () const{
    return adflex(-value, -deriv);
  }
  adflex operator* (const adflex &other) const{
    int newn = std::max((int)deriv.size(),(int)other.deriv.size());
    Vector d1 = value * resize(other.deriv,newn);
    Vector d2 = resize(deriv,newn) * other.value;
    d1 += d2;
    return adflex(value * other.value,
		  d1);
  }
  adflex operator/ (const adflex &other) const{
    int newn = std::max((int)deriv.size(),(int)other.deriv.size());
    Type res = value / other.value;
    Vector d1 = resize(deriv,newn);
    Vector d2 = res * resize(other.deriv,newn);
    d1 -= d2;
 
    return adflex(res,
		  d1 / other.value );
  }
  /* Comparison operators */
#define COMPARISON_OPERATOR(OP)			\
  template<class other>				\
  bool operator OP (const other &x) const{	\
    return (value OP x);			\
  }
  COMPARISON_OPERATOR(<)
  COMPARISON_OPERATOR(>)
  COMPARISON_OPERATOR(<=)
  COMPARISON_OPERATOR(>=)
  COMPARISON_OPERATOR(==)
  COMPARISON_OPERATOR(!=)
#undef COMPARISON_OPERATOR
  /* Combine ad with other types (constants) */
  adflex operator+ (const double &x) const{
    return adflex(value + x, deriv);
  }
  adflex operator- (const double &x) const{
    return adflex(value - x, deriv);
  }
  adflex operator* (const double &x) const{
    return adflex(value * x, deriv * x);
  }
  adflex operator/ (const double &x) const{
    return adflex(value / x, deriv / x);
  }
  /* Note: 'this' and 'other' may point to the same object */
  adflex& operator+=(const adflex &other){
    value += other.value;
    int newn = std::max((int)deriv.size(),(int)other.deriv.size());
    deriv = resize(deriv,newn);
    deriv += resize(other.deriv,newn);
    return *this;
  }
  adflex& operator-=(const adflex &other){
    value -= other.value;
    int newn = std::max((int)deriv.size(),(int)other.deriv.size());
    deriv = resize(deriv,newn);
    deriv -= resize(other.deriv,newn);
    return *this;
  }
  adflex& operator*=(const adflex &other){
    if (this != &other) {
      int newn = std::max((int)deriv.size(),(int)other.deriv.size());
      deriv = resize(deriv,newn);
      deriv *= other.value;
      deriv += resize(other.deriv,newn) * value;
      value *= other.value;
    } else {
      deriv *= value * 2.;
      value *= value;
    }
    return *this;
  }
  adflex& operator/=(const adflex &other){
    int newn = std::max((int)deriv.size(),(int)other.deriv.size());
    deriv = resize(deriv,newn);
    value /= other.value;
    deriv -= resize(other.deriv,newn) * value;
    deriv /= other.value;
    return *this;
  }
};
/* Binary operators where a constant is first argument */
template<class T>
adflex<T> operator+ (const double &x, const adflex<T> &y) {
  return y + x;
}
template<class T>
adflex<T> operator- (const double &x, const adflex<T> &y) {
  return -(y - x);
}
template<class T>
adflex<T> operator* (const double &x, const adflex<T> &y) {
  return y * x;
}
template<class T>
adflex<T> operator/ (const double &x, const adflex<T> &y) {
  T value = x / y.value;
  return adflex<T>(value, T(-value / y.value) * y.deriv);
}
/* Unary operators with trivial derivatives */
#define UNARY_MATH_ZERO_DERIV(F)		\
  template<class T>			\
  double F (const adflex<T> &x){		\
    return F(x.value);				\
  }
using ::floor; using ::ceil;
using ::trunc; using ::round;
UNARY_MATH_ZERO_DERIV(floor)
UNARY_MATH_ZERO_DERIV(ceil)
UNARY_MATH_ZERO_DERIV(trunc)
UNARY_MATH_ZERO_DERIV(round)
template<class T>
double sign(const T &x){return (x > 0) - (x < 0);}
bool isfinite(const double &x)CSKIP( {return std::isfinite(x);} )
  template<class T, class V>
  bool isfinite(const adflex<T> &x){return isfinite(x.value);}
#undef UNARY_MATH_ZERO_DERIV
/* Unary operators with non-trivial derivatives */
#define UNARY_MATH_DERIVATIVE(F,DF)			\
  template<class T>				\
  adflex<T> F (const adflex<T> &x){		\
    return adflex<T>(F (x.value),			\
			T(DF(x.value)) * x.deriv);	\
  }
using ::exp;  using ::log;
using ::sin;  using ::cos;  using ::tan;
using ::sinh; using ::cosh; using ::tanh;
using ::sqrt; using ::fabs;
template<class T> T D_tan(const T &x) {
  T y = cos(x); return 1. / (y * y);
}
template<class T> T D_tanh(const T &x) {
  T y = cosh(x); return 1. / (y * y);
}
UNARY_MATH_DERIVATIVE(exp, exp)
UNARY_MATH_DERIVATIVE(log, 1.0/)
UNARY_MATH_DERIVATIVE(sin, cos)
UNARY_MATH_DERIVATIVE(cos, -sin)
UNARY_MATH_DERIVATIVE(tan, D_tan)
UNARY_MATH_DERIVATIVE(sinh, cosh)
UNARY_MATH_DERIVATIVE(cosh, sinh)
UNARY_MATH_DERIVATIVE(tanh, D_tanh)
UNARY_MATH_DERIVATIVE(sqrt, 0.5/sqrt)
UNARY_MATH_DERIVATIVE(fabs, sign)
using ::expm1; using ::log1p;
UNARY_MATH_DERIVATIVE(expm1, exp)
template<class T> T D_log1p(const T &x) {return 1. / (x + 1.);}
UNARY_MATH_DERIVATIVE(log1p, D_log1p)
/* asin, acos, atan */
using ::asin; using ::acos; using ::atan;
template<class T> T D_asin(const T &x) {
  return 1. / sqrt(1. - x * x);
}
template<class T> T D_acos(const T &x) {
  return -1. / sqrt(1. - x * x);
}
template<class T> T D_atan(const T &x) {
  return 1. / (1. + x * x);
}
UNARY_MATH_DERIVATIVE(asin, D_asin)
UNARY_MATH_DERIVATIVE(acos, D_acos)
UNARY_MATH_DERIVATIVE(atan, D_atan)
#undef UNARY_MATH_DERIVATIVE
/* A few more ... */
template<class T>
adflex<T> pow (const adflex<T> &x, const adflex<T> &y){
  return exp(y * log(x));
}
using ::pow;
template<class T>
adflex<T> pow (const adflex<T> &x, const double &y){
  return adflex<T> (pow(x.value, y), // Note: x.value could be 0
		       T( y * pow(x.value, y - 1.) ) * x.deriv);
}
/* Comparison operators where a constant is first argument */
#define COMPARISON_OPERATOR_FLIP(OP1, OP2)			\
  template<class T>					\
  bool operator OP1 (const double &x, const adflex<T> &y) {	\
    return y OP2 x;						\
  }
COMPARISON_OPERATOR_FLIP(<,>)
COMPARISON_OPERATOR_FLIP(<=,>=)
COMPARISON_OPERATOR_FLIP(>,<)
COMPARISON_OPERATOR_FLIP(>=,<=)
COMPARISON_OPERATOR_FLIP(==,==)
COMPARISON_OPERATOR_FLIP(!=,!=)
#undef COMPARISON_OPERATOR_FLIP
/* Utility: Return the value of a tiny_adflex type */
//double asDouble(double x) CSKIP( {return x;} )
template<class T>
double asDouble (const adflex<T> &x){
  return asDouble(x.value);
}
/* Utility: Return the max absolute value of all members of a
   tiny_adflex type */
double max_fabs(double x) CSKIP( {return fabs(x);} )
  template<class T>
  double max_fabs (const adflex<T> &x){
  double ans = max_fabs(x.value);
  for(int i=0; i<x.deriv.size(); i++) {
    double tmp = max_fabs(x.deriv[i]);
    ans = (tmp > ans ? tmp : ans);
  }
  return ans;
}
/* R-specific derivatives (rely on Rmath)*/
#ifdef R_RCONFIG_H
extern "C" {
  /* See 'R-API: entry points to C-code' (Writing R-extensions) */
  double	Rf_lgammafn(double);
  double	Rf_psigamma(double, double);
}
template<int deriv>
double lgamma(const double &x) {
  return Rf_psigamma(x, deriv-1);
}
template<>
double lgamma<0>(const double &x) CSKIP( {return Rf_lgammafn(x);} )
  double lgamma(const double &x) CSKIP( {return lgamma<0>(x);} )
  template<int deriv, class T>
  adflex<T> lgamma (const adflex<T> &x){
  return adflex<T> (lgamma< deriv >(x.value),
		       T(lgamma< deriv + 1 >(x.value)) * x.deriv);
}
template<class T>
adflex<T> lgamma (const adflex<T> &x){
  return lgamma<0>(x);
}
#endif
/* Print method */
template<class T>
std::ostream &operator<<(std::ostream &os, const adflex<T> &x) {
  os << "{";
  os << " value=" << x.value;
  os << " deriv=" << x.deriv;
  os << "}";
  return os;
}

/* Interface to higher order derivatives. Example:

   typedef tiny_ad::variable<3, 2> Float; // Track 3rd order derivs wrt. 2 parameters
   Float a (1.23, 0);                     // Let a = 1.23 have parameter index 0
   Float b (2.34, 1);                     // Let b = 2.34 have parameter index 1
   Float y = sin(a + b);                  // Run the algorithm
   y.getDeriv();                          // Get all 3rd order derivatives
*/
#define VARIABLEFLEX(order, scalartype) variableflex<order, scalartype>
template<int order, class Double=double>
struct variableflex : adflex< VARIABLEFLEX(order-1, Double)> {
  typedef adflex< VARIABLEFLEX(order-1, Double) > Base;
			  typedef variableflex<order-1, Double> Type;
			  variableflex() { /* Do not zero-initialize */ }
			  variableflex(Base x) : Base(x) {}
			  variableflex(double x) : Base(x) {}
			  variableflex(double x, int id) : Base(x) {
			    setid(id);
			  }
			  template<class Constant>
			  variableflex(Constant x) {
			    Base::value = x; //Base::deriv.setZero();
			  }
			  template<class Constant>
			  variableflex(Constant x, int id) {
			    Base::value = x;
			    Base::deriv = resize(Base::deriv, std::max((int)Base::deriv.size(),id));
			    setid(id);
			  }
			  void setid(int i0, int count = 0){
			    Base::deriv = resize(Base::deriv, std::max((int)Base::deriv.size(),i0+1));
			    this->value.setid(i0, count);
			    this->deriv[i0].setid(i0, count + 1);
			  }
			  vector<Double> getDeriv(){
			    int nvar = Base::deriv.size();
			    int result_size = pow(nvar,order);
			    int stride = result_size / nvar;
			    vector<Double> ans(result_size);
			    ans.setZero();
			    for(int i=0; i<nvar; i++){
			      ans.segment(i * stride, stride) = resize(this->deriv[i].getDeriv(),stride);
			    }
			    return ans;
			  }
			  };
#undef VARIABLE
  template<class Double>
  struct variableflex<1, Double> : adflex<Double>{
    typedef adflex<Double> Base;
    variableflex<1, Double>() { /* Do not zero-initialize */ }
    variableflex<1, Double>(Base x) : Base(x) {}
    variableflex<1, Double>(double x) : Base(x) {}
    variableflex<1, Double>(double x, int id) : Base(x) {
      setid(id);
    }
    template<class Constant>
    variableflex<1, Double>(Constant x) {
      Base::value = x;// Base::deriv.setZero();
    }
    template<class Constant>
    variableflex<1, Double>(Constant x, int id) {
      Base::value = x;
      Base::deriv = resize(Base::deriv, std::max((int)Base::deriv.size(),id+1));			   
      setid(id);
    }
    void setid(int i0, int count = 0){
      Base::deriv = resize(Base::deriv, std::max((int)Base::deriv.size(),i0+1));
      if(count == 0)
	this->deriv[i0] = 1.0;
      if(count == 1)
	this->value = 1.0;
    }
    vector<Double> getDeriv(){
      return this->deriv;
    }
  };








#define TMB_BIND_ATOMIC_FLEX_PART(NAME,CALL, LASTVAR)		\
TMB_ATOMIC_VECTOR_FUNCTION(					\
 NAME,								\
    (size_t)							\
  pow((double)							\
      LASTVAR + 1,						\
      CppAD::Integer(tx[tx.size()-1]))				\
  ,								\
 int order = CppAD::Integer(tx[tx.size()-1]);			\
 int nvar = LASTVAR + 1;					\
 atomic::tiny_vec_ref<double> tyref(&ty[0], ty.size());		\
 if(order==0) {							\
   typedef double Float;					\
   CppAD::vector<Float> x(tx);					\
   ty[0] = CALL;						\
 } else if (order==1) {						\
   typedef variableflex<1> Float;				\
   CppAD::vector<Float> x(tx.size());				\
   for(int i = 0; i < (int)tx.size() - 1; ++i){			\
     x[i] = tx[i];						\
     x[i].deriv = resize(x[i].deriv,nvar);			\
     if(i < nvar)						\
       x[i].setid(i);						\
   }								\
   x[x.size()-1] = order;					\
   tyref = CALL.getDeriv();					\
 } else if (order==2) {						\
   typedef variableflex<2> Float;				\
   CppAD::vector<Float> x(tx.size());				\
   for(int i = 0; i < (int)tx.size() - 1; ++i){			\
     x[i] = tx[i];						\
     x[i].deriv = resize(x[i].deriv,nvar);			\
     if(i < nvar)						\
       x[i].setid(i);						\
   }								\
   x[x.size()-1] = order;					\
   tyref = CALL.getDeriv();					\
 } else if (order==3) {						\
   typedef variableflex<3> Float;				\
   CppAD::vector<Float> x(tx.size());				\
   for(int i = 0; i < (int)tx.size() - 1; ++i){			\
     x[i] = tx[i];						\
     x[i].deriv = resize(x[i].deriv,nvar);			\
     if(i < nvar)						\
       x[i].setid(i);						\
   }								\
   x[x.size()-1] = order;					\
   tyref = CALL.getDeriv();					\
 } else {							\
   Rf_error("Order not implemented");				\
 }								\
 ,								\
  int nvar =  LASTVAR + 1;					\
 CppAD::vector<Type> tx_(tx);					\
 tx_[tx.size() - 1] = tx_[tx.size() - 1] + Type(1.0);		\
 vector<Type> tmp = NAME(tx_);					\
 matrix<Type> m = tmp.matrix();					\
 m.resize(nvar, m.size() / nvar);				\
 vector<Type> w = py;						\
 vector<Type> px_ = m * w.matrix();				\
 px = px_;							\
 px[tx.size() - 1] = 0;						\
 )
 



#define TMB_BIND_ATOMIC_FLEX(NAME,CALL)				\
  TMB_BIND_ATOMIC_FLEX_PART(NAME,CALL,tx.size()-2)


#endif

