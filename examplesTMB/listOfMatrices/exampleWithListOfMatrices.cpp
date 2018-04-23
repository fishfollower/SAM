#include <TMB.hpp>
#include <iostream>
using namespace Eigen;
using namespace tmbutils;

using std::cout;
using std::endl;


template<class Type>
struct my_list : vector<matrix<Type> > {

  my_list(SEXP x){ // Constructor
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asMatrix<Type>(sm);
    }
  }
};



template<class Type>
Type objective_function<Type>::operator() ()
{

  // x = list of sparse matrices
  DATA_STRUCT(x, my_list);

  PARAMETER(test);

  std::cout<<"First matrix"<<endl;
  cout<<x(0)<<endl;
  std::cout<<"Second matrix"<<endl;
  cout<<x(1)<<endl;


  Type nll = 0;

  return(nll);

}
