#ifndef TMB_H
#define TMB_H

/* #ifdef TMB_PRECOMPILE */
/* #undef WITH_LIBTMB */
/* #undef TMB_PRECOMPILE */
/* #undef CSKIP */
/* #undef IF_TMB_PRECOMPILE */
/* #undef TMB_EXTERN */
/* // Redefine */
/* #undef  WITH_LIBTMB */
/* #define TMB_PRECOMPILE */
/* #define CSKIP(...) __VA_ARGS__ */
/* #define IF_TMB_PRECOMPILE(...) __VA_ARGS__ */
/* #else */
/* /\** \file */
/*     \brief Include this file to extract declarations only */
/* *\/ */
/* #undef WITH_LIBTMB */
/* #undef TMB_PRECOMPILE */
/* #undef CSKIP */
/* #undef IF_TMB_PRECOMPILE */
/* #undef TMB_EXTERN */
/* // Redefine */
/* #define WITH_LIBTMB */
/* #undef  TMB_PRECOMPILE */
/* #define CSKIP(...) ; */
/* #define IF_TMB_PRECOMPILE(...) */
/* #endif */

#ifndef TMB_PRECOMPILE
#ifdef TMB_ALREADY_PRECOMPILED
#undef WITH_LIBTMB
#undef TMB_PRECOMPILE
#undef CSKIP
#undef IF_TMB_PRECOMPILE
#undef TMB_EXTERN
// Redefine
#define WITH_LIBTMB
#include <tmb_enable_header_only.hpp>
#undef  TMB_PRECOMPILE
#define CSKIP(...) ;
#define IF_TMB_PRECOMPILE(...)
#define TMB_EXTERN extern
#endif
#endif

#include <TMB.hpp>


struct Rint {
  int i;
  inline Rint () { }
  inline Rint (double x) : i(IntegerFromReal(x)) { }
  inline operator int() const { return i; }
  // From src/main/coerce.c
  inline int IntegerFromReal(double x) {
    using std::isnan;
    if (ISNAN(x))
      return NA_INTEGER;
    else if (x >= INT_MAX+1. || x <= INT_MIN ) {
      return NA_INTEGER;
    }
    return (int) x;
  }
};

template<>
inline vector<int> asVector(SEXP x) {
  return asVector<Rint>(x).cast<int>();
}

template<>
inline matrix<int> asMatrix(SEXP x) {
  return asMatrix<Rint>(x).cast<int>();
}

namespace tmbutils {
template<>
inline array<int> asArray(SEXP x) {
  array<Rint> tmp = asArray<Rint>(x);
  return array<int>(tmp.vectorcopy.cast<int>(), tmp.dim);
}
}

#ifdef DATA_INTEGER
#undef DATA_INTEGER
#endif
#define DATA_INTEGER(name) int name(asVector<Rint>(     \
getListElement(TMB_OBJECTIVE_PTR -> data,               \
#name, &isNumericScalar))[0]);



#ifdef WITH_LIBTMB
#define TMB_SPEC(...) extern template __VA_ARGS__;
#else
#define TMB_SPEC(...) template __VA_ARGS__;
#endif


/* Specific from Eigen */
using Eigen::SparseMatrix;
using Eigen::Triplet;
using Eigen::LLT;
using Eigen::Dynamic;

TMB_SPEC(class Eigen::SparseMatrix<double>);
TMB_SPEC(class Eigen::SparseMatrix<TMBad::ad_aug>);
TMB_SPEC(class Eigen::Triplet<double>);
TMB_SPEC(class Eigen::Triplet<TMBad::ad_aug>);
/* Specific from Eigen */
TMB_SPEC(class Eigen::LLT< Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >);
TMB_SPEC(class Eigen::LLT< Matrix<TMBad::ad_aug, Eigen::Dynamic, Eigen::Dynamic> >);


/* Specific from TMB */
using density::MVNORM_t;
using density::GMRF_t;
using density::SCALE_t;
using density::UNSTRUCTURED_CORR_t;
using Eigen::SparseMatrix;
using tmbutils::invertSparseMatrix;

TMB_SPEC(matrix<double> tmbutils::invertSparseMatrix(Eigen::SparseMatrix<double> A));
TMB_SPEC(matrix<TMBad::ad_aug> tmbutils::invertSparseMatrix(Eigen::SparseMatrix<TMBad::ad_aug> A));

#ifndef TMB_ALREADY_PRECOMPILED
TMB_SPEC(class density::MVNORM_t<double >);
TMB_SPEC(class density::MVNORM_t<TMBad::ad_aug >);

/* TMB_SPEC(class density::GMRF_t<double >); */
/* TMB_SPEC(class density::GMRF_t<TMBad::ad_aug >); */
#endif

/* TMB_SPEC(class density::SCALE_t<GMRF_t<double > >); */
/* TMB_SPEC(class density::SCALE_t<GMRF_t<TMBad::ad_aug > >); */

TMB_SPEC(class density::UNSTRUCTURED_CORR_t<double>);
TMB_SPEC(class density::UNSTRUCTURED_CORR_t<TMBad::ad_aug>);

TMB_SPEC(struct tmbutils::vector<matrix<double> >);
TMB_SPEC(struct tmbutils::vector<matrix<TMBad::ad_aug> >);
TMB_SPEC(struct tmbutils::vector<int>);
TMB_SPEC(struct tmbutils::vector<double>);
TMB_SPEC(struct tmbutils::vector<TMBad::ad_aug>);

TMB_SPEC(struct tmbutils::array<int>);
TMB_SPEC(struct tmbutils::array<double>);
TMB_SPEC(struct tmbutils::array<TMBad::ad_aug>);

TMB_SPEC(struct tmbutils::matrix<int>);
TMB_SPEC(struct tmbutils::matrix<double>);
TMB_SPEC(struct tmbutils::matrix<TMBad::ad_aug>);

TMB_SPEC(struct data_indicator<vector<double>,double>);
TMB_SPEC(struct data_indicator<vector<TMBad::ad_aug>,TMBad::ad_aug>);

TMB_SPEC(double dnorm(double,double,double,int));
TMB_SPEC(TMBad::ad_aug dnorm(TMBad::ad_aug,TMBad::ad_aug,TMBad::ad_aug,int));



#endif
