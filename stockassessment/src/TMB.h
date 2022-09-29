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

#include <TMB.hpp>


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

TMB_SPEC(class density::MVNORM_t<double >);
TMB_SPEC(class density::MVNORM_t<TMBad::ad_aug >);

TMB_SPEC(class density::GMRF_t<double >);
TMB_SPEC(class density::GMRF_t<TMBad::ad_aug >);

TMB_SPEC(class density::SCALE_t<GMRF_t<double > >);
TMB_SPEC(class density::SCALE_t<GMRF_t<TMBad::ad_aug > >);

TMB_SPEC(class UNSTRUCTURED_CORR_t<double>);
TMB_SPEC(class UNSTRUCTURED_CORR_t<TMBad::ad_aug>);

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
