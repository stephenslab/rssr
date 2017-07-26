#ifndef RSSR_TBB_H
#define RSSR_TBB_H

#include <RcppEigen.h>
#include "rssr_types.h"



Rcpp::DataFrame grid_search_rss_varbvsr_tls(
    const  Rcpp::NumericMatrix &SiRiS,
    const Rcpp::NumericVector &sigma_beta,
    const Rcpp::NumericVector &logodds,
    const Rcpp::NumericVector &betahat,
    const Rcpp::NumericVector &se,
    const Rcpp::NumericVector &alpha0,
    const Rcpp::NumericVector &mu0,
    const Rcpp::NumericVector &SiRiSr0,
    const double tolerance,
    const int itermax,
    Rcpp::LogicalVector lnz_tol,
    const int n,const int grainsize);


Rcpp::DataFrame grid_search_rss_varbvsr_norm_tls(const  Rcpp::NumericMatrix &SiRiS,
						 const Rcpp::NumericVector &sigma_beta,
						 const Rcpp::NumericVector &betahat,
						 const Rcpp::NumericVector &se,
						 const Rcpp::NumericVector &mu0,
						 const Rcpp::NumericVector &SiRiSr0,
						 const double tolerance,
						 const int itermax,
						 Rcpp::LogicalVector lnz_tol,
						 const int n,const int grainsize);


Rcpp::DataFrame grid_search_rss_varbvsr_sparse(
    const  Eigen::Map<Eigen::SparseMatrix<double> > SiRiS,
    const Rcpp::NumericVector &sigma_beta,
    const Rcpp::NumericVector &logodds,
    const Rcpp::NumericVector &betahat,
    const Rcpp::NumericVector &se,
    const Rcpp::NumericVector &alpha0,
    const Rcpp::NumericVector &mu0,
    const Rcpp::NumericVector &SiRiSr0,
    const double tolerance,
    const int itermax,
    Rcpp::LogicalVector lnz_tol,
    const int n,const int grainsize);

#endif
